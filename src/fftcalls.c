#include "ransomfft.h"
#include "stdint.h"
#include "string.h"
#include <errno.h>
#include <pthread.h>

// Following gives the same as FFTW's fftwf_alignment_of when
// BYTE_COUNT = 16, which is what we need for SSE.
// 0 means that it is aligned on BYTE_COUNT boundaries
#define is_aligned(POINTER, BYTE_COUNT) \
    ((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT)

void read_wisdom(void)
{
    FILE *wisdomfile;
    static char wisdomfilenm[120];

    /* First try to import the system wisdom if available */
    fftwf_import_system_wisdom();
    sprintf(wisdomfilenm, "%s/lib/fftw_wisdom.txt", getenv("PRESTO"));
    wisdomfile = fopen(wisdomfilenm, "r");
    if (wisdomfile == NULL) {
        printf("Warning:  Couldn't open '%s'\n"
               "          You should run 'makewisdom'.  See $PRESTO/INSTALL.\n",
               wisdomfilenm);
    } else {
        if (!fftwf_import_wisdom_from_file(wisdomfile))
            printf("Warning:  '%s' is not up-to-date.\n"
                   "          You should run 'makewisdom'.  See $PRESTO/INSTALL.\n",
                   wisdomfilenm);
        fclose(wisdomfile);
    }
    // The following resets errno if one of the wisdom files was not found
    errno = 0;
}


void fftwcallsimple(fcomplex * data, long nn, int isign)
/* Simple FFTW calling function for testing */
{
    static int firsttime = 1;
    fftwf_plan plan;
    if (firsttime) {
        read_wisdom();
        firsttime = 0;
    }
    // Note: We need to use FFTW_ESTIMATE since other
    // plan-making destroys the input and output arrays
    plan = fftwf_plan_dft_1d(nn, (fftwf_complex *) data,
                             (fftwf_complex *) data, isign, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
}

// 定义缓存大小的宏
#define CACHE_SIZE 256

// 全局互斥锁
pthread_mutex_t fftw_lock = PTHREAD_MUTEX_INITIALIZER;

// 定义变量
fftwf_plan plancache_forward[CACHE_SIZE] = {NULL};
fftwf_plan plancache_inverse[CACHE_SIZE] = {NULL};
int aligncache[CACHE_SIZE] = {-99};
int firsttime = 1;
int lastslot = 0, lastused[CACHE_SIZE] = {0};
long nncache[CACHE_SIZE] = {0};
unsigned int planflag;
fcomplex *datacopy;

void fftwcall(fcomplex *indata, long nn, int isign) {
    fftwf_complex *dataptr = (fftwf_complex *)indata;
    int ii, indata_align, slot, incache = 0, oldestplan = 0;

    indata_align = is_aligned(indata, 16);

    if (firsttime) {
        pthread_mutex_lock(&fftw_lock);
        read_wisdom();
        pthread_mutex_unlock(&fftw_lock);
        firsttime = 0;
    }

    ii = 0;
    slot = lastslot;
    while (ii < CACHE_SIZE) {
        if (nn == nncache[slot] && indata_align == aligncache[slot]) {
            incache = 1;
            break;
        }
        slot = (slot + 1) % CACHE_SIZE;
        ii++;
    }

    if (!incache) {
        pthread_mutex_lock(&fftw_lock);
        for (ii = CACHE_SIZE - 1; ii >= 0; ii--) {
            if (lastused[ii] > lastused[oldestplan]) {
                oldestplan = ii;
            }
        }

        if (plancache_forward[oldestplan])
            fftwf_destroy_plan(plancache_forward[oldestplan]);
        if (plancache_inverse[oldestplan])
            fftwf_destroy_plan(plancache_inverse[oldestplan]);

        planflag = (nn > 16384) ? FFTW_ESTIMATE : FFTW_MEASURE;
        datacopy = gen_cvect(nn);
        memcpy(datacopy, dataptr, nn * sizeof(fcomplex));

        plancache_forward[oldestplan] = fftwf_plan_dft_1d(nn, dataptr, dataptr, -1, planflag);
        plancache_inverse[oldestplan] = fftwf_plan_dft_1d(nn, dataptr, dataptr, +1, planflag);

        memcpy(dataptr, datacopy, nn * sizeof(fcomplex));
        vect_free(datacopy);
        nncache[oldestplan] = nn;
        aligncache[oldestplan] = indata_align;
        lastused[oldestplan] = 0;
        slot = oldestplan;
        pthread_mutex_unlock(&fftw_lock);
    }

    lastused[slot] = 0;
    for (int j = 1; j < CACHE_SIZE; j++) {
        lastused[(slot + j) % CACHE_SIZE]++;
    }
    lastslot = slot;

    if (isign == -1) {
        fftwf_execute_dft(plancache_forward[slot], dataptr, dataptr);
    } else {
        fftwf_execute_dft(plancache_inverse[slot], dataptr, dataptr);
    }
}