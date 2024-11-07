#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <glib.h>
#include "cuda_runtime.h"

#define LINE_SIZE 256

// Define a structure representing (ws_len, zs_len, fftlen)
typedef struct
{
    int x;
    int y;
    int z;
} Triple;

// Hash function
guint triple_hash(gconstpointer key)
{
    const Triple *t = key;
    return g_int_hash(&t->x) ^ g_int_hash(&t->y) ^ g_int_hash(&t->z);
}

// Comparison function
gboolean triple_equal(gconstpointer a, gconstpointer b)
{
    const Triple *t1 = a;
    const Triple *t2 = b;
    return t1->x == t1->x && t1->y == t1->y && t1->z == t1->z;
}

// Global hash table variable
GHashTable *cufftPlan_cache;

// Function to initialize the hash table
void initialize_cache()
{
    cufftPlan_cache = g_hash_table_new(triple_hash, triple_equal);
}

void extractDoubleValue(const char *str, double *value)
{
    char *equalSign = strchr(str, '=');
    if (equalSign != NULL)
    {
        *value = atof(equalSign + 1);
    }
    else
    {
        printf("未找到 '=' 号。\n");
    }
}

/* Stepsize in Fourier Freq */
#define ACCEL_NUMBETWEEN 2
/* Stepsize in Fourier Freq */
#define ACCEL_DR 0.5
/* Reciprocal of ACCEL_DR */
#define ACCEL_RDR 2
/* Stepsize in Fourier F-dot */
#define ACCEL_DZ 2
/* Reciprocal of ACCEL_DZ */
#define ACCEL_RDZ 0.5
/* Stepsize in Fourier F-dot-dot */
#define ACCEL_DW 20
/* Reciprocal of ACCEL_DW */
#define ACCEL_RDW 0.05

/*  Number of bins on each side of the central frequency  */
/*    to sum for Fourier interpolation (low accuracy)     */
#define NUMFINTBINS 16

#define DBLCORRECT 1e-14

#define NEAREST_INT(x) (int)(x < 0 ? x - 0.5 : x + 0.5)

/* Return 2**n */
#define index_to_twon(n) (1 << n)

// #define MAX(a, b) ((a) > (b) ? (a) : (b))

/*  Constants used in the interpolation routines */
typedef enum
{
    LOWACC,
    HIGHACC
} presto_interp_acc;

typedef struct subharminfo
{
    int numharm;      /* The number of sub-harmonics */
    int harmnum;      /* The sub-harmonic number (fundamental = numharm) */
    int zmax;         /* The maximum Fourier f-dot for this harmonic */
    int wmax;         /* The maximum Fourier f-dot-dot for this harmonic */
    int numkern_zdim; /* Number of kernels calculated in the z dimension */
    int numkern_wdim; /* Number of kernels calculated in the w dimension */
    int numkern;      /* Total number of kernels in the vector */
    int fftlen;
} subharminfo;

// 存储多个 subharmonic 的结构体
typedef struct
{
    int harm_fract;
    int subharmonic_wlo;
    unsigned short *subharmonic_zinds;
    unsigned short *subharmonic_rinds;
    float *subharmonic_powers;
    int subharmonic_numzs;
    int subharmonic_numrs;
} SubharmonicMap;

FILE *chkfopen(char *path, const char *mode)
{
    FILE *file;

    if ((file = fopen(path, mode)) == NULL)
    {
        perror("\nError in chkfopen()");
        printf("   path = '%s'\n", path);
        exit(-1);
    }
    return (file);
}

long long chkfilelen(FILE *file, size_t size)
{
    int filenum, rt;
    struct stat buf;

    filenum = fileno(file);
    rt = fstat(filenum, &buf);
    if (rt == -1)
    {
        perror("\nError in chkfilelen()");
        printf("\n");
        exit(-1);
    }
    return (long long)(buf.st_size / size);
}

static inline double calc_required_r(double harm_fract, double rfull)
{
    return rint(ACCEL_RDR * rfull * harm_fract) * ACCEL_DR;
}

static inline int calc_required_z(double harm_fract, double zfull)
{
    return NEAREST_INT(ACCEL_RDZ * zfull * harm_fract) * ACCEL_DZ;
}

static inline int calc_required_w(double harm_fract, double wfull)
{
    return NEAREST_INT(ACCEL_RDW * wfull * harm_fract) * ACCEL_DW;
}

int z_resp_halfwidth(double z, presto_interp_acc accuracy)
/*  Return the approximate kernel half width in FFT bins required    */
/*  to achieve a fairly high accuracy correlation based correction   */
/*  or interpolation for a Fourier signal with constant f-dot. (i.e  */
/*  a constant frequency derivative)                                 */
/*  Arguments:                                                       */
/*    'z' is the Fourier Frequency derivative (# of bins the signal  */
/*       smears over during the observation).                        */
/*    'accuracy' is either LOWACC or HIGHACC.                        */
/*  Notes:                                                           */
/*    The result must be multiplied by 2*'numbetween' to get the     */
/*    length of the array required to hold such a kernel.            */
{
    int m = (int)(0.5 * 1.1 * fabs(z));
    if (accuracy == HIGHACC)
    {
        m += NUMFINTBINS * 3;
    }
    else
    {
        m += NUMFINTBINS;
    }
    return m;
}

int w_resp_halfwidth(double z, double w, presto_interp_acc accuracy)
/*  Return the approximate kernel half width in FFT bins required    */
/*  to achieve a fairly high accuracy correlation based correction   */
/*  or interpolation for a Fourier signal with an f-dot that (i.e    */
/*  varies linearly in time -- a constant f-dotdot)                  */
/*  Arguments:                                                       */
/*    'z' is the average Fourier Frequency derivative (# of bins     */
/*       the signal smears over during the observation).             */
/*    'w' is the Fourier Frequency 2nd derivative (change in the     */
/*       Fourier f-dot during the observation).                      */
/*    'accuracy' is either LOWACC or HIGHACC.                        */
/*  Notes:                                                           */
/*    The result must be multiplied by 2*'numbetween' to get the     */
/*    length of the array required to hold such a kernel.            */
{
    if (fabs(w) < 1.0e-7)
        return z_resp_halfwidth(z, accuracy);
    double r0 = 0.5 * (w / 6.0 - z); // Starting deviation from r_avg
    double r1 = 0.5 * (w / 6.0 + z); // Ending deviation from r_avg
    // We need to know the maximum deviation from r_avg
    double maxdev = fabs(r0) > fabs(r1) ? fabs(r0) : fabs(r1);
    // If the extrema of the parabola is within 0 < u < 1, then
    // it will be a new freq minimum or maximum
    double u_ext = 0.5 - z / w;
    if (u_ext > 0.0 && u_ext < 1.0)
    {
        double z0 = z - w / 2.0; // Starting z
        // Value of r at the extremum
        double r_ext = 0.5 * w * u_ext * u_ext + z0 * u_ext + r0;
        maxdev = fabs(r_ext) > maxdev ? fabs(r_ext) : maxdev;
    }
    if (accuracy == HIGHACC)
    {
        return (int)(1.1 * maxdev) + NUMFINTBINS * 3;
    }
    else
    {
        return (int)(1.1 * maxdev) + NUMFINTBINS;
    }
}

long long next2_to_n(long long x)
/* Return the first value of 2^n >= x */
{
    long long i = 1;

    while (i < x)
        i <<= 1;
    return i;
}

int fftlen_from_kernwidth(int kernwidth)
/* return the length of the optimal FFT to use for correlations with
 * some kernel width kernwidth.  This assumes FFTW. */
{
    // The following nummbers were determined using FFTW 3.3.7 on an
    // AVX-enabled processor.  Metric used was max throughput of good
    // correlated data.
    if (kernwidth < 6)
        return 128;
    else if (kernwidth < 52)
        return 256;
    else if (kernwidth < 67)
        return 512;
    else if (kernwidth < 378)
        return 1024;
    else if (kernwidth < 664)
        return 2048;
    else if (kernwidth < 1672)
        return 4096;
    else if (kernwidth < 3015)
        return 10240;
    else if (kernwidth < 3554)
        return 15360;
    else if (kernwidth < 6000)
        return 25600;
    else
        return next2_to_n(kernwidth * 5);
}

int next_good_fftlen(int N)
/* Return one of the shortest, yet best performing, FFT lengths larger
 * than N.  This assumes FFTW. */
{
    int fftlens[17] = {128, 192, 256, 384, 512, 768, 1024, 1280, 2048, 4096,
                       5120, 7680, 10240, 12288, 15360, 16384, 25600};
    int ii = 0;
    if (N <= fftlens[0])
        return fftlens[0];
    if (N > fftlens[16])
        return next2_to_n(N);
    while (N > fftlens[ii])
        ii++;
    return fftlens[ii];
}

/* Return x such that 2**x = n */
static inline int twon_to_index(int n)
{
    int x = 0;

    while (n > 1)
    {
        n >>= 1;
        x++;
    }
    return x;
}

static int calc_fftlen(int numharm, int harmnum, int max_zfull, int max_wfull, int obs_corr_uselen)
/* The fft length needed to properly process a subharmonic */
{
    int bins_needed, end_effects;
    double harm_fract;

    harm_fract = (double)harmnum / (double)numharm;
    bins_needed = (int)ceil(obs_corr_uselen * harm_fract) + 2;
    end_effects = 2 * ACCEL_NUMBETWEEN *
                  w_resp_halfwidth(calc_required_z(harm_fract, max_zfull),
                                   calc_required_w(harm_fract, max_wfull), LOWACC);
    return next_good_fftlen(bins_needed + end_effects);
}

static void init_subharminfo(int numharm, int harmnum, int zmax, int wmax, subharminfo *shi, int obs_fftlen, int obs_corr_uselen)
{
    int ii, jj, fftlen;
    double harm_fract;

    harm_fract = (double)harmnum / (double)numharm;
    shi->numharm = numharm;
    shi->harmnum = harmnum;
    shi->zmax = calc_required_z(harm_fract, zmax);
    shi->wmax = calc_required_w(harm_fract, wmax);
    if (numharm == 1 && harmnum == 1)
        fftlen = obs_fftlen;
    else
        fftlen = calc_fftlen(numharm, harmnum, zmax, wmax, obs_corr_uselen);
    shi->fftlen = fftlen;
    shi->numkern_zdim = (shi->zmax / ACCEL_DZ) * 2 + 1;
    shi->numkern_wdim = (shi->wmax / ACCEL_DW) * 2 + 1;
    shi->numkern = shi->numkern_zdim * shi->numkern_wdim;
}

typedef float rawtype_part;
typedef struct FCOMPLEX
{
    rawtype_part r, i;
} fcomplex;

subharminfo **create_subharminfos(int zmax, int wmax, int obs_fftlen, int obs_numharmstages, int obs_corr_uselen)
{
    double kern_ram_use = 0;
    int ii, jj, harmtosum, fftlen;
    subharminfo **shis;

    shis = (subharminfo **)malloc(obs_numharmstages * sizeof(subharminfo *));
    /* Prep the fundamental (actually, the highest harmonic) */
    shis[0] = (subharminfo *)malloc(2 * sizeof(subharminfo));
    init_subharminfo(1, 1, (int)zmax, (int)wmax, &shis[0][0], obs_fftlen, obs_corr_uselen);
    kern_ram_use += shis[0][0].numkern * fftlen * sizeof(fcomplex); // in Bytes
                                                                    /* Prep the sub-harmonics if needed */
    for (ii = 1; ii < obs_numharmstages; ii++)
    {
        harmtosum = index_to_twon(ii);
        shis[ii] = (subharminfo *)malloc(harmtosum * sizeof(subharminfo));
        for (jj = 1; jj < harmtosum; jj += 2)
        {
            init_subharminfo(harmtosum, jj, (int)zmax,
                             (int)wmax, &shis[ii][jj - 1], obs_fftlen, obs_corr_uselen);
        }
    }
    return shis;
}

long long get_fkern_size(subharminfo **shis, int obs_numharmstages)
{
    long long kern_ram_use = 0;
    int ii, jj, harmtosum;
    kern_ram_use += shis[0][0].numkern * shis[0][0].fftlen * sizeof(fcomplex); // in Bytes
    for (ii = 1; ii < obs_numharmstages; ii++)
    {
        harmtosum = index_to_twon(ii);
        for (jj = 1; jj < harmtosum; jj += 2)
        {
            kern_ram_use += shis[ii][jj - 1].numkern * shis[ii][jj - 1].fftlen * sizeof(fcomplex); // in Bytes
        }
    }
    return kern_ram_use;
}

long long get_init_constant_device_size(int obs_numharmstages, int wmax)
{
    int single_loop = 1 << (obs_numharmstages - 1);
    int shi_numkern_wdim = ((calc_required_w(1, wmax)) / ACCEL_DW) * 2 + 1;
    int fundamental_numws = shi_numkern_wdim;
    int subw_size = single_loop * fundamental_numws;
    int numharmstages_size = obs_numharmstages;
    long long init_constant_device_size = subw_size * sizeof(int) + numharmstages_size * sizeof(float) + numharmstages_size * sizeof(int) + numharmstages_size * sizeof(double);
    return init_constant_device_size;
}

long long get_full_tmpdat_array_size(subharminfo **subharminfs, int batchsize, int obs_numharmstages)
{
    int ws_len_max = 0, zs_len_max = 0, fft_len_max = 0;
    subharminfo *shii = &subharminfs[0][0];
    ws_len_max = shii->numkern_wdim;
    zs_len_max = shii->numkern_zdim;
    fft_len_max = shii->fftlen;
    Triple t1 = {shii->numkern_wdim, shii->numkern_zdim, shii->fftlen};
    if (!g_hash_table_contains(cufftPlan_cache, &t1))
    {
        g_hash_table_insert(cufftPlan_cache, g_memdup(&t1, sizeof(Triple)), NULL);
    }

    int stage, harmtosum, harm;

    for (stage = 1; stage < obs_numharmstages; stage++)
    {
        harmtosum = 1 << stage;
        for (harm = 1; harm < harmtosum; harm += 2)
        {
            subharminfo *shi = &subharminfs[stage][harm - 1];
            fft_len_max = MAX(shi->fftlen, fft_len_max);
            zs_len_max = MAX(shi->numkern_zdim, zs_len_max);
            ws_len_max = MAX(shi->numkern_wdim, ws_len_max);
            Triple t2 = {shi->numkern_wdim, shi->numkern_zdim, shi->fftlen};
            if (!g_hash_table_contains(cufftPlan_cache, &t2))
            {
                g_hash_table_insert(cufftPlan_cache, g_memdup(&t2, sizeof(Triple)), NULL);
            }
        }
    }

    long long array_size = batchsize * sizeof(fcomplex) * fft_len_max * ws_len_max * zs_len_max; // 最大不超过 batch_size_max
    return array_size;
}

double get_rlo(char *filename)
{
    // Copy the original filename to a newly allocated string
    char *infname = malloc(strlen(filename) + 1);
    strcpy(infname, filename);

    // Find the position of the last '.' character
    char *dot = strrchr(infname, '.');
    if (dot != NULL)
    {
        // Replace the extension with 'inf'
        strcpy(dot + 1, "inf");
    }
    else
    {
        // If no '.' character is found, append ".inf" directly to the string
        strcat(infname, ".inf");
    }

    FILE *file;
    char line[LINE_SIZE];
    double value10 = 0.0, value11 = 0.0;
    int lineNumber = 0;

    file = fopen(infname, "r");
    if (file == NULL)
    {
        printf("cannot open file: %s\n", infname);
        return 1;
    }

    while (fgets(line, sizeof(line), file) != NULL)
    {
        lineNumber++;
        if (lineNumber == 10)
        {
            extractDoubleValue(line, &value10);
        }
        if (lineNumber == 11)
        {
            extractDoubleValue(line, &value11);
        }
    }

    fclose(file);

    double rlo = floor(value10 * value11);
    return rlo;
}

long long single_powers_size(int numharm, int harmnum, double fullrlo, double fullrhi, subharminfo *shi, int obs_corr_uselen)
{
    double drlo, drhi, harm_fract;
    /* Calculate and get the required amplitudes */
    harm_fract = (double)harmnum / (double)numharm;
    drlo = calc_required_r(harm_fract, fullrlo);
    drhi = calc_required_r(harm_fract, fullrhi);
    int numrs = (int)((ceil(drhi) - floor(drlo)) * ACCEL_RDR + DBLCORRECT) + 1;
    if (numharm == 1 && harmnum == 1)
    {
        numrs = obs_corr_uselen;
    }
    else
    {
        if (numrs % ACCEL_RDR)
            numrs = (numrs / ACCEL_RDR + 1) * ACCEL_RDR;
    }
    int numzs = shi->numkern_zdim;
    int numws = shi->numkern_wdim;
    long long single_size = (numrs * numzs * numws) * sizeof(float);
    return single_size;
}

long long get_powers_size(char *filename, subharminfo **subharminfs, int obs_corr_uselen, long long obs_highestbin, int obs_numharmstages)
{
    double rlo = get_rlo(filename);
    int rstep = obs_corr_uselen * ACCEL_DR;
    double startr = rlo;
    double lastr = 0;
    double nextr = 0;
    int stage, harmtosum, harm;
    long long max_powers_size = 0;
    while (startr + rstep < obs_highestbin)
    {
        long long loop_powers_size = 0;
        nextr = startr + rstep;
        lastr = nextr - ACCEL_DR;
        long long fundamental_size = single_powers_size(1, 1, startr, lastr, &subharminfs[0][0], obs_corr_uselen);
        loop_powers_size += fundamental_size;
        if (obs_numharmstages > 1)
        {
            for (stage = 1; stage < obs_numharmstages; stage++)
            {
                harmtosum = 1 << stage;
                for (harm = 1; harm < harmtosum; harm += 2)
                {
                    subharminfo *shi = &subharminfs[stage][harm - 1];
                    long long subharmonic_size = single_powers_size(harmtosum, harm, startr, lastr, &subharminfs[stage][harm - 1], obs_corr_uselen);
                    loop_powers_size += subharmonic_size;
                }
            }
        }
        startr = nextr;
        max_powers_size = MAX(max_powers_size, loop_powers_size);
    }
    return max_powers_size;
}

long long get_cufftPlan_size()
{
    long long cufftPlan_size = 0;

    // Manually iterate over the keys of the hash table
    GHashTableIter iter;
    gpointer key, value;
    g_hash_table_iter_init(&iter, cufftPlan_cache);
    while (g_hash_table_iter_next(&iter, &key, &value))
    {
        Triple *t = (Triple *)key;
        double sum = (t->x * t->y * t->z) * sizeof(fcomplex);
        cufftPlan_size += sum;
    }

    // Release the hash table
    g_hash_table_destroy(cufftPlan_cache);
    cufftPlan_size = cufftPlan_size;
    printf("Total cufftPlan size:%.3f MB\n", cufftPlan_size / (1024.0 * 1024.0));
    return cufftPlan_size;
}

void accelsearch_maxsize(int arg1, int arg2, int arg3, int arg4, char* arg5){
    subharminfo **subharminfs;
    int zmax = arg1;
    int wmax = arg2;
    int numharm = arg3;
    int batchsize = arg4;
    char *filename = arg5;

    // 初始化全局哈希表
    initialize_cache();
    FILE *fftfile = chkfopen(filename, "rb");
    long long total_size = 0;

    // size_t before_mem, after_mem, total_mem;
    // cudaMemGetInfo(&before_mem, &total_mem);
    // long long cuda_init_size = total_mem - before_mem;
    long long cuda_init_size = 299892736; //286 MB
    // printf("GPU init mem:%.3f MB\n", cuda_init_size / (1024.0 * 1024.0));
    // printf("GPU init mem:%lld\n", cuda_init_size);
    total_size += cuda_init_size;

    long long obs_numbins = chkfilelen(fftfile, sizeof(fcomplex));
    long long obs_highestbin = obs_numbins - 1;
    int obs_maxkernlen = 2 * ACCEL_NUMBETWEEN * w_resp_halfwidth(zmax, wmax, LOWACC);
    int obs_fftlen = fftlen_from_kernwidth(obs_maxkernlen);
    if (obs_fftlen < 2048)
        obs_fftlen = 2048; // This gives slightly better speed empirically
    int obs_corr_uselen = obs_fftlen - obs_maxkernlen;
    if (obs_corr_uselen % ACCEL_RDR)
        obs_corr_uselen = obs_corr_uselen / ACCEL_RDR * ACCEL_RDR;
    int obs_numharmstages = twon_to_index(numharm) + 1;

    subharminfs = create_subharminfos(zmax, wmax, obs_fftlen, obs_numharmstages, obs_corr_uselen);

    long long kern_ram_use = get_fkern_size(subharminfs, obs_numharmstages);
    total_size += kern_ram_use;
    long long init_constant_device_size = get_init_constant_device_size(obs_numharmstages, wmax);
    total_size += init_constant_device_size;
    long long full_tmpdat_array_size = get_full_tmpdat_array_size(subharminfs, batchsize, obs_numharmstages);
    total_size += full_tmpdat_array_size;

    long long subharmonics_add_size = (1 << (obs_numharmstages - 1)) * batchsize * sizeof(SubharmonicMap);
    total_size += subharmonics_add_size;

    long long search_results_size;
    if (batchsize > 8)
        search_results_size = 16 * full_tmpdat_array_size / batchsize;
    else
        search_results_size = 8 * full_tmpdat_array_size / batchsize;
    total_size += search_results_size;

    long long max_powers_size = batchsize * get_powers_size(filename, subharminfs, obs_corr_uselen, obs_highestbin, obs_numharmstages);

    // long long cufftPlan_size = get_cufftPlan_size();
    // total_size += cufftPlan_size;

    total_size += max_powers_size;
    total_size += (size_t)(obs_corr_uselen * sizeof(unsigned short) * batchsize * 2) * 15 + (size_t)(sizeof(fcomplex) * 4096 * batchsize) * 16 + batchsize * sizeof(int) * 2;
    total_size = total_size*1.1;
    double td = total_size / (1024.0 * 1024.0);
    printf("%.3f\n", td);
}
int main(int argc, char *argv[])
{
    // cudaFree(0);
    accelsearch_maxsize(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), argv[5]);
    return 0;
}