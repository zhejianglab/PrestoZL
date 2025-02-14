/*
 * Copyright (c) 2024 Zhejiang Lab
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "accel.h"
#include <sys/time.h>
#include <math.h>
#include <assert.h>
/*#undef USEMMAP*/

#ifdef USEMMAP
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#include "cuda_runtime.h"
#include "cuda_helper.h"

// Use OpenMP
#ifdef _OPENMP
#include <omp.h>
extern void set_openmp_numthreads(int numthreads);
#endif

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define NEAREST_INT(x) (int)(x < 0 ? x - 0.5 : x + 0.5)

void free_ffdotpows_cu_batch(ffdotpows_cu *ffd_array, int batch_size,
                             cudaStream_t sub_stream);

// Structures to store multiple subharmonics
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

subharminfo **create_subharminfos_cu(accelobs *obs);

typedef struct
{
    long long index;
    float pow;
    float sig;
} SearchValue;

void fuse_add_search_batch(ffdotpows_cu *fundamentals,
                           SubharmonicMap *subhmap,
                           int stages,
                           int fundamental_num,
                           cudaStream_t stream,
                           SearchValue *search_results,
                           unsigned long long int *search_nums,
                           long long pre_size,
                           int proper_batch_size,
                           int max_searchnum,
                           int *too_large);

GSList *insert_to_cands(
    int fundamental_numrs,
    int fundamental_numzs,
    int fundamental_numws,
    long long *fundamental_rlos,
    int fundamental_zlo,
    int fundamental_wlo,
    int proper_batch_size,
    double *numindeps,
    GSList *cands,
    SearchValue *search_results,
    unsigned long long int *search_num,
    long long single_batch_size,
    int numharmstages,
    cudaStream_t main_stream,
    cudaStream_t sub_stream);

void sort_search_results(SearchValue *search_results, unsigned long long int search_num);
void clear_cache();
void subharm_fderivs_vol_cu_batch(
    ffdotpows_cu *ffdot_array,
    int numharm,
    int harmnum,
    double *fullrlo_array,
    double *fullrhi_array,
    subharminfo *shi,
    accelobs *obs,
    cudaStream_t stream,
    fcomplex *full_tmpdat_array,
    fcomplex *full_tmpout_array,
    int batch_size,
    fcomplex *fkern,
    int inds_idx);

void init_inds_array(int size);

fcomplex *fkern_host_to_dev(subharminfo **subharminfs, int numharmstages, int **offset_array);
void init_constant_device(int *subw_host, int subw_size, float *powcuts_host, int *numharms_host, double *numindeps_host, int numharmstages_size);

extern float calc_median_powers(fcomplex *amplitudes, int numamps);
extern void zapbirds(double lobin, double hibin, FILE *fftfile, fcomplex *fft);

long long timeInMilliseconds(void);

static inline int readenv_atoi(char *env)
{
    char *p;
    if ((p = getenv(env)))
        return (atoi(p));
    else
        return (0);
}

static void print_percent_complete(int current, int number, char *what, int reset)
{
    static int newper = 0, oldper = -1;

    if (reset)
    {
        oldper = -1;
        newper = 0;
    }
    else
    {
        newper = (int)(current / (float)(number) * 100.0);
        if (newper < 0)
            newper = 0;
        if (newper > 100)
            newper = 100;
        if (newper > oldper)
        {
            printf("\rAmount of %s complete = %3d%%", what, newper);
            fflush(stdout);
            oldper = newper;
        }
    }
}

// structures for Map key
typedef struct
{
    int harmtosum;
    int harm;
} MapKey;

typedef struct
{
    float startr;
    int harmtosum;
    int harm;
} StartrHarmTuple;

// structures for Map value
typedef struct
{
    // startr batch
    double *startr_array;
    double *lastr_array;
    // shi
    subharminfo *shi;
    // StartrHarmTuple idx batch
    StartrHarmTuple *tuple_array;
    int count;
} MapValue;

// structures for Map Entry
typedef struct
{
    MapKey key;
    MapValue value;
    ffdotpows_cu *subharmonics_batch; // Compute the array pointer for the resulting batch, and store it here for convenient release later
} MapEntry;

// Compare whether the two keys are equal
int compareKeys(MapKey *key1, MapKey *key2)
{
    return key1->harmtosum == key2->harmtosum &&
           key1->harm == key2->harm;
}

// search map
MapEntry *getMap(MapEntry *map, int map_size, MapKey key)
{
    for (int i = 0; i < map_size; i++)
    {
        if (compareKeys(&key, &map[i].key))
        {
            return &map[i];
        }
    }
    return; // If no matching key is found, return NULL
}

void freeMap(MapEntry *map, int *map_size)
{
    for (int i = 0; i < *map_size; i++)
    {
        free(map[i].value.startr_array);
        free(map[i].value.lastr_array);
        free(map[i].value.tuple_array);
    }
}

// insert or update map
void insertOrUpdateMap(MapEntry *map, int *map_size, MapKey key, float startr, float lastr, subharminfo *shi, int harmtosum, int harm, int max_map_size)
{
    // Search the map for a matching key
    for (int i = 0; i < *map_size; i++)
    {
        if (compareKeys(&key, &map[i].key))
        {
            // if a matching key is found, update the entry
            map[i].value.startr_array[map[i].value.count] = startr;
            map[i].value.lastr_array[map[i].value.count] = lastr;
            map[i].value.tuple_array[map[i].value.count] = (StartrHarmTuple){startr, harmtosum, harm};
            map[i].value.count += 1;
            return;
        }
    }
    // if no matching key is found, add a new entry
    if (*map_size < max_map_size)
    {
        map[*map_size].key = key;
        map[*map_size].value.startr_array[0] = startr;
        map[*map_size].value.lastr_array[0] = lastr;
        map[*map_size].value.shi = shi;
        map[*map_size].value.tuple_array[0] = (StartrHarmTuple){startr, harmtosum, harm};
        map[*map_size].value.count += 1;
        (*map_size)++;
    }
    else
    {
        // handle the case of map overflow
        printf("Map overflow\n");
    }
}

#define FLOAT_EPSILON 1e-6 // this value can be adjusted as needed

// print the map to a file
void printMapToFile(const char *filename, MapEntry map[], int map_size)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        perror("Error opening file");
        return;
    }

    for (int i = 0; i < map_size; i++)
    {
        fprintf(file, "MapEntry %d:\n", i);

        fprintf(file, "  Key: harmtosum=%d, harm=%d, shi=%p\n",
                map[i].key.harmtosum, map[i].key.harm, (void *)map[i].value.shi);

        fprintf(file, "  Values:\n");
        for (int j = 0; j < map[i].value.count; j++)
        {
            fprintf(file, "    startr: %f, lastr: %f, tuple: (startr: %f, harmtosum: %d, harm: %d)\n",
                    map[i].value.startr_array[j], map[i].value.lastr_array[j],
                    map[i].value.tuple_array[j].startr, map[i].value.tuple_array[j].harmtosum, map[i].value.tuple_array[j].harm);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

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

int compare_ffdotpows_cu(ffdotpows_cu *a, ffdotpows_cu *b)
{
    if (a->rlo != b->rlo || a->zlo != b->zlo || a->wlo != b->wlo ||
        a->numrs != b->numrs || a->numzs != b->numzs || a->numws != b->numws)
    {
        return 1; // Basic type fields not match
    }

    for (int i = 0; i < a->numrs; ++i)
    {
        if (a->rinds[i] != b->rinds[i])
        {
            return 1; // rinds not match
        }
    }

    for (int i = 0; i < a->numzs; ++i)
    {
        if (a->zinds[i] != b->zinds[i])
        {
            return 1; // zinds not match
        }
    }

    int totalPowers = a->numrs * a->numzs * a->numws;
    float *a_powers_host = (float *)malloc(totalPowers * sizeof(float));
    float *b_powers_host = (float *)malloc(totalPowers * sizeof(float));

    if (a_powers_host == NULL || b_powers_host == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        free(a_powers_host);
        free(b_powers_host);
        return 1;
    }

    // Copy from GPU memory to Host memory
    cudaMemcpy(a_powers_host, a->powers, totalPowers * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(b_powers_host, b->powers, totalPowers * sizeof(float), cudaMemcpyDeviceToHost);

    for (int i = 0; i < totalPowers; ++i)
    {
        if (a_powers_host[i] != b_powers_host[i])
        {
            free(a_powers_host);
            free(b_powers_host);
            return 1; // powers not match
        }
    }

    free(a_powers_host);
    free(b_powers_host);
    return 0; // All fields match
}

void accelsearch_CPU1(int argc, char *argv[], subharminfo ***subharminfs_ptr, accelobs *obs_ptr, infodata *idata_ptr, Cmdline **cmd_ptr)
{
    int ii;
    subharminfo **subharminfs;
    accelobs obs;
    infodata idata;
    Cmdline *cmd;

    // cudaFree(0);
    /* Call usage() if we have no command line arguments */

    if (argc == 1)
    {
        Program = argv[0];
        printf("\n");
        usage();
        exit(1);
    }

    /* Parse the command line using the excellent program Clig */

    cmd = parseCmdline(argc, argv);
    *cmd_ptr = cmd;

#ifdef DEBUG
    showOptionValues();
#endif

    /* Create the accelobs structure */
    create_accelobs(&obs, &idata, cmd, 1);
    *obs_ptr = obs;
    *idata_ptr = idata;

    /* Zap birdies if requested and if in memory */
    if (cmd->zaplistP && !obs.mmap_file && obs.fft)
    {
        int numbirds;
        double *bird_lobins, *bird_hibins, hibin;

        /* Read the Standard bird list */
        numbirds = get_birdies(cmd->zaplist, obs.T, cmd->baryv,
                               &bird_lobins, &bird_hibins);

        /* Zap the birdies */
        printf("Zapping them using a barycentric velocity of %.5gc.\n\n",
               cmd->baryv);
        hibin = obs.N / 2;
        for (ii = 0; ii < numbirds; ii++)
        {
            if (bird_lobins[ii] >= hibin)
                break;
            if (bird_hibins[ii] >= hibin)
                bird_hibins[ii] = hibin - 1;
            zapbirds(bird_lobins[ii], bird_hibins[ii], NULL, obs.fft);
        }

        vect_free(bird_lobins);
        vect_free(bird_hibins);
    }

    /* Generate the correlation kernels */

    printf("\nGenerating correlation kernels:\n");
    subharminfs = create_subharminfos(&obs);
    *subharminfs_ptr = subharminfs;

    if (cmd->ncpus > 1)
    {
#ifdef _OPENMP
        set_openmp_numthreads(cmd->ncpus);
#endif
    }
    else
    {
#ifdef _OPENMP
        omp_set_num_threads(1); // Explicitly turn off OpenMP
#endif
        // printf("Starting the search.\n\n");
    }

    /* Don't use the *.txtcand files on short in-memory searches */
    if (!obs.dat_input)
    {
        printf("  Working candidates in a test format are in '%s'.\n\n",
               obs.workfilenm);
    }

    /* Function pointers to make code a bit cleaner */
    void (*fund_to_ffdot)() = NULL;
    void (*inmem_add_subharm)() = NULL;
    if (obs.inmem)
    {
        if (cmd->otheroptP)
        {
            fund_to_ffdot = &fund_to_ffdotplane_trans;
            inmem_add_subharm = &inmem_add_ffdotpows_trans;
        }
        else
        {
            fund_to_ffdot = &fund_to_ffdotplane;
            inmem_add_subharm = &inmem_add_ffdotpows;
        }
    }
    else
    {
        if (cmd->otheroptP)
        {
            printf("otheroptP is not supported\n");
            exit(-1);
        }
    }
}

int accelsearch_GPU(accelobs obs, subharminfo **subharminfs, GSList **cands_ptr, Cmdline *cmd)
{
    int ii, rstep;
    GSList *cands = NULL;
    /* The step-size of blocks to walk through the input data */
    rstep = obs.corr_uselen * ACCEL_DR;
    // Create cuda stream
    cudaStream_t main_stream, sub_stream;
    CUDA_CHECK(cudaStreamCreate(&main_stream));
    CUDA_CHECK(cudaStreamCreate(&sub_stream));

    // Init subw_host，powcuts_host，numharms_host，numindeps_host
    int single_loop = 1 << (obs.numharmstages - 1);
    int fundamental_wlo = NEAREST_INT(ACCEL_RDW * obs.wlo) * ACCEL_DW;
    int fundamental_zlo = NEAREST_INT(ACCEL_RDZ * obs.zlo) * ACCEL_DZ;
    int fundamental_numws = subharminfs[0][0].numkern_wdim;
    int fundamental_numzs = subharminfs[0][0].numkern_zdim;
    int fundamental_numrs = obs.corr_uselen;
    long long fundamental_size = fundamental_numrs * fundamental_numzs * fundamental_numws;
    int *subw_host = malloc(single_loop * fundamental_numws * sizeof(int *));
    float *powcuts_host = (float *)malloc(obs.numharmstages * sizeof(float *));
    int *numharms_host = (int *)malloc(obs.numharmstages * sizeof(int *));
    double *numindeps_host = (double *)malloc(obs.numharmstages * sizeof(double *));
    if (obs.numharmstages > 1)
    {
        int id = 0;
        for (int ss = 0; ss < obs.numharmstages; ss++)
        {
            int hs = 1;
            if (ss > 0)
            {
                hs = 1 << ss;
            }
            powcuts_host[ss] = obs.powcut[twon_to_index(hs)];
            numharms_host[ss] = hs;
            numindeps_host[ss] = (double)obs.numindep[twon_to_index(hs)];
            for (int hh = 1; hh < hs; hh += 2)
            {
                double harm_fract = (double)hh / (double)hs;
                for (int ii = 0; ii < fundamental_numws; ii++)
                {
                    int ww = fundamental_wlo + ii * ACCEL_DW;
                    int subw = NEAREST_INT(ACCEL_RDW * ww * harm_fract) * ACCEL_DW;
                    subw_host[id * fundamental_numws + ii] = subw;
                }
                id++;
            }
        }
    }
    // Copy to the GPU. Since it occupies less space than the GPU constant memory size, place it in the GPU constant memory
    init_constant_device(subw_host, single_loop * fundamental_numws, powcuts_host, numharms_host, numindeps_host, obs.numharmstages);

    // Saving shi->kern.data from CPU to GPU
    fcomplex *fkern_gpu; // shi->kern.data
    /*When saving different keys, store the starting address of kern.data in fkern_gpu.
    Since the length of each kern varies, it is necessary to record the offset*/
    int **offset_array;

    offset_array = (int **)malloc(obs.numharmstages * sizeof(int *));
    offset_array[0] = (int *)malloc(1 * sizeof(int));
    int jj;
    for (ii = 1; ii < obs.numharmstages; ii++)
    {
        jj = 1 << ii;
        offset_array[ii] = (int *)malloc(jj * sizeof(int));
    }

    fkern_gpu = fkern_host_to_dev(subharminfs, obs.numharmstages, offset_array);
    init_inds_array(single_loop);

    /* Start the main search loop */
    {
        int max_map_size = (1 << obs.numharmstages);
        MapEntry map[max_map_size];
        double startr, lastr, nextr;

        int map_array_size = (int)(((double)obs.highestbin - obs.rlo) / ((double)obs.corr_uselen * ACCEL_DR)) + 1;
        for (int i = 0; i < max_map_size; i++)
        {
            map[i].value.startr_array = malloc(map_array_size * sizeof(double));
            map[i].value.lastr_array = malloc(map_array_size * sizeof(double));
            map[i].value.tuple_array = malloc(map_array_size * sizeof(StartrHarmTuple));
            map[i].value.count = 0;
        }

        /* Populate the saved F-Fdot plane at low freqs for in-memory
         * searches of harmonics that are below obs.rlo */
        if (obs.inmem)
        {
            printf("inmem is not supported\n");
            exit(-1);
        }

        int stage, harmtosum, harm;
        int map_size = 0;
        // Start collecting batch data
        startr = obs.rlo;
        lastr = 0;
        nextr = 0;
        int fundamental_cnt = 0;
        while (startr + rstep < obs.highestbin)
        {
            fundamental_cnt++;
            nextr = startr + rstep;
            lastr = nextr - ACCEL_DR;
            if (obs.numharmstages > 1)
            {
                for (stage = 1; stage < obs.numharmstages; stage++)
                {
                    harmtosum = 1 << stage;
                    for (harm = 1; harm < harmtosum; harm += 2)
                    {
                        subharminfo *shi = &subharminfs[stage][harm - 1];
                        MapKey key = {harmtosum, harm};
                        assert(map_size < max_map_size);
                        insertOrUpdateMap(map, &map_size, key, startr, lastr, shi, harmtosum, harm, max_map_size);
                    }
                }
            }
            startr = nextr;
        }

        // Batch processing, map_size=15
        int batch_size_max = 0, ws_len_max = 0, zs_len_max = 0, fft_len_max = 0;
        subharminfo *shii = &subharminfs[0][0];
        ws_len_max = shii->numkern_wdim;
        zs_len_max = shii->numkern_zdim;
        fft_len_max = shii->kern[0][0].fftlen;

        for (int i = 0; i < map_size; i++)
        {
            fft_len_max = MAX(map[i].value.shi->kern[0][0].fftlen, fft_len_max);
            zs_len_max = MAX(map[i].value.shi->numkern_zdim, zs_len_max);
            ws_len_max = MAX(map[i].value.shi->numkern_wdim, ws_len_max);
            batch_size_max = MAX(map[i].value.count, batch_size_max);
        }
        int proper_batch_size = MIN(cmd->batchsize, batch_size_max);
        long long array_size = proper_batch_size * sizeof(fcomplex) * fft_len_max * ws_len_max * zs_len_max; // Maximum does not exceed batch_size_max
        // printf("batch_size_max:%d, fft_len_max:%d, zs_len_max:%d, ws_len_max:%d, proper_batch_size: %d, array_size: %.2fMB\n",
        //        batch_size_max,
        //        fft_len_max,
        //        zs_len_max,
        //        ws_len_max,
        //        proper_batch_size,
        //        array_size * 1.0 / 1024 / 1024);
        fcomplex *full_tmpdat_array;
        cudaMallocAsync(&full_tmpdat_array, array_size, sub_stream); // Allocate an array of pointers on the device

        SubharmonicMap *subharmonics_add;
        CUDA_CHECK(cudaMallocAsync(&subharmonics_add, (1 << (obs.numharmstages - 1)) * proper_batch_size * sizeof(SubharmonicMap), main_stream));
        SubharmonicMap *subharmonics_add_host = (SubharmonicMap *)malloc((1 << (obs.numharmstages - 1)) * proper_batch_size * sizeof(SubharmonicMap));

        SearchValue *search_results_base;
        int search_num_base;
        if (proper_batch_size > 8)
            search_num_base = 8;
        else
            search_num_base = 4;
        long long single_search_size = sizeof(SearchValue) * fft_len_max * ws_len_max * zs_len_max;
        int search_num_array = 1;
        long long search_results_size = search_num_array * search_num_base * single_search_size;
        int max_searchnum = search_results_size / (sizeof(SearchValue)) * 0.75;
        int *too_large = (int *)malloc(sizeof(int));
        *too_large = 0;
        int *d_too_large;
        cudaMallocAsync(&d_too_large, sizeof(int), sub_stream);
        cudaMemset(d_too_large, 0, sizeof(int));
        cudaMallocAsync(&search_results_base, search_results_size, sub_stream);
        SearchValue *search_results = search_results_base;

        unsigned long long int *search_nums;
        CUDA_CHECK(cudaMalloc((void **)&search_nums, sizeof(unsigned long long int)));
        CUDA_CHECK(cudaMemsetAsync(search_nums, 0, sizeof(unsigned long long int), sub_stream));

        /* Reset indices if needed and search for real */
        if (obs.numharmstages > 1)
        { /* Search the subharmonics */
            int tk = 0, num_expand = 0;
            int batch_size = map[0].value.count;
            int total_batch = (batch_size + proper_batch_size - 1) / proper_batch_size;
            long long single_batch_size = obs.numharmstages * proper_batch_size * fundamental_size;
            int current_batch = -1;
            long long *fundamental_rlos = (long long *)malloc(batch_size * sizeof(long long *));
            SearchValue *search_results_host;
            long long total_search_num = 0;
            long long search_results_host_size = 0;

            for (int j = 0; j < batch_size; j += proper_batch_size)
            {
                current_batch++;
                int current_batch_size = proper_batch_size < batch_size - j ? proper_batch_size : batch_size - j;
                ffdotpows_cu *fundamentals = malloc(current_batch_size * sizeof(ffdotpows_cu));
                double *startr_array = map[0].value.startr_array;
                double *lastr_array = map[0].value.lastr_array;

                subharm_fderivs_vol_cu_batch(
                    fundamentals,
                    1,
                    1,
                    startr_array + j,
                    lastr_array + j,
                    &subharminfs[0][0],
                    &obs,
                    main_stream,
                    full_tmpdat_array,
                    full_tmpdat_array,
                    current_batch_size,
                    fkern_gpu,
                    0);

                for (int ib = 0; ib < current_batch_size; ib++)
                {
                    fundamental_rlos[j + ib] = fundamentals[ib].rlo;
                }

                int kk = 0, num_expand = 0;
                for (stage = 1; stage < obs.numharmstages; stage++)
                {
                    harmtosum = 1 << stage;
                    num_expand += harmtosum / 2;
                    for (harm = 1; harm < harmtosum; harm += 2)
                    {
                        // prepare batch
                        MapKey k = {harmtosum, harm};
                        MapEntry *map_entry = getMap(map, map_size, k);
                        if (map_entry == NULL)
                        {
                            fprintf(stderr, "map_entry should not be null\n");
                            exit(EXIT_FAILURE);
                        }

                        ffdotpows_cu *subharmonics_batch = malloc(batch_size * sizeof(ffdotpows_cu));
                        subharminfo *shi_local = map_entry->value.shi;

                        double *sub_startr_array = map_entry[0].value.startr_array + j;
                        double *sub_lastr_array = map_entry[0].value.lastr_array + j;
                        int harmtosum_local = map_entry[0].key.harmtosum;
                        int harm_local = map_entry[0].key.harm;

                        subharm_fderivs_vol_cu_batch(
                            subharmonics_batch,
                            harmtosum_local,
                            harm_local,
                            sub_startr_array,
                            sub_lastr_array,
                            shi_local,
                            &obs,
                            main_stream,
                            full_tmpdat_array,
                            full_tmpdat_array,
                            current_batch_size,
                            fkern_gpu + offset_array[stage][harm - 1],
                            kk);

                        for (int i = 0; i < current_batch_size; i++)
                        {
                            subharmonics_add_host[kk * current_batch_size + i].harm_fract = kk;
                            subharmonics_add_host[kk * current_batch_size + i].subharmonic_numrs = subharmonics_batch[i].numrs;
                            subharmonics_add_host[kk * current_batch_size + i].subharmonic_numzs = subharmonics_batch[i].numzs;
                            subharmonics_add_host[kk * current_batch_size + i].subharmonic_powers = subharmonics_batch[i].powers;
                            subharmonics_add_host[kk * current_batch_size + i].subharmonic_rinds = subharmonics_batch[i].rinds;
                            subharmonics_add_host[kk * current_batch_size + i].subharmonic_zinds = subharmonics_batch[i].zinds;
                            subharmonics_add_host[kk * current_batch_size + i].subharmonic_wlo = subharmonics_batch[i].wlo;
                        }
                        kk++;
                        free(subharmonics_batch);
                    }
                }

                CUDA_CHECK(cudaMemcpyAsync(subharmonics_add, subharmonics_add_host, (1 << (obs.numharmstages - 1)) * proper_batch_size * sizeof(SubharmonicMap), cudaMemcpyHostToDevice, main_stream));
                fuse_add_search_batch(fundamentals, subharmonics_add, (obs.numharmstages - 1), current_batch_size, main_stream, search_results, search_nums, (long long)(current_batch * single_batch_size), proper_batch_size, max_searchnum, d_too_large);

                // copy results to host
                cudaStreamSynchronize(main_stream);
                cudaMemcpy(too_large, d_too_large, sizeof(int), cudaMemcpyDeviceToHost);
                if (too_large[0])
                {
                    free_subharmonic_cu_batch(subharmonics_add_host, current_batch_size, num_expand, sub_stream);
                    free_ffdotpows_cu_batch(fundamentals, current_batch_size, sub_stream);
                    cudaFreeAsync(full_tmpdat_array, sub_stream);
                    cudaFreeAsync(search_results, sub_stream);
                    cudaFreeAsync(subharmonics_add, sub_stream);
                    cudaFreeAsync(search_nums, sub_stream);
                    cudaFreeAsync(d_too_large, sub_stream);
                    // clean cufftPlan
                    clear_cache();
                    cudaFree(fkern_gpu);
                    CUDA_CHECK(cudaStreamDestroy(main_stream));
                    CUDA_CHECK(cudaStreamDestroy(sub_stream));

                    freeMap(map, &max_map_size);
                    free(subharmonics_add_host);
                    free(fundamental_rlos);
                    free(subw_host);
                    free(powcuts_host);
                    free(numharms_host);
                    free(numindeps_host);
                    free(offset_array[0]);
                    for (ii = 1; ii < obs.numharmstages; ii++)
                    {
                        int jj = 1 << ii;
                        free(offset_array[ii]);
                    }
                    free(offset_array);
                    free_subharminfos(&obs, subharminfs);
                    subharminfs = NULL;
                    *cands_ptr = NULL;
                    free(too_large);
                    return 1;
                }

                unsigned long long int *search_nums_current = (unsigned long long int *)malloc(sizeof(unsigned long long int));
                CUDA_CHECK(cudaMemcpyAsync(search_nums_current, search_nums, sizeof(unsigned long long int), cudaMemcpyDeviceToHost, main_stream));
                unsigned long long int search_num = search_nums_current[0];
                if (search_num > search_num_array * search_num_base * single_search_size / (sizeof(SearchValue) * 2))
                {
                    // Too much data in search_results, transferring back to CPU
                    // Sort the data first because indices in later batches are always greater than earlier ones, so partial sorting won't affect overall order.
                    sort_search_results(search_results, search_num); 
                    if (total_search_num == 0)                       // search_results_host memory has not been allocated yet
                    {
                        search_results_host_size = search_num_array * search_num_base * single_search_size;
                        search_results_host = (SearchValue *)malloc(search_results_host_size);
                        cudaStreamSynchronize(main_stream);
                        CUDA_CHECK(cudaMemcpy(search_results_host, search_results, search_num * sizeof(SearchValue), cudaMemcpyDeviceToHost));
                        total_search_num = search_num;
                        CUDA_CHECK(cudaMemsetAsync(search_nums, 0, sizeof(unsigned long long int), sub_stream));
                    }
                    else if ((search_num + total_search_num) * sizeof(SearchValue) < search_results_host_size) // Previously allocated memory for search_results_host is sufficient
                    {
                        cudaStreamSynchronize(main_stream);
                        CUDA_CHECK(cudaMemcpy(search_results_host + total_search_num, search_results, search_num * sizeof(SearchValue), cudaMemcpyDeviceToHost));
                        total_search_num += search_num;
                        CUDA_CHECK(cudaMemsetAsync(search_nums, 0, sizeof(unsigned long long int), sub_stream));
                    }
                    else // Previously allocated memory for search_results_host is insufficient to store search_num + total_search_num elements
                    {
                        // 1. Allocate larger contiguous memory
                        search_results_host_size += search_num_array * search_num_base * single_search_size;
                        SearchValue *search_results_host_append;
                        search_results_host_append = (SearchValue *)malloc(search_results_host_size);
                        // 2. Copy existing data in search_results_host to the new memory
                        memcpy(search_results_host_append, search_results_host, total_search_num * sizeof(SearchValue));
                        // 3. Transfer search_num elements from GPU to the new memory
                        cudaStreamSynchronize(main_stream);
                        CUDA_CHECK(cudaMemcpy(search_results_host_append + total_search_num, search_results, search_num * sizeof(SearchValue), cudaMemcpyDeviceToHost));
                        // 4. Free the old memory
                        free(search_results_host);
                        // 5. Update the pointer
                        search_results_host = search_results_host_append;
                        // 6. Update total_search_num and search_nums
                        total_search_num += search_num;
                        CUDA_CHECK(cudaMemsetAsync(search_nums, 0, sizeof(unsigned long long int), sub_stream));
                    }
                }
                free_subharmonic_cu_batch(subharmonics_add_host, current_batch_size, num_expand, sub_stream);
                free_ffdotpows_cu_batch(fundamentals, current_batch_size, sub_stream);
            }
            cudaStreamSynchronize(main_stream);
            unsigned long long int *search_nums_host = (unsigned long long int *)malloc(sizeof(unsigned long long int));
            CUDA_CHECK(cudaMemcpy(search_nums_host, search_nums, sizeof(unsigned long long int), cudaMemcpyDeviceToHost));
            unsigned long long int search_num = search_nums_host[0];
            
            if (search_num > 0)
            {
                sort_search_results(search_results, search_num);
                if (total_search_num == 0) // search_results_host memory has not been allocated yet
                {
                search_results_host = (SearchValue *)malloc(search_num * sizeof(SearchValue));
                CUDA_CHECK(cudaMemcpy(search_results_host, search_results, search_num * sizeof(SearchValue), cudaMemcpyDeviceToHost));
                total_search_num = search_num;
                }
                else if ((search_num + total_search_num) * sizeof(SearchValue) < search_results_host_size) // Previously allocated memory for search_results_host is sufficient
                {
                    CUDA_CHECK(cudaMemcpy(search_results_host + total_search_num, search_results, search_num * sizeof(SearchValue), cudaMemcpyDeviceToHost));
                    total_search_num += search_num;
                }
                else // Previously allocated memory for search_results_host is insufficient to store search_num + total_search_num elements
                {
                    // 1. Allocate larger contiguous memory
                    SearchValue *search_results_host_append;
                    search_results_host_append = (SearchValue *)malloc((search_num + total_search_num) * sizeof(SearchValue));
                    // 2. Copy existing data in search_results_host to the new memory
                    memcpy(search_results_host_append, search_results_host, total_search_num * sizeof(SearchValue));
                    // 3. Transfer search_num elements from GPU to the new memory
                    CUDA_CHECK(cudaMemcpy(search_results_host_append + total_search_num, search_results, search_num * sizeof(SearchValue), cudaMemcpyDeviceToHost));
                    // 4. Free the old memory
                    free(search_results_host);
                    // 5. Update the pointer
                    search_results_host = search_results_host_append;
                    // 6. Update total_search_num
                    total_search_num += search_num;
                }
            }
            cudaFreeAsync(full_tmpdat_array, sub_stream);
            cudaFreeAsync(search_results, sub_stream);
            cudaFreeAsync(subharmonics_add, sub_stream);
            cudaFreeAsync(search_nums, sub_stream);
            cudaFreeAsync(d_too_large, sub_stream);
            // clean cufftPlan
            clear_cache();
            cudaFree(fkern_gpu);
            CUDA_CHECK(cudaStreamDestroy(main_stream));
            CUDA_CHECK(cudaStreamDestroy(sub_stream));

            cands = insert_to_cands(fundamental_numrs, fundamental_numzs, fundamental_numws, fundamental_rlos, fundamental_zlo, fundamental_wlo, proper_batch_size, numindeps_host, cands, search_results_host, total_search_num, single_batch_size, obs.numharmstages, main_stream, sub_stream);

            freeMap(map, &max_map_size);
            free(subharmonics_add_host);
            free(fundamental_rlos);
            free(subw_host);
            free(powcuts_host);
            free(numharms_host);
            free(numindeps_host);
            free(offset_array[0]);
            for (ii = 1; ii < obs.numharmstages; ii++)
            {
                int jj = 1 << ii;
                free(offset_array[ii]);
            }
            free(offset_array);
            free(too_large);
        }
    }
    free_subharminfos(&obs, subharminfs);
    subharminfs = NULL;
    *cands_ptr = cands;
    return 0;
}

void accelsearch_CPU2(GSList **cands, accelobs *obs, infodata *idata, Cmdline *cmd)
{
    /* Candidate list trimming and optimization */
    int ii;
    int numcands = g_slist_length(*cands);
    printf("numcands after while:%d\n", numcands);
    GSList *listptr;
    accelcand *cand;
    fourierprops *props;

    if (numcands)
    {

        /* Sort the candidates according to the optimized sigmas */
        *cands = sort_accelcands(*cands);

        /* Eliminate (most of) the harmonically related candidates */
        if ((cmd->numharm > 1) && !(cmd->noharmremoveP))
            eliminate_harmonics(*cands, &numcands);

        /* Now optimize each candidate and its harmonics */
        // print_percent_complete(0, 0, NULL, 1);
        listptr = *cands;
        for (ii = 0; ii < numcands; ii++)
        {
            // print_percent_complete(ii, numcands, "optimization", 0);
            cand = (accelcand *)(listptr->data);
            optimize_accelcand(cand, obs);
            listptr = listptr->next;
        }
        // print_percent_complete(ii, numcands, "optimization", 0);

        /* Calculate the properties of the fundamentals */
        props = (fourierprops *)malloc(sizeof(fourierprops) * numcands);
        listptr = *cands;
        for (ii = 0; ii < numcands; ii++)
        {
            cand = (accelcand *)(listptr->data);
            /* In case the fundamental harmonic is not significant,  */
            /* send the originally determined r and z from the       */
            /* harmonic sum in the search.  Note that the derivs are */
            /* not used for the computations with the fundamental.   */
            calc_props(cand->derivs[0], cand->r, cand->z, cand->w, props + ii);
            /* Override the error estimates based on power */
            props[ii].rerr = (float)(ACCEL_DR) / cand->numharm;
            props[ii].zerr = (float)(ACCEL_DZ) / cand->numharm;
            props[ii].werr = (float)(ACCEL_DW) / cand->numharm;
            listptr = listptr->next;
        }

        /* Write the fundamentals to the output text file */
        output_fundamentals(props, *cands, obs, idata);

        /* Write the harmonics to the output text file */
        output_harmonics(*cands, obs, idata);

        /* Write the fundamental fourierprops to the cand file */
        obs->workfile = chkfopen(obs->candnm, "wb");
        chkfwrite(props, sizeof(fourierprops), numcands, obs->workfile);
        fclose(obs->workfile);
        free(props);
        printf("\n\n");
    }
    else
    {
        printf("No candidates above sigma = %.2f were found.\n\n", obs->sigma);
    }

    printf("Final candidates in binary format are in '%s'.\n", obs->candnm);
    printf("Final Candidates in a text format are in '%s'.\n\n", obs->accelnm);

    free_accelobs(obs);
    g_slist_foreach(*cands, free_accelcand, NULL);
    g_slist_free(*cands);
    obs = NULL;
    idata = NULL;
    cands = NULL;
}

int main(int argc, char *argv[])
{
    subharminfo **subharminfs;
    accelobs obs;
    infodata idata;
    GSList *cands = NULL;
    Cmdline *cmd;

    accelsearch_CPU1(argc, argv, &subharminfs, &obs, &idata, &cmd);
    int too_large = accelsearch_GPU(obs, subharminfs, &cands, cmd);
    accelsearch_CPU2(&cands, &obs, &idata, cmd);
    return 0;
}
