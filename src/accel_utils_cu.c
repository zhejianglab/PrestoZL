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
#include "cuda_runtime.h"
#include "cuda_helper.h"
#include <sys/time.h>
#include <assert.h>

#define NEAREST_INT(x) (int)(x < 0 ? x - 0.5 : x + 0.5)

unsigned short **inds_array;

void init_inds_array(int size)
{
    inds_array = (unsigned short **)malloc(size * sizeof(unsigned short *));
    for (int i = 0; i < size; i++)
    {
        inds_array[i] = NULL;
    }
}

void free_inds_array()
{
    free(inds_array);
}

typedef struct
{
    long long index;
    float pow;
    float sig;
} SearchValue;

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

// Comparison function for comparing the index field of two SearchValue structures
int compare(const void *a, const void *b)
{
    const SearchValue *valueA = (const SearchValue *)a;
    const SearchValue *valueB = (const SearchValue *)b;

    if (valueA->index < valueB->index)
    {
        return -1;
    }
    else if (valueA->index > valueB->index)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/* Return 2**n */
#define index_to_twon(n) (1 << n)

float ***gen_f3Darr_cu(long nhgts, long nrows, long ncols, cudaStream_t stream);

kernel **gen_kernmatrix_cu(int numz, int numw);

long long timeInMilliseconds(void)
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    return (((long long)tv.tv_sec) * 1000) + (tv.tv_usec / 1000);
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

static inline double calc_required_r(double harm_fract, double rfull)
/* Calculate the 'r' you need for subharmonic  */
/* harm_fract = harmnum / numharm if the       */
/* 'r' at the fundamental harmonic is 'rfull'. */
{
    return rint(ACCEL_RDR * rfull * harm_fract) * ACCEL_DR;
}

static inline int calc_required_z(double harm_fract, double zfull)
/* Calculate the 'z' you need for subharmonic  */
/* harm_fract = harmnum / numharm if the       */
/* 'z' at the fundamental harmonic is 'zfull'. */
{
    return NEAREST_INT(ACCEL_RDZ * zfull * harm_fract) * ACCEL_DZ;
}

static inline int calc_required_w(double harm_fract, double wfull)
/* Calculate the maximum 'w' needed for the given subharmonic  */
/* harm_fract = harmnum / numharm if the       */
/* 'w' at the fundamental harmonic is 'wfull'. */
{
    return NEAREST_INT(ACCEL_RDW * wfull * harm_fract) * ACCEL_DW;
}

static inline int index_from_r(double r, double lor)
/* Return an index for a Fourier Freq given an array that */
/* has stepsize ACCEL_DR and low freq 'lor'.              */
{
    return (int)((r - lor) * ACCEL_RDR + DBLCORRECT);
}

static inline int index_from_z(double z, double loz)
/* Return an index for a Fourier Fdot given an array that */
/* has stepsize ACCEL_DZ and low freq dot 'loz'.              */
{
    return (int)((z - loz) * ACCEL_RDZ + DBLCORRECT);
}

static inline int index_from_w(double w, double low)
/* Return an index for a Fourier Fdotdot given an array that */
/* has stepsize ACCEL_DW and low freq dotdot 'low'.              */
{
    return (int)((w - low) * ACCEL_RDW + DBLCORRECT);
}

static fcomplex *gen_cvect_cu(long length)
{
    fcomplex *v;

    CUDA_CHECK(cudaMallocManaged(&v, (size_t)(sizeof(fcomplex) * length), cudaMemAttachGlobal));

    return v;
}

static int calc_fftlen(int numharm, int harmnum, int max_zfull, int max_wfull, accelobs *obs)
/* The fft length needed to properly process a subharmonic */
{
    int bins_needed, end_effects;
    double harm_fract;

    harm_fract = (double)harmnum / (double)numharm;
    bins_needed = (int)ceil(obs->corr_uselen * harm_fract) + 2;
    end_effects = 2 * ACCEL_NUMBETWEEN *
                  w_resp_halfwidth(calc_required_z(harm_fract, max_zfull),
                                   calc_required_w(harm_fract, max_wfull), LOWACC);
    return next_good_fftlen(bins_needed + end_effects);
}

static void init_kernel(int z, int w, int fftlen, kernel *kern)
{
    int numkern;
    fcomplex *tempkern;

    kern->z = z;
    kern->w = w;
    kern->fftlen = fftlen;
    kern->numbetween = ACCEL_NUMBETWEEN;
    kern->kern_half_width = w_resp_halfwidth((double)z, (double)w, LOWACC);
    numkern = 2 * kern->numbetween * kern->kern_half_width;
    kern->numgoodbins = kern->fftlen - numkern;
    kern->data = gen_cvect_cu(kern->fftlen);
    tempkern = gen_w_response(0.0, kern->numbetween, kern->z, kern->w, numkern);
    place_complex_kernel(tempkern, numkern, kern->data, kern->fftlen);
    vect_free(tempkern);
    COMPLEXFFT(kern->data, kern->fftlen, -1);

    // pre fetch kern->data
    CUDA_CHECK(cudaMemPrefetchAsync(kern->data, (size_t)(sizeof(fcomplex) * kern->fftlen), 0, NULL));
}

static void init_subharminfo_cu(int numharm, int harmnum, int zmax, int wmax, subharminfo *shi, accelobs *obs)
/* Note:  'zmax' is the overall maximum 'z' in the search while
          'wmax' is the overall maximum 'w' in the search       */
{
    int ii, jj, fftlen;
    double harm_fract;

    harm_fract = (double)harmnum / (double)numharm;
    shi->numharm = numharm;
    shi->harmnum = harmnum;
    shi->zmax = calc_required_z(harm_fract, zmax);
    shi->wmax = calc_required_w(harm_fract, wmax);
    if (numharm > 1)
    {
        CUDA_CHECK(cudaMallocManaged(&(shi->rinds), (size_t)(obs->corr_uselen * sizeof(unsigned short)), cudaMemAttachGlobal));
        CUDA_CHECK(cudaMallocManaged(&(shi->zinds), (size_t)(obs->corr_uselen * sizeof(unsigned short)), cudaMemAttachGlobal));
    }
    if (numharm == 1 && harmnum == 1)
        fftlen = obs->fftlen;
    else
        fftlen = calc_fftlen(numharm, harmnum, zmax, wmax, obs);
    shi->numkern_zdim = (shi->zmax / ACCEL_DZ) * 2 + 1;
    shi->numkern_wdim = (shi->wmax / ACCEL_DW) * 2 + 1;
    shi->numkern = shi->numkern_zdim * shi->numkern_wdim;
    /* Allocate 2D array of kernels, with dimensions being z and w */
    shi->kern = gen_kernmatrix_cu(shi->numkern_zdim, shi->numkern_wdim);

    /* Actually append kernels to each array element */
    for (ii = 0; ii < shi->numkern_wdim; ii++)
    {
        for (jj = 0; jj < shi->numkern_zdim; jj++)
        {
            init_kernel(-shi->zmax + jj * ACCEL_DZ,
                        -shi->wmax + ii * ACCEL_DW, fftlen, &shi->kern[ii][jj]);
        }
    }
}

float *gen_f3Darr_flat(long nhgts, long nrows, long ncols, cudaStream_t main_stream)
{
    float *c;

    CUDA_CHECK(cudaMallocAsync(&c, (size_t)(nhgts * nrows * ncols * sizeof(float)), main_stream));

    return c;
}

void free_ffdotpows_cu(ffdotpows *ffd,
                       cudaStream_t sub_stream)
{
    CUDA_CHECK(cudaFreeAsync(ffd->powers, sub_stream));
    free(ffd);
}

void free_ffdotpows_cu_batch(ffdotpows_cu *ffd_array, int batch_size,
                             cudaStream_t sub_stream)
{
    CUDA_CHECK(cudaFreeAsync(ffd_array[0].powers, sub_stream));
    free(ffd_array);
}

void free_subharmonic_cu_batch(SubharmonicMap *ffd_array, int batch_size, int num_expand,
                               cudaStream_t sub_stream)
{
    for (int i = 0; i < num_expand; i++)
    {
        SubharmonicMap *ffd = &ffd_array[i * batch_size];
        CUDA_CHECK(cudaFreeAsync(ffd->subharmonic_powers, sub_stream));
        CUDA_CHECK(cudaFreeAsync(inds_array[i], sub_stream));
    }
}

void free_subharminfo_cu(subharminfo *shi)
{
    int ii, jj;

    for (ii = 0; ii < shi->numkern_wdim; ii++)
    {
        for (jj = 0; jj < shi->numkern_zdim; jj++)
        {
            CUDA_CHECK(cudaFree((&shi->kern[ii][jj])->data));
        }
    }

    if (shi->numharm > 1)
    {
        CUDA_CHECK(cudaFree(shi->rinds));
        CUDA_CHECK(cudaFree(shi->zinds));
    }

    CUDA_CHECK(cudaFree(shi->kern[0]));
    CUDA_CHECK(cudaFree(shi->kern));
}

void free_subharminfos_cu(accelobs *obs, subharminfo **shis)
{
    int ii, jj, harmtosum;

    /* Free the sub-harmonics */
    if (!obs->inmem)
    {
        for (ii = 1; ii < obs->numharmstages; ii++)
        {
            harmtosum = index_to_twon(ii);
            for (jj = 1; jj < harmtosum; jj += 2)
            {
                free_subharminfo_cu(&shis[ii][jj - 1]);
            }
            free(shis[ii]);
        }
    }

    /* Free the fundamental */
    free_subharminfo_cu(&shis[0][0]);
    free(shis[0]);
    /* Free the container */
    free(shis);
}

subharminfo **create_subharminfos_cu(accelobs *obs)
{
    double kern_ram_use = 0;
    int ii, jj, harmtosum, fftlen;
    subharminfo **shis;

    shis = (subharminfo **)malloc(obs->numharmstages * sizeof(subharminfo *));
    /* Prep the fundamental (actually, the highest harmonic) */
    shis[0] = (subharminfo *)malloc(2 * sizeof(subharminfo));
    init_subharminfo_cu(1, 1, (int)obs->zhi, (int)obs->whi, &shis[0][0], obs);
    fftlen = obs->fftlen;
    kern_ram_use += shis[0][0].numkern * fftlen * sizeof(fcomplex); // in Bytes
    if (obs->numw)
        printf("  Harm  1/1 : %5d kernels, %4d < z < %-4d and %5d < w < %-5d (%5d pt FFTs)\n",
               shis[0][0].numkern, -shis[0][0].zmax, shis[0][0].zmax,
               -shis[0][0].wmax, shis[0][0].wmax, fftlen);
    else
        printf("  Harm  1/1 : %5d kernels, %4d < z < %-4d (%d pt FFTs)\n",
               shis[0][0].numkern, -shis[0][0].zmax, shis[0][0].zmax, fftlen);
    /* Prep the sub-harmonics if needed */
    if (!obs->inmem)
    {
        for (ii = 1; ii < obs->numharmstages; ii++)
        {
            harmtosum = index_to_twon(ii);
            shis[ii] = (subharminfo *)malloc(harmtosum * sizeof(subharminfo));
            for (jj = 1; jj < harmtosum; jj += 2)
            {
                init_subharminfo_cu(harmtosum, jj, (int)obs->zhi,
                                    (int)obs->whi, &shis[ii][jj - 1], obs);
                fftlen = calc_fftlen(harmtosum, jj, (int)obs->zhi, (int)obs->whi, obs);
                kern_ram_use += shis[ii][jj - 1].numkern * fftlen * sizeof(fcomplex); // in Bytes
                if (obs->numw)
                    printf("  Harm %2d/%-2d: %5d kernels, %4d < z < %-4d and %5d < w < %-5d (%5d pt FFTs)\n",
                           jj, harmtosum, shis[ii][jj - 1].numkern,
                           -shis[ii][jj - 1].zmax, shis[ii][jj - 1].zmax,
                           -shis[ii][jj - 1].wmax, shis[ii][jj - 1].wmax, fftlen);
                else
                    printf("  Harm %2d/%-2d: %5d kernels, %4d < z < %-4d (%d pt FFTs)\n",
                           jj, harmtosum, shis[ii][jj - 1].numkern,
                           -shis[ii][jj - 1].zmax, shis[ii][jj - 1].zmax, fftlen);
            }
        }
    }
    printf("Total RAM used by correlation kernels:  %.3f GB\n", kern_ram_use / (1 << 30));
    return shis;
}

void deep_copy_ffdotpows_cpu2cu(ffdotpows_cu *ffdot, subharminfo *shi, int corr_uselen, unsigned short *inds, int b, int batch_size, cudaStream_t stream)
{
    // alloc and copy rinds
    ffdot->rinds = &inds[b * corr_uselen];
    CUDA_CHECK(cudaMemcpyAsync(ffdot->rinds, shi->rinds, corr_uselen * sizeof(unsigned short), cudaMemcpyHostToDevice, stream));
    // alloc and copy zinds
    ffdot->zinds = &inds[batch_size * corr_uselen + b * corr_uselen];
    CUDA_CHECK(cudaMemcpyAsync(ffdot->zinds, shi->zinds, corr_uselen * sizeof(unsigned short), cudaMemcpyHostToDevice, stream));
}

fcomplex *fkern_host_to_dev(subharminfo **subharminfs, int numharmstages, int **offset_array)
{
    int fkern_size = 0;
    int harm, harmtosum, stage;
    int numkern, fftlen;
    int offset_base = 0, offset_tmp;
    fcomplex *fkern_cpu;
    fcomplex *fkern_gpu;

    for (stage = 0; stage < numharmstages; stage++)
    {
        harmtosum = (stage == 0 ? 2 : 1 << stage);
        for (harm = 1; harm < harmtosum; harm += 2)
        {
            numkern = subharminfs[stage][harm - 1].numkern;
            fftlen = subharminfs[stage][harm - 1].kern[0][0].fftlen;
            offset_array[stage][harm - 1] = fkern_size;
            fkern_size = fkern_size + numkern * fftlen;
        }
    }

    fkern_cpu = (fcomplex *)malloc(sizeof(fcomplex) * fkern_size);
    cudaMalloc((void **)&fkern_gpu, sizeof(fcomplex) * fkern_size);

    for (stage = 0; stage < numharmstages; stage++)
    {
        harmtosum = (stage == 0 ? 2 : 1 << stage);
        for (harm = 1; harm < harmtosum; harm += 2)
        {

            offset_base = offset_array[stage][harm - 1];

            fftlen = subharminfs[stage][harm - 1].kern[0][0].fftlen;
            numkern = subharminfs[stage][harm - 1].numkern;

            for (int ii = 0; ii < numkern; ii++)
            {
                for (int jj = 0; jj < fftlen; jj++)
                {
                    offset_tmp = offset_base + fftlen * ii + jj;
                    fkern_cpu[offset_tmp].r = subharminfs[stage][harm - 1].kern[0][ii].data[jj].r;
                    fkern_cpu[offset_tmp].i = subharminfs[stage][harm - 1].kern[0][ii].data[jj].i;
                }
            }
        }
    }
    cudaMemcpy(fkern_gpu, fkern_cpu, sizeof(fcomplex) * fkern_size, cudaMemcpyHostToDevice);

    free(fkern_cpu);

    return fkern_gpu;
}

// batch create subharmonics
void subharm_fderivs_vol_cu_batch(
    ffdotpows_cu *ffdot_array,
    int numharm,
    int harmnum,
    double *fullrlo_array, // store batch fullrlo
    double *fullrhi_array, // store batch fullrhi
    subharminfo *shi,      // store batch subharminfo
    accelobs *obs,         // store batch accelobs
    cudaStream_t stream,
    fcomplex *full_tmpdat_array,
    fcomplex *full_tmpout_array,
    int batch_size,
    fcomplex *fkern,
    int inds_idx) // size of batch
{
    // local variables needed
    int ii, numdata, fftlen, binoffset;
    long long lobin;
    float powargr, powargi;
    double drlo, drhi, harm_fract;
    fcomplex *data, *pdata;

    // prepare pdata_dev
    fcomplex *pdata_dev;
    unsigned short *inds;
    fftlen = shi->kern[0][0].fftlen;
    CUDA_CHECK(cudaMallocAsync(&pdata_dev, (size_t)(sizeof(fcomplex) * fftlen * batch_size), stream));
    if (!(numharm == 1 && harmnum == 1))
    {
        CUDA_CHECK(cudaMallocAsync(&(inds_array[inds_idx]), (size_t)(obs->corr_uselen * sizeof(unsigned short) * batch_size * 2), stream));
    }
    // loop through each batch
    for (int b = 0; b < batch_size; b++)
    {
        double fullrlo = fullrlo_array[b];
        double fullrhi = fullrhi_array[b];

        ffdotpows_cu *ffdot = &ffdot_array[b];
        /* Calculate and get the required amplitudes */
        harm_fract = (double)harmnum / (double)numharm;
        drlo = calc_required_r(harm_fract, fullrlo);
        drhi = calc_required_r(harm_fract, fullrhi);
        ffdot->rlo = (long long)floor(drlo);
        ffdot->zlo = calc_required_z(harm_fract, obs->zlo);
        ffdot->wlo = calc_required_w(harm_fract, obs->wlo);

        /* Initialize the lookup indices */
        if (numharm > 1 && !obs->inmem)
        {
            double rr, subr;
            for (ii = 0; ii < obs->corr_uselen; ii++)
            {
                rr = fullrlo + ii * ACCEL_DR;
                subr = calc_required_r(harm_fract, rr);
                shi->rinds[ii] = index_from_r(subr, ffdot->rlo);
            }
            double zz, subz;
            for (ii = 0; ii < obs->numz; ii++)
            {
                zz = obs->zlo + ii * ACCEL_DZ;
                subz = calc_required_z(harm_fract, zz);
                shi->zinds[ii] = index_from_z(subz, ffdot->zlo);
            }
        }

        // The +1 below is important!
        ffdot->numrs = (int)((ceil(drhi) - floor(drlo)) * ACCEL_RDR + DBLCORRECT) + 1;
        if (numharm == 1 && harmnum == 1)
        {
            ffdot->numrs = obs->corr_uselen;
        }
        else
        {
            if (ffdot->numrs % ACCEL_RDR)
                ffdot->numrs = (ffdot->numrs / ACCEL_RDR + 1) * ACCEL_RDR;
        }
        ffdot->numzs = shi->numkern_zdim;
        ffdot->numws = shi->numkern_wdim;
        if (numharm == 1 && harmnum == 1)
        {
            ffdot->rinds = shi->rinds;
            ffdot->zinds = shi->zinds;
        }
        else
        {
            deep_copy_ffdotpows_cpu2cu(ffdot, shi, obs->corr_uselen, (inds_array[inds_idx]), b, batch_size, stream);
        }
        /* Determine the largest kernel halfwidth needed to analyze the current subharmonic */
        /* Verified numerically that, as long as we have symmetric z's and w's, */
        /* shi->kern[0][0].kern_half_width is the maximal halfwidth over the range of w's and z's */
        binoffset = shi->kern[0][0].kern_half_width;
        lobin = ffdot->rlo - binoffset;
        numdata = fftlen / ACCEL_NUMBETWEEN;
        data = get_fourier_amplitudes(lobin, numdata, obs);
        // Create the output power array

        // Normalize the Fourier amplitudes
        if (obs->nph > 0.0)
        {
            //  Use freq 0 normalization if requested (i.e. photons)
            double norm = 1.0 / sqrt(obs->nph);
            for (ii = 0; ii < numdata; ii++)
            {
                data[ii].r *= norm;
                data[ii].i *= norm;
            }
        }
        else if (obs->norm_type == 0)
        {
            // default block median normalization
            float *powers;
            double norm;
            powers = gen_fvect(numdata);
            for (ii = 0; ii < numdata; ii++)
                powers[ii] = POWER(data[ii].r, data[ii].i);
            norm = 1.0 / sqrt(median(powers, numdata) / log(2.0));
            vect_free(powers);
            for (ii = 0; ii < numdata; ii++)
            {
                data[ii].r *= norm;
                data[ii].i *= norm;
            }
        }
        else
        {
            // optional running double-tophat local-power normalization
            float *powers, *loc_powers;
            powers = gen_fvect(numdata);
            for (ii = 0; ii < numdata; ii++)
            {
                powers[ii] = POWER(data[ii].r, data[ii].i);
            }
            loc_powers = corr_loc_pow(powers, numdata);
            for (ii = 0; ii < numdata; ii++)
            {
                float norm = invsqrtf(loc_powers[ii]);
                data[ii].r *= norm;
                data[ii].i *= norm;
            }
            vect_free(powers);
            vect_free(loc_powers);
        }

        // Prep, spread, and FFT the data
        pdata = gen_cvect(fftlen);

        spread_no_pad(data, fftlen / ACCEL_NUMBETWEEN, pdata, fftlen, ACCEL_NUMBETWEEN);

        // Note COMPLEXFFT is not thread-safe because of wisdom caching
        COMPLEXFFT(pdata, fftlen, -1);

        // copy pdata to pdata_dev, host to device
        CUDA_CHECK(cudaMemcpyAsync(&pdata_dev[b * fftlen], pdata, (size_t)(sizeof(fcomplex) * fftlen), cudaMemcpyHostToDevice, stream));

        vect_free(data);
        vect_free(pdata);
    }
    int powers_len = 0;
    for (int b = 0; b < batch_size; b++)
    {
        ffdotpows_cu *ffdot = &ffdot_array[b];
        powers_len += ffdot->numws * ffdot->numzs * ffdot->numrs;
    }

    float *powers_dev_batch;
    CUDA_CHECK(cudaMallocAsync(&powers_dev_batch, (size_t)(powers_len * sizeof(float)), stream));
    // TODO Free powers_dev_batch
    int idx = 0;
    int *idx_array = (int *)malloc(batch_size * sizeof(int));
    for (int b = 0; b < batch_size; b++)
    {
        ffdotpows_cu *ffdot = &ffdot_array[b];
        ffdot->powers = &powers_dev_batch[idx];
        idx_array[b] = idx;
        idx += ffdot->numws * ffdot->numzs * ffdot->numrs;
    }

    do_fft_batch(fftlen, binoffset, ffdot_array, shi, pdata_dev, idx_array, full_tmpdat_array, full_tmpout_array, batch_size, fkern, stream);
    CUDA_CHECK(cudaFreeAsync(pdata_dev, stream));
    free(idx_array);
}

static accelcand *create_accelcand(float power, float sigma,
                                   int numharm, double r, double z, double w)
{
    accelcand *obj;

    obj = (accelcand *)malloc(sizeof(accelcand));
    obj->power = power;
    obj->sigma = sigma;
    obj->numharm = numharm;
    obj->r = r;
    obj->z = z;
    obj->w = w;
    obj->pows = NULL;
    obj->hirs = NULL;
    obj->hizs = NULL;
    obj->hiws = NULL;
    obj->derivs = NULL;
    return obj;
}

static GSList *insert_new_accelcand(GSList *list, float power, float sigma,
                                    int numharm, double rr, double zz, double ww, int *added)
/* Checks the current list to see if there is already */
/* a candidate within ACCEL_CLOSEST_R bins.  If not,  */
/* it adds it to the list in increasing freq order.   */
{
    GSList *tmp_list = list, *prev_list = NULL, *new_list;
    double prev_diff_r = ACCEL_CLOSEST_R + 1.0, next_diff_r;

    *added = 0;
    if (!list)
    {
        new_list = g_slist_alloc();
        new_list->data =
            (gpointer *)create_accelcand(power, sigma, numharm, rr, zz, ww);
        *added = 1;
        return new_list;
    }

    /* Find the correct position in the list for the candidate */

    while ((tmp_list->next) && (((accelcand *)(tmp_list->data))->r < rr))
    {
        prev_list = tmp_list;
        tmp_list = tmp_list->next;
    }
    next_diff_r = fabs(rr - ((accelcand *)(tmp_list->data))->r);
    if (prev_list)
        prev_diff_r = fabs(rr - ((accelcand *)(prev_list->data))->r);

    /* Similar candidate(s) is(are) present */

    if (prev_diff_r < ACCEL_CLOSEST_R)
    {
        /* Overwrite the prev cand */
        if (((accelcand *)(prev_list->data))->sigma < sigma)
        {
            free_accelcand(prev_list->data, NULL);
            prev_list->data = (gpointer *)create_accelcand(power, sigma,
                                                           numharm, rr, zz, ww);
            *added = 1;
        }
        if (next_diff_r < ACCEL_CLOSEST_R)
        {
            if (((accelcand *)(tmp_list->data))->sigma < sigma)
            {
                free_accelcand(tmp_list->data, NULL);
                if (*added)
                {
                    /* Remove the next cand */
                    list = g_slist_remove_link(list, tmp_list);
                    g_slist_free_1(tmp_list);
                }
                else
                {
                    /* Overwrite the next cand */
                    tmp_list->data = (gpointer *)create_accelcand(power, sigma,
                                                                  numharm, rr, zz, ww);
                    *added = 1;
                }
            }
        }
    }
    else if (next_diff_r < ACCEL_CLOSEST_R)
    {
        /* Overwrite the next cand */
        if (((accelcand *)(tmp_list->data))->sigma < sigma)
        {
            free_accelcand(tmp_list->data, NULL);
            tmp_list->data = (gpointer *)create_accelcand(power, sigma,
                                                          numharm, rr, zz, ww);
            *added = 1;
        }
    }
    else
    { /* This is a new candidate */
        new_list = g_slist_alloc();
        new_list->data =
            (gpointer *)create_accelcand(power, sigma, numharm, rr, zz, ww);
        *added = 1;
        if (!tmp_list->next &&
            (((accelcand *)(tmp_list->data))->r < (rr - ACCEL_CLOSEST_R)))
        {
            tmp_list->next = new_list;
            return list;
        }
        if (prev_list)
        {
            prev_list->next = new_list;
            new_list->next = tmp_list;
        }
        else
        {
            new_list->next = list;
            return new_list;
        }
    }
    return list;
}

static GSList *insert_new_accelcand_last(GSList *list, float power, float sigma,
                                         int numharm, double rr, double zz, double ww, int *added, GSList **tmp_list_ptr, GSList **prev_list_ptr)
/* Checks the current list to see if there is already */
/* a candidate within ACCEL_CLOSEST_R bins.  If not,  */
/* it adds it to the list in increasing freq order.   */
/* Record last insert position, search new cand from  */
/* head or last position.                             */
{
    GSList *new_list;
    double prev_diff_r = ACCEL_CLOSEST_R + 1.0, next_diff_r;

    *added = 0;
    if (!list)
    {
        *tmp_list_ptr = g_slist_alloc();
        (*tmp_list_ptr)->data =
            (gpointer *)create_accelcand(power, sigma, numharm, rr, zz, ww);
        *added = 1;
        return *tmp_list_ptr;
    }

    /* Initialize tmp_list and prev_list with the values from the pointers */
    GSList *tmp_list = *tmp_list_ptr;
    GSList *prev_list = *prev_list_ptr;

    /* Find the correct position in the list for the candidate */

    while ((tmp_list->next) && (((accelcand *)(tmp_list->data))->r < rr))
    {
        prev_list = tmp_list;
        tmp_list = tmp_list->next;
    }
    next_diff_r = fabs(rr - ((accelcand *)(tmp_list->data))->r);
    if (prev_list)
    {
        prev_diff_r = fabs(rr - ((accelcand *)(prev_list->data))->r);
    }

    /* Similar candidate(s) is(are) present */

    if (prev_diff_r < ACCEL_CLOSEST_R)
    {
        /* Overwrite the prev cand */
        if (((accelcand *)(prev_list->data))->sigma < sigma)
        {
            free_accelcand(prev_list->data, NULL);
            prev_list->data = (gpointer *)create_accelcand(power, sigma,
                                                           numharm, rr, zz, ww);
            *added = 1;
        }
        if (next_diff_r < ACCEL_CLOSEST_R)
        {
            if (((accelcand *)(tmp_list->data))->sigma < sigma)
            {
                free_accelcand(tmp_list->data, NULL);
                if (*added)
                {
                    /* Remove the next cand */
                    list = g_slist_remove_link(list, tmp_list);
                    g_slist_free_1(tmp_list);
                    tmp_list = prev_list->next;
                    if (!tmp_list)
                    {
                        tmp_list = list;
                        prev_list = NULL;
                    }
                }
                else
                {
                    /* Overwrite the next cand */
                    tmp_list->data = (gpointer *)create_accelcand(power, sigma,
                                                                  numharm, rr, zz, ww);
                    *added = 1;
                }
            }
        }
    }
    else if (next_diff_r < ACCEL_CLOSEST_R)
    {
        /* Overwrite the next cand */
        if (((accelcand *)(tmp_list->data))->sigma < sigma)
        {
            free_accelcand(tmp_list->data, NULL);
            tmp_list->data = (gpointer *)create_accelcand(power, sigma,
                                                          numharm, rr, zz, ww);
            *added = 1;
        }
    }
    else
    { /* This is a new candidate */
        new_list = g_slist_alloc();
        new_list->data =
            (gpointer *)create_accelcand(power, sigma, numharm, rr, zz, ww);
        *added = 1;
        if (!tmp_list->next &&
            (((accelcand *)(tmp_list->data))->r < (rr - ACCEL_CLOSEST_R)))
        {
            tmp_list->next = new_list;
            *tmp_list_ptr = tmp_list;
            *prev_list_ptr = prev_list;
            return list;
        }
        if (prev_list)
        {
            prev_list->next = new_list;
            new_list->next = tmp_list;
            *tmp_list_ptr = new_list;
            *prev_list_ptr = prev_list;
            return list;
        }
        else
        {
            new_list->next = list;
            *tmp_list_ptr = list;
            *prev_list_ptr = new_list;
            return new_list;
        }
    }
    *tmp_list_ptr = tmp_list;
    *prev_list_ptr = prev_list;
    return list;
}

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
    cudaStream_t sub_stream)
{
    int ii, jj, kk;

    long long fundamental_size = fundamental_numrs * fundamental_numzs * fundamental_numws;

    if (search_num <= 0)
        return cands;

    // record last insert position
    GSList *tmp_list = cands;
    GSList *prev_list = NULL;
    double last_rr = 0;

    for (int hitnum = 0; hitnum < search_num; hitnum++)
    {
        int current_batch = search_results[hitnum].index / single_batch_size;
        long long in_batch_id = search_results[hitnum].index % single_batch_size;
        int b = in_batch_id / (numharmstages * fundamental_size);
        long long act_id = in_batch_id % (numharmstages * fundamental_size);
        int stage = act_id / fundamental_size;
        int nh = 1 << stage;
        long long numindep = numindeps[stage];
        int fund_id = current_batch * proper_batch_size + b;
        int index = act_id % fundamental_size;
        ii = (int)(index / (fundamental_numzs * fundamental_numrs));
        jj = (int)((index / fundamental_numrs) % fundamental_numzs);
        kk = (int)(index % fundamental_numrs);

        float pow, sig;
        double rr, zz, ww;
        int added = 0;

        pow = search_results[hitnum].pow;
        sig = search_results[hitnum].sig;
        rr = (fundamental_rlos[fund_id] + kk * (double)ACCEL_DR) / (double)nh;
        zz = (fundamental_zlo + jj * (double)ACCEL_DZ) / (double)nh;
        ww = (fundamental_wlo + ii * (double)ACCEL_DW) / (double)nh;

        /*If the current rr >= the previous rr, continue searching from the last position;
        otherwise, start searching again from the beginning of cands list.*/
        if (rr <= last_rr || tmp_list == prev_list)
        {
            tmp_list = cands;
            prev_list = NULL;
        }

        {
            cands = insert_new_accelcand_last(cands, pow, sig, nh,
                                              rr, zz, ww, &added, &tmp_list, &prev_list);
        }
        last_rr = rr;
    }

    free(search_results);

    return cands;
}
