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

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

#include "accel.h"
#include "device_launch_parameters.h"
#include "cufft.h"
#include "cuda_helper.h"
#include <math.h>

__constant__ int subw_device[12006];
__constant__ float powcuts_device[6];
__constant__ int numharms_device[6];
__constant__ double numindeps_device[6];

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

// Structures to store search results
typedef struct
{
	long long index;
	float pow;
	float sig;
} SearchValue;

extern "C" extern float ***gen_f3Darr_cu(long nhgts, long nrows, long ncols, cudaStream_t stream);

extern "C" extern void do_fft_batch(int fftlen, int binoffset, ffdotpows_cu *ffdot_array, subharminfo *shi, fcomplex *pdata_array, int *idx_array, fcomplex *full_tmpdat_array, fcomplex *full_tmpout_array, int batch_size, fcomplex *fkern, cudaStream_t stream);

extern "C" extern void fuse_add_search_batch(ffdotpows_cu *fundamentals, SubharmonicMap *subhmap, int stages, int fundamental_num, cudaStream_t stream, SearchValue *search_results, unsigned long long int *search_nums, long long pre_size, int proper_batch_size, int max_searchnum, int *too_large);

extern "C" extern kernel **gen_kernmatrix_cu(int numz, int numw);

extern "C" extern void sort_search_results(SearchValue *search_results, unsigned long long int search_num);

extern "C" extern void init_constant_device(int *subw_host, int subw_size, float *powcuts_host, int *numharms_host, double *numindeps_host, int numharmstages_size);

extern "C" extern void clear_cache();

void init_constant_device(int *subw_host, int subw_size, float *powcuts_host, int *numharms_host, double *numindeps_host, int numharmstages_size)
{
	CUDA_CHECK(cudaMemcpyToSymbol(subw_device, subw_host, subw_size * sizeof(int), 0, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpyToSymbol(powcuts_device, powcuts_host, numharmstages_size * sizeof(float), 0, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpyToSymbol(numharms_device, numharms_host, numharmstages_size * sizeof(int), 0, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpyToSymbol(numindeps_device, numindeps_host, numharmstages_size * sizeof(double), 0, cudaMemcpyHostToDevice));
}

#define NEAREST_INT(x) (int)(x < 0 ? x - 0.5 : x + 0.5)

#define CUDA_THREAD_BLOCK_NUM 64

#define CHECK_CUFFT_ERRORS(call)                                                         \
	{                                                                                    \
		cufftResult_t err;                                                               \
		if ((err = (call)) != CUFFT_SUCCESS)                                             \
		{                                                                                \
			fprintf(stderr, "cuFFT error %d:%s at %s:%d\n", err, _cudaGetErrorEnum(err), \
					__FILE__, __LINE__);                                                 \
			exit(1);                                                                     \
		}                                                                                \
	}

#define MAX_EMBED_DIMS 3

typedef struct
{
	int rank;
	int n[1];
	int inembed[MAX_EMBED_DIMS];
	int istride;
	int idist;
	int onembed[MAX_EMBED_DIMS];
	int ostride;
	int odist;
	int batch;
} CufftParams;

#define MAX_CACHE_SIZE 3000

typedef struct
{
	CufftParams key;
	cufftHandle value;
} CufftCacheEntry;

CufftCacheEntry cache[MAX_CACHE_SIZE];
int cacheIndex = 0;
int cacheSize = 0;

bool find_in_cache(const CufftParams *params, cufftHandle *cu_plan)
{
	for (int i = 0; i < cacheSize; ++i)
	{
		if (memcmp(params, &cache[i].key, sizeof(CufftParams)) == 0)
		{
			*cu_plan = cache[i].value;
			return true;
		}
	}
	return false;
}

void add_to_cache(const CufftParams *params, cufftHandle cu_plan)
{
	// Replace the earliest entry
	cache[cacheIndex].key = *params;
	cache[cacheIndex].value = cu_plan;

	// Update cache index
	cacheIndex = (cacheIndex + 1) % MAX_CACHE_SIZE;
	if (cacheSize < MAX_CACHE_SIZE)
	{
		cacheSize += 1;
	}
}

void clear_cache()
{
	for (int i = 0; i < cacheSize; ++i)
	{
		cufftDestroy(cache[i].value);
	}
	cacheIndex = 0;
	cacheSize = 0;
}

static const char *_cudaGetErrorEnum(cufftResult error)
{
	switch (error)
	{
	case CUFFT_SUCCESS:
		return "CUFFT_SUCCESS";

	case CUFFT_INVALID_PLAN:
		return "CUFFT_INVALID_PLAN";

	case CUFFT_ALLOC_FAILED:
		return "CUFFT_ALLOC_FAILED";

	case CUFFT_INVALID_TYPE:
		return "CUFFT_INVALID_TYPE";

	case CUFFT_INVALID_VALUE:
		return "CUFFT_INVALID_VALUE";

	case CUFFT_INTERNAL_ERROR:
		return "CUFFT_INTERNAL_ERROR";

	case CUFFT_EXEC_FAILED:
		return "CUFFT_EXEC_FAILED";

	case CUFFT_SETUP_FAILED:
		return "CUFFT_SETUP_FAILED";

	case CUFFT_INVALID_SIZE:
		return "CUFFT_INVALID_SIZE";

	case CUFFT_UNALIGNED_DATA:
		return "CUFFT_UNALIGNED_DATA";
	}

	return "<unknown>";
}

__device__ __forceinline__ static size_t calc_index_dev(int ws, int zs, int ws_len, int zs_len, int single_element_size)
{
	size_t index = (zs * ws_len + ws) * single_element_size;
	return index;
}

__global__ void powers_init(long nhgts, long nrows, long ncols, float ***c)
{
	long i, j;
	for (j = 1; j < nrows; j++)
	{
		c[0][j] = c[0][j - 1] + ncols;
	}

	for (i = 1; i < nhgts; i++)
	{
		c[i] = c[i - 1] + nrows;
		c[i][0] = c[i - 1][0] + nrows * ncols;
		for (j = 1; j < nrows; j++)
		{
			c[i][j] = c[i][j - 1] + ncols;
		}
	}
}

float ***gen_f3Darr_cu(long nhgts, long nrows, long ncols, cudaStream_t stream)
{
	float ***c;

	CUDA_CHECK(cudaMallocManaged(&c, (size_t)(nhgts * sizeof(float **))));

	CUDA_CHECK(cudaMallocManaged(&c[0], (size_t)((nhgts * nrows) * sizeof(float *))));

	CUDA_CHECK(cudaMallocManaged(&c[0][0], (size_t)((nhgts * nrows * ncols) * sizeof(float))));

	powers_init<<<1, 1, 0, stream>>>(nhgts, nrows, ncols, c);

	CUDA_CHECK(cudaGetLastError());

	return c;
}

__device__ __forceinline__ int calc_required_w_cu(double harm_fract, double wfull)
/* Calculate the maximum 'w' needed for the given subharmonic  */
/* harm_fract = harmnum / numharm if the       */
/* 'w' at the fundamental harmonic is 'wfull'. */
{
	return NEAREST_INT(ACCEL_RDW * wfull * harm_fract) * ACCEL_DW;
}

__device__ __forceinline__ int index_from_w_cu(double w, double low)
/* Return an index for a Fourier Fdotdot given an array that */
/* has stepsize ACCEL_DW and low freq dotdot 'low'.              */
{
	return (int)((w - low) * ACCEL_RDW + DBLCORRECT);
}

__device__ double devlpl(double a[], double *x)
{
	static double devlpl, term;
	static int i;
	term = a[4];
	for (i = 3; i >= 0; i--)
		term = a[i] + term * *x;
	devlpl = term;
	return devlpl;
}

__device__ double stvaln(double *p)
{
	double xden[5] = {
		0.993484626060e-1, 0.588581570495e0, 0.531103462366e0, 0.103537752850e0,
		0.38560700634e-2};
	double xnum[5] = {
		-0.322232431088e0, -1.000000000000e0, -0.342242088547e0, -0.204231210245e-1,
		-0.453642210148e-4};
	// int K1 = 5;
	double sign, y, z;

	if (*p <= 0.5e0)
	{
		sign = -1.0e0;
		z = *p;
	}
	else
	{
		sign = 1.0e0;
		z = 1.0e0 - *p;
	}

	y = __sqrtf(-(2.0e0 * __logf(z)));
	double result = y + devlpl(xnum, &y) / devlpl(xden, &y);
	return sign * result;
}

__device__ double fifdint(double a)
/* a     -     number to be truncated */
{
	long temp;
	temp = (long)(a);
	return (double)(temp);
}

__device__ void cumnor(double *arg, double *result, double *ccum)
{
	static double a[5] = {
		2.2352520354606839287e00, 1.6102823106855587881e02, 1.0676894854603709582e03,
		1.8154981253343561249e04, 6.5682337918207449113e-2};
	static double b[4] = {
		4.7202581904688241870e01, 9.7609855173777669322e02, 1.0260932208618978205e04,
		4.5507789335026729956e04};
	static double c[9] = {
		3.9894151208813466764e-1, 8.8831497943883759412e00, 9.3506656132177855979e01,
		5.9727027639480026226e02, 2.4945375852903726711e03, 6.8481904505362823326e03,
		1.1602651437647350124e04, 9.8427148383839780218e03, 1.0765576773720192317e-8};
	static double d[8] = {
		2.2266688044328115691e01, 2.3538790178262499861e02, 1.5193775994075548050e03,
		6.4855582982667607550e03, 1.8615571640885098091e04, 3.4900952721145977266e04,
		3.8912003286093271411e04, 1.9685429676859990727e04};
	static double half = 0.5e0;
	static double p[6] = {
		2.1589853405795699e-1, 1.274011611602473639e-1, 2.2235277870649807e-2,
		1.421619193227893466e-3, 2.9112874951168792e-5, 2.307344176494017303e-2};
	static double one = 1.0e0;
	static double q[5] = {
		1.28426009614491121e00, 4.68238212480865118e-1, 6.59881378689285515e-2,
		3.78239633202758244e-3, 7.29751555083966205e-5};
	static double sixten = 1.60e0;
	static double sqrpi = 3.9894228040143267794e-1;
	static double thrsh = 0.66291e0;
	static double root32 = 5.656854248e0;
	static double zero = 0.0e0;
	// static int K1 = 1;
	// static int K2 = 2;
	static int i;
	static double del, eps, temp, x, xden, xnum, y, xsq, min;
	/*
	------------------------------------------------------------------
	  Machine dependent constants
	------------------------------------------------------------------
	*/
	eps = __powf(2, (double)(-52)) * 0.5e0;
	double binv = 1.0 / 2;
	min = __powf(2, (double)(-1019)) * binv * binv * binv;
	x = *arg;
	y = fabs(x);
	if (y <= thrsh)
	{
		/*
		------------------------------------------------------------------
		  Evaluate  anorm  for  |X| <= 0.66291
		------------------------------------------------------------------
		*/
		xsq = zero;
		if (y > eps)
			xsq = x * x;
		xnum = a[4] * xsq;
		xden = xsq;
		for (i = 0; i < 3; i++)
		{
			xnum = (xnum + a[i]) * xsq;
			xden = (xden + b[i]) * xsq;
		}
		*result = x * (xnum + a[3]) / (xden + b[3]);
		temp = *result;
		*result = half + temp;
		*ccum = half - temp;
	}
	/*
	------------------------------------------------------------------
	  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
	------------------------------------------------------------------
	*/
	else if (y <= root32)
	{
		xnum = c[8] * y;
		xden = y;
		for (i = 0; i < 7; i++)
		{
			xnum = (xnum + c[i]) * y;
			xden = (xden + d[i]) * y;
		}
		*result = (xnum + c[7]) / (xden + d[7]);
		xsq = fifdint(y * sixten) / sixten;
		del = (y - xsq) * (y + xsq);
		*result = __expf(-(xsq * xsq * half)) * __expf(-(del * half)) * *result;
		*ccum = one - *result;
		if (x > zero)
		{
			temp = *result;
			*result = *ccum;
			*ccum = temp;
		}
	}
	/*
	------------------------------------------------------------------
	  Evaluate  anorm  for |X| > sqrt(32)
	------------------------------------------------------------------
	*/
	else
	{
		*result = zero;
		xsq = one / (x * x);
		xnum = p[5] * xsq;
		xden = xsq;
		for (i = 0; i < 4; i++)
		{
			xnum = (xnum + p[i]) * xsq;
			xden = (xden + q[i]) * xsq;
		}
		*result = xsq * (xnum + p[4]) / (xden + q[4]);
		*result = (sqrpi - *result) / y;
		xsq = fifdint(x * sixten) / sixten;
		del = (x - xsq) * (x + xsq);
		*result = __expf(-(xsq * xsq * half)) * __expf(-(del * half)) * *result;
		*ccum = one - *result;
		if (x > zero)
		{
			temp = *result;
			*result = *ccum;
			*ccum = temp;
		}
	}
	if (*result < min)
		*result = 0.0e0;
	if (*ccum < min)
		*ccum = 0.0e0;
}

__device__ double dinvnr(double *p, double *q)
{
	const int maxit = 100;
	const double eps = 1.0e-13;
	const double r2pi = 0.3989422804014326e0;
	const double nhalf = -0.5e0;
#define dennor(x) (r2pi * __expf(nhalf * (x) * (x)))

	double pp = (*p <= *q) ? *p : *q; // Find the minimum value between p and q
	double strtx = stvaln(&pp);
	double xcur = strtx;
	double cum, ccum, dx;
	double result;

	// NEWTON ITERATIONS
	for (int i = 1; i <= maxit; i++)
	{
		cumnor(&xcur, &cum, &ccum);
		dx = (cum - pp) / dennor(xcur);
		xcur -= dx;
		if (fabs(dx / xcur) < eps)
		{
			result = xcur; // If the NEWTON iteration succeeds
			break;
		}

		if (i == maxit)
		{ // If the loop ends and the accuracy requirement is still not met, the NEWTON iteration is considered to have failed
			result = strtx;
		}
	}

	// In the case where p>q, the sign of the result needs to be adjusted
	if (*p > *q)
	{
		result = -result;
	}

	return result;
#undef dennor
}

__device__ void cdfnor2(double *p, double *q, double *x)
{
	double pq = *p + *q;

	double z = dinvnr(p, q);
	*x = z;
}

__device__ double extended_equiv_gaussian_sigma(double logp)
{
	double t, num, denom;

	t = sqrt(-2.0 * logp);
	num = 2.515517 + t * (0.802853 + t * 0.010328);
	denom = 1.0 + t * (1.432788 + t * (0.189269 + t * 0.001308));
	return t - num / denom;
}

__device__ double equivalent_gaussian_sigma2(double logp)
/* Return the approximate significance in Gaussian sigmas */
/* corresponding to a natural log probability logp        */
{
	double x;
	if (logp < -600.0)
	{
		x = extended_equiv_gaussian_sigma(logp);
	}
	else
	{
		// int status;
		double p, q;
		q = __expf(logp);
		p = 1.0 - q;
		x = dinvnr(&p, &q);
	}
	if (x < 0.0)
		return 0.0;
	else
		return x;
}

__device__ float log_asymtotic_incomplete_gamma(float a, float z)
/*
  log_asymtotic_incomplete_gamma(double a, double z):
	  Return the natural log of the incomplete gamma function in
		  its asymtotic limit as z->infty.  This is from Abramowitz
		  and Stegun eqn 6.5.32.
*/
{
	float x = 1.0, newxpart = 1.0, term = 1.0;
	int ii = 1;

	while (fabs(newxpart) > 1e-15)
	{
		term *= (a - ii);
		newxpart = term / __powf(z, ii);
		x += newxpart;
		ii += 1;
	}
	return (a - 1.0) * __logf(z) - z + __logf(x);
}

__device__ float log_asymtotic_gamma(float z)
/*
  log_asymtotic_gamma(double z):
	  Return the natural log of the gamma function in its asymtotic limit
		  as z->infty.  This is from Abramowitz and Stegun eqn 6.1.41.
*/
{
	float x, y;

	x = (z - 0.5) * __logf(z) - z + 0.91893853320467267;
	y = 1.0 / (z * z);
	x += (((-5.9523809523809529e-4 * y + 7.9365079365079365079365e-4) * y - 2.7777777777777777777778e-3) * y + 8.3333333333333333333333e-2) / z;
	return x;
}

__device__ double fifdmax1(double a, double b)
/* a     -      first number */
/* b     -      second number */
{
	if (a < b)
		return b;
	else
		return a;
}

__device__ double Xgamm(double *a)
{
	double d = .41893853320467274178e0;
	double pi = 3.1415926535898e0;
	double r1 = .820756370353826e-03;
	double r2 = -.595156336428591e-03;
	double r3 = .793650663183693e-03;
	double r4 = -.277777777770481e-02;
	double r5 = .833333333333333e-01;
	double p[7] = {
		.539637273585445e-03, .261939260042690e-02, .204493667594920e-01,
		.730981088720487e-01, .279648642639792e+00, .553413866010467e+00, 1.0e0};
	double q[7] = {
		-.832979206704073e-03, .470059485860584e-02, .225211131035340e-01,
		-.170458969313360e+00, -.567902761974940e-01, .113062953091122e+01, 1.0e0};
	double result, bot, g, lnx, s, t, top, w, x, z;

	x = *a;
	if (fabs(*a) < 15.0e0)
	{
		t = 1.0e0;
		for (int j = 1; j <= (int)(*a) - 1; ++j)
		{
			x -= 1.0e0;
			t *= x;
		}

		top = p[0];
		bot = q[0];
		for (int i = 1; i < 7; ++i)
		{
			top = p[i] + x * top;
			bot = q[i] + x * bot;
		}
		result = (*a < 1.0e0) ? (top / bot) / t : t * (top / bot);
	}
	else
	{
		if (fabs(*a) >= 1.e3)
			return 0.0;
		t = 1.0e0 / (x * x);
		g = ((((r1 * t + r2) * t + r3) * t + r4) * t + r5) / x;
		lnx = __logf(x);
		z = x;
		g = d + g + (z - 0.5e0) * (lnx - 1.e0);
		result = __expf(g);
		if (*a < 0.0e0)
			result = 1.0e0 / (result * s) / x;
	}
	return result;
}

__device__ long fifidint(double a)
/* a - number to be truncated */
{
	return (long)(a);
}

__device__ double gam1(double *a)
{
	double s1 = .273076135303957e+00;
	double s2 = .559398236957378e-01;
	double p[7] = {
		.577215664901533e+00, -.409078193005776e+00, -.230975380857675e+00,
		.597275330452234e-01, .766968181649490e-02, -.514889771323592e-02,
		.589597428611429e-03};
	double q[5] = {
		.100000000000000e+01, .427569613095214e+00, .158451672430138e+00,
		.261132021441447e-01, .423244297896961e-02};

	double gam1, bot, t, top, w;

	double x = *a;
	x = x - 1.0;
	if (x == 0)
	{
		gam1 = 0.0e0;
	}
	else
	{
		top = (((((p[6] * x + p[5]) * x + p[4]) * x + p[3]) * x + p[2]) * x + p[1]) * x + p[0];
		bot = (((q[4] * x + q[3]) * x + q[2]) * x + q[1]) * x + 1.0e0;
		w = top / bot;
		gam1 = x / *a * (w - 0.5e0 - 0.5e0);
	}

	return gam1;
}

__device__ void gratio2(double *a, double *x, double *ans, double *qans)
{
	static int i, max, m, n;
	static double wk[20];
	static double e, acc, e0, x0, t1, u, r, amn, t, sum, c;

	e = __powf(2, (double)(-52));

	if (*a * *x == 0.0e0)
	{
		*ans = 1.0e0;
		*qans = 0.0e0;
		return;
	}
	acc = fifdmax1(5.e-15, e);
	e0 = .25e-3;
	x0 = 31.0e0;

	if (*a >= 20.0e0 && (*a > *x || *x >= x0))
	{
		t1 = *a * __logf(*x) - *x;
		r = __expf(t1) / Xgamm(a);
	}
	else
	{
		i = (int)(*a);
		sum = exp(-*x);
		t = sum;
		n = 1;
		c = 0.0e0;
		while (n != i)
		{
			n += 1;
			c += 1.0e0;
			t = *x * t / c;
			sum += t;
		}
		*qans = sum;
		*ans = 0.5e0 + (0.5e0 - *qans);
		return;
	}
	if (r != 0.0e0)
	{
		amn = *a - 1.0e0;
		t = amn / *x;
		wk[0] = t;
		for (n = 2; n <= 20; n++)
		{
			amn -= 1.0e0;
			t *= (amn / *x);
			if (fabs(t) <= 1.e-3)
				break;
			wk[n - 1] = t;
		}
		sum = t;
		while (fabs(t) > acc)
		{
			amn -= 1.0e0;
			t *= (amn / *x);
			sum += t;
		}
		max = n - 1;
		for (m = 1; m <= max; m++)
		{
			n -= 1;
			sum += wk[n - 1];
		}
		*qans = r / *x * (1.0e0 + sum);
		*ans = 0.5e0 + (0.5e0 - *qans);
	}
}

__device__ void cdfchi2(double *p, double *q, double *x, double *df)
{
	static double a, xx;
	a = *df * 0.5e0;
	xx = *x * 0.5e0;
	if (xx <= 0.0e0)
	{
		*p = 0.0e0;
		*q = 1.0e0;
	}
	else
	{
		// x is greater than 0, call the gratio routine
		gratio2(&a, &xx, p, q);
	}
}

__device__ float chi2_logp(float chi2, float dof)
/* Return the natural log probability corresponding to a chi^2 value */
/* of chi2 given dof degrees of freedom. */
{
	float logp;

	if (chi2 / dof > 15.0 || (dof > 150 && chi2 / dof > 6.0))
	{
		logp = log_asymtotic_incomplete_gamma(0.5 * dof, 0.5 * chi2) -
			   log_asymtotic_gamma(0.5 * dof);
	}
	else
	{
		double p, q, df = dof, x = chi2;

		cdfchi2(&p, &q, &x, &df);
		logp = log(q);
	}
	return logp;
}

__device__ float candidate_sigma_cu(float power, int numsum, float numtrials)
{
	float sigma;
	float logp, chi2, dof;

	if (power <= 0.0)
	{
		sigma = 0.0;
	}
	else
	{
		// Get the natural log probability
		chi2 = 2.0 * power;
		dof = 2.0 * numsum;
		logp = chi2_logp(chi2, dof);

		// Correct for numtrials
		logp += __logf(numtrials);

		// Convert to sigma
		sigma = equivalent_gaussian_sigma2(logp);
	}

	return sigma;
}

__global__ void fuse_add_search_batch_kernel(
	int fundamental_numrs,
	int fundamental_numzs,
	int fundamental_numws,
	float *fundamental_powers_flat,
	SubharmonicMap *subhmap,
	int stages,
	int fundamental_num,
	long long fundamental_size,
	SearchValue *search_results,
	unsigned long long int *search_nums,
	long long pre_size,
	int proper_batch_size,
	int max_searchnum,
	int *too_large)
{
	int f = blockIdx.x / fundamental_numws;
	int ii = blockIdx.x % fundamental_numws;
	int jj = blockIdx.y;
	int kk = blockIdx.z * blockDim.x + threadIdx.x;

	if (kk < fundamental_numrs && *too_large == 0)
	{
		int fundamental_index = matrix_3d_index(ii, jj, kk, fundamental_numzs, fundamental_numrs);
		float tmp = fundamental_powers_flat[f * fundamental_size + fundamental_index];

		int stage = 0;
		if (tmp > powcuts_device[stage])
		{
			unsigned long long int index = atomicAdd(&search_nums[0], 1ULL);
			if (index >= max_searchnum)
			{
				*too_large = 1;
				return;
			}
			search_results[index].index = (long long)(pre_size) + (long long)(stage * fundamental_size + f * (stages + 1) * fundamental_size + fundamental_index);
			search_results[index].pow = tmp;
			float sig = candidate_sigma_cu(tmp, numharms_device[stage], numindeps_device[stage]);
			search_results[index].sig = sig;
		}

		int pre = 0;
		for (stage = 1; stage <= stages; stage++)
		{
			int harmtosum = 1 << (stage - 1);
			for (int b = 0; b < harmtosum; b++)
			{
				SubharmonicMap subh = subhmap[(b + pre) * fundamental_num + f];
				int harm_fract = subh.harm_fract;
				float *subharmonic_powers_flat = subh.subharmonic_powers;
				int subharmonic_wlo = subh.subharmonic_wlo;
				unsigned short *subharmonic_zinds = subh.subharmonic_zinds;
				unsigned short *subharmonic_rinds = subh.subharmonic_rinds;
				int subw = subw_device[harm_fract * fundamental_numws + ii];
				int wind = ((subw - subharmonic_wlo) / ACCEL_DW);

				int zind = subharmonic_zinds[jj];
				int rind = subharmonic_rinds[kk];
				int subharmonic_index = matrix_3d_index(wind, zind, rind, subh.subharmonic_numzs, subh.subharmonic_numrs);
				tmp += subharmonic_powers_flat[subharmonic_index];
			}
			pre += harmtosum;

			if (tmp > powcuts_device[stage])
			{
				unsigned long long int index = atomicAdd(&search_nums[0], 1ULL);
				if (index >= max_searchnum)
				{
					*too_large = 1;
					return;
				}
				else
				{
					search_results[index].index = (long long)(pre_size) + (long long)(stage * fundamental_size + f * (stages + 1) * fundamental_size + fundamental_index);
					search_results[index].pow = tmp;
					float sig = candidate_sigma_cu(tmp, numharms_device[stage], numindeps_device[stage]);
					search_results[index].sig = sig;
				}
			}
		}
	}
}

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
						   int *too_large)
{
	int threads = 128;
	ffdotpows_cu *fundamental = &fundamentals[0];
	dim3 gridDim(fundamental->numws * fundamental_num, fundamental->numzs, (fundamental->numrs + threads - 1) / threads);

	long long fundamental_size = fundamental->numws * fundamental->numzs * fundamental->numrs;

	fuse_add_search_batch_kernel<<<gridDim, threads, 0, stream>>>(
		fundamental->numrs,
		fundamental->numzs,
		fundamental->numws,
		fundamental->powers,
		subhmap,
		stages,
		fundamental_num,
		fundamental_size,
		search_results,
		search_nums,
		pre_size,
		proper_batch_size,
		max_searchnum,
		too_large);

	CUDA_CHECK(cudaGetLastError());
}

kernel **gen_kernmatrix_cu(int numz, int numw)
{
	kernel **kerns;

	CUDA_CHECK(cudaMallocManaged(&kerns, (size_t)numw * sizeof(kernel *)));

	if (!kerns)
	{
		perror("\nError in 1st malloc() in gen_kernmatrix()");
		printf("\n");
		exit(-1);
	}

	CUDA_CHECK(cudaMallocManaged(&kerns[0], (size_t)((numz * numw) * sizeof(kernel))));

	if (!kerns[0])
	{
		perror("\nError in 2nd malloc() in gen_kernmatrix()");
		printf("\n");
		exit(-1);
	}
	for (int ii = 1; ii < numw; ii++)
	{
		kerns[ii] = kerns[ii - 1] + numz;
	}

	return kerns;
}

__global__ void pre_fft_kernel_batch_float4(fcomplex *pdata_array, fcomplex *full_tmpdat_array, fcomplex *fkern_gpu,
											int batch_size, int ws_len, int zs_len, int fftlen)
{
	int b = blockIdx.x / ws_len; // batch index
	int ii = blockIdx.x % ws_len;
	int jj = blockIdx.y;
	int kk = blockIdx.z * blockDim.x + threadIdx.x;

	if (kk < fftlen / 2) // Since float4 is used, the effective length is half of the original
	{
		float4 *pdata = (float4 *)(&pdata_array[b * fftlen]);
		float4 *fdata = (float4 *)(&full_tmpdat_array[b * (fftlen * ws_len * zs_len) + calc_index_dev(ii, jj, ws_len, zs_len, fftlen)]);
		float4 *fkern = (float4 *)(&fkern_gpu[(ii * zs_len + jj) * fftlen]);

		// Read the float4 in pdata and fkern, each contains two complex numbers
		float4 p = pdata[kk];
		float4 k = fkern[kk];

		// Perform complex multiplication
		float2 p1 = make_float2(p.x, p.y);
		float2 k1 = make_float2(k.x, k.y);
		float2 p2 = make_float2(p.z, p.w);
		float2 k2 = make_float2(k.z, k.w);

		// First complex multiplication
		float2 res1 = make_float2(p1.x * k1.x + p1.y * k1.y, p1.y * k1.x - p1.x * k1.y);
		// Second complex multiplication
		float2 res2 = make_float2(p2.x * k2.x + p2.y * k2.y, p2.y * k2.x - p2.x * k2.y);

		// Store results into fdata
		fdata[kk] = make_float4(res1.x, res1.y, res2.x, res2.y);
	}
}

__global__ void after_fft_kernel_batch(float *powers, fcomplex *full_tmpout_array, int offset, float norm, int ws_len, int zs_len, int *rs_len_array, int fftlen, int *idx_array, int max_rs_len)
{
	int b = blockIdx.x / ws_len;
	int ii = blockIdx.x % ws_len;
	int linear_index = blockIdx.y * blockDim.x + threadIdx.x;
	int jj = linear_index / (max_rs_len);
	int kk = linear_index % (max_rs_len);

	int rs_len = rs_len_array[b];

	if (jj < zs_len && kk < rs_len)
	{
		fcomplex *full_tmpout = &full_tmpout_array[b * (fftlen * ws_len * zs_len)];
		fcomplex *output = &full_tmpout[calc_index_dev(ii, jj, ws_len, zs_len, fftlen)];

		// Turn the good parts of the result into powers and store
		// them in the output matrix
		float *fdata = (float *)output;
		const int ind = 2 * (kk + offset);
		int index = matrix_3d_index(ii, jj, kk, zs_len, rs_len);
		powers[idx_array[b] + index] = (fdata[ind] * fdata[ind] +
										fdata[ind + 1] * fdata[ind + 1]) *
									   norm;
	}
}

// Define a comparison function object for comparing the index fields of two SearchValue structures
struct CompareSearchValue
{
	__host__ __device__ bool operator()(const SearchValue &a, const SearchValue &b) const
	{
		return a.index < b.index;
	}
};

void sort_search_results(SearchValue *search_results, unsigned long long int search_num)
{
	// Wrap the raw pointer with thrust::device_ptr so that it can be used in thrust algorithms
	thrust::device_ptr<SearchValue> dev_ptr(search_results);

	// Use thrust sorting with the custom comparator
	thrust::sort(thrust::device, dev_ptr, dev_ptr + search_num, CompareSearchValue());
}

void do_fft_batch(int fftlen, int binoffset, ffdotpows_cu *ffdot_array, subharminfo *shi, fcomplex *pdata_array, int *idx_array,
				  fcomplex *full_tmpdat_array, fcomplex *full_tmpout_array, int batch_size, fcomplex *fkern, cudaStream_t stream)
{
	int ws_len_global = ffdot_array[0].numws;
	int zs_len_global = ffdot_array[0].numzs;

	cufftHandle *cu_plan_array = (cufftHandle *)malloc(batch_size * sizeof(cufftHandle));

	// 1. handle operations before the pre_fft_kernel call, prepare cu_plan
	for (int b = 0; b < batch_size; b++)
	{
		ffdotpows_cu *ffdot = &ffdot_array[b];
		const float norm = 1.0 / (fftlen * fftlen);
		const int offset = binoffset * ACCEL_NUMBETWEEN;
		int ws_len = ffdot->numws;
		int zs_len = ffdot->numzs;
		int rs_len = ffdot->numrs;

		size_t full_size = ws_len * zs_len;

		cufftHandle cu_plan;
		int rank = 1;
		int n[1] = {fftlen};
		int istride = 1, idist = fftlen;
		int ostride = 1, odist = fftlen;
		int inembed[3] = {ws_len, zs_len, fftlen};
		int onembed[3] = {ws_len, zs_len, fftlen};

		CufftParams params;
		params.rank = rank;
		params.n[0] = fftlen;
		params.inembed[0] = ws_len;
		params.inembed[1] = zs_len;
		params.inembed[2] = fftlen;
		params.onembed[0] = ws_len;
		params.onembed[1] = zs_len;
		params.onembed[2] = fftlen;
		params.istride = istride;
		params.ostride = ostride;
		params.idist = fftlen;
		params.odist = fftlen;
		params.batch = ws_len * zs_len;

		if (!find_in_cache(&params, &cu_plan))
		{
			CHECK_CUFFT_ERRORS(cufftPlanMany(&cu_plan, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, ws_len * zs_len));
			add_to_cache(&params, cu_plan);
		}
		cu_plan_array[b] = cu_plan;
		CHECK_CUFFT_ERRORS(cufftSetStream(cu_plan, stream));
	}

	// 2. run pre_fft_kernel
	int threads_pre = 512;
	dim3 blocks_pre(batch_size * ws_len_global, zs_len_global, (fftlen / 2 + threads_pre - 1) / threads_pre);

	pre_fft_kernel_batch_float4<<<blocks_pre, threads_pre, 0, stream>>>(
		pdata_array, full_tmpdat_array, fkern, batch_size, ws_len_global, zs_len_global, fftlen);

	CUDA_CHECK(cudaGetLastError());

	// 3. run cufftExecC2C
	for (int b = 0; b < batch_size; b++)
	{
		fcomplex *full_tmpdat = &full_tmpdat_array[b * (fftlen * ws_len_global * zs_len_global)];
		fcomplex *full_tmpout = &full_tmpout_array[b * (fftlen * ws_len_global * zs_len_global)];

		CHECK_CUFFT_ERRORS(cufftExecC2C(cu_plan_array[b], (cufftComplex *)full_tmpdat, (cufftComplex *)full_tmpout, CUFFT_INVERSE));
	}

	// 4. prepare to run after_fft_kernel
	int max_rs_len = 0;
	int *rs_len_array = (int *)malloc(batch_size * sizeof(int));
	for (int b = 0; b < batch_size; b++)
	{
		ffdotpows_cu *ffdot = &ffdot_array[b];
		int rs_len = ffdot->numrs;
		if (rs_len > max_rs_len)
		{
			max_rs_len = rs_len;
		}
		rs_len_array[b] = rs_len;
	}
	int *idx_array_device, *rs_len_array_device;

	// Allocate memory for the array on the device
	cudaMallocAsync((void **)&idx_array_device, batch_size * sizeof(int), stream);
	cudaMallocAsync((void **)&rs_len_array_device, batch_size * sizeof(int), stream);

	// Copy data from host to device
	cudaMemcpyAsync(idx_array_device, idx_array, batch_size * sizeof(int), cudaMemcpyHostToDevice, stream);
	cudaMemcpyAsync(rs_len_array_device, rs_len_array, batch_size * sizeof(int), cudaMemcpyHostToDevice, stream);

	// 5. run after_fft_kernel
	const int offset = binoffset * ACCEL_NUMBETWEEN;
	const float norm = 1.0 / (fftlen * fftlen);
	int threads_after = 512;
	dim3 blocks_after(ws_len_global * batch_size, (zs_len_global * max_rs_len + threads_after - 1) / threads_after, 1);

	after_fft_kernel_batch<<<blocks_after, threads_after, 0, stream>>>(ffdot_array[0].powers, full_tmpout_array, offset, norm, ws_len_global, zs_len_global, rs_len_array_device, fftlen, idx_array_device, max_rs_len);

	// Release the allocated resources
	free(cu_plan_array);
	free(rs_len_array);
	cudaFreeAsync(idx_array_device, stream);
	cudaFreeAsync(rs_len_array_device, stream);
}
