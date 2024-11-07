#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include "cuda_runtime.h"
#include "backend_common.h"
#include "presto.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))

void prepsubband_maxsize(int nsub, double lodm, double dmstep, int numdms, int downsamp, char *fitsfilename){
    int spectra_per_subint = 1024;
    struct spectra_info s;
    spectra_info_set_defaults(&s);
    s.filenames = &fitsfilename;
    s.num_files = 1;
    identify_psrdatatype(&s, 1);
    if (s.datatype != PSRFITS)
    {
        printf("prepsubband_maxsize only support *.fits\n");
    }
    read_rawdata_files(&s);

    double *dms = gen_dvect(numdms);
    for (int ii = 0; ii < numdms; ii++)
    {
        dms[ii] = lodm + ii * dmstep;
    }
    double maxdm = dms[numdms - 1];

    double BW_ddelay = delay_from_dm(maxdm, s.lo_freq) - delay_from_dm(maxdm, s.lo_freq + (s.num_channels - 1) * s.df);
    int blocksperread = ((int)(BW_ddelay / s.dt) / spectra_per_subint + 1);

    long long total_size = 0;

    long long cuda_init_size = 299892736; // 286 MB
    total_size += cuda_init_size;

    int gpu_worklen = s.spectra_per_subint * blocksperread;
    gpu_worklen /= downsamp;
    long long cudaMalloc_size = sizeof(float) * gpu_worklen * downsamp * nsub * 4 + sizeof(float) * gpu_worklen * nsub * 2 + sizeof(float) * gpu_worklen * numdms + sizeof(int) * numdms * nsub;
    total_size += cudaMalloc_size;

    double td = total_size / (1024.0 * 1024.0);
    printf("%.3f\n", td);
}

int main(int argc, char *argv[])
{
    int nsub = atoi(argv[1]);
    double lodm = atof(argv[2]);
    double dmstep = atof(argv[3]);
    int numdms = atoi(argv[4]);
    int downsamp = atoi(argv[5]);
    char *fitsfilename = argv[6];
    prepsubband_maxsize(nsub, lodm, dmstep, numdms, downsamp, fitsfilename);
    return 0;
}