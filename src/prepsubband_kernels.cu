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
#include "cuda_helper.h"

extern "C" void DownsampDataArray(float *currentdsdata_gpu, float *currentdata_gpu, int worklen, int bs, int nsub);
extern "C" void getDedisp(float *currentdsdata_gpu, float *lastdsdata_gpu, float *outdata_gpu, int worklen, int nsub, int numdms, int *offsets);

__global__ void DownsampDataArray_kernel(float *currentdsdata_gpu, float *currentdata_gpu, int nsub, int worklen, int bs)
{
    int ThrPerBlk = blockDim.x;
    int MYbid = blockIdx.x;
    int MYtid = threadIdx.x;
    int MYgtid = ThrPerBlk * MYbid + MYtid;
    if (MYgtid >= worklen * nsub)
        return;

    int j;
    float val = 0;
    int jj = (int)(MYgtid / worklen);
    int ii = MYgtid - jj * worklen;

    int index = jj + ii * bs * nsub;
    int dsindex = jj + ii * nsub;
    for (j = 0; j < bs; j++)
    {
        val += (currentdata_gpu[index]);
        index += nsub;
    }
    currentdsdata_gpu[dsindex] = val / bs;
}

__global__ void getDedispGPU(float *currentdsdata_gpu, float *lastdsdata_gpu, float *outdata_gpu, int worklen, int nsub, int numdms, int *offsets)
{
    int i, j, b;
    int ThrPerBlk = blockDim.x;
    int MYbid = blockIdx.x;
    int MYtid = threadIdx.x;
    int MYgtid = ThrPerBlk * MYbid + MYtid;
    if (MYgtid >= worklen * numdms)
        return;

    float val = 0;
    int dmi = (int)(MYgtid / worklen); // dms
    int xi = MYgtid - dmi * worklen;   // worklen
    for (i = 0; i < nsub; i++)
    {
        int offset = offsets[dmi * nsub + i];
        int jj = i + (xi + offset) * nsub;
        if (xi < (worklen - offset))
        {
            val += lastdsdata_gpu[jj];
        }
        else
        {
            val += currentdsdata_gpu[i + (xi - (worklen - offset)) * nsub];
        }
    }
    outdata_gpu[MYgtid] = val;
}

void DownsampDataArray(float *currentdsdata_gpu, float *currentdata_gpu, int worklen, int bs, int nsub)
{
    int BlkPerRow = (worklen * nsub - 1 + 1024) / 1024;
    DownsampDataArray_kernel<<<BlkPerRow, 1024>>>(currentdsdata_gpu, currentdata_gpu, nsub, worklen, bs);
}

void getDedisp(float *currentdsdata_gpu, float *lastdsdata_gpu, float *outdata_gpu, int worklen, int nsub, int numdms, int *offsets)
{
    int Blocksize = 1024;
    int BlkPerRow = (worklen * numdms - 1 + Blocksize) / Blocksize;
    getDedispGPU<<<BlkPerRow, Blocksize>>>(currentdsdata_gpu, lastdsdata_gpu, outdata_gpu, worklen, nsub, numdms, offsets);
}