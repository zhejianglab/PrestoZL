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
