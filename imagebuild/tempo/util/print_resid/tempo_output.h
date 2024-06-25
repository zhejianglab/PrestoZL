#ifndef _TEMPO_OUTPUT_H
#define _TEMPO_OUTPUT_H
#include <stdio.h>
struct residual {
    double toa;
    double res_phase;
    double res_sec;
    double ophase;
    double rf_bary;
    double weight;
    double err_us;
    double prefit_sec;
    double aux; // This is red resids (s) in GLS mode, ddm in standard mode
};
int read_fortran_record(FILE *f, char *buf, int max);
int read_resid(FILE *f, struct residual **r, int *n);
#endif
