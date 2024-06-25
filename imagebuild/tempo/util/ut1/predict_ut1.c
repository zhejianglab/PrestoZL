/* predict_ut1.c
 * P. Demorest, 2010/04
 * Extend ut1 time series using linear prediction.
 * Note, this requires GSL to be available.
 * Build: gcc -o predict_ut1 predict_ut1.c `gsl-config --libs`
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#define UT1_TAI "ut1_tai.small"

int main(int argc, char *argv[]) {

    // Read input file
    FILE *f = fopen(UT1_TAI, "r");
    if (f==NULL) { 
        fprintf(stderr, "Error opening %s\n", UT1_TAI);
        exit(1);
    }
    char line[256];
    double *mjd0, *tdiff0, *terr0;
    int npts=1; 
    int i=0, j;
    mjd0 = (double *)malloc(sizeof(double) * npts);
    tdiff0 = (double *)malloc(sizeof(double) * npts);
    terr0 = (double *)malloc(sizeof(double) * npts);
    while (fgets(line, 256, f)!=NULL) {
        if (line[0]=='#') continue;
        if (sscanf(line, "%lf %lf %lf", 
                    &mjd0[i], &tdiff0[i], &terr0[i])!=3) continue;
        npts++;
        i++;
        mjd0 = (double *)realloc(mjd0, sizeof(double) * npts);
        tdiff0 = (double *)realloc(tdiff0, sizeof(double) * npts);
        terr0 = (double *)realloc(terr0, sizeof(double) * npts);
    }
    npts--; // Final read always failed
    //printf("Read %d pts\n", npts);
    fclose(f);

    // Subtract off overall linear trend
    double td_mean=0.0, td_slope=0.0;
    double cov00, cov01, cov11, sumsq;
    gsl_fit_linear(mjd0, 1, tdiff0, 1, npts, &td_mean, &td_slope, 
            &cov00, &cov01, &cov11, &sumsq);

    // Convert to GSL vectors, subtracting off line
    gsl_vector *mjd = gsl_vector_alloc(npts);
    gsl_vector *tdiff = gsl_vector_alloc(npts);
    gsl_vector *terr = gsl_vector_alloc(npts);
    for (i=0; i<npts; i++) {
        gsl_vector_set(mjd, i, mjd0[i]);
        double tdiff_linear = td_mean + td_slope*mjd0[i];
        gsl_vector_set(tdiff, i, tdiff0[i] - tdiff_linear);
        gsl_vector_set(terr, i, terr0[i]);
        //printf("%+.10f %+e %+e\n", mjd0[i], tdiff0[i], tdiff_linear);
    }
    free(mjd0);
    free(tdiff0);
    free(terr0);

    // Compute autocorr func
    // TODO: use errors as weights?
    int nlags = 512;
    gsl_vector *acf = gsl_vector_alloc(nlags);
    int il;
    for (il=0; il<nlags; il++) {
        gsl_vector_set(acf, il, 0.0);
        for (i=0; i<npts-il; i++) {
            double t0 = gsl_vector_get(tdiff, i);
            double t1 = gsl_vector_get(tdiff, i+il);
            double tmp = gsl_vector_get(acf, il);
            tmp += t0*t1;
            gsl_vector_set(acf, il, tmp);
        }
        double tmp = gsl_vector_get(acf, il);
        tmp /= (double)(npts);
        gsl_vector_set(acf, il, tmp);
    }

    // Form corr matrix, simple LU-solve
    int np = nlags - 1;
    gsl_matrix *R = gsl_matrix_alloc(np, np);
    gsl_vector *acf1 = gsl_vector_alloc(np);
    for (i=0; i<np; i++) {
        gsl_vector_set(acf1, i, gsl_vector_get(acf, i+1));
        for (j=0; j<np; j++) {
            gsl_matrix_set(R, i, j, gsl_vector_get(acf, abs(i-j)));
        }
    }
    gsl_vector *a = gsl_vector_alloc(np);
    gsl_permutation *p = gsl_permutation_alloc(np);
    int ss;
    gsl_linalg_LU_decomp(R, p, &ss);
    gsl_linalg_LU_solve(R, p, acf1, a);

    //for (i=0; i<np; i++) {
    //    printf("# %d %+e\n", i+1, gsl_vector_get(a,i));
    //}

    // Predict ahead by some number of points
    int next = 50;
    gsl_vector *mjd_ext = gsl_vector_alloc(npts + next);
    gsl_vector *tdiff_ext = gsl_vector_alloc(npts + next);
    gsl_vector *terr_ext = gsl_vector_alloc(npts + next);
    for (i=0; i<npts; i++) {
        gsl_vector_set(mjd_ext, i, gsl_vector_get(mjd, i));
        gsl_vector_set(tdiff_ext, i, gsl_vector_get(tdiff, i));
        gsl_vector_set(terr_ext, i, gsl_vector_get(terr, i));
    }
    for (i=npts; i<npts+next; i++) {
        gsl_vector_set(mjd_ext, i, gsl_vector_get(mjd_ext,i-1) + 1.0);
        double td_ext = 0.0;
        for (j=0; j<np; j++) {
            double sign = 1.0;
            //if (j>0) { sign = -1.0; }
            //else { sign = 1.0; }
            td_ext += sign*gsl_vector_get(a,j)*gsl_vector_get(tdiff_ext,i-j-1);
        }
        gsl_vector_set(tdiff_ext, i, 2.0*gsl_vector_get(tdiff_ext,i-1)-td_ext);
        //gsl_vector_set(tdiff_ext, i, td_ext);
        gsl_vector_set(terr_ext, i, 0.0);
    }

    for (i=npts; i<npts+next; i++) {
        double tdiff_linear = td_mean + td_slope*gsl_vector_get(mjd_ext, i);
        //double tdiff_linear=0.0;
        printf("%.2f %+.4f %.4f\n", 
                gsl_vector_get(mjd_ext, i),
                gsl_vector_get(tdiff_ext, i) + tdiff_linear,
                gsl_vector_get(terr_ext, i)
              );
    }

}
