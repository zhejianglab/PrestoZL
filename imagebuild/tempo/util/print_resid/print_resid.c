/* print_resid.c 
 *
 * Read a resid2.tmp file and do stuff with it
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include "tempo_output.h"

#define MAX_OUTS 10

void usage() {
    printf(
            "Usage: print_resid [output_options] resid_file\n"
            "Options determine what is printed, in what order (max %d):\n"
            "  -m, --mjd        MJD of barycentric TOA\n"
            "  -p, --res_phase  Residual in phase (turns)\n"
            "  -t, --res_sec    Residual in time (sec)\n"
            "  -r, --res_us     Residual in time (us)\n"
            "  -o, --ophase     Orbital phase (turns)\n"
            "  -f, --rf         Barycentric RF (MHz)\n"
            "  -w, --weight     Weight of point in fit\n"
            "  -e, --err        TOA uncertainty (us)\n"
            "  -i, --prefit_sec Pre-fit residual (sec)\n"
            "  -d, --aux        Delta DM (non-GLS) or red residual (s, GLS)\n"
            "  -I, --info       Labels from info.tmp\n"
            "Other options:\n"
            "  -W, --white      Print whitened residuals (GLS only)\n"
            "  -s, --stats      Print stats at beginning\n"
            "  -b, --bands      Print blank lines between bands\n"
            "  -z, --zap        Don't print zero-weighted points\n"
            "  -h, --help       Print this message\n"
            "Calling with no args is equivalent to:\n"
            "  print_resid -mfreoI resid2.tmp\n", 
            MAX_OUTS
          );
}

/* Calc weighted or not rms */
double resid_rms(int npts, struct residual *r, int weight) {
    int i;
    double rsum=0.0, wsum=0.0, w;
    for (i=0; i<npts; i++) {
        w = weight ? r[i].weight : 1.0; 
        rsum += w*r[i].res_sec*r[i].res_sec;
        wsum += w;
    }
    rsum /= wsum;
    rsum = sqrt(rsum);
    rsum *= 1e6;
    return(rsum);
}

/* Daily avg rms - residuals should be already ~sorted by date */
double resid_rms_daily(int npts, struct residual *r, int weight) {
    int i;
    double rsum=0.0, wsum=0.0, drsum=0.0, dwsum=0.0, w;
    w = weight ? r[0].weight : 1.0; 
    drsum = w*r[0].res_sec;
    dwsum = w;
    for (i=1; i<npts; i++) {
        if (fabs(r[i].toa-r[i-1].toa)>1.0) { 
            drsum /= dwsum;
            rsum += dwsum*drsum*drsum;
            wsum += dwsum;
            drsum=0.0;
            dwsum=0.0;
        } 
        w = weight ? r[i].weight : 1.0; 
        drsum += w*r[i].res_sec;
        dwsum += w;
    }
    drsum /= dwsum;
    rsum += dwsum*drsum*drsum;
    wsum += dwsum;
    rsum /= wsum;
    rsum = sqrt(rsum);
    rsum *= 1e6;
    return(rsum);
}

int main(int argc, char *argv[]) {
    int rv;

    /* Parse opts */
    static struct option long_opts[] = {
        {"help",       0, NULL, 'h'},
        {"stats",      0, NULL, 's'},
        {"bands",      0, NULL, 'b'},
        {"zap",        0, NULL, 'z'},
        {"white",      0, NULL, 'W'},
        {"mjd",        0, NULL, 'm'},
        {"res_phase",  0, NULL, 'p'},
        {"res_sec",    0, NULL, 't'},
        {"res_us",     0, NULL, 'r'},
        {"ophase",     0, NULL, 'o'},
        {"rf",         0, NULL, 'f'},
        {"weight",     0, NULL, 'w'},
        {"err",        0, NULL, 'e'},
        {"prefit_sec", 0, NULL, 'i'},
        {"aux",        0, NULL, 'd'},
        {"info",       0, NULL, 'I'},
        {0,0,0,0}
    };
    int opt, opti;
    char outputs[MAX_OUTS];
    int nout=0;
    outputs[0]='\0';
    int do_stat=0, do_band=0, zap_zero_wt=0, do_info=0, do_white=0;
    while ((opt=getopt_long(argc,argv,"hsbzWmptrofweidI",long_opts,&opti))!=-1) {
        switch (opt) {
            case 'm':
            case 'p':
            case 't':
            case 'r':
            case 'o':
            case 'f':
            case 'w':
            case 'e':
            case 'i':
            case 'd':
            case 'I':
                outputs[nout] = opt;
                nout++;
                if (nout==MAX_OUTS) { 
                    fprintf(stderr, "Too many options.\n");
                    usage();
                    exit(1);
                }
                outputs[nout] = '\0';
                if (opt=='I') do_info=1;
                break;
            case 's':
                do_stat=1;
                break;
            case 'b':
                do_band=1;
                break;
            case 'z':
                zap_zero_wt=1;
                break;
            case 'W':
                do_white=1;
                break;
            case 'h':
            default:
                usage();
                exit(0);
                break;

        }
    }

    /* Fill default options if none given */
    if (nout==0) { 
      /* Check if info.tmp exists */
      FILE *ftmp = fopen("info.tmp","r");
      if (ftmp) {
        sprintf(outputs, "mfreoI"); 
        nout=6; 
        do_info=1; 
        fclose(ftmp);
      } else {
        sprintf(outputs, "mfreo"); 
        nout=5; 
      }
    }

    /* Use default filename if none given */
    char fname[256];
    if (optind>=argc) {
        sprintf(fname, "resid2.tmp");
    } else {
        sprintf(fname, "%s", argv[optind]);
    }

    FILE *rf=NULL;
    struct residual *r=NULL;
    int npts;
    rf = fopen(fname, "r");
    if (rf==NULL) {
        fprintf(stderr, "Error opening %s\n", fname);
        if (optind>=argc) { usage(); }
        exit(1);
    }
    rv = read_resid(rf, &r, &npts);
    if (rv<0) { 
        fprintf(stderr, "Error parsing %s\n", fname);
        fclose(rf);
        exit(0);
    }
    fclose(rf);

    int i,j;
    char **info_lines = NULL;
    const int info_len = 80;
    if (do_info) {
      info_lines = (char**)malloc(sizeof(char*)*npts);
      for (i=0; i<npts; i++) 
        info_lines[i] = (char*)malloc(sizeof(char)*info_len);
      FILE *info_file = fopen("info.tmp","r");
      if (info_file==NULL) {
        fprintf(stderr, "Error opening info.tmp\n");
        exit(1);
      }
      for (i=0; i<npts; i++) {
        rv = fscanf(info_file, "%s\n", info_lines[i]);
        if (rv!=1) {
          fprintf(stderr, "Error reading from info.tmp (i=%d npts=%d rv=%d)\n",
              i, npts, rv);
          exit(1);
        }
      }
      fclose(info_file);
    }

    /* Whiten residuals if needed.  This only works if the fit was done
     * in GLS mode.  There is no way to check this from the resid2.tmp
     * file so we need to trust that people know what they are doing :) 
     */
    if (do_white) {
        for (i=0; i<npts; i++) {
          double p = r[i].res_sec / r[i].res_phase;
          r[i].res_sec -= r[i].aux;
          r[i].res_phase -= r[i].aux/p;
        }
    }


    /* Print some summary statistics at top */
    if (do_stat) {

        /* Basic info */
        printf("# File: %s\n", fname);
        printf("# Ntoa: %d\n", npts);

        /* Compute anything else and put it here */
        printf("#  RMS(raw): %.4f\n", resid_rms(npts, r, 0));
        printf("# WRMS(raw): %.4f\n", resid_rms(npts, r, 1));
        printf("#  RMS(day): %.4f\n", resid_rms_daily(npts, r, 0));
        printf("# WRMS(day): %.4f\n", resid_rms_daily(npts, r, 1));

        /* Last line, describes outputs */
        printf("# Data: ");
        for (i=0; i<nout; i++) {
            switch (outputs[i]) {
                case 'm':
                    printf("MJD");
                    break;
                case 'p':
                    printf("Resid(turns)");
                    break;
                case 't':
                    printf("Resid(s)");
                    break;
                case 'r':
                    printf("Resid(us)");
                    break;
                case 'o':
                    printf("Ophase(turns)");
                    break;
                case 'f':
                    printf("RF(MHz)");
                    break;
                case 'w':
                    printf("Weight");
                    break;
                case 'e':
                    printf("Err(us)");
                    break;
                case 'i':
                    printf("Prefit(s)");
                    break;
                case 'd':
                    printf("DDM");
                    break;
                case 'I':
                    printf("Info");
                    break;
            }
            if (i==nout-1) { printf("\n"); } else { printf(" "); }
        }
    }

    /* Print outputs */
    float last_rf=0.0;
    for (i=0; i<npts; i++) {
        if ((fabs(r[i].rf_bary-last_rf)>50.0) && (i>0) && do_band) 
            printf("\n\n"); 
        if (zap_zero_wt && r[i].weight==0.0) continue; 
        for (j=0; j<nout; j++) {
            switch (outputs[j]) {
                case 'm':
                    printf("%15.9f", r[i].toa);
                    break;
                case 'p':
                    printf("%+.8e", r[i].res_phase);
                    break;
                case 't':
                    printf("%+.8e", r[i].res_sec);
                    break;
                case 'r':
                    printf("%+.8e", r[i].res_sec*1e6);
                    break;
                case 'o':
                    printf("%.8f", r[i].ophase);
                    break;
                case 'f':
                    printf("%9.4f", r[i].rf_bary);
                    break;
                case 'w':
                    printf("%.4e", r[i].weight);
                    break;
                case 'e':
                    printf("%6.3e", r[i].err_us);
                    break;
                case 'i':
                    printf("%+.8e", r[i].prefit_sec);
                    break;
                case 'd':
                    printf("%+.8e", r[i].aux);
                    break;
                case 'I':
                    printf("%s", info_lines[i]);
                    break;
            }
            if (j==nout-1) { printf("\n"); }
            else { printf(" "); }
        }
        last_rf = r[i].rf_bary;
    }

    exit(0);
}

