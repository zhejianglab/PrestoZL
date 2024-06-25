/* tempo_output.c
 *
 * Functions for dealing with tempo's binary output files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "tempo_output.h"

/* Read a fortran style "record" from file f.
 * First 4 bytes are a long int giving record length (N), 
 * followed by N bytes of data, followed by the same 4 byte
 * long N.
 *
 * TODO: Some versions of fortran use 8-byte record markers.
 *       Would be good to determine this automatically...
 *
 * Errors if first/last N's are not equal, or if N
 * is greater than max
 */
int read_fortran_record(FILE *f, char *buf, int max) {
    uint32_t n, m;
    int rv;

    rv = fread(&n, sizeof(uint32_t), 1, f);
    if (rv!=1) { 
        if (feof(f)) { return(0); }
        fprintf(stderr, "read_fortran_record: Error reading file\n");
        return(rv); 
    }
    if (n>max) { 
        fprintf(stderr, "read_fortran_record: max bytes exceeded.\n");
        return(-1);
    }
    rv = fread(buf, 1, n, f);
    if (rv!=n) { 
        fprintf(stderr, "read_fortran_record: Expected %d bytes, read %d\n",
                n, rv);
        return(rv); }
    rv = fread(&m, sizeof(uint32_t), 1, f);
    if (m!=n) { 
        fprintf(stderr, "read_fortran_record: invalid record format\n");
        return(-2);
    }
    return(n);
}

/* Read a open resid2.tmp file, alloc memory, return results
 * via residual structure.  So far no checking of contents, 
 * it is just mapped direct to the residual struct.
 */
int read_resid(FILE *f, struct residual **r, int *n) {
    int rv;
    char buf[80];
    *n=0;
    if (sizeof(struct residual)!=72) { 
        fprintf(stderr, "read_resid:  sizeof(struct residual)=%ld\n",
                sizeof(struct residual));
        return(-1);
    }
    rv = read_fortran_record(f, buf, 72);
    while (rv==72) {
        *r = (struct residual *)realloc(*r, sizeof(struct residual) * (*n+1));
        memcpy(&((*r)[*n]), buf, 72);
        (*n)++;
        rv = read_fortran_record(f, buf, 72);
    }
    if (rv==0) { return(0); } else { return(-1); }
}
