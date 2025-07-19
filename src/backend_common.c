#include <sys/types.h>
#include <pwd.h>
#include <ctype.h>
#include <sys/time.h>
#include "backend_common.h"
#include "misc_utils.h"
#include "fftw3.h"

static long long currentspectra = 0;
static int using_MPI = 0;

#define SWAP(a,b) tmpswap=(a);(a)=(b);(b)=tmpswap;

extern int clip_times(float *rawdata, int ptsperblk, int numchan,
                      float clip_sigma, float *good_chan_levels);
extern void float_dedisp(float *data, float *lastdata,
                         int numpts, int numchan,
                         int *delays, float approx_mean, float *result);
extern void dedisp_subbands(float *data, float *lastdata,
                            int numpts, int numchan,
                            int *delays, int numsubbands, float *result);
extern void dedisp_subbands_cache(unsigned char *data, float *data_scl, float *data_offs, unsigned char *lastdata, float *lastdata_scl, float *lastdata_offs,
                     int numpts, int numchan,
                     int *delays, int numsubbands, float *result);
extern short transpose_float(float *a, int nx, int ny, unsigned char *move,
                             int move_size);
extern double DATEOBS_to_MJD(char *dateobs, int *mjd_day, double *mjd_fracday);
extern void read_filterbank_files(struct spectra_info *s);
extern void read_PSRFITS_files(struct spectra_info *s);
extern void init_static_values(struct spectra_info *s);
extern fftwf_plan plan_transpose(int rows, int cols, float *in, float *out);
extern int *ranges_to_ivect(char *str, int minval, int maxval, int *numvals);

void psrdatatype_description(char *outstr, psrdatatype ptype)
{
    if (ptype == SIGPROCFB)
        strcpy(outstr, "SIGPROC filterbank");
    else if (ptype == PSRFITS)
        strcpy(outstr, "PSRFITS");
    else if (ptype == SCAMP)
        strcpy(outstr, "SCAMP 1-bit filterbank");
    else if (ptype == BPP)
        strcpy(outstr, "GBT BCPM");
    else if (ptype == WAPP)
        strcpy(outstr, "GBT/Caltech Spigot");
    else if (ptype == SPIGOT)
        strcpy(outstr, "WAPP");
    else if (ptype == SUBBAND)
        strcpy(outstr, "PRESTO subband");
    else if (ptype == DAT)
        strcpy(outstr, "PRESTO time series of floats");
    else if (ptype == SDAT)
        strcpy(outstr, "PRESTO time series of shorts");
    else if (ptype == EVENTS)
        strcpy(outstr, "Event list");
    else
        strcpy(outstr, "Unknown");
    return;
}

void set_using_MPI(void)
{
    using_MPI = 1;
}

void close_rawfiles(struct spectra_info *s)
{
    int ii;
    if (s->datatype == PSRFITS) {
        int status = 0;
        for (ii = 0; ii < s->num_files; ii++)
            fits_close_file(s->fitsfiles[ii], &status);
        free(s->fitsfiles);
    } else {
        for (ii = 0; ii < s->num_files; ii++)
            fclose(s->files[ii]);
        free(s->files);
    }
}

void read_rawdata_files(struct spectra_info *s)
{
    if (s->datatype == SIGPROCFB)
        read_filterbank_files(s);
    else if (s->datatype == PSRFITS)
        read_PSRFITS_files(s);
    else if (s->datatype == SCAMP)
        exit(1);
    else if (s->datatype == BPP)
        exit(1);
    else if (s->datatype == WAPP)
        exit(1);
    else if (s->datatype == SPIGOT)
        exit(1);
    return;
}

void add_padding(float *fdata, float *padding, int numchan, int numtopad)
{
    int ii;
    for (ii = 0; ii < numtopad; ii++)
        memcpy(fdata + ii * numchan, padding, numchan * sizeof(float));
}


void identify_psrdatatype(struct spectra_info *s, int output)
{
    char *root, *suffix, ctmp[40];

    /* Split the filename into a rootname and a suffix */
    if (split_root_suffix(s->filenames[0], &root, &suffix) == 0) {
        fprintf(stderr,
                "Error!:  The input filename (%s) must have a suffix!",
                s->filenames[0]);
        exit(-1);
    } else {
        if (strcmp(suffix, "dat") == 0)
            s->datatype = DAT;
        else if (strcmp(suffix, "sdat") == 0)
            s->datatype = SDAT;
        else if (strncmp(suffix, "sub0", 4) == 0)
            s->datatype = SUBBAND;
        else if (strcmp(suffix, "events") == 0)
            s->datatype = EVENTS;
        else if (strcmp(suffix, "bcpm1") == 0 || strcmp(suffix, "bcpm2") == 0)
            s->datatype = BPP;
        else if (strcmp(suffix, "fil") == 0 || strcmp(suffix, "fb") == 0)
            s->datatype = SIGPROCFB;
        else if ((strcmp(suffix, "fits") == 0) || (strcmp(suffix, "sf") == 0)) {
            if (strstr(root, "spigot_5") != NULL)
                s->datatype = SPIGOT;
//            else if (is_PSRFITS(s->filenames[0])) s->datatype = PSRFITS;
            else
                s->datatype = PSRFITS;
        } else if (strcmp(suffix, "pkmb") == 0)
            s->datatype = SCAMP;
        else if (isdigit(suffix[0]) && isdigit(suffix[1]) && isdigit(suffix[2]))
            s->datatype = WAPP;
        else
            s->datatype = UNSET;
    }
    psrdatatype_description(ctmp, s->datatype);
    // if (output)
    //     printf("Assuming the data are %s format...\n", ctmp);
    free(root);
    free(suffix);
}


void spectra_info_set_defaults(struct spectra_info *s)
{
    strcpy(s->telescope, "unset");
    strcpy(s->observer, "unset");
    strcpy(s->source, "unset");
    strcpy(s->frontend, "unset");
    strcpy(s->backend, "unset");
    strcpy(s->project_id, "unset");
    strcpy(s->date_obs, "unset");
    strcpy(s->ra_str, "unset");
    strcpy(s->dec_str, "unset");
    strcpy(s->poln_type, "unset");
    strcpy(s->poln_order, "unset");
    s->datatype = UNSET;
    s->N = 0;
    s->T = 0.0;
    s->dt = 0.0;
    s->fctr = 0.0;
    s->lo_freq = 0.0;
    s->hi_freq = 0.0;
    s->orig_df = 0.0;
    s->chan_dm = 0.0;
    s->df = 0.0;
    s->BW = 0.0;
    s->ra2000 = 0.0;
    s->dec2000 = 0.0;
    s->azimuth = 0.0;
    s->zenith_ang = 0.0;
    s->beam_FWHM = 0.0;
    s->time_per_subint = 0.0;
    s->scan_number = 0;
    s->tracking = 1;
    s->orig_num_chan = 0;
    s->num_channels = 0;
    s->num_polns = 0;
    s->num_beams = 1;
    s->beamnum = 0;
    s->summed_polns = 1;
    s->FITS_typecode = 0;
    s->bits_per_sample = 0;
    s->bytes_per_spectra = 0;
    s->samples_per_spectra = 0;
    s->bytes_per_subint = 0;
    s->spectra_per_subint = 0;
    s->samples_per_subint = 0;
    s->min_spect_per_read = 0;
    s->num_files = 0;
    s->offs_sub_col = 0;
    s->dat_wts_col = 0;
    s->dat_offs_col = 0;
    s->dat_scl_col = 0;
    s->data_col = 0;
    s->apply_scale = 0;
    s->apply_offset = 0;
    s->apply_weight = 0;
    s->apply_flipband = 0;
    s->remove_zerodm = 0;
    s->use_poln = 0;
    s->flip_bytes = 0;
    s->num_ignorechans = 0;
    s->zero_offset = 0.0;
    s->clip_sigma = 0.0;
    s->start_MJD = NULL;
    s->files = NULL;
    s->fitsfiles = NULL;
    s->padvals = NULL;
    s->header_offset = NULL;
    s->start_subint = NULL;
    s->num_subint = NULL;
    s->ignorechans = NULL;
    s->ignorechans_str = NULL;
    s->start_spec = NULL;
    s->num_spec = NULL;
    s->num_pad = NULL;
};


void print_spectra_info(struct spectra_info *s)
// Output a spectra_info structure in human readable form
{
    char ctmp[40];

    psrdatatype_description(ctmp, s->datatype);
    printf("From the %s file '%s':\n", ctmp, s->filenames[0]);
    if (strcmp(s->telescope, "unset") != 0)
        printf("                  Telescope = %s\n", s->telescope);
    if (strcmp(s->observer, "unset") != 0)
        printf("                   Observer = %s\n", s->observer);
    if (strcmp(s->source, "unset") != 0)
        printf("                Source Name = %s\n", s->source);
    if (strcmp(s->frontend, "unset") != 0)
        printf("                   Frontend = %s\n", s->frontend);
    if (strcmp(s->backend, "unset") != 0)
        printf("                    Backend = %s\n", s->backend);
    if (strcmp(s->project_id, "unset") != 0)
        printf("                 Project ID = %s\n", s->project_id);
    if (strcmp(s->date_obs, "unset") != 0)
        printf("            Obs Date String = %s\n", s->date_obs);
    if (s->datatype == PSRFITS) {
        int itmp;
        double dtmp;
        DATEOBS_to_MJD(s->date_obs, &itmp, &dtmp);
        sprintf(ctmp, "%.14f", dtmp);
        printf("  MJD start time (DATE-OBS) = %5i.%14s\n", itmp, ctmp + 2);
        printf("     MJD start time (STT_*) = %19.14Lf\n", s->start_MJD[0]);
    } else {
        printf("             MJD start time = %19.14Lf\n", s->start_MJD[0]);
    }
    printf("                   RA J2000 = %s\n", s->ra_str);
    printf("             RA J2000 (deg) = %-17.15g\n", s->ra2000);
    printf("                  Dec J2000 = %s\n", s->dec_str);
    printf("            Dec J2000 (deg) = %-17.15g\n", s->dec2000);
    printf("                  Tracking? = %s\n", s->tracking ? "True" : "False");
    printf("              Azimuth (deg) = %-.7g\n", s->azimuth);
    printf("           Zenith Ang (deg) = %-.7g\n", s->zenith_ang);
    if (strcmp(s->poln_type, "unset") != 0)
        printf("          Polarization type = %s\n", s->poln_type);
    if (s->num_polns >= 2 && !s->summed_polns)
        printf("            Number of polns = %d\n", s->num_polns);
    else if (s->summed_polns)
        printf("            Number of polns = 2 (summed)\n");
    else
        printf("            Number of polns = 1\n");
    if (strcmp(s->poln_order, "unset") != 0)
        printf("         Polarization order = %s\n", s->poln_order);
    printf("           Sample time (us) = %-17.15g\n", s->dt * 1e6);
    printf("         Central freq (MHz) = %-17.15g\n", s->fctr);
    printf("          Low channel (MHz) = %-17.15g\n", s->lo_freq);
    printf("         High channel (MHz) = %-17.15g\n", s->hi_freq);
    printf("        Channel width (MHz) = %-17.15g\n", s->df);
    printf("         Number of channels = %d\n", s->num_channels);
    if (s->chan_dm != 0.0) {
        printf("   Orig Channel width (MHz) = %-17.15g\n", s->orig_df);
        printf("    Orig Number of channels = %d\n", s->orig_num_chan);
        printf("    DM used for chan dedisp = %-17.15g\n", s->chan_dm);
    }
    printf("      Total Bandwidth (MHz) = %-17.15g\n", s->BW);
    if (s->num_beams > 0)
        printf("                       Beam = %d of %d\n", s->beamnum, s->num_beams);
    printf("            Beam FWHM (deg) = %.3f\n", s->beam_FWHM);
    printf("         Spectra per subint = %d\n", s->spectra_per_subint);
    if (s->start_subint != NULL && s->datatype != SIGPROCFB)
        printf("            Starting subint = %d\n", s->start_subint[0]);
    if (s->num_subint != NULL && s->datatype != SIGPROCFB)
        printf("           Subints per file = %d\n", s->num_subint[0]);
    printf("           Spectra per file = %lld\n", s->num_spec[0]);
    printf("      Time per subint (sec) = %-.12g\n", s->time_per_subint);
    printf("        Time per file (sec) = %-.12g\n", s->num_spec[0] * s->dt);
    printf("            bits per sample = %d\n", s->bits_per_sample);
    if (s->bits_per_sample==8)
        printf("          Are bytes signed? = %s\n",
               s->signedints ? "True" : "False");
    {
        int fact = 1;
        if ((s->datatype == PSRFITS) && (s->bits_per_sample < 8))
            fact = 8 / s->bits_per_sample;
        printf("          bytes per spectra = %d\n", s->bytes_per_spectra / fact);
        printf("        samples per spectra = %d\n", s->samples_per_spectra);
        printf("           bytes per subint = %d\n", s->bytes_per_subint / fact);
        printf("         samples per subint = %d\n", s->samples_per_subint);
    }
    printf("                zero offset = %-17.15g\n", s->zero_offset);
    printf("           Invert the band? = %s\n",
           (s->apply_flipband > 0) ? "True" : "False");
    if (s->header_offset != NULL)
        printf("       bytes in file header = %d\n", s->header_offset[0]);
    if (s->datatype == PSRFITS) {
        int ii, numhdus, hdutype, status = 0;
        char comment[120];
        printf("  PSRFITS Specific info:\n");
        fits_get_num_hdus(s->fitsfiles[0], &numhdus, &status);
        printf("                       HDUs = primary, ");
        for (ii = 2; ii < numhdus + 1; ii++) {
            fits_movabs_hdu(s->fitsfiles[0], ii, &hdutype, &status);
            fits_read_key(s->fitsfiles[0], TSTRING, "EXTNAME", ctmp, comment,
                          &status);
            printf("%s%s", ctmp, (ii < numhdus) ? ", " : "\n");
        }
        printf("              FITS typecode = %d\n", s->FITS_typecode);
        printf("                DATA column = %d\n", s->data_col);
        printf("             Apply scaling? = %s\n",
               s->apply_scale ? "True" : "False");
        printf("             Apply offsets? = %s\n",
               s->apply_offset ? "True" : "False");
        printf("             Apply weights? = %s\n",
               s->apply_weight ? "True" : "False");
    }
}


void print_spectra_info_summary(struct spectra_info *s)
// Print the basic details of the files that are being processed
{
    int ii, nn;
    printf("    Number of files = %d\n", s->num_files);
    if (s->num_polns >= 2 && !s->summed_polns)
        printf("       Num of polns = %d\n", s->num_polns);
    else if (s->summed_polns)
        printf("       Num of polns = 2 (summed)\n");
    else
        printf("       Num of polns = 1\n");
    printf("  Center freq (MHz) = %.8g\n", s->fctr);
    printf("    Num of channels = %d\n", s->num_channels);
    printf("    Sample time (s) = %-14.14g\n", s->dt);
    printf("     Spectra/subint = %d\n", s->spectra_per_subint);
    printf("   Total points (N) = %lld\n", s->N);
    printf("     Total time (s) = %-14.14g\n", s->T);
    printf("     Clipping sigma = %.3f\n", s->clip_sigma);
    if (s->zero_offset != 0.0)
        printf("        zero offset = %-17.15g\n", s->zero_offset);
    printf("   Invert the band? = %s\n", (s->apply_flipband > 0) ? "True" : "False");
    printf("          Byteswap? = %s\n", s->flip_bytes ? "True" : "False");
    printf("     Remove zeroDM? = %s\n", s->remove_zerodm ? "True" : "False");
    if (s->datatype == PSRFITS) {
        printf("     Apply scaling? = %s\n", s->apply_scale ? "True" : "False");
        printf("     Apply offsets? = %s\n", s->apply_offset ? "True" : "False");
        printf("     Apply weights? = %s\n", s->apply_weight ? "True" : "False");
    }
    if (s->num_ignorechans)
    printf("  Ignoring channels = %s\n", s->ignorechans_str);
    printf("\nFile  Start Spec   Samples     Padding        Start MJD\n");
    printf("----  ----------  ----------  ----------  --------------------\n");
    if (s->datatype == SUBBAND || s->datatype == DAT ||
        s->datatype == EVENTS || s->datatype == SDAT)
        nn = 1;
    else
        nn = s->num_files;
    for (ii = 0; ii < nn; ii++)
        printf("%-4d  %10lld  %10lld  %10lld  %19.14Lf\n", ii + 1,
               s->start_spec[ii], s->num_spec[ii], s->num_pad[ii], s->start_MJD[ii]);
    printf("\n");
}


void spectra_info_to_inf(struct spectra_info *s, infodata * idata)
// Convert a spectra_info structure into an infodata structure
{
    int ii, index = 2;
    char ctmp[100];
    struct passwd *pwd;

    strcpy(idata->object, s->source);
    hours2hms(s->ra2000 / 15.0, &(idata->ra_h), &(idata->ra_m), &(idata->ra_s));
    deg2dms(s->dec2000, &(idata->dec_d), &(idata->dec_m), &(idata->dec_s));
    strcpy(idata->telescope, s->telescope);
    strcpy(idata->instrument, s->backend);
    idata->num_chan = s->num_channels;
    idata->dt = s->dt;
    // DATEOBS_to_MJD(s->date_obs, &(idata->mjd_i), &(idata->mjd_f));
    idata->mjd_i = (int) (s->start_MJD[0]);
    idata->mjd_f = s->start_MJD[0] - idata->mjd_i;
    idata->N = s->N;
    idata->freqband = s->BW;
    idata->chan_wid = s->df;
    idata->freq = s->lo_freq;
    idata->fov = s->beam_FWHM * 3600.0; // in arcsec
    idata->bary = 0;
    idata->numonoff = 0;
    strcpy(idata->band, "Radio");
    pwd = getpwuid(geteuid());
    strcpy(idata->analyzer, pwd->pw_name);
    strcpy(idata->observer, s->observer);
    if (s->summed_polns)
        sprintf(ctmp,
                "2 polns were summed.  Samples have %d bits.", s->bits_per_sample);
    else
        sprintf(ctmp, "%d polns were not summed.  Samples have %d bits.",
                s->num_polns, s->bits_per_sample);
    sprintf(idata->notes, "Project ID %s, Date: %s.\n    %s\n",
            s->project_id, s->date_obs, ctmp);
    // If we have more than one input file, see it we need onoff bins
    if (s->num_files == 1 && s->num_pad[0] == 0) {
        idata->numonoff = 0;
        return;
    }
    /* Determine the topocentric onoff bins */
    idata->numonoff = 1;
    idata->onoff[0] = 0.0;
    idata->onoff[1] = s->num_spec[0] - 1.0;
    for (ii = 1; ii < s->num_files; ii++) {
        if (s->num_pad[ii - 1]) {
            idata->onoff[index] = idata->onoff[index - 1] + s->num_pad[ii - 1];
            idata->onoff[index + 1] = idata->onoff[index] + s->num_spec[ii];
            idata->numonoff++;
            index += 2;
        } else {
            idata->onoff[index - 1] += s->num_spec[ii];
        }
    }
    if (s->num_pad[s->num_files - 1]) {
        idata->onoff[index] = idata->onoff[index - 1] + s->num_pad[s->num_files - 1];
        idata->onoff[index + 1] = idata->onoff[index];
        idata->numonoff++;
    }
}


long long offset_to_spectra(long long specnum, struct spectra_info *s)
// This routine offsets into the raw data files to the spectra
// 'specnum'.  It returns the current spectra number.
{
    long long retval;
    retval = s->offset_to_spectra(specnum, s);
    currentspectra = retval;
    return retval;
}


void set_currentspectra(long long specnum)
// Set the static variable currentspectra to specnum.  This is used
// 'by mpiprepsubband for the non-master CPUs to know what spectra they
// are getting for masking.
{
    currentspectra = specnum;
}


int read_rawblocks(float *fdata, int numsubints, struct spectra_info *s,
                   int *padding)
// This routine reads numsubints rawdata blocks from raw radio pulsar
// data. The floating-point filterbank data is returned in rawdata
// which must have a size of numsubints * s->samples_per_subint.  The
// number of blocks read is returned.  If padding is returned as 1,
// then padding was added and statistics should not be calculated.
{
    long long ii, loffset;
    int retval = 0, gotblock = 0, pad = 0, numpad = 0, numvals;
    static float *rawdata = NULL;
    static int firsttime = 1;

    numvals = s->spectra_per_subint * s->num_channels;
    if (firsttime) {
        // Needs to be twice as large for buffering if adding observations together
        rawdata = gen_fvect(2 * numvals);
        firsttime = 0;
    }
    *padding = 0;
    for (ii = 0; ii < numsubints; ii++) {
        gotblock = s->get_rawblock(rawdata, s, &pad);
        if (gotblock == 0)
            break;
        retval += gotblock;
        loffset = ii * numvals;
        memcpy(fdata + loffset, rawdata, numvals * sizeof(float));
        if (pad)
            numpad++;
    }
    if (gotblock == 0) {        // Now fill the rest of the data with padding
        for (; ii < numsubints; ii++) {
            long long jj, loffset2;
            loffset = ii * numvals;
            for (jj = 0; jj < s->spectra_per_subint; jj++) {
                loffset2 = loffset + jj * s->num_channels;
                memcpy(fdata + loffset2, s->padvals,
                       s->num_channels * sizeof(float));
            }
        }
        numpad++;
    }

    /* Return padding 'true' if any block was padding */
    if (numpad)
        *padding = 1;
    return retval;
}


int read_psrdata(float *fdata, int numspect, struct spectra_info *s,
                 int *delays, int *padding,
                 int *maskchans, int *nummasked, mask * obsmask)
// This routine reads numspect from the raw pulsar data defined in
// "s". Time delays and a mask are applied to each channel.  It
// returns the # of points read if successful, 0 otherwise.  If
// padding is returned as 1, then padding was added and statistics
// should not be calculated.  maskchans is an array of length numchans
// contains a list of the number of channels that were masked.  The #
// of channels masked is returned in nummasked.  obsmask is the mask
// structure to use for masking.
{
    int numread = 0;
    double starttime = 0.0;
    long long ii, jj, templen, loffset;
    static float *tmpswap, *rawdata1, *rawdata2;
    static float *currentdata, *lastdata;
    static int firsttime = 1, numsubints = 1, allocd = 0, mask = 0;
    static double duration = 0.0;

    *nummasked = 0;
    if (firsttime) {
        if (numspect % s->spectra_per_subint) {
            fprintf(stderr,
                    "Error!:  numspect %d must be a multiple of %d in read_psrdata()!\n",
                    numspect, s->spectra_per_subint);
            exit(-1);
        } else
            numsubints = numspect / s->spectra_per_subint;
        if (obsmask->numchan)
            mask = 1;
        // The following can overflow a regular int
        templen = ((long long) numsubints) * s->spectra_per_subint * s->num_channels;
        if (templen > 1000000000L) {
            printf("\nWARNING:  Trying to allocate %.2f GB of RAM in read_psrdata()!!\n",
                   templen*8L/1e9);
            printf("    This will possibly fail.  Is the dispersive delay across the\n");
            printf("    band longer than (or comparable to) the duration of the input file??\n");
        }
        rawdata1 = gen_fvect(templen);
        rawdata2 = gen_fvect(templen);
        allocd = 1;
        duration = numsubints * s->time_per_subint;
        currentdata = rawdata1;
        lastdata = rawdata2;
    }

    /* Read, convert and de-disperse */
    if (allocd) {
        while (1) {
            starttime = currentspectra * s->dt;
            numread = read_rawblocks(currentdata, numsubints, s, padding);
            if (mask)
                *nummasked = check_mask(starttime, duration, obsmask, maskchans);
            currentspectra += numread * s->spectra_per_subint;

            /* Clip nasty RFI if requested and we're not masking all the channels */
            if ((s->clip_sigma > 0.0) && !(mask && (*nummasked == -1)))
                clip_times(currentdata, numspect, s->num_channels, s->clip_sigma,
                           s->padvals);

            if (mask) {
                if (*nummasked == -1) { /* If all channels are masked */
                    for (ii = 0; ii < numspect; ii++) {
                        loffset = ii * s->num_channels;
                        memcpy(currentdata + loffset,
                               s->padvals, s->num_channels * sizeof(float));
                    }
                } else if (*nummasked > 0) {    /* Only some of the channels are masked */
                    int channum;
                    for (ii = 0; ii < numspect; ii++) {
                        loffset = ii * s->num_channels;
                        for (jj = 0; jj < *nummasked; jj++) {
                            channum = maskchans[jj];
                            currentdata[loffset + channum] = s->padvals[channum];
                        }
                    }
                }
            }

            if (s->num_ignorechans) { // These are channels we explicitly zero
                int channum;
                for (ii = 0; ii < numspect; ii++) {
                    loffset = ii * s->num_channels;
                    for (jj = 0; jj < s->num_ignorechans; jj++) {
                        channum = s->ignorechans[jj];
                        currentdata[loffset + channum] = 0.0;
                    }
                }
            }

            if (!firsttime)
                float_dedisp(currentdata, lastdata, numspect, s->num_channels,
                             delays, 0.0, fdata);

            SWAP(currentdata, lastdata);
            if (firsttime)
                firsttime = 0;
            else {
                if (numread != numsubints) {
                    vect_free(rawdata1);
                    vect_free(rawdata2);
                    allocd = 0;
                }
                break;
            }
        }
        return numsubints * s->spectra_per_subint;
    } else {
        return 0;
    }
}


void get_channel(float chandat[], int channum, int numsubints, float rawdata[],
                 struct spectra_info *s)
// Return the values for channel 'channum' in 'chandat' of a block of
// 'numsubints' floating-point spectra data stored in 'rawdata'.
// 'rawdata' should have been initialized and then filled using
// read_rawblocks(), and 'chandat' must have at least 'numsubints' *
// 's->spectra_per_subint' spaces.  Channel 0 is assumed to be the
// lowest freq channel.
{
    long long ii, jj, numspec = numsubints * s->spectra_per_subint;

    if (channum > s->num_channels || channum < 0) {
        fprintf(stderr, "Error!: channum = %d is out of range in get_channel()!\n",
                channum);
        exit(-1);
    }
    /* Check to see if we are explicitly zeroing this channel */
    if (s->num_ignorechans) {
        for (ii = 0; ii < s->num_ignorechans; ii++) {
            if (channum==s->ignorechans[ii]) { // zero it
                for (jj = 0; jj < numspec; jj++)
                    chandat[jj] = 0.0;
                return;
            }
        }
    }
    /* Else select the correct channel */
    for (ii = 0, jj = channum; ii < numspec; ii++, jj += s->num_channels)
        chandat[ii] = rawdata[jj];
}

int prep_subbands(float *fdata, float *rawdata, int *delays, int numsubbands,
                  struct spectra_info *s, int transpose,
                  int *maskchans, int *nummasked, mask * obsmask)
// This routine preps a block of raw spectra for subbanding.  It uses
// dispersion delays in 'delays' to de-disperse the data into
// 'numsubbands' subbands.  It stores the resulting data in vector
// 'fdata' of length 'numsubbands' * 's->spectra_per_subint'.  The low
// freq subband is stored first, then the next highest subband etc,
// with 's->spectra_per_subint' floating points per subband. It
// returns the # of points read if succesful, 0 otherwise.
// 'maskchans' is an array of length numchans which contains a list of
// the number of channels that were masked.  The # of channels masked
// is returned in 'nummasked'.  'obsmask' is the mask structure to use
// for masking.  If 'transpose'==0, the data will be kept in time
// order instead of arranged by subband as above.
{
    int ii, jj, offset;
    double starttime = 0.0;
    static float *tmpswap, *rawdata1, *rawdata2;
    static float *currentdata, *lastdata;
    static int firsttime = 1, mask = 0;
    static fftwf_plan tplan1, tplan2;

    *nummasked = 0;
    if (firsttime) {
        if (obsmask->numchan)
            mask = 1;
        rawdata1 = gen_fvect(s->spectra_per_subint * s->num_channels);
        rawdata2 = gen_fvect(s->spectra_per_subint * s->num_channels);
        currentdata = rawdata1;
        lastdata = rawdata2;
        // Make plans to do fast transposes using FFTW
        tplan1 = plan_transpose(s->spectra_per_subint, s->num_channels,
                                currentdata, currentdata);
        tplan2 = plan_transpose(numsubbands, s->spectra_per_subint, fdata, fdata);
    }

    /* Read and de-disperse */
    memcpy(currentdata, rawdata,
           s->spectra_per_subint * s->num_channels * sizeof(float));
    starttime = currentspectra * s->dt; // or -1 subint?
    if (mask)
        *nummasked = check_mask(starttime, s->time_per_subint, obsmask, maskchans);

    /* Clip nasty RFI if requested and we're not masking all the channels */
    if ((s->clip_sigma > 0.0) && !(mask && (*nummasked == -1)))
        clip_times(currentdata, s->spectra_per_subint, s->num_channels,
                   s->clip_sigma, s->padvals);

    if (mask) {
        if (*nummasked == -1) { /* If all channels are masked */
            for (ii = 0; ii < s->spectra_per_subint; ii++)
                memcpy(currentdata + ii * s->num_channels,
                       s->padvals, s->num_channels * sizeof(float));
        } else if (*nummasked > 0) {    /* Only some of the channels are masked */
            int channum;
            for (ii = 0; ii < s->spectra_per_subint; ii++) {
                offset = ii * s->num_channels;
                for (jj = 0; jj < *nummasked; jj++) {
                    channum = maskchans[jj];
                    currentdata[offset + channum] = s->padvals[channum];
                }
            }
        }
    }

    if (s->num_ignorechans) { // These are channels we explicitly zero
        int channum;
        for (ii = 0; ii < s->spectra_per_subint; ii++) {
            offset = ii * s->num_channels;
            for (jj = 0; jj < s->num_ignorechans; jj++) {
                channum = s->ignorechans[jj];
                currentdata[offset + channum] = 0.0;
            }
        }
    }

    // In mpiprepsubband, the nodes do not call read_subbands() where
    // currentspectra gets incremented.
    if (using_MPI)
        currentspectra += s->spectra_per_subint;

    // Now transpose the raw block of data so that the times in each
    // channel are the most rapidly varying index
    fftwf_execute_r2r(tplan1, currentdata, currentdata);

    if (firsttime) {
        SWAP(currentdata, lastdata);
        firsttime = 0;
        return 0;
    } else {
        dedisp_subbands(currentdata, lastdata, s->spectra_per_subint,
                        s->num_channels, delays, numsubbands, fdata);
        SWAP(currentdata, lastdata);
        // Transpose the resulting data into spectra as a function of time
        if (transpose == 0)
            fftwf_execute_r2r(tplan2, fdata, fdata);
        return s->spectra_per_subint;
    }
}


int read_subbands(float *fdata, int *delays, int numsubbands,
                  struct spectra_info *s, int transpose, int *padding,
                  int *maskchans, int *nummasked, mask * obsmask)
// This routine reads a spectral block/subint from the input raw data
// files. The routine uses dispersion delays in 'delays' to
// de-disperse the data into 'numsubbands' subbands.  It stores the
// resulting data in vector 'fdata' of length 'numsubbands' *
// 's->spectra_per_subint'.  The low freq subband is stored first,
// then the next highest subband etc, with 's->spectra_per_subint'
// floating points per subband.  It returns the # of points read if
// successful, 0 otherwise. If padding is returned as 1, then padding
// was added and statistics should not be calculated.  'maskchans' is
// an array of length numchans which contains a list of the number of
// channels that were masked.  The # of channels masked is returned in
// 'nummasked'.  'obsmask' is the mask structure to use for masking.
// If 'transpose'==0, the data will be kept in time order instead of
// arranged by subband as above.
{
    static int firsttime = 1;
    static float *frawdata;

    if (firsttime) {
        // Check to make sure there isn't more dispersion across a
        // subband than time in a block of data
        if (delays[0] > s->spectra_per_subint) {
            perror("\nError: there is more dispersion across a subband than time\n"
                   "in a block of data.  Increase spectra_per_subint if possible.");
            exit(-1);
        }
        // Needs to be twice as large for buffering if adding observations together
        frawdata = gen_fvect(2 * s->num_channels * s->spectra_per_subint);
        if (!s->get_rawblock(frawdata, s, padding)) {
            perror("Error: problem reading the raw data file in read_subbands()");
            exit(-1);
        }
        if (0 != prep_subbands(fdata, frawdata, delays, numsubbands, s,
                               transpose, maskchans, nummasked, obsmask)) {
            perror("Error: problem initializing prep_subbands() in read_subbands()");
            exit(-1);
        }
        firsttime = 0;
    }
    if (!s->get_rawblock(frawdata, s, padding)) {
        return 0;
    }
    if (prep_subbands(fdata, frawdata, delays, numsubbands, s, transpose,
                      maskchans, nummasked, obsmask) == s->spectra_per_subint) {
        currentspectra += s->spectra_per_subint;
        return s->spectra_per_subint;
    } else {
        return 0;
    }
}

int read_subbands_log(float *fdata, int *delays, int numsubbands,
                  struct spectra_info *s, int transpose, int *padding,
                  int *maskchans, int *nummasked, mask * obsmask, long long *data_size, long *total_microseconds)
// This routine reads a spectral block/subint from the input raw data
// files. The routine uses dispersion delays in 'delays' to
// de-disperse the data into 'numsubbands' subbands.  It stores the
// resulting data in vector 'fdata' of length 'numsubbands' *
// 's->spectra_per_subint'.  The low freq subband is stored first,
// then the next highest subband etc, with 's->spectra_per_subint'
// floating points per subband.  It returns the # of points read if
// successful, 0 otherwise. If padding is returned as 1, then padding
// was added and statistics should not be calculated.  'maskchans' is
// an array of length numchans which contains a list of the number of
// channels that were masked.  The # of channels masked is returned in
// 'nummasked'.  'obsmask' is the mask structure to use for masking.
// If 'transpose'==0, the data will be kept in time order instead of
// arranged by subband as above.
{
    static int firsttime = 1;
    static float *frawdata;

    if (firsttime) {
        // Check to make sure there isn't more dispersion across a
        // subband than time in a block of data
        if (delays[0] > s->spectra_per_subint) {
            perror("\nError: there is more dispersion across a subband than time\n"
                   "in a block of data.  Increase spectra_per_subint if possible.");
            exit(-1);
        }
        // Needs to be twice as large for buffering if adding observations together
        frawdata = gen_fvect(2 * s->num_channels * s->spectra_per_subint);
        if (!s->get_rawblock_log(frawdata, s, padding, data_size, total_microseconds)) {
            perror("Error: problem reading the raw data file in read_subbands()");
            exit(-1);
        }
        if (0 != prep_subbands(fdata, frawdata, delays, numsubbands, s,
                               transpose, maskchans, nummasked, obsmask)) {
            perror("Error: problem initializing prep_subbands() in read_subbands()");
            exit(-1);
        }
        firsttime = 0;
    }
    if (!s->get_rawblock_log(frawdata, s, padding, data_size, total_microseconds)) {
        return 0;
    }
    if (prep_subbands(fdata, frawdata, delays, numsubbands, s, transpose,
                      maskchans, nummasked, obsmask) == s->spectra_per_subint) {
        currentspectra += s->spectra_per_subint;
        return s->spectra_per_subint;
    } else {
        return 0;
    }
}


void flip_band(float *fdata, struct spectra_info *s)
// Flip the bandpass
{
    float ftmp;
    int ii, jj, looffs, hioffs;
    for (ii = 0; ii < s->spectra_per_subint; ii++) {
        looffs = ii * s->num_channels;
        hioffs = looffs + s->num_channels - 1;
        for (jj = 0; jj < s->num_channels / 2; jj++, looffs++, hioffs--) {
            ftmp = fdata[looffs];
            fdata[looffs] = fdata[hioffs];
            fdata[hioffs] = ftmp;
        }
    }
}


int *get_ignorechans(char *ignorechans_str, int minchan, int maxchan,
                     int *num_ignorechans, char **filestr)
// Parse an ignorechans string (or file it was from) and return it as
// a vector of length num_ignorechans.  If the channel string came
// from a file, return the channel string in the file in filestr.
{
    int ii;
    char *parsestr, *sptr;
    FILE *file;
    long long filelen;
    
    if ((file = fopen(ignorechans_str, "r")) != NULL) {
        // If so, see how big it is in bytes
        filelen = chkfilelen(file, 1);
        // if not too big, assume that this is a string to parse
        if (filelen < 100000L) {
            *filestr = (char *) malloc(filelen+1);
            // Now read lines of the file, until the first character
            // is not a comment marker or a newline
            do {
                sptr = fgets(*filestr, filelen, file);
                // Remove newline if needed
                if ((*filestr)[strlen(*filestr)-1] == '\n')
                    (*filestr)[strlen(*filestr)-1] = '\0';
                if (sptr != NULL &&
                    sptr[0] != '\n' &&
                    sptr[0] != '#' &&
                    0 != (ii = strlen(sptr))) { // This is a good line
                    // Copy the line read into parsestr
                    parsestr = (char *) malloc(strlen(*filestr)+1);
                    strcpy(parsestr, *filestr);
                    fclose(file);
                    break;
                } else {
                    if (feof(file)) {
                        sprintf(*filestr,
                                "Error:  end-of-file while looking for range string in get_ignorechans()\n");
                        perror(*filestr);
                        exit(EXIT_FAILURE);
                    }
                }
            } while (1);
            
        } else {
            parsestr = (char *) malloc(1000);
            sprintf(parsestr,
                    "Error:  '%s' is a file, but too big to parse in get_ignorechans()\n",
                    ignorechans_str);
            perror(parsestr);
            exit(EXIT_FAILURE);
        }
    } else {
        // Input string name is not a file, so we will parse it directly
        parsestr = ignorechans_str;
        *filestr = NULL; // Not being used
    }
    return ranges_to_ivect(parsestr, minchan, maxchan, num_ignorechans);
}

void SclData(float *intputDATA, int x, int y, unsigned char *DATA,float *sclArray, float *offsArray, int thread)
{
  int i;

#ifdef _OPENMP
// #pragma omp parallel for default(shared)
#pragma omp parallel for num_threads(thread) shared(intputDATA, DATA, sclArray, offsArray)
#endif
  for(i=0; i<y; i++)
  {
    int j;
    float min, max;
    float value;
    min = max = 0.0f;
    for(j=0; j<x; j++)
    {
      value = intputDATA[i*x+j];
      min = (min > value) ? value : min;
      max = (max < value) ? value : max;
    }

    float range;
    float scale_bk, offset_bk;

    range = max - min;
    if (max == 0.0f && min == 0.0f)
    {
        scale_bk = 1.0;
        offset_bk = 0.0;
    }
    else 
    {
        if (range==0.0f)
            scale_bk = max / 255.0f;
        else
          scale_bk = range / 255.0f;
        offset_bk = min;
    }
    for(j=0; j<x; j++)
    {
        DATA[i*x+j] = (unsigned char) ((intputDATA[i*x+j] - offset_bk) / scale_bk + 0.5f);
    }
    sclArray[i] = scale_bk;
    offsArray[i] = offset_bk;
  }
}

int write_prep_subbands_cache(float *fdata, float *rawdata, int *delays, int numsubbands,
                  struct spectra_info *s, int transpose,
                  int *maskchans, int *nummasked, mask * obsmask, int thread)
// This routine preps a block of raw spectra for subbanding.  It uses
// dispersion delays in 'delays' to de-disperse the data into
// 'numsubbands' subbands.  It stores the resulting data in vector
// 'fdata' of length 'numsubbands' * 's->spectra_per_subint'.  The low
// freq subband is stored first, then the next highest subband etc,
// with 's->spectra_per_subint' floating points per subband. It
// returns the # of points read if succesful, 0 otherwise.
// 'maskchans' is an array of length numchans which contains a list of
// the number of channels that were masked.  The # of channels masked
// is returned in 'nummasked'.  'obsmask' is the mask structure to use
// for masking.  If 'transpose'==0, the data will be kept in time
// order instead of arranged by subband as above.
{
    int ii, jj, offset;
    double starttime = 0.0;
    static float *tmpswap, *rawdata1, *rawdata2;
    static float *currentdata, *lastdata;
    
    static float *currentdata_scl;
    static float *currentdata_offs;
    static unsigned char *currentdata_data;

    static int firsttime = 1, mask = 0;
    static fftwf_plan tplan1;

    *nummasked = 0;
    if (firsttime) {
        if (obsmask->numchan)
            mask = 1;
        rawdata1 = gen_fvect(s->spectra_per_subint * s->num_channels);
        rawdata2 = gen_fvect(s->spectra_per_subint * s->num_channels);
        
        currentdata_scl = malloc(s->num_channels * sizeof(float));
        currentdata_offs = malloc(s->num_channels * sizeof(float));
        currentdata_data = malloc(s->spectra_per_subint * s->num_channels * sizeof(unsigned char));

        currentdata = rawdata1;
        lastdata = rawdata2;
        // Make plans to do fast transposes using FFTW
        tplan1 = plan_transpose(s->spectra_per_subint, s->num_channels,
                                currentdata, currentdata);
        // tplan2 = plan_transpose(numsubbands, s->spectra_per_subint, fdata, fdata);
    }

    /* Read and de-disperse */
    memcpy(currentdata, rawdata,
           s->spectra_per_subint * s->num_channels * sizeof(float));
    starttime = currentspectra * s->dt; // or -1 subint?
    if (mask)
        *nummasked = check_mask(starttime, s->time_per_subint, obsmask, maskchans);

    /* Clip nasty RFI if requested and we're not masking all the channels */
    if ((s->clip_sigma > 0.0) && !(mask && (*nummasked == -1)))
        clip_times(currentdata, s->spectra_per_subint, s->num_channels,
                   s->clip_sigma, s->padvals);

    if (mask) {
        if (*nummasked == -1) { /* If all channels are masked */
            for (ii = 0; ii < s->spectra_per_subint; ii++)
                memcpy(currentdata + ii * s->num_channels,
                       s->padvals, s->num_channels * sizeof(float));
            // fftwf_execute_r2r(tplan1, currentdata, currentdata);
        } else if (*nummasked > 0) {    /* Only some of the channels are masked */
            int channum;
            for (ii = 0; ii < s->spectra_per_subint; ii++) {
                offset = ii * s->num_channels;
                for (jj = 0; jj < *nummasked; jj++) {
                    channum = maskchans[jj];
                    currentdata[offset + channum] = s->padvals[channum];
                }
            }
        }
    }

    if (s->num_ignorechans) { // These are channels we explicitly zero
        int channum;
        for (ii = 0; ii < s->spectra_per_subint; ii++) {
            offset = ii * s->num_channels;
            for (jj = 0; jj < s->num_ignorechans; jj++) {
                channum = s->ignorechans[jj];
                currentdata[offset + channum] = 0.0;
            }
        }
    }

    // In mpiprepsubband, the nodes do not call read_subbands() where
    // currentspectra gets incremented.
    if (using_MPI)
        currentspectra += s->spectra_per_subint;

    // Now transpose the raw block of data so that the times in each
    // channel are the most rapidly varying index
    fftwf_execute_r2r(tplan1, currentdata, currentdata);

    if (firsttime) {
    
        /* Write cache file*/
        int prepvalue = 0;
        SclData(currentdata, s->spectra_per_subint, s->num_channels, currentdata_data, currentdata_scl, currentdata_offs, thread);
        fwrite(currentdata_scl, sizeof(float), s->num_channels, s->cacheFile);
        fwrite(currentdata_offs, sizeof(float), s->num_channels, s->cacheFile);
        fwrite(currentdata_data, sizeof(unsigned char), s->spectra_per_subint * s->num_channels, s->cacheFile);
        
        firsttime = 0;
        SWAP(currentdata, lastdata);
        return 0;
    } else {

        /* Write cache file*/
        SclData(currentdata, s->spectra_per_subint, s->num_channels, currentdata_data, currentdata_scl, currentdata_offs, thread);
        fwrite(currentdata_scl, sizeof(float), s->num_channels, s->cacheFile);
        fwrite(currentdata_offs, sizeof(float), s->num_channels, s->cacheFile);
        fwrite(currentdata_data, sizeof(unsigned char), s->spectra_per_subint * s->num_channels, s->cacheFile);
        
        SWAP(currentdata, lastdata);
        // Transpose the resulting data into spectra as a function of time
        // if (transpose == 0)
        //     fftwf_execute_r2r(tplan2, fdata, fdata);
        return s->spectra_per_subint;
    }
}

int write_subbands_cache(float *fdata, int *delays, int numsubbands,
                  struct spectra_info *s, int transpose, int *padding,
                  int *maskchans, int *nummasked, mask * obsmask, int thread)
// This routine reads a spectral block/subint from the input raw data
// files. The routine uses dispersion delays in 'delays' to
// de-disperse the data into 'numsubbands' subbands.  It stores the
// resulting data in vector 'fdata' of length 'numsubbands' *
// 's->spectra_per_subint'.  The low freq subband is stored first,
// then the next highest subband etc, with 's->spectra_per_subint'
// floating points per subband.  It returns the # of points read if
// successful, 0 otherwise. If padding is returned as 1, then padding
// was added and statistics should not be calculated.  'maskchans' is
// an array of length numchans which contains a list of the number of
// channels that were masked.  The # of channels masked is returned in
// 'nummasked'.  'obsmask' is the mask structure to use for masking.
// If 'transpose'==0, the data will be kept in time order instead of
// arranged by subband as above.
{
    static int firsttime = 1;
    static float *frawdata;

    if (firsttime) {
        // Check to make sure there isn't more dispersion across a
        // subband than time in a block of data
        if (delays[0] > s->spectra_per_subint) {
            perror("\nError: there is more dispersion across a subband than time\n"
                   "in a block of data.  Increase spectra_per_subint if possible.");
            exit(-1);
        }
        // Needs to be twice as large for buffering if adding observations together
        frawdata = gen_fvect(2 * s->num_channels * s->spectra_per_subint);
        if (!s->get_rawblock(frawdata, s, padding)) {
            perror("Error: problem reading the raw data file in read_subbands()");
            int ifget_rawblock = 0;
            fwrite(padding, sizeof(int), 1, s->cacheFile); 
            fwrite(&ifget_rawblock, sizeof(int), 1, s->cacheFile);
            exit(-1);
        }
        else
        {
            int ifget_rawblock = 1;
            fwrite(padding, sizeof(int), 1, s->cacheFile); 
            fwrite(&ifget_rawblock, sizeof(int), 1, s->cacheFile);
        }

        if (0 != write_prep_subbands_cache(fdata, frawdata, delays, numsubbands, s,
                               transpose, maskchans, nummasked, obsmask, thread)) {
            perror("Error: problem initializing write_prep_subbands_cache() in read_subbands()");
            exit(-1);
        }
        firsttime = 0;
    }
    
    if (!s->get_rawblock(frawdata, s, padding)) {
        int ifget_rawblock = 0;
        fwrite(padding, sizeof(int), 1, s->cacheFile);
        fwrite(&ifget_rawblock, sizeof(int), 1, s->cacheFile);
        return 0;
    }
    else
    {
        int ifget_rawblock = 1;
        fwrite(padding, sizeof(int), 1, s->cacheFile); 
        fwrite(&ifget_rawblock, sizeof(int), 1, s->cacheFile);
    }

    if (write_prep_subbands_cache(fdata, frawdata, delays, numsubbands, s, transpose,
                      maskchans, nummasked, obsmask, thread) == s->spectra_per_subint) {
        currentspectra += s->spectra_per_subint;
        return s->spectra_per_subint;
    } else {
        return 0;
    }
}

int write_subbands_cache_log(float *fdata, int *delays, int numsubbands,
                  struct spectra_info *s, int transpose, int *padding,
                  int *maskchans, int *nummasked, mask * obsmask, int thread, long long *data_size, long *total_microseconds)
// This routine reads a spectral block/subint from the input raw data
// files. The routine uses dispersion delays in 'delays' to
// de-disperse the data into 'numsubbands' subbands.  It stores the
// resulting data in vector 'fdata' of length 'numsubbands' *
// 's->spectra_per_subint'.  The low freq subband is stored first,
// then the next highest subband etc, with 's->spectra_per_subint'
// floating points per subband.  It returns the # of points read if
// successful, 0 otherwise. If padding is returned as 1, then padding
// was added and statistics should not be calculated.  'maskchans' is
// an array of length numchans which contains a list of the number of
// channels that were masked.  The # of channels masked is returned in
// 'nummasked'.  'obsmask' is the mask structure to use for masking.
// If 'transpose'==0, the data will be kept in time order instead of
// arranged by subband as above.
{
    static int firsttime = 1;
    static float *frawdata;

    if (firsttime) {
        // Check to make sure there isn't more dispersion across a
        // subband than time in a block of data
        if (delays[0] > s->spectra_per_subint) {
            perror("\nError: there is more dispersion across a subband than time\n"
                   "in a block of data.  Increase spectra_per_subint if possible.");
            exit(-1);
        }
        // Needs to be twice as large for buffering if adding observations together
        frawdata = gen_fvect(2 * s->num_channels * s->spectra_per_subint);
        if (!s->get_rawblock_log(frawdata, s, padding, data_size, total_microseconds)) {
            perror("Error: problem reading the raw data file in read_subbands()");
            int ifget_rawblock = 0;
            fwrite(padding, sizeof(int), 1, s->cacheFile); 
            fwrite(&ifget_rawblock, sizeof(int), 1, s->cacheFile);
            exit(-1);
        }
        else
        {
            int ifget_rawblock = 1;
            fwrite(padding, sizeof(int), 1, s->cacheFile); 
            fwrite(&ifget_rawblock, sizeof(int), 1, s->cacheFile);
        }

        if (0 != write_prep_subbands_cache(fdata, frawdata, delays, numsubbands, s,
                               transpose, maskchans, nummasked, obsmask, thread)) {
            perror("Error: problem initializing write_prep_subbands_cache() in read_subbands()");
            exit(-1);
        }
        firsttime = 0;
    }
    
    if (!s->get_rawblock_log(frawdata, s, padding, data_size, total_microseconds)) {
        int ifget_rawblock = 0;
        fwrite(padding, sizeof(int), 1, s->cacheFile);
        fwrite(&ifget_rawblock, sizeof(int), 1, s->cacheFile);
        return 0;
    }
    else
    {
        int ifget_rawblock = 1;
        fwrite(padding, sizeof(int), 1, s->cacheFile); 
        fwrite(&ifget_rawblock, sizeof(int), 1, s->cacheFile);
    }

    if (write_prep_subbands_cache(fdata, frawdata, delays, numsubbands, s, transpose,
                      maskchans, nummasked, obsmask, thread) == s->spectra_per_subint) {
        currentspectra += s->spectra_per_subint;
        return s->spectra_per_subint;
    } else {
        return 0;
    }
}

// 保存 spectra_info 到文件
void serialize_spectra_info(struct spectra_info *s, const char *filename) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        perror("Failed to open file for writing");
        exit(EXIT_FAILURE);
    }

    // 保存标量字段
    fwrite(&(s->N), sizeof(long long), 1, file);
    fwrite(&(s->T), sizeof(double), 1, file);
    fwrite(&(s->dt), sizeof(double), 1, file);
    fwrite(&(s->fctr), sizeof(double), 1, file);
    fwrite(&(s->lo_freq), sizeof(double), 1, file);
    fwrite(&(s->hi_freq), sizeof(double), 1, file);
    fwrite(&(s->orig_df), sizeof(double), 1, file);
    fwrite(&(s->chan_dm), sizeof(double), 1, file);
    fwrite(&(s->df), sizeof(double), 1, file);
    fwrite(&(s->BW), sizeof(double), 1, file);
    fwrite(&(s->ra2000), sizeof(double), 1, file);
    fwrite(&(s->dec2000), sizeof(double), 1, file);
    fwrite(&(s->azimuth), sizeof(double), 1, file);
    fwrite(&(s->zenith_ang), sizeof(double), 1, file);
    fwrite(&(s->beam_FWHM), sizeof(double), 1, file);
    fwrite(&(s->time_per_subint), sizeof(double), 1, file);
    fwrite(&(s->datatype), sizeof(psrdatatype), 1, file);
    fwrite(&(s->scan_number), sizeof(int), 1, file);
    fwrite(&(s->tracking), sizeof(int), 1, file);
    fwrite(&(s->orig_num_chan), sizeof(int), 1, file);
    fwrite(&(s->num_channels), sizeof(int), 1, file);
    fwrite(&(s->num_polns), sizeof(int), 1, file);
    fwrite(&(s->num_beams), sizeof(int), 1, file);
    fwrite(&(s->beamnum), sizeof(int), 1, file);
    fwrite(&(s->summed_polns), sizeof(int), 1, file);
    fwrite(&(s->FITS_typecode), sizeof(int), 1, file);
    fwrite(&(s->bits_per_sample), sizeof(int), 1, file);
    fwrite(&(s->bytes_per_spectra), sizeof(int), 1, file);
    fwrite(&(s->samples_per_spectra), sizeof(int), 1, file);
    fwrite(&(s->bytes_per_subint), sizeof(int), 1, file);
    fwrite(&(s->spectra_per_subint), sizeof(int), 1, file);
    fwrite(&(s->samples_per_subint), sizeof(int), 1, file);
    fwrite(&(s->min_spect_per_read), sizeof(int), 1, file);
    fwrite(&(s->num_files), sizeof(int), 1, file);
    fwrite(&(s->offs_sub_col), sizeof(int), 1, file);
    fwrite(&(s->dat_wts_col), sizeof(int), 1, file);
    fwrite(&(s->dat_offs_col), sizeof(int), 1, file);
    fwrite(&(s->dat_scl_col), sizeof(int), 1, file);
    fwrite(&(s->data_col), sizeof(int), 1, file);
    fwrite(&(s->apply_scale), sizeof(int), 1, file);
    fwrite(&(s->apply_offset), sizeof(int), 1, file);
    fwrite(&(s->apply_weight), sizeof(int), 1, file);
    fwrite(&(s->apply_flipband), sizeof(int), 1, file);
    fwrite(&(s->signedints), sizeof(int), 1, file);
    fwrite(&(s->remove_zerodm), sizeof(int), 1, file);
    fwrite(&(s->use_poln), sizeof(int), 1, file);
    fwrite(&(s->flip_bytes), sizeof(int), 1, file);
    fwrite(&(s->num_ignorechans), sizeof(int), 1, file);
    fwrite(&(s->zero_offset), sizeof(float), 1, file);
    fwrite(&(s->clip_sigma), sizeof(float), 1, file);

    // 保存字符串字段
    fwrite(s->telescope, sizeof(char), 40, file);
    fwrite(s->observer, sizeof(char), 100, file);
    fwrite(s->source, sizeof(char), 100, file);
    fwrite(s->frontend, sizeof(char), 100, file);
    fwrite(s->backend, sizeof(char), 100, file);
    fwrite(s->project_id, sizeof(char), 40, file);
    fwrite(s->date_obs, sizeof(char), 40, file);
    fwrite(s->ra_str, sizeof(char), 40, file);
    fwrite(s->dec_str, sizeof(char), 40, file);
    fwrite(s->poln_type, sizeof(char), 40, file);
    fwrite(s->poln_order, sizeof(char), 40, file);
    fwrite(s->cacheFileName, sizeof(char), 1024, file);

    // 保存动态数组
    for (int i = 0; i < s->num_files; i++) {
        fwrite(&(s->start_MJD[i]), sizeof(long double), 1, file);
        fwrite(&(s->start_spec[i]), sizeof(long long), 1, file);
        fwrite(&(s->num_spec[i]), sizeof(long long), 1, file);
        fwrite(&(s->num_pad[i]), sizeof(long long), 1, file);
        fwrite(&(s->start_subint[i]), sizeof(int), 1, file);
        fwrite(&(s->num_subint[i]), sizeof(int), 1, file);
    }

    // 保存文件名数组
    for (int i = 0; i < s->num_files; i++) {
        size_t len = strlen(s->filenames[i]) + 1; // 包括 '\0'
        fwrite(&len, sizeof(size_t), 1, file);
        fwrite(s->filenames[i], sizeof(char), len, file);
    }

    // 保存 padvals 数组
    fwrite(s->padvals, sizeof(float), s->num_channels, file);

    fclose(file);
}

// 从文件加载 spectra_info
void deserialize_spectra_info(struct spectra_info *s, const char *filename) {
    FILE *file = fopen(filename, "rb");
    printf("spec file name: %s\n", filename);
    if (!file) {
        perror("Failed to open file for reading");
        exit(EXIT_FAILURE);
    }

    // 读取标量字段
    fread(&(s->N), sizeof(long long), 1, file);
    fread(&(s->T), sizeof(double), 1, file);
    fread(&(s->dt), sizeof(double), 1, file);
    fread(&(s->fctr), sizeof(double), 1, file);
    fread(&(s->lo_freq), sizeof(double), 1, file);
    fread(&(s->hi_freq), sizeof(double), 1, file);
    fread(&(s->orig_df), sizeof(double), 1, file);
    fread(&(s->chan_dm), sizeof(double), 1, file);
    fread(&(s->df), sizeof(double), 1, file);
    fread(&(s->BW), sizeof(double), 1, file);
    fread(&(s->ra2000), sizeof(double), 1, file);
    fread(&(s->dec2000), sizeof(double), 1, file);
    fread(&(s->azimuth), sizeof(double), 1, file);
    fread(&(s->zenith_ang), sizeof(double), 1, file);
    fread(&(s->beam_FWHM), sizeof(double), 1, file);
    fread(&(s->time_per_subint), sizeof(double), 1, file);
    fread(&(s->datatype), sizeof(psrdatatype), 1, file);
    fread(&(s->scan_number), sizeof(int), 1, file);
    fread(&(s->tracking), sizeof(int), 1, file);
    fread(&(s->orig_num_chan), sizeof(int), 1, file);
    fread(&(s->num_channels), sizeof(int), 1, file);
    fread(&(s->num_polns), sizeof(int), 1, file);
    fread(&(s->num_beams), sizeof(int), 1, file);
    fread(&(s->beamnum), sizeof(int), 1, file);
    fread(&(s->summed_polns), sizeof(int), 1, file);
    fread(&(s->FITS_typecode), sizeof(int), 1, file);
    fread(&(s->bits_per_sample), sizeof(int), 1, file);
    fread(&(s->bytes_per_spectra), sizeof(int), 1, file);
    fread(&(s->samples_per_spectra), sizeof(int), 1, file);
    fread(&(s->bytes_per_subint), sizeof(int), 1, file);
    fread(&(s->spectra_per_subint), sizeof(int), 1, file);
    fread(&(s->samples_per_subint), sizeof(int), 1, file);
    fread(&(s->min_spect_per_read), sizeof(int), 1, file);
    fread(&(s->num_files), sizeof(int), 1, file);
    fread(&(s->offs_sub_col), sizeof(int), 1, file);
    fread(&(s->dat_wts_col), sizeof(int), 1, file);
    fread(&(s->dat_offs_col), sizeof(int), 1, file);
    fread(&(s->dat_scl_col), sizeof(int), 1, file);
    fread(&(s->data_col), sizeof(int), 1, file);
    fread(&(s->apply_scale), sizeof(int), 1, file);
    fread(&(s->apply_offset), sizeof(int), 1, file);
    fread(&(s->apply_weight), sizeof(int), 1, file);
    fread(&(s->apply_flipband), sizeof(int), 1, file);
    fread(&(s->signedints), sizeof(int), 1, file);
    fread(&(s->remove_zerodm), sizeof(int), 1, file);
    fread(&(s->use_poln), sizeof(int), 1, file);
    fread(&(s->flip_bytes), sizeof(int), 1, file);
    fread(&(s->num_ignorechans), sizeof(int), 1, file);
    fread(&(s->zero_offset), sizeof(float), 1, file);
    fread(&(s->clip_sigma), sizeof(float), 1, file);

    // 读取字符串字段
    fread(s->telescope, sizeof(char), 40, file);
    fread(s->observer, sizeof(char), 100, file);
    fread(s->source, sizeof(char), 100, file);
    fread(s->frontend, sizeof(char), 100, file);
    fread(s->backend, sizeof(char), 100, file);
    fread(s->project_id, sizeof(char), 40, file);
    fread(s->date_obs, sizeof(char), 40, file);
    fread(s->ra_str, sizeof(char), 40, file);
    fread(s->dec_str, sizeof(char), 40, file);
    fread(s->poln_type, sizeof(char), 40, file);
    fread(s->poln_order, sizeof(char), 40, file);
    fread(s->cacheFileName, sizeof(char), 1024, file);

    // 初始化动态数组
    s->start_MJD = (long double *)malloc(sizeof(long double) * s->num_files);
    s->start_spec = (long long *)malloc(sizeof(long long) * s->num_files);
    s->num_spec = (long long *)malloc(sizeof(long long) * s->num_files);
    s->num_pad = (long long *)malloc(sizeof(long long) * s->num_files);
    s->start_subint = (int *)malloc(sizeof(int) * s->num_files);
    s->num_subint = (int *)malloc(sizeof(int) * s->num_files);

    for (int i = 0; i < s->num_files; i++) {
        fread(&(s->start_MJD[i]), sizeof(long double), 1, file);
        fread(&(s->start_spec[i]), sizeof(long long), 1, file);
        fread(&(s->num_spec[i]), sizeof(long long), 1, file);
        fread(&(s->num_pad[i]), sizeof(long long), 1, file);
        fread(&(s->start_subint[i]), sizeof(int), 1, file);
        fread(&(s->num_subint[i]), sizeof(int), 1, file);
    }

    // 读取文件名数组
    // s->filenames = (char **)malloc(sizeof(char *) * s->num_files);
    // for (int i = 0; i < s->num_files; i++) {
    //     size_t len;
    //     fread(&len, sizeof(size_t), 1, file);
    //     s->filenames[i] = (char *)malloc(len);
    //     fread(s->filenames[i], sizeof(char), len, file);
    // }

    // 初始化 padvals
    s->padvals = (float *)malloc(sizeof(float) * s->num_channels);
    fread(s->padvals, sizeof(float), s->num_channels, file);

    // 初始化 cacheFile
    // s->cacheFile = fopen(s->cacheFileName, "wb+");

    init_static_values(s);

    fclose(file);
}

int compare_spectra_info(const struct spectra_info *s1, const struct spectra_info *s2) {
    // 比较标量字段
    if (strcmp(s1->telescope, s2->telescope) != 0) return 0;
    if (strcmp(s1->observer, s2->observer) != 0) return 0;
    if (strcmp(s1->source, s2->source) != 0) return 0;
    if (strcmp(s1->frontend, s2->frontend) != 0) return 0;
    if (strcmp(s1->backend, s2->backend) != 0) return 0;
    if (strcmp(s1->project_id, s2->project_id) != 0) return 0;
    if (strcmp(s1->date_obs, s2->date_obs) != 0) return 0;
    if (strcmp(s1->ra_str, s2->ra_str) != 0) return 0;
    if (strcmp(s1->dec_str, s2->dec_str) != 0) return 0;
    if (strcmp(s1->poln_type, s2->poln_type) != 0) return 0;
    if (strcmp(s1->poln_order, s2->poln_order) != 0) return 0;

    if (s1->N != s2->N) return 0;
    if (fabs(s1->T - s2->T) > 1e-6) return 0; // 浮点数比较
    if (fabs(s1->dt - s2->dt) > 1e-6) return 0;
    if (fabs(s1->fctr - s2->fctr) > 1e-6) return 0;
    if (fabs(s1->lo_freq - s2->lo_freq) > 1e-6) return 0;
    if (fabs(s1->hi_freq - s2->hi_freq) > 1e-6) return 0;
    if (fabs(s1->orig_df - s2->orig_df) > 1e-6) return 0;
    if (fabs(s1->chan_dm - s2->chan_dm) > 1e-6) return 0;
    if (fabs(s1->df - s2->df) > 1e-6) return 0;
    if (fabs(s1->BW - s2->BW) > 1e-6) return 0;
    if (fabs(s1->ra2000 - s2->ra2000) > 1e-6) return 0;
    if (fabs(s1->dec2000 - s2->dec2000) > 1e-6) return 0;
    if (fabs(s1->azimuth - s2->azimuth) > 1e-6) return 0;
    if (fabs(s1->zenith_ang - s2->zenith_ang) > 1e-6) return 0;
    if (fabs(s1->beam_FWHM - s2->beam_FWHM) > 1e-6) return 0;
    if (fabs(s1->time_per_subint - s2->time_per_subint) > 1e-6) return 0;

    if (s1->datatype != s2->datatype) return 0;
    if (s1->scan_number != s2->scan_number) return 0;
    if (s1->tracking != s2->tracking) return 0;
    if (s1->orig_num_chan != s2->orig_num_chan) return 0;
    if (s1->num_channels != s2->num_channels) return 0;
    if (s1->num_polns != s2->num_polns) return 0;
    if (s1->num_beams != s2->num_beams) return 0;
    if (s1->beamnum != s2->beamnum) return 0;
    if (s1->summed_polns != s2->summed_polns) return 0;
    if (s1->FITS_typecode != s2->FITS_typecode) return 0;
    if (s1->bits_per_sample != s2->bits_per_sample) return 0;
    if (s1->bytes_per_spectra != s2->bytes_per_spectra) return 0;
    if (s1->samples_per_spectra != s2->samples_per_spectra) return 0;
    if (s1->bytes_per_subint != s2->bytes_per_subint) return 0;
    if (s1->spectra_per_subint != s2->spectra_per_subint) return 0;
    if (s1->samples_per_subint != s2->samples_per_subint) return 0;
    if (s1->min_spect_per_read != s2->min_spect_per_read) return 0;
    if (s1->num_files != s2->num_files) return 0;

    // 比较指针字段
    // 比较 start_MJD 数组
    if (s1->start_MJD != NULL && s2->start_MJD != NULL) {
        for (int i = 0; i < s1->num_files; i++) {
            if (fabs(s1->start_MJD[i] - s2->start_MJD[i]) > 1e-6) return 0;
        }
    }

    // 比较 filenames 数组
    for (int i = 0; i < s1->num_files; i++) {
        if (strcmp(s1->filenames[i], s2->filenames[i]) != 0) return 0;
    }

    // 比较 ignorechans 数组
    if (s1->num_ignorechans == s2->num_ignorechans && s1->ignorechans != NULL && s2->ignorechans != NULL) {
        for (int i = 0; i < s1->num_ignorechans; i++) {
            if (s1->ignorechans[i] != s2->ignorechans[i]) return 0;
        }
    }

    // 比较 ignorechans_str
    if (s1->ignorechans_str != NULL && s2->ignorechans_str != NULL) {
        if (strcmp(s1->ignorechans_str, s2->ignorechans_str) != 0) return 0;
    }

    // 比较 start_spec 数组
    if (s1->start_spec != NULL && s2->start_spec != NULL) {
        for (int i = 0; i < s1->num_files; i++) {
            if (s1->start_spec[i] != s2->start_spec[i]) return 0;
        }
    }

    // 比较 num_spec 数组
    if (s1->num_spec != NULL && s2->num_spec != NULL) {
        for (int i = 0; i < s1->num_files; i++) {
            if (s1->num_spec[i] != s2->num_spec[i]) return 0;
        }
    }

    // 比较 num_pad 数组
    if (s1->num_pad != NULL && s2->num_pad != NULL) {
        for (int i = 0; i < s1->num_files; i++) {
            if (s1->num_pad[i] != s2->num_pad[i]) return 0;
        }
    }

    // 如果所有字段都一致
    return 1;
}

int read_prep_subbands_cache(float *fdata, float *rawdata, int *delays, int numsubbands,
                  struct spectra_info *s, int transpose,
                  int *maskchans, int *nummasked, mask * obsmask)
// This routine preps a block of raw spectra for subbanding.  It uses
// dispersion delays in 'delays' to de-disperse the data into
// 'numsubbands' subbands.  It stores the resulting data in vector
// 'fdata' of length 'numsubbands' * 's->spectra_per_subint'.  The low
// freq subband is stored first, then the next highest subband etc,
// with 's->spectra_per_subint' floating points per subband. It
// returns the # of points read if succesful, 0 otherwise.
// 'maskchans' is an array of length numchans which contains a list of
// the number of channels that were masked.  The # of channels masked
// is returned in 'nummasked'.  'obsmask' is the mask structure to use
// for masking.  If 'transpose'==0, the data will be kept in time
// order instead of arranged by subband as above.
{
    int ii, jj, offset;
    double starttime = 0.0;
    static float *tmpswap;// *rawdata1, *rawdata2;
    // static float *currentdata, *lastdata;

    static float *currentdata_scl, *lastdata_scl;
    static float *currentdata_offs, *lastdata_offs;
    static unsigned char *currentdata_data, *lastdata_data;

    static int firsttime = 1, mask = 0;
    static fftwf_plan tplan1, tplan2;

    *nummasked = 0;
    if (firsttime) {
        if (obsmask->numchan)
            mask = 1;
            
        currentdata_scl = malloc(s->num_channels * sizeof(float));
        currentdata_offs = malloc(s->num_channels * sizeof(float));
        currentdata_data = malloc(s->spectra_per_subint * s->num_channels * sizeof(unsigned char));

        lastdata_scl = malloc(s->num_channels * sizeof(float));
        lastdata_offs = malloc(s->num_channels * sizeof(float));
        lastdata_data = malloc(s->spectra_per_subint * s->num_channels * sizeof(unsigned char));

        tplan2 = plan_transpose(numsubbands, s->spectra_per_subint, fdata, fdata);
    }

    /* Read and de-disperse */
    fread(currentdata_scl, sizeof(float), s->num_channels, s->cacheFile);
    fread(currentdata_offs, sizeof(float), s->num_channels, s->cacheFile);
    fread(currentdata_data, sizeof(unsigned char), s->spectra_per_subint * s->num_channels, s->cacheFile);
    
        
    starttime = currentspectra * s->dt; // or -1 subint?

    // In mpiprepsubband, the nodes do not call read_subbands() where
    // currentspectra gets incremented.
    if (using_MPI)
        currentspectra += s->spectra_per_subint;

    if (firsttime) {
        SWAP(currentdata_scl, lastdata_scl);
        SWAP(currentdata_offs, lastdata_offs);
        SWAP(currentdata_data, lastdata_data);
        
        firsttime = 0;
        return 0;
    } else {
        dedisp_subbands_cache(currentdata_data, currentdata_scl, currentdata_offs, lastdata_data, lastdata_scl, lastdata_offs, s->spectra_per_subint,
                        s->num_channels, delays, numsubbands, fdata);
        // SWAP(currentdata, lastdata);
        SWAP(currentdata_scl, lastdata_scl);
        SWAP(currentdata_offs, lastdata_offs);
        SWAP(currentdata_data, lastdata_data);

        // Transpose the resulting data into spectra as a function of time
        if (transpose == 0)
            fftwf_execute_r2r(tplan2, fdata, fdata);
        return s->spectra_per_subint;
    }
}

int read_subbands_cache(float *fdata, int *delays, int numsubbands,
                  struct spectra_info *s, int transpose, int *padding,
                  int *maskchans, int *nummasked, mask * obsmask)
// This routine reads a spectral block/subint from the input raw data
// files. The routine uses dispersion delays in 'delays' to
// de-disperse the data into 'numsubbands' subbands.  It stores the
// resulting data in vector 'fdata' of length 'numsubbands' *
// 's->spectra_per_subint'.  The low freq subband is stored first,
// then the next highest subband etc, with 's->spectra_per_subint'
// floating points per subband.  It returns the # of points read if
// successful, 0 otherwise. If padding is returned as 1, then padding
// was added and statistics should not be calculated.  'maskchans' is
// an array of length numchans which contains a list of the number of
// channels that were masked.  The # of channels masked is returned in
// 'nummasked'.  'obsmask' is the mask structure to use for masking.
// If 'transpose'==0, the data will be kept in time order instead of
// arranged by subband as above.
{
    static int firsttime = 1;
    static float *frawdata;
    int ifget_rawblock;

    if (firsttime) {
        // Check to make sure there isn't more dispersion across a
        // subband than time in a block of data
        if (delays[0] > s->spectra_per_subint) {
            perror("\nError: there is more dispersion across a subband than time\n"
                   "in a block of data.  Increase spectra_per_subint if possible.");
            exit(-1);
        }
        // Needs to be twice as large for buffering if adding observations together
        // frawdata = gen_fvect(2 * s->num_channels * s->spectra_per_subint);

        fread(padding, sizeof(int), 1, s->cacheFile);
        fread(&ifget_rawblock, sizeof(int), 1, s->cacheFile);

        if (!ifget_rawblock) {
            perror("Error: problem reading the raw data file in read_subbands()");
            exit(-1);
        }
        if (0 != read_prep_subbands_cache(fdata, frawdata, delays, numsubbands, s,
                               transpose, maskchans, nummasked, obsmask)) {
            perror("Error: problem initializing read_prep_subbands_cache() in read_subbands()");
            exit(-1);
        }
        firsttime = 0;
    }

    fread(padding, sizeof(int), 1, s->cacheFile);
    fread(&ifget_rawblock, sizeof(int), 1, s->cacheFile);


    if (!ifget_rawblock) {
        return 0;
    }
    if (read_prep_subbands_cache(fdata, frawdata, delays, numsubbands, s, transpose,
                      maskchans, nummasked, obsmask) == s->spectra_per_subint) {
        currentspectra += s->spectra_per_subint;
        return s->spectra_per_subint;
    } else {
        return 0;
    }
}

int read_prep_subbands_cache_log(float *fdata, float *rawdata, int *delays, int numsubbands,
                  struct spectra_info *s, int transpose,
                  int *maskchans, int *nummasked, mask * obsmask, long long *data_size, long *total_microseconds)
// This routine preps a block of raw spectra for subbanding.  It uses
// dispersion delays in 'delays' to de-disperse the data into
// 'numsubbands' subbands.  It stores the resulting data in vector
// 'fdata' of length 'numsubbands' * 's->spectra_per_subint'.  The low
// freq subband is stored first, then the next highest subband etc,
// with 's->spectra_per_subint' floating points per subband. It
// returns the # of points read if succesful, 0 otherwise.
// 'maskchans' is an array of length numchans which contains a list of
// the number of channels that were masked.  The # of channels masked
// is returned in 'nummasked'.  'obsmask' is the mask structure to use
// for masking.  If 'transpose'==0, the data will be kept in time
// order instead of arranged by subband as above.
{
    int ii, jj, offset;
    double starttime = 0.0;
    static float *tmpswap;// *rawdata1, *rawdata2;
    // static float *currentdata, *lastdata;

    static float *currentdata_scl, *lastdata_scl;
    static float *currentdata_offs, *lastdata_offs;
    static unsigned char *currentdata_data, *lastdata_data;

    static int firsttime = 1, mask = 0;
    static fftwf_plan tplan1, tplan2;

    *nummasked = 0;
    if (firsttime) {
        if (obsmask->numchan)
            mask = 1;
            
        currentdata_scl = malloc(s->num_channels * sizeof(float));
        currentdata_offs = malloc(s->num_channels * sizeof(float));
        currentdata_data = malloc(s->spectra_per_subint * s->num_channels * sizeof(unsigned char));

        lastdata_scl = malloc(s->num_channels * sizeof(float));
        lastdata_offs = malloc(s->num_channels * sizeof(float));
        lastdata_data = malloc(s->spectra_per_subint * s->num_channels * sizeof(unsigned char));

        tplan2 = plan_transpose(numsubbands, s->spectra_per_subint, fdata, fdata);
    }

    /* Read and de-disperse */    
    struct timeval start, end; // 定义两个时间结构体
    // 获取开始时间
    gettimeofday(&start, NULL);
    fread(currentdata_scl, sizeof(float), s->num_channels, s->cacheFile);
    fread(currentdata_offs, sizeof(float), s->num_channels, s->cacheFile);
    fread(currentdata_data, sizeof(unsigned char), s->spectra_per_subint * s->num_channels, s->cacheFile);    
    // 获取结束时间
    gettimeofday(&end, NULL);

    // 计算执行时间（单位：微秒）
    long seconds = end.tv_sec - start.tv_sec;            // 秒部分的差值
    long microseconds = end.tv_usec - start.tv_usec;     // 微秒部分的差值
    *total_microseconds += seconds * 1000000 + microseconds; // 转换为总微秒
    *data_size += s->num_channels * sizeof(float) * 2 + s->spectra_per_subint * s->num_channels * sizeof(unsigned char);
        
    starttime = currentspectra * s->dt; // or -1 subint?

    // In mpiprepsubband, the nodes do not call read_subbands() where
    // currentspectra gets incremented.
    if (using_MPI)
        currentspectra += s->spectra_per_subint;

    if (firsttime) {
        SWAP(currentdata_scl, lastdata_scl);
        SWAP(currentdata_offs, lastdata_offs);
        SWAP(currentdata_data, lastdata_data);
        
        firsttime = 0;
        return 0;
    } else {
        dedisp_subbands_cache(currentdata_data, currentdata_scl, currentdata_offs, lastdata_data, lastdata_scl, lastdata_offs, s->spectra_per_subint,
                        s->num_channels, delays, numsubbands, fdata);
        // SWAP(currentdata, lastdata);
        SWAP(currentdata_scl, lastdata_scl);
        SWAP(currentdata_offs, lastdata_offs);
        SWAP(currentdata_data, lastdata_data);

        // Transpose the resulting data into spectra as a function of time
        if (transpose == 0)
            fftwf_execute_r2r(tplan2, fdata, fdata);
        return s->spectra_per_subint;
    }
}

int read_subbands_cache_log(float *fdata, int *delays, int numsubbands,
                  struct spectra_info *s, int transpose, int *padding,
                  int *maskchans, int *nummasked, mask * obsmask, long long *data_size, long *total_microseconds)
// This routine reads a spectral block/subint from the input raw data
// files. The routine uses dispersion delays in 'delays' to
// de-disperse the data into 'numsubbands' subbands.  It stores the
// resulting data in vector 'fdata' of length 'numsubbands' *
// 's->spectra_per_subint'.  The low freq subband is stored first,
// then the next highest subband etc, with 's->spectra_per_subint'
// floating points per subband.  It returns the # of points read if
// successful, 0 otherwise. If padding is returned as 1, then padding
// was added and statistics should not be calculated.  'maskchans' is
// an array of length numchans which contains a list of the number of
// channels that were masked.  The # of channels masked is returned in
// 'nummasked'.  'obsmask' is the mask structure to use for masking.
// If 'transpose'==0, the data will be kept in time order instead of
// arranged by subband as above.
{
    static int firsttime = 1;
    static float *frawdata;
    int ifget_rawblock;

    if (firsttime) {
        // Check to make sure there isn't more dispersion across a
        // subband than time in a block of data
        if (delays[0] > s->spectra_per_subint) {
            perror("\nError: there is more dispersion across a subband than time\n"
                   "in a block of data.  Increase spectra_per_subint if possible.");
            exit(-1);
        }
        // Needs to be twice as large for buffering if adding observations together
        // frawdata = gen_fvect(2 * s->num_channels * s->spectra_per_subint);

        fread(padding, sizeof(int), 1, s->cacheFile);
        fread(&ifget_rawblock, sizeof(int), 1, s->cacheFile);

        if (!ifget_rawblock) {
            perror("Error: problem reading the raw data file in read_subbands()");
            exit(-1);
        }
        if (0 != read_prep_subbands_cache_log(fdata, frawdata, delays, numsubbands, s,
                               transpose, maskchans, nummasked, obsmask, data_size, total_microseconds)) {
            perror("Error: problem initializing read_prep_subbands_cache() in read_subbands()");
            exit(-1);
        }
        firsttime = 0;
    }

    fread(padding, sizeof(int), 1, s->cacheFile);
    fread(&ifget_rawblock, sizeof(int), 1, s->cacheFile);


    if (!ifget_rawblock) {
        return 0;
    }
    if (read_prep_subbands_cache_log(fdata, frawdata, delays, numsubbands, s, transpose,
                      maskchans, nummasked, obsmask, data_size, total_microseconds) == s->spectra_per_subint) {
        currentspectra += s->spectra_per_subint;
        return s->spectra_per_subint;
    } else {
        return 0;
    }
}
