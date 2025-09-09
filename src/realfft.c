/*     Real-Valued Data FFT Program        */
/*          by Scott Ransom                */
/*            Version 3.0                  */

#include <time.h>
#include <sys/times.h>
#include "clk_tck.h"
#include "misc_utils.h"
#include "chkio.h"
#include "ransomfft.h"
#include "vectors.h"
#include "realfft_cmd.h"

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/*  This program calculates the FFT of a file containing    */
/*  a number of single-precision floats representing        */
/*  real numbers.  (i.e. a normal time series)              */
/*  The data is assumed to be located in the directory,     */
/*  specified in the input filename(s).                     */
/*  Input filename(s) must include '.dat' or '.fft'         */
/*        suffixes.  The output file(s) will have the       */
/*        appropriate other suffix.                         */
/*  Do not end paths in '/'.                                */
/*  Scratch file(s) are the same size as the input file(s). */
/*  If '-inv' is specified, the file to be transformed      */
/*        should end in '.fft'.  Otherwise, it should end   */
/*        in '.dat'.                                        */

static char* log_timestamp() {
    static char buffer[20]; // 注意使用 static 以确保返回的字符串在调用后仍然有效
    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", timeinfo);

    return buffer;
}

int realfft_main(int argc, char *argv[]) 
{
    FILE *datfile, *tmpfile = NULL, *outfile;
    char *datdir = NULL, *datfilenm = NULL, *tmpfilenm = NULL, *outfilenm = NULL;
    float *data;
    int isign = -1, numfiles;
    long long numdata = 0, filelen, maxfilelen = 0;
    struct tms runtimes;
    double ttim, stim, utim, tott;
    Cmdline *cmd;
    char datsuffix[] = "dat";
    char outsuffix[] = "fft";
    char tmpsuffix[] = "tmp";

    /* Call usage() if we have no command line arguments */

    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(1);
    }

    /* Parse the command line using the excellent program Clig */

    cmd = parseCmdline(argc, argv);
    printf("[%s] Realfft Start: %s\n", log_timestamp(), cmd->full_cmd_line);

#ifdef DEBUG
    showOptionValues();
#endif

    tott = times(&runtimes) / (double) CLK_TCK;
    // printf("\n");
    // printf("   Real-Valued Data FFT Program v3.0\n");
    // printf("        by Scott M. Ransom\n\n");

    /* Get our file information */
    long long data_size = 0;
    long total_microseconds = 0;
    struct timeval start, end; // 定义两个时间结构体

    numfiles = cmd->argc;
    for (int fi = 0 ; fi < numfiles ; fi++) {
        printf("[%s] %4d:  Processing data in '%s'\n", log_timestamp(), fi + 1, cmd->argv[fi]);
        {
            int hassuffix = 0, filenmlen;
            char *filenm, *root, *suffix;

            split_path_file(cmd->argv[fi], &datdir, &filenm);
            hassuffix = split_root_suffix(filenm, &root, &suffix);
            if (hassuffix) {
                if (strcmp(suffix, "fft") == 0) {
                    isign = 1;
                    strcpy(datsuffix, "fft");
                    strcpy(outsuffix, "dat");
                }
                free(suffix);
            }
            free(filenm);
            filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
            datfilenm = (char *) calloc(filenmlen, 1);
            sprintf(datfilenm, "%s/%s.%s", datdir, root, datsuffix);
            if (cmd->tmpdirP) {
                filenmlen = strlen(cmd->tmpdir) + 1 + strlen(root) + 5;
                tmpfilenm = (char *) calloc(filenmlen, 1);
                sprintf(tmpfilenm, "%s/%s.%s", cmd->tmpdir, root, tmpsuffix);
            } else {
                filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
                tmpfilenm = (char *) calloc(filenmlen, 1);
                sprintf(tmpfilenm, "%s/%s.%s", datdir, root, tmpsuffix);
            }
            if (cmd->outdirP) {
                filenmlen = strlen(cmd->outdir) + 1 + strlen(root) + 5;
                outfilenm = (char *) calloc(filenmlen, 1);
                sprintf(outfilenm, "%s/%s.%s", cmd->outdir, root, outsuffix);
            } else {
                filenmlen = strlen(datdir) + 1 + strlen(root) + 5;
                outfilenm = (char *) calloc(filenmlen, 1);
                sprintf(outfilenm, "%s/%s.%s", datdir, root, outsuffix);
            }
            free(root);
            free(datdir);
        }

        /* Force a forward or inverse transform.   */
        /* Note that the suffixes do _not_ change! */

        if (cmd->forwardP)
            isign = -1;
        if (cmd->inverseP)
            isign = 1;
        if (cmd->diskfftP && cmd->memfftP) {
            printf("\n[%s] You cannot take both an in- and out-of-core FFT!\n\n", log_timestamp());
            exit(1);
        }

        /* Open and check data files */

        datfile = chkfopen(datfilenm, "r");
        filelen = chkfilelen(datfile, 1);
        if (filelen > maxfilelen)
            maxfilelen = filelen;
        numdata = filelen / sizeof(float);
        if (isign == -1) {
            if (filelen % sizeof(float)) {
                printf("\n[%s] Input file does not contain the correct number of\n", log_timestamp());
                printf("   bytes for it to be floating point data!  Skipping: %s.\n\n", cmd->full_cmd_line);
                continue;
            }
            printf("\n[%s] Data OK.  There are %lld floats.\n\n", log_timestamp(), numdata);
        } else {
            if (filelen % sizeof(fcomplex)) {
                printf("\n[%s] Input file does not contain the correct number of\n", log_timestamp());
                printf("   bytes for it to be single precision complex data!  Skipping: %s.\n\n", cmd->full_cmd_line);
                continue;
            }
            printf("\n[%s] Data OK.  There are %lld complex points.\n\n", log_timestamp(), numdata / 2);
        }
        printf("[%s] %4d:   Result will be in '%s'\n", log_timestamp(), fi + 1, outfilenm);

        /*  Start the transform sequence  */

        if ((numdata > MAXREALFFT || cmd->diskfftP) && !cmd->memfftP) {

            /*  Perform Two-Pass, Out-of-Core, FFT  */

            if (isign == -1) {
                printf("\n[%s] Performing out-of-core two-pass forward FFT on data.\n", log_timestamp());
            } else {
                printf("\n[%s] Performing out-of-core two-pass inverse FFT on data.\n", log_timestamp());
            }

            /* Copy the input files if we want to keep them */

            if (!cmd->deleteP) {
                int status, slen, suf;
                char *cmd, *root, *suffix;

                slen = strlen(datfilenm);
                cmd = (char *) calloc(2 * slen + 10, 1);
                suf = split_root_suffix(datfilenm, &root, &suffix);
                sprintf(cmd, "cp %s %s.bak", datfilenm, root);
                if ((status = (system(cmd))) == -1 || status == 127) {
                    perror("\nSystem call (cp) failed");
                    printf("\n");
                    exit(1);
                }
                if (suf)
                    free(suffix);
                free(root);
                free(cmd);
            }

            /* Close the input files and re-open them in write mode */

            fclose(datfile);
            datfile = chkfopen(datfilenm, "r+");
            tmpfile = chkfopen(tmpfilenm, "w+");
            if (isign == 1) {
                realfft_scratch_inv(datfile, tmpfile, numdata);
            } else {
                realfft_scratch_fwd(datfile, tmpfile, numdata);
            }

            /* Remove the scratch files */

            fclose(tmpfile);
            remove(tmpfilenm);

            /* Change the output filename to the correct suffix and   */
            /* rename the back-up data files if needed.               */

            {
                int slen, suf;
                char *file1, *file2, *root, *suffix;

                slen = strlen(datfilenm);
                file1 = (char *) calloc(slen + 5, 1);
                file2 = (char *) calloc(slen + 5, 1);
                suf = split_root_suffix(datfilenm, &root, &suffix);
                sprintf(file1, "%s.%s", root, datsuffix);
                sprintf(file2, "%s.%s", root, outsuffix);
                rename(file1, file2);
                if (!cmd->deleteP) {
                    sprintf(file1, "%s.bak", root);
                    rename(file1, datfilenm);
                }
                if (suf)
                    free(suffix);
                free(root);
                free(file1);
                free(file2);
            }

        } else {

            /* Perform standard FFT for real functions  */

            outfile = chkfopen(outfilenm, "w");
            if (isign == -1) {
                printf("\n[%s] Performing in-core forward FFT on data:\n", log_timestamp());
            } else {
                printf("\n[%s] Performing in-core inverse FFT on data:\n", log_timestamp());
            }
            printf("[%s]    Reading.\n", log_timestamp());
            data = gen_fvect(numdata);
            if (cmd->IOlogP){
                gettimeofday(&start, NULL);
                chkfread(data, sizeof(float), numdata, datfile);
                gettimeofday(&end, NULL);
                long seconds = end.tv_sec - start.tv_sec;            // 秒部分的差值
                long microseconds = end.tv_usec - start.tv_usec;     // 微秒部分的差值
                total_microseconds += seconds * 1000000 + microseconds; // 转换为总微秒
                data_size += numdata * sizeof(float);
            }else{
                chkfread(data, sizeof(float), numdata, datfile);
            }
            printf("[%s]    Transforming.\n", log_timestamp());
            realfft(data, numdata, isign);
            /* fftwcall((fcomplex *)data, numdata/2, isign); */
            /* tablesixstepfft((fcomplex *)data, numdata/2, isign); */
            printf("[%s]    Writing.\n", log_timestamp());
            chkfwrite(data, sizeof(float), numdata, outfile);
            fclose(outfile);

            /* Delete the input files if requested */

            if (cmd->deleteP)
                remove(datfilenm);
        }

        /* Close our input files */

        fclose(datfile);
    }

    /* Output the timing information */

    printf("[%s] Realfft Finished:%s.\n\n", log_timestamp(), cmd->full_cmd_line);
    printf("Timing summary:\n");
    tott = times(&runtimes) / (double) CLK_TCK - tott;
    utim = runtimes.tms_utime / (double) CLK_TCK;
    stim = runtimes.tms_stime / (double) CLK_TCK;
    ttim = utim + stim;
    printf("  CPU usage: %.3f sec total (%.3f sec user, %.3f sec system)\n",
           ttim, utim, stim);
    printf("  Total time elapsed:  %.3f sec (%.3f sec/file)\n\n", tott, tott / numfiles);

    if(cmd->IOlogP){
        if ((numdata > MAXREALFFT || cmd->diskfftP) && !cmd->memfftP){
            printf("[%s] IOlog: numdata too large and use out-core forward FFT: %s\n", log_timestamp(), cmd->full_cmd_line);
        }else{
            printf("[%s] IOlog: %s read %.3f MB data, use %.3f ms, %.3f GB/s\n", log_timestamp(), cmd->full_cmd_line, (double)data_size/(1024.0*1024.0), (double)total_microseconds/(1000), ((double)data_size/(1024.0*1024.0*1024.0))/((double)total_microseconds/(1000000)));
        }
    }

    /*
       fftw_print_max_memory_usage();
       fftw_check_memory_leaks();
     */

    free(datfilenm);
    free(tmpfilenm);
    free(outfilenm);
    return 0;
}

/* 定义共享库接口 */
int call_realfft(int argc, char *argv[]) {
    return realfft_main(argc, argv);
}

/* 保留原来的 main 函数 */
int main(int argc, char *argv[]) {
    return realfft_main(argc, argv);
}