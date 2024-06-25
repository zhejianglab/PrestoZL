/*  $Id$   */

#include <stdio.h> 
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/stat.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <fcntl.h>

#ifdef HAVE_STDINT_H
#include <stdint.h>
typedef int64_t offs_type;
#else
typedef long long offs_type;
#endif

/*  FORTRAN:  fd = close(filedes)      */
int close_(filedes)
int *filedes;
{
return(close(*filedes));
}

/*  FORTRAN:  fd = open(filnam,flags,mode)  */
int open_(filnam,flags,mode)
char filnam[];
int *flags;
int *mode;
{
  return(open(filnam,*flags,*mode));
}

/*  FORTRAN:  nflag = isetflag() */
/*  Set flag for writing to a file.  Specifically:
      - read and write to file
      - create file if it doesn't exist
      - truncate any existing file 
    This used to be hard-wired in Fortran code, but different libraries
    use different flag values.  
*/
int isetflag_()
{
  return((int)(O_TRUNC | O_CREAT | O_RDWR));
}

/* FORTRAN:  fd = creat(filnam,mode) */
int creat_(filnam,mode)
char filnam[];
int *mode;
{
  return(creat(filnam,*mode));
}
/* FORTRAN:  nread = read(fd,buf,n) */
int read_(fd,buf,n)
int *fd,*n;
char buf[];
{
  return(read(*fd,buf,*n));
}
/* FORTRAN:  nwrt = write(fd,buf,n) */
int write_(fd,buf,n)
int *fd,*n;
char buf[];
{
  return(write(*fd,buf,*n));
}
/* FORTRAN: ns = lseek(fd,offset,origin) */
int lseek_(fd,offset,origin)
int *fd,*offset,*origin;
{
  return(lseek(*fd,*offset,*origin));
}
/* times(2) */
int times_(buf)
int buf[];
{
  return (times((struct tms *)buf));
}
/* returns 32 bit random numbers to FORTRAN program */
int iran_(xsubi)
unsigned short xsubi[];
{
  return (jrand48(xsubi));
}
int exit_(n)
int *n;
{
  printf("\n\n");
  exit(*n);
}

void cprnt_(string,n)
char string[];
int *n;
{
  int i;

  for(i=0;i< *n;i++)
    printf("%s",string);
  fflush(stdout);
  return;
}

void byterev_(i,n)
int *i,*n;
{
int j,k;

for(j=0;j<*n;j++) {
  k  = (*i<<24) | ((*i>>24)&255);
  k  = k | ((*i&65280)<<8);
  *i = k | ((*i>>8)&65280);
  i++;
}
}

void dbyterev_(i,n)
int *i,*n;
{
int j,k,l;

for(j=0;j<*n;j++) {
  k = (*i<<24) | ((*i>>24)&255);
  k = k | ((*i&65280)<<8);
  k = k | ((*i>>8)&65280);
  i++;
  l = (*i<<24) | ((*i>>24)&255);
  l = l | ((*i&65280)<<8);
  l = l | ((*i>>8)&65280);
  *(i-1) = l;
  *i = k;
  i++;
}
}

void wswap_(i,n)
int *i,*n;
{
int j,k;

for(j=0;j<*n;j++) {
  k = *i;
  *i = *(i+1);
  *(i+1) = k;
  i+=2;
}
}

/*

double hrtime_()
{
  int tv[2],tz[2];
  gettimeofday(tv,tz);
  return(tv[0]+1.e-6*tv[1]);
}

*/

int time_()
{
  return(time(0));
}

/*
int lstat_(char *filename, int *buf)
{
  return(lstat(filename,buf));
}
*/

void usleep_(long int *n)
{
  usleep(*n);
}




/* mallocx and freex are variants on code suggested
   by Richard Dodgson (U.Tasmania)  

   sample use:  if a real*8 array of size N is needed:
      integer ipointer, aoff
      real*8 a(1)
      ipointer = mallocx(a(1),n,8,aoff)

   ... we will compare its use to array b(1:n), which is allocated
       in the normal way:
      real*8 b(n)

   ... now use a(1+aoff) to a(n+aoff) as array references, i.e.,
       add aoff to the index of all references to a(); for example,
       to set the array to 1..n, use:

      do i = 1, n
         a(i+aoff) = i
         b(i) = i
      end do

   ...   when sending the array to a subroutine, use a(1+aoff) in
         the subroutine call [normally, one would write a(1) or a]
      call somesub(a(1+aoff))
      call somesub(b(1))   ...or...
      call somesub(b)

      call freex(ipointer)

   Finally, note that linux complains if the fortran code calling
   mallocx or freex uses different types of arguments in different calls.  
   To eliminate annoying warning messages, use mallocxi and freexi
   for integers and mallocxd and freexd for doubles.
*/



/* use global storage to retain "original" address of arrays and
   "malloc'd" address of arrays  */

#define MAXADDR 10
int naddr=0;
void *addr1[MAXADDR], *addr2[MAXADDR];



void mallocx_(void *a, int *nelem, int *size, offs_type *idx) {

  offs_type n;   /* change to 64 bit */

  if (naddr==MAXADDR) {
    printf ("Error: can't allocate more than %d arrays\n", MAXADDR); exit (1); }

  /* Following prevents use of mmap() in gcc/glibc.  This is necessary
     because use of mmap on alphas may produce pointers separated by
     more than 2^31, which can't be handled by fortran INT*4     */
#ifdef M_MMAP_MAX
  mallopt(M_MMAP_MAX,0);
#endif

  addr1[naddr] = a;                          /* store original address */
  addr2[naddr] = malloc((*nelem)*(*size));   /* allocate the memory    */
  if (addr2[naddr]==(void *)NULL)
    { printf ("Error: out of memory\n"); exit(1); }
       /* char* cast in next line forces calculation to be done in bytes */

                                             /* calculate array offset */
  n = (offs_type)((char *)addr2[naddr]-(char *)addr1[naddr])/(*size);
  // Don't need this with 64-bit addresses:
  //    /* catch cases when n can't be cast to int in 64-bit systems */
  //if (n>2147483647 || n<-2147483647)    
  //  { printf ("Error: array address more than 2^31 from original location\n");
  //    exit(1); }
  *idx = n;

  naddr++;
}

void freex_(void *a) {
  int i;
  int flag;
  flag = 1;
  for (i=0; i<naddr; i++) {
    if (addr1[i]==a) {
      free(addr2[i]);
      naddr--;    /* shorten the array */
      addr1[i] = addr1[naddr];
      addr2[i] = addr2[naddr];
      flag = 0;
      break;
    }
  }
  if (flag)
    printf ("Warning: can't free array\n"); 
}




void freexi_(int *a) {
  freex_((void *)a);
}

void freexd_(double *a) {
  freex_((void *)a);
}


void mallocxi_(int *a, int *nelem, int *size, offs_type *idx) {
  mallocx_((void *)a, nelem, size, idx);
}

void mallocxd_(double *a, int *nelem, int *size, offs_type *idx) {
  mallocx_((void *)a, nelem, size, idx);
}

