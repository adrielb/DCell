#ifndef _MISC_H
#define _MISC_H
#include <stdio.h>

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif
#define ABS(a) ((a)>=0?(a):-(a))
#define SQR(x) ((x)*(x))
#define SIGN(x) ((x)==0?0:(x)<0?-1:1)
#define SIGN2(x,y) ((x)*SIGN(y))
#define IMPLIES(a,b) (!(a)||(b))
/*static __inline__ double sqr(double x) { return x*x; }*/
#define MACH_EPS 1e-16  /* something that really should be in math.h */


/* FA accesses a 1D FORTRAN array using FORTRAN semantics, ie, from 1
 * If you compiled your FORTRAN code using -r8, then integers take up
 * 64 bits, although half of them are not used. (Repeat after me:
 * "I love FORTRAN, honest I do!")
 */

#ifdef sgi
#define R8 1    /* 1 for no, 2 for yes if you compiled f77 -r8 */
#else
#define R8 1
#endif

#define IFA(name, index) (*((name) + R8*((index)-1)))   /* integer array */
#define RFA(name, index) (*((name) + (index)-1))        /* Real array */

#define elsif else if
#define until(x) while(!(x))

#if 0
extern int fprintf(const FILE *, const char *, ...);
extern int printf(const char *, ...);
extern char *puts(char[]), fputs(char[], FILE*);
extern int fflush(FILE*);
extern int fscanf(const FILE *, const char *, ...);
extern int scanf(const char *, ...);
#endif

/* enum _boolEnum { false=0, true=1, maybe=2 };*/
/* TODO: CONFLICTING WITH GLOBAL ARRAY Boolean, renamed to boolean, but need to clear out all code except for heap */
typedef unsigned char boolean;

#define false (boolean)0
#define true  (boolean)1
#define maybe (boolean)2


typedef union _voidInt {
    int i;
    long l;
    long long ll;
    float f;
    double d;
    void *v;
} foint;

/* The comparison function type: used by heaps, binary trees and sorts.
*/
typedef int (*pCmpFcn)(foint, foint);   

/* Copy a foint.  In all instances, you are expected to know what the
** foint actually is, and return a copy of it.  If a FointCopy function
** pointer is ever NULL, the code will do a shallow copy.
*/
typedef foint (*pFointCopyFcn)(foint);   
typedef void (*pFointFreeFcn)(foint);

/* this is is the general error return value for most abstract data types */
extern const foint ABSTRACT_ERROR;

/*
** The following three functions print out warnings, apologies (ie,
** fatal limitations in the program that could be worked around with
** more effort on the part of the library creator), and fatal errors.
** The messages are always sent to stderr.  In addition, it is assumed
** that these messages should also *always* be printed to the terminal,
** so if neither stdout nor stderr are a terminal, then we open /dev/tty
** and put it there anyway.  None of these touch stdout.
*/
extern void Warning(const char *fmt, ...);
extern void Apology(const char *fmt, ...);
extern void Fatal(const char *fmt, ...);  /* generates an assertion failure */

extern void *Malloc(int);
extern void *Calloc(int, int);
extern void *Realloc(void *ptr, int newSize);
void *Memdup(void *v, int n);
#define Free free


/* return current user time used in seconds.  Returns -1 if error. */
double uTime(void);

char *Int2BitString(char word[33], unsigned i);

/*
** Be careful!  These may not work as you think if i is not unsigned.
** Rotate s bits of i by n positions and mask out anything above bit s.
*/
#if __sun__ /* gcc 2.6.3 bug */
unsigned long long _NotZero;
#define RotLeft(type,i,s,n) ((_NotZero = ~(type)0), ((((i) << (n)) | ((i) >> (s-(n)))) & (_NotZero >> (8*sizeof(i)-(s)))))
#define RotRight(type,i,s,n) ((_NotZero = ~(type)0), ((((i) >> (n)) | ((i) << (s-(n)))) & (_NotZero >> (8*sizeof(i)-(s)))))
#else
#define RotLeft(type,i,s,n) ((((i) << (n)) | ((i) >> (s-(n)))) & (~(type)0 >> (8*sizeof(i)-(s))))
#define RotRight(type,i,s,n) ((((i) >> (n)) | ((i) << (s-(n)))) & (~(type)0 >> (8*sizeof(i)-(s))))
#endif

void PrintArray(FILE *fp, int n, int *array);

long long IIntPow(int, int);
double IntPow(double, int);
int Log2(int);	/* integer part of log to the base 2 */

/* Greatest Common Divisor */
int gcd(int a, int b);
int IsPrime(long long n);
int PrimeFactors(int n, int count[]);

#endif  /* _MISC_H */