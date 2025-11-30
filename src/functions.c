/* This file contains source code adapted from the wavelet R package
   (version 4.7.3), written by Guy Nason and others.
   https://cran.r-project.org/package=wavethresh */

#include <R.h> 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

/* For boundary condition handling */
#define PERIODIC        1
#define SYMMETRIC       2

/* For the type of wavelet decomposition */
#define WAVELET     1   /* The standard decomposition */
#define STATION     2   /* The stationary decomposition */

/* Threshold types */
#define HARD    1
#define SOFT    2

/*
 * ACCESSC handles negative accesses, as well as those that exceed the number
 * of elements
 */

#define ACCESS(image, size, i, j)       *(image + (i)*(size) + (j))
#define ACCESSC(c, firstC, lengthC, ix, bc) *(c+reflect(((ix)-(firstC)),(lengthC),(bc)))
#define ACCESSD(l, i)   *(Data + (*LengthData*(l)) + (i))
#define POINTD(l,i) (Data + (*LengthData*(l)) + (i))
#define POINTC(l,i) (Carray +(*LengthData*(l)) + (i))

/*
 * The next three are exclusively for the stationary wavelet packet algorithm
 * WPST
 */
#define NPKTS(level, nlev)  (1 << (2*(nlev-level)))
#define PKTLENGTH(level)    (1 << level)

#define ACCWPST(a, level, avixstart, pkix, i) *((a) + *(avixstart+(level))+(pkix)*PKTLENGTH(level)+i)

/* Optimiser parameters */

#define R   0.61803399  /* The golden ratio for bisection searches */
#define Cons    (1.0-R)     /* For bisection searches          */

/* These next 3 are for the ipndacw code */
#define ACCESSW(w,j,k)  *(*(w+j)+k)
#define max(a,b)        ((a) > (b) ? (a) : (b))
#define min(a,b)        ((a) > (b) ? (b) : (a))

#define CEIL(i) ( ((i)>0) ? ( ((i)+1)/2):((i)/2) )

/* Note this is different to what is in wavepackst.c */
#define AVBPOINTD(w, l,i) (w + (nlevels*(i)) + (l))

struct complex  {
    double *realval;
    double *imagval;
    };

/* FUNCTION DECLARATIONS
 * --------------------- */
int reflect(int n, int lengthC, int bc);

void wavedecomp(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error);

void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc);

void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);

/* FUNCTION DEFINITIONS
 * -------------------- */
void destroycomplex(struct complex *a)
{
free((void *)a->realval);
free((void *)a->imagval);
free((void *)a);
}

/*
 * A C version of getpacket
 *
 * Warning. The argument list is the same as for the WaveThresh Splus version
 * except that nlevels here should be one more!
 */

double *getpacket(double *wst, int nlevels, int level, int index, int *error)
/* int nlevels::    This looks like it should be nlevels+1 for some reason */
{
register int i;
double *packet;
int PacketLength;

PacketLength = 1 << level;

if ((packet = (double *)malloc((unsigned)PacketLength*sizeof(double)))==NULL){
    *error = 3;
    return(NULL);
    }


for(i=0; i< PacketLength; ++i)
    *(packet+i) = *AVBPOINTD(wst, level, (index*PacketLength+i));

return(packet);
}
/* Rotate a vector */

/* Vector: a_1, a_2, a_3, ..., a_{n-1}, a_n

    becomes

       a_2, a_3, a_4, ..., a_n, a_1

   rotateback() does the opposite

*/

void rotater(double *book, int length)
{
register int i;
double tmp;

tmp = *book;

for(i=0; i<length-1; ++i)
        *(book+i) = *(book+i+1);

*(book+length-1) = tmp;
}

void rotateback(double *book, int length)
{
register int i;
double tmp;

tmp = *(book+length-1);

for(i= length-1; i>0; --i)
    *(book+i) = *(book+i-1);

*book = tmp;
}

/*
 * CONVOLVE -   Do filter H filter convolution with boundary
 */

void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc)
/*---------------------
 * Argument description
 *---------------------
double *c_in::    	Input data                       
int LengthCin::   	Length of this array                 
int firstCin::    	The first C value                    
double *H::   		Filter                       
int LengthH::     	Length of filter                 
double *c_out::   	Output data                      
int firstCout::   	First index of C array               
int lastCout::    	Last index of C array                
int type::		Type of wavelet decomposition            
int step_factor::	For stationary wavelets only             
int bc::      		Method of boundary correction PERIODIC, SYMMETRIC    
 *---------------------*/
{
double sum;
register int k;
register int count_out;
register int m;
register int cfactor;   /* This determines what sort of dilation we do  */
            /* and depends on the type argument     */
int reflect(int n, int lengthC, int bc);

count_out = 0;

switch(type)    {

    case WAVELET:   /*  Ordinary wavelets   */
            cfactor = 2;    /* Pick every other coefficient */
            break;

    case STATION:   /* Stationary wavelets  */
            cfactor = 1;    /* Pick every coefficient   */
            break;


    default:    /* This is an error, one of the above must have */
            /* been picked */
            /* However, this must be tested in a previous   */
            /* routine.                 */
            cfactor=0;       /* MAN: added for total cover: shouldn't happen */
            break;
        }

for(k=firstCout; k<=lastCout; ++k)  {
    sum = 0.0;


    for(m=0; m<LengthH; ++m)    {


        sum += *(H+m) * ACCESSC(c_in, firstCin, LengthCin,
            ((step_factor*m)+(cfactor*k)),bc);
        }

    *(c_out + count_out) = sum;
    ++count_out;
    }
}

void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc)
/*---------------------
 * Argument description
 *---------------------
double *c_in::    	Input data                       
int LengthCin::   	Length of this array                 
int firstCin::    	The first C value                    
double *H::   		Filter                       
int LengthH::     	Length of filter                 
double *d_out::   	Output data                      
int firstDout::   	First index of C array               
int lastDout::    	Last index of C array                
int type::		Type of wavelet decomposition            
int step_factor::	For stationary wavelets only             
int bc::      		Method of boundary correction PERIODIC, SYMMETRIC    
 *---------------------*/
{
double sum;
double tmp;
register int k;
register int count_out;
register int m;
register int cfactor;

int reflect(int n, int lengthC, int bc);

count_out = 0;

switch(type)    {

    case WAVELET:   /*  Ordinary wavelets   */
            cfactor = 2;    /* Pick every other coefficient */
            break;

    case STATION:   /* Stationary wavelets  */
            cfactor = 1;    /* Pick every coefficient   */
            break;


    default:    /* This is an error, one of the above must have */
            /* been picked */
            /* However, this must be tested in a previous   */
            /* routine.                 */
            cfactor=0;       /* MAN: added for total cover: shouldn't happen */
            break;
        }

for(k=firstDout; k<=lastDout; ++k)  {
    sum = 0.0;


    for(m=0; m<LengthH; ++m)    {

        tmp = ACCESSC(c_in, firstCin, LengthCin,
                (cfactor*k+(step_factor*(1-m))),bc);
        
        if (m&1)    /* odd */
            sum += *(H+m) *  tmp;
        else
            sum -= *(H+m) *  tmp;
        
        }

    *(d_out + count_out) = sum;
    ++count_out;
    }
}

/* Works out reflection, as REFLECT, but reports access errors */
int reflect(int n, int lengthC, int bc)
{

if ((n >= 0) && (n < lengthC))
    return(n);
else if (n<0)   {
    if (bc==PERIODIC)   {
        /*
        n = lengthC+n;
        */
        n = n%lengthC + lengthC*((n%lengthC)!=0);
        if (n < 0)      {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
            REprintf("reflect: left info from right\n");
	    error("This should not happen. Stopping.\n");
            }
        else
            return(n);
        }

    else if (bc==SYMMETRIC) {
        n = -1-n;
        if (n >= lengthC)       {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
	    error("This should not happen. Stopping.\n");
            }
        else
            return(n);
        }

    else    {
        REprintf("reflect: Unknown boundary correction");
        REprintf("value of %d\n", bc);
        error("This should not happen. Stopping.\n");
        }

    }
else    {
    if (bc==PERIODIC)   {
        /*
        Rprintf("periodic extension, was %d (%d) now ",n,lengthC);
        n = n - lengthC; 
        */
        n %= lengthC;
        /*
        Rprintf("%d\n", n);
        */
        if (n >= lengthC)   {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
            REprintf("reflect: right info from left\n");
	    error("This should not happen. Stopping.\n");
            }
        else
            return(n);
        }
    else if (bc==SYMMETRIC) {
        n = 2*lengthC - n - 1;
        if (n<0)        {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
	    error("This should not happen. Stopping.\n");
            }
        else
            return(n);
        }
    else    {
        REprintf("reflect: Unknown boundary correction\n");
	error("This should not happen. Stopping.\n");
        }


    }
/* Safety */
REprintf("reflect: SHOULD NOT HAVE REACHED THIS POINT\n");
error("This should not happen. Stopping.\n");
return(0); /* for lint only */
}

void wavedecomp(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error)
/*---------------------
 * Argument description
 *---------------------
 double *C::       Input data, and the subsequent smoothed data 
 double *D::       The wavelet coefficients                     
 double *H::       The smoothing filter H                       
 int *LengthH::    Length of smoothing filter                   
 int *levels::     The number of levels in this decomposition   
 int *firstC::     The first possible C coef at a given level   
 int *lastC::      The last possible C coef at a given level    
 int *offsetC::    Offset from C[0] for certain level's coeffs  
 int *firstD::     The first possible D coef at a given level   
 int *lastD::      The last possible D coef at a given level    
 int *offsetD::    Offset from D[0] for certain level's coeffs  
 int *type::       The type of wavelet decomposition        
 int *bc::         Method of boundary correction        
 int *error::      Error code                                   
 *---------------------*/
{
register int next_level,at_level;
register int step_factor;   /* Controls width of filter for station */
register int verbose;   /* Controls message printing, passed in error var*/

void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc);
void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);

if (*error == 1l)   /* Error switches on verbosity */
    verbose = 1;
else
    verbose = 0;

switch(*bc) {

    case PERIODIC:  /* Periodic boundary conditions */
        if (verbose) Rprintf("Periodic boundary method\n");
        break;

    case SYMMETRIC: /* Symmetric boundary conditions */
        if (verbose) Rprintf("Symmetric boundary method\n");
        break;

    default:    /* The bc must be one of the above */
        Rprintf("Unknown boundary correction method\n");
        *error = 1;
        return;
    }

switch(*type)   {

    case WAVELET:   /* Standard wavelets */
        if (verbose) Rprintf("Standard wavelet decomposition\n");
        break;

    case STATION:   /* Stationary wavelets */
        if (verbose) Rprintf("Stationary wavelet decomposition\n");
        break;

    default:    /* The type must be of one the above */
        if (verbose) Rprintf("Unknown decomposition type\n");
        *error = 2;
        return;
    }
        
if (verbose) Rprintf("Decomposing into level: ");

*error = 0;

step_factor = 1;    /* This variable should *always* be 1 for standard
             * wavelets. It should start at 1 for stationary
             * wavelets and multiply itself by 2 each stage
             */

for(next_level = *levels - 1; next_level >= 0; --next_level)    {

    if (verbose)
        Rprintf("%d ", next_level);

    at_level = next_level + 1;

/* For stationary wavelets we need to define a step factor.
 * This widens the span of the filter. At the top level (*levels->*levels-1)
 * it is one, as usual. Then for the next step it becomes 2, then 4 etc.
 */

    convolveC( (C+*(offsetC+at_level)),
        (int)(*(lastC+ at_level) - *(firstC+at_level)+1),
        (int)(*(firstC+at_level)),
        H,
        (int)*LengthH,
        (C+*(offsetC+next_level)),
        (int)(*(firstC+next_level)),
        (int)(*(lastC+next_level)) , (int)*type,
        step_factor, (int)*bc);

    convolveD( (C+*(offsetC+at_level)),
                (int)(*(lastC+ at_level) - *(firstC+at_level)+1),
                (int)(*(firstC+at_level)),
                H,
                (int)*LengthH,
        (D+*(offsetD+next_level)),
        (int)(*(firstD+next_level)),
        (int)(*(lastD+next_level)), (int)*type,
        step_factor, (int)*bc );

    if (*type == STATION)
        step_factor *= 2;   /* Any half decent compiler should
                     * know what to do here ! */
    }
if (verbose)
    Rprintf("\n");
return;
}

/*
 * Functions to do complex arithmetic
 *
 */

/*
 * Addition: a+ib + c+id = a+c +i(b+d) = e + i f
 */

void comadd(double a, double b, double c, double d, double *e, double *f)
{
*e = a+c;
*f = b+d;
}

/*
 * Subtraction: a+ib - c+id = a+c -i(b+d) = e + i f
 */

void comsub(double a, double b, double c, double d, double *e, double *f)
{
*e = a-c;
*f = b-d; 
}

/*
 * Multiplication: (a+ib)(c+id) = ac-bd +i(bc+ad) = e + i f
 */

void commul(double a, double b, double c, double d, double *e, double *f)
{
*e = (a*c - b*d);
*f = (b*c + a*d);
}


/*
 * Division: (a+ib)(c+id) = (ac+bd +i(bc-ad))/(c^2+d^2) = e + i f
 */

void comdiv(double a, double b, double c, double d, double *e, double *f)
{
double tmp;

tmp = c*c + d*d;

*e = (a*c + b*d)/tmp;
*f = (b*c - a*d)/tmp;
}

/*
 * This routine is identical to the convolve.c routine except it
 * does it for complex wavelets.

 * COMCONC  -   Do filter H filter convolution with boundary
 */


void comconC(double *c_inR, double *c_inI,
	int LengthCin, int firstCin,
	double *HR, double *HI, int LengthH,
	double *c_outR, double *c_outI,
	int LengthCout, int firstCout, int lastCout,
	int type, int step_factor, int bc)
/*---------------------
 * Argument description
 *---------------------
double *c_inR::   Input data (real)                    
double *c_inI::   Input data (imaginary)               
int LengthCin::   Length of this array                 
int firstCin::     <-- MAN: added since missing...     
double *HR::      Lowpass Filter                   
double *HI::      Lowpass Filter                   
int LengthH::     Length of filter                 
double *c_outR::  Output data (real)                   
double *c_outI::  Output data (imaginary)              
int LengthCout::  Length of above array                
int firstCout::   First index of C array               
int lastCout::    Last index of C array                
int type::        Type of wavelet decomposition            
int step_factor:: For stationary wavelets only             
int bc::          Method of boundary correction PERIODIC, SYMMETRIC    
 *---------------------*/
{
double sumR,sumI;
double a,b,c,d,e,f;
register int k;
register int count_out;
register int m;
register int cfactor;   /* This determines what sort of dilation we do  */
            /* and depends on the type argument     */

count_out = 0;

switch(type)    {

    case WAVELET:   /*  Ordinary wavelets   */
            cfactor = 2;    /* Pick every other coefficient */
            break;

    case STATION:   /* Stationary wavelets  */
            cfactor = 1;    /* Pick every coefficient   */
            break;


    default:    /* This is an error, one of the above must have */
            /* been picked */
            /* However, this must be tested in a previous   */
            /* routine.                 */
            cfactor=0;       /* MAN: added for total cover: shouldn't happen */
            break;
        }

for(k=firstCout; k<=lastCout; ++k)  {
    sumR = 0.0;
    sumI = 0.0;


    for(m=0; m<LengthH; ++m)    {

        a = *(HR + m);  /* real part */
        b = *(HI + m);  /* imaginary part */

        c = ACCESSC(c_inR, firstCin, LengthCin,
            ((step_factor*m)+(cfactor*k)),bc);
        d = ACCESSC(c_inI, firstCin, LengthCin,
            ((step_factor*m)+(cfactor*k)),bc);

        commul(a,b,c,d,&e, &f);

        sumR += e;
        sumI += f;
        }

    *(c_outR + count_out) = sumR;
    *(c_outI + count_out) = sumI;
    ++count_out;
    }
}

void comconD(double *c_inR, double *c_inI,
	int LengthCin, int firstCin,
	double *GR, double *GI, int LengthH,
	double *d_outR, double *d_outI,
	int LengthDout, int firstDout, int lastDout,
	int type, int step_factor, int bc)
/*---------------------
 * Argument description
 *---------------------
double *c_inR::   Input data (real)                    
double *c_inI::   Input data (imaginary)               
int LengthCin::   Length of this array                 
int firstCin::     <-- MAN: added since missing...     
double *GR::      Lowpass Filter                   
double *GI::      Lowpass Filter                   
int LengthH::     Length of filter                 
double *d_outR::  Output data (real)                   
double *d_outI::  Output data (imaginary)              
int LengthDout::  Length of above array                
int firstDout::   First index of C array               
int lastDout::    Last index of C array                
int type::        Type of wavelet decomposition            
int step_factor:: For stationary wavelets only             
int bc::          Method of boundary correction PERIODIC, SYMMETRIC    
 *---------------------*/
{
double sumR, sumI;
double a,b,c,d,e,f;
register int k;
register int count_out;
register int m;
register int cfactor;

count_out = 0;

switch(type)    {

    case WAVELET:   /*  Ordinary wavelets   */
            cfactor = 2;    /* Pick every other coefficient */
            break;

    case STATION:   /* Stationary wavelets  */
            cfactor = 1;    /* Pick every coefficient   */
            break;


    default:    /* This is an error, one of the above must have */
            /* been picked */
            /* However, this must be tested in a previous   */
            /* routine.                 */
            cfactor=0;       /* MAN: added for total cover: shouldn't happen */
            break;
        }

for(k=firstDout; k<=lastDout; ++k)  {
    sumR = 0.0;
    sumI = 0.0;


    for(m=0; m<LengthH; ++m)    {

        a = *(GR+m);

        b = *(GI+m);

        c = ACCESSC(c_inR, firstCin, LengthCin,
            ((step_factor*m)+(cfactor*k)),bc);
        d = ACCESSC(c_inI, firstCin, LengthCin,
            ((step_factor*m)+(cfactor*k)),bc);
        /*
        c = ACCESSC(c_inR, firstCin, LengthCin,
                (cfactor*k+(step_factor*(1-m))),bc);
        d = ACCESSC(c_inI, firstCin, LengthCin,
                (cfactor*k+(step_factor*(1-m))),bc);

        Rprintf("%d: (%lf, %lf)* (%lf, %lf)\n", a,b,c,d);
        */
        commul(a,b,c,d,&e,&f);

        sumR += e;
        sumI += f;
        }

    *(d_outR + count_out) = sumR;
    *(d_outI + count_out) = sumI;
    ++count_out;
    }
}

/*
 * Complex version of wavelet transform
 */

void comwd(double *CR, double *CI, int *LengthC,
	   double *DR, double *DI, int *LengthD,
	double *HR, double *HI, double *GR, double *GI, int *LengthH,
	int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error)
/*---------------------
 * Argument description
 *---------------------
double *CR::             Input data, and the subsequent smoothed data 
double *CI::             Input data, and the subsequent smoothed data 
int *LengthC::           Length of C array                            
double *DR::             The wavelet coefficients                     
double *DI::             The wavelet coefficients                     
int *LengthD::           Length of D array                            
double *HR::             The smoothing filter H                       
double *HI::             The smoothing filter H                       
double *GR::             The highpass filter H                       
double *GI::             The highpass filter H                       
int *LengthH::           Length of smoothing filter                   
int *levels::            The number of levels in this decomposition   
int *firstC::            The first possible C coef at a given level   
int *lastC::             The last possible C coef at a given level    
int *offsetC::           Offset from C[0] for certain level's coeffs  
int *firstD::            The first possible D coef at a given level   
int *lastD::             The last possible D coef at a given level    
int *offsetD::           Offset from D[0] for certain level's coeffs  
int *type::       	 The type of wavelet decomposition        
int *bc::        	 Method of boundary correction        
int *error::             Error code                                   
 *---------------------*/
{
register int next_level,at_level;
register int step_factor;   /* Controls width of filter for station */
register int verbose;   /* Controls message printing, passed in error var*/

if (*error == 1)   /* Error switches on verbosity */
    verbose = 1;
else
    verbose = 0;

switch(*bc) {

    case PERIODIC:  /* Periodic boundary conditions */
        if (verbose) Rprintf("Periodic boundary method\n");
        break;

    case SYMMETRIC: /* Symmetric boundary conditions */
        if (verbose) Rprintf("Symmetric boundary method\n");
        break;

    default:    /* The bc must be one of the above */
        Rprintf("Unknown boundary correction method\n");
        *error = 1;
        return;
        break;
    }

switch(*type)   {

    case WAVELET:   /* Standard wavelets */
        if (verbose) Rprintf("Standard wavelet decomposition\n");
        break;

    case STATION:   /* Stationary wavelets */
        if (verbose) Rprintf("Stationary wavelet decomposition\n");
        break;

    default:    /* The type must be of one the above */
        if (verbose) Rprintf("Unknown decomposition type\n");
        *error = 2;
        return;
        break;
    }
        
if (verbose) Rprintf("Decomposing into level: ");

*error = 0;

step_factor = 1;    /* This variable should *always* be 1 for standard
             * wavelets. It should start at 1 for stationary
             * wavelets and multiply itself by 2 each stage
             */

for(next_level = *levels - 1; next_level >= 0; --next_level)    {

    if (verbose)
        Rprintf("%d ", next_level);

    at_level = next_level + 1;

/* For stationary wavelets we need to define a step factor.
 * This widens the span of the filter. At the top level (*levels->*levels-1)
 * it is one, as usual. Then for the next step it becomes 2, then 4 etc.
 */

    comconC( (CR+*(offsetC+at_level)),
           (CI+*(offsetC+at_level)),
        (int)(*(lastC+ at_level) - *(firstC+at_level)+1),
        (int)(*(firstC+at_level)),
        HR, HI,
        (int)*LengthH,
        (CR+*(offsetC+next_level)),
        (CI+*(offsetC+next_level)),
        (int)(*(lastC+next_level) - *(firstC+next_level)+1),
        (int)(*(firstC+next_level)),
        (int)(*(lastC+next_level)) , (int)*type,
        step_factor, (int)*bc);

    comconD( (CR+*(offsetC+at_level)),
           (CI+*(offsetC+at_level)),
                (int)(*(lastC+ at_level) - *(firstC+at_level)+1),
                (int)(*(firstC+at_level)),
                GR, GI,
                (int)*LengthH,
        (DR+*(offsetD+next_level)),
        (DI+*(offsetD+next_level)),
        (int)(*(lastD+next_level) - *(lastD+next_level)+1),
        (int)(*(firstD+next_level)),
        (int)(*(lastD+next_level)), (int)*type,
        step_factor, (int)*bc );

    if (*type == STATION)
        step_factor *= 2;   /* Any half decent compiler should
                     * know what to do here ! */
    }
if (verbose)
    Rprintf("\n");
return;
}

/*
 * COMCBR: Does the reconstruction convolution
 */

void comcbr(double *c_inR, double *c_inI,
	int LengthCin, int firstCin, int lastCin,
	double *d_inR, double *d_inI,
	int LengthDin, int firstDin, int lastDin,
	double *HR, double *HI, double *GR, double *GI, int LengthH,
	double *c_outR, double *c_outI, int LengthCout, int firstCout,
	int lastCout, int type, int bc)
{
register int n,k;
register int cfactor;
double sumCR, sumCI, sumDR, sumDI;
double a,b,c,d,e,f;

switch(type)    {

    case WAVELET:   /* Standard wavelets */
        cfactor = 2;
        break;

    case STATION:   /* Stationary wavelets */
        cfactor = 1;
        break;

    default:    /* This should never happen */
        cfactor=0;       /* MAN: added for total cover: shouldn't happen */
        break;
    }


/* Compute each of the output C */

for(n=firstCout; n<=lastCout; ++n)  {

    /* We want  n+1-LengthH <= 2*k to start off */


    k = CEIL(n+1-LengthH);

    sumCR = 0.0;
    sumCI = 0.0;
    sumDR = 0.0;
    sumDI = 0.0;

    while( cfactor*k <= n ) {

        a = *(HR + n - cfactor*k);
        b = *(HI + n - cfactor*k);

        c = ACCESSC(c_inR, firstCin, LengthCin, k, bc);
        d = ACCESSC(c_inI, firstCin, LengthCin, k, bc);

        commul(a,b,c,d, &e, &f);

        sumCR += e;
        sumCI += f;

        /* Now D part */

        a = *(GR + n - cfactor*k);
        b = *(GI + n - cfactor*k);

        c = ACCESSC(d_inR, firstDin, LengthDin, k, bc);
        d = ACCESSC(d_inI, firstDin, LengthDin, k, bc);

        commul(a,b,c,d, &e, &f);

        sumDR += e;
        sumDI += f;

        ++k;
        }

    sumCR += sumDR;
    sumCI += sumDI;

    ACCESSC(c_outR, firstCout, LengthCout, n, bc) = sumCR;
    ACCESSC(c_outI, firstCout, LengthCout, n, bc) = sumCI;
    }

}

/* comAB Do the basis averaging for complex WST*/

/*
 * Error codes
 *
 * 1,2  -   Memory error in creating clR, clI
 * 3,4  -   Memory error in creating crR, crI
 * 3    -   Memory error in creating packet (getpacket)
 */


struct complex *comAB(double *wstR, double *wstI, double *wstCR, double *wstCI,
    int nlevels, int level, int ix1, int ix2,
    double *HR, double *HI, double *GR, double *GI, int LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *wstR::    Wavelet coefficients, non-dec, real            
double *wstI::    Wavelet coefficients, non-dec, imag            
double *wstCR::   Father wav. coeffs, non-dec, real            
double *wstCI::   Father wav. coeffs, non-dec, imag            
int nlevels::     The original length of the data         
int level::       The level to reconstruct             
int ix1::         The "left" packet index              
int ix2::         The "right" packet index             
double *HR,*HI::  Smoothing filter                 
double *GR,*GI::  Detail filter                    
int LengthH::     The length of the filter             
int *error::      Error code                       
 *---------------------*/
{
register int i;
double *clR, *clI;
double *crR, *crI;
struct complex *genericC;
struct complex *answer;
double *genCR, *genCI;  /* Generic Cs for when we need real and imag */
double *genDR, *genDI;  /* Generic Cs for when we need real and imag */
int LengthC;
int LengthCin;

void comcbr(double *c_inR, double *c_inI,
	int LengthCin, int firstCin, int lastCin,
	double *d_inR, double *d_inI,
	int LengthDin, int firstDin, int lastDin,
	double *HR, double *HI, double *GR, double *GI, int LengthH,
	double *c_outR, double *c_outI, int LengthCout, int firstCout,
	int lastCout, int type, int bc);
double *getpacket(double *wst, int nlevels, int level, int index, int *error);
struct complex *comAB(double *wstR, double *wstI, double *wstCR, double *wstCI,
    int nlevels, int level, int ix1, int ix2,
    double *HR, double *HI, double *GR, double *GI, int LengthH, int *error);
void rotateback(double *book, int length);
void destroycomplex(struct complex *a);

*error = 0;

/*
 * Now we must create cl and cr. These will contain the reconstructions
 * from the left and right packets respectively. The length of these
 * vectors depends upon the level we're at.
 */

LengthC = 1 << (level+1);
LengthCin = 1 << level;

/*
 * Create cl and cr: real and imaginary
 */

if ((clR = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 1;
    return(NULL);
    }

if ((clI = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 2;
    return(NULL);
    }

if ((crR = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 3;
    return(NULL);
    }

if ((crI = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 4;
    return(NULL);
    }

/*
 * What we do next depends on the level.
 *
 * If level is zero then we've recursed all the way down to the bottom of
 * the tree. And we can reconstruct the 2-vectors one-up-the-tree by using
 * good old conbar().
 *
 * If the level is not zero then we construct at that stage using conbar()
 * but to obtain the Cs we recurse. 
 */

if (level != 0) {

    /* Get C's at this level by asking the next level down. */

    genericC = comAB(wstR, wstI, wstCR, wstCI,
            nlevels, level-1, 2*ix1, 2*ix1+1,
            HR, HI, GR, GI, LengthH, error);

    if (*error != 0)
        return(NULL); 

    /* Get D's straight from the wst matrix */

    genDR = getpacket(wstR, nlevels, level, ix1, error);
    genDI = getpacket(wstI, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    comcbr(genericC->realval, genericC->imagval, LengthCin, 0, LengthCin-1, 
           genDR, genDI, LengthCin, 0, LengthCin-1,
           HR, HI, GR, GI, LengthH,
           clR, clI, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    destroycomplex(genericC);
    free((void *)genDR);
    free((void *)genDI);

    /* Now do the RHS */
    
    genericC = comAB(wstR, wstI, wstCR, wstCI, nlevels, level-1,
        2*ix2, 2*ix2+1,
        HR, HI, GR, GI, LengthH, error);

    if (*error != 0)
        return(NULL); 

    /* Get D's straight from the wst matrix */

    genDR = getpacket(wstR, nlevels, level, ix2, error);
    genDI = getpacket(wstI, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    comcbr(genericC->realval, genericC->imagval, LengthCin, 0, LengthCin-1,
           genDR, genDI, LengthCin, 0, LengthCin-1,
           HR, HI, GR, GI, LengthH,
           crR, crI, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    /* Rotate the RHS back */

    rotateback(crR, LengthC);
    rotateback(crI, LengthC);

    /* Can get rid of generics now */

    destroycomplex(genericC);
    free((void *)genDR);
    free((void *)genDI);
    }

else    {
    /* Have to really do it! */

    genCR = getpacket(wstCR, nlevels, level, ix1, error);
    genCI = getpacket(wstCI, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    genDR = getpacket(wstR, nlevels, level, ix1, error);
    genDI = getpacket(wstI, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    comcbr(genCR, genCI, LengthCin, 0, LengthCin-1, 
           genDR, genDI, LengthCin, 0, LengthCin-1,
           HR, HI, GR, GI, LengthH,
           clR, clI, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    free((void *)genCR);
    free((void *)genCI);
    free((void *)genDR);
    free((void *)genDI);

    genCR = getpacket(wstCR, nlevels, level, ix2, error);
    genCI = getpacket(wstCI, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    genDR = getpacket(wstR, nlevels, level, ix2, error);
    genDI = getpacket(wstI, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    comcbr(genCR, genCI, LengthCin, 0, LengthCin-1,
           genDR, genDI, LengthCin, 0, LengthCin-1,
           HR, HI, GR, GI, LengthH,
           crR, crI, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    /* Rotate the RHS back */

    rotateback(crR, LengthC);
    rotateback(crI, LengthC);

    free((void *)genCR);
    free((void *)genCI);
    free((void *)genDR);
    free((void *)genDI);
    }

for(i=0; i<LengthC; ++i)    {
    *(clR+i) = ((double)0.5)*( *(clR+i) + *(crR+i) );
    *(clI+i) = ((double)0.5)*( *(clI+i) + *(crI+i) );
    }

if ((answer=(struct complex *)malloc((unsigned)sizeof(struct complex)))==NULL) {
    *error = 5l;
    return(NULL);
    }

answer->realval = clR;
answer->imagval = clI;

return(answer);
}

void comAB_WRAP(double *wstR, double *wstI, double *wstCR, double *wstCI,
    int *LengthData, int *level,
    double *HR, double *HI, double *GR, double *GI, int *LengthH,
    double *answerR, double *answerI, int *error)
/*---------------------
 * Argument description
 *---------------------
double *wstR::        Wavelet coefficients - real          
double *wstI::        Wavelet coefficients - imag          
double *wstCR::       Father coeffs - real               
double *wstCI::       Father coeffs - imag             
int *LengthData::
int *level::
double *HR, *HI::     Smoothing filter, real and imag      
double *GR, *GI::     Detail filter, real and imag         
int *LengthH::
double *answerR, *answerI::   Real and imag of answer      
int *error::
 *---------------------*/
{
register int i;
int nlevels;
struct complex *acopy;
struct complex *comAB(double *wstR, double *wstI, double *wstCR, double *wstCI,
    int nlevels, int level, int ix1, int ix2,
    double *HR, double *HI, double *GR, double *GI, int LengthH, int *error);
void destroycomplex(struct complex *a);

nlevels = 2 + (int)*level;

acopy =  comAB(wstR, wstI, wstCR, wstCI, nlevels, (int)*level, 0, 1,
        HR, HI, GR, GI, (int)*LengthH, error);

for(i=0; i< (int)*LengthData; ++i)  {
    *(answerR+i) = *(acopy->realval+i);
    *(answerI+i) = *(acopy->imagval+i);
    }

destroycomplex(acopy);
}

/*
 * CONBAR: Does the reconstruction convolution
 */

#define CEIL(i) ( ((i)>0) ? ( ((i)+1)/2):((i)/2) )

void conbar(double *c_in, int LengthCin, int firstCin,
	double *d_in, int LengthDin, int firstDin,
	double *H, int LengthH,
	double *c_out, int LengthCout, int firstCout, int lastCout,
	int type, int bc)
{
register int n,k;
register int cfactor;
double sumC, sumD;

int reflect(int n, int lengthC, int bc);

switch(type)    {

    case WAVELET:   /* Standard wavelets */
        cfactor = 2;
        break;

    case STATION:   /* Stationary wavelets */
        cfactor = 1;
        break;

    default:    /* This should never happen */
        cfactor=0;       /* MAN: added for total cover: shouldn't happen */
        break;
    }


/* Compute each of the output C */

for(n=firstCout; n<=lastCout; ++n)  {

    /* We want  n+1-LengthH <= 2*k to start off */


    k = CEIL(n+1-LengthH);

    sumC = 0.0;

    while( cfactor*k <= n ) {

        sumC += *(H + n - cfactor*k)*ACCESSC(c_in, firstCin, LengthCin,
                    k, bc);

        ++k;
        }

    /* Now do D part */

    k = CEIL(n-1);

    sumD = 0.0;

    while( cfactor*k <= (LengthH +n -2) )   {

        sumD += *(H+1+cfactor*k-n) * ACCESSC(d_in, firstDin, LengthDin,
                    k, bc);

        ++k;

        }

    if (n & 1)      /* n odd */
        sumC -= sumD;
    else
        sumC += sumD;

    ACCESSC(c_out, firstCout, LengthCout, n, bc) = sumC;
    }

}
/*
 * CONBARL: Wrapper called by SPlus conbar() to call C conbar.
 */

void conbarL(double *c_in, int *LengthCin, int *firstCin,
	double *d_in, int *LengthDin, int *firstDin,
	double *H, int *LengthH,
	double *c_out, int *LengthCout, int *firstCout, int *lastCout,
	int *type, int *bc)
{
int LLengthCin;
int LfirstCin;
int LLengthDin;
int LfirstDin;
int LLengthH;
int LLengthCout;
int LfirstCout;
int LlastCout;
int Ltype;
int Lbc;
void conbar(double *c_in, int LengthCin, int firstCin,
	double *d_in, int LengthDin, int firstDin,
	double *H, int LengthH,
	double *c_out, int LengthCout, int firstCout, int lastCout,
	int type, int bc);

LLengthCin = (int)*LengthCin;
LfirstCin = (int)*firstCin;
LLengthDin = (int)*LengthDin;
LfirstDin = (int)*firstDin;
LLengthH = (int)*LengthH;
LLengthCout = (int)*LengthCout;
LfirstCout = (int)*firstCout;
LlastCout = (int)*lastCout;
Ltype = (int)*type;
Lbc = (int)*bc;


conbar(c_in, LLengthCin, LfirstCin,
       d_in, LLengthDin, LfirstDin,
       H, LLengthH,
       c_out, LLengthCout, LfirstCout, LlastCout, Ltype, Lbc);
}

/* AV_BASIS Do the basis averaging */

/*
 * Error codes
 *
 * 1    -   Memory error in creating cl
 * 2    -   Memory error in creating cr
 * 3    -   Memory error in creating packet (getpacket)
 */


double *av_basis(double *wst, double *wstC, int nlevels, int level,
	int ix1, int ix2, double *H, int LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *wst::	The stationary wavelet decomposition         
double *wstC::  The stationary wavelet decomposition         
int nlevels::  	The original length of the data         
int level::     The level to reconstruct             
int ix1::       The "left" packet index              
int ix2::       The "right" packet index             
double *H::     The filter                       
int LengthH::   The length of the filter             
int *error::    Error code                       
 *---------------------*/
{
register int i;
double *cl;
double *cr;
double *genericC;
double *genericD;
void conbar(double *c_in, int LengthCin, int firstCin,
	double *d_in, int LengthDin, int firstDin,
	double *H, int LengthH,
	double *c_out, int LengthCout, int firstCout, int lastCout,
	int type, int bc);
double *getpacket(double *wst, int nlevels, int level, int index, int *error);
double *av_basis(double *wst, double *wstC, int nlevels, int level,
	int ix1, int ix2, double *H, int LengthH, int *error);
int LengthC;
int LengthCin;

void rotateback(double *book, int length);

*error = 0;

/*
 * Now we must create cl and cr. These will contain the reconstructions
 * from the left and right packets respectively. The length of these
 * vectors depends upon the level we're at.
 */

LengthC = 1 << (level+1);
LengthCin = 1 << level;

/*
 * Create cl and cr
 */

if ((cl = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 1;
    return(NULL);
    }

if ((cr = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 2;
    return(NULL);
    }
    
/*
 * What we do next depends on the level.
 *
 * If level is zero then we've recursed all the way down to the bottom of
 * the tree. And we can reconstruct the 2-vectors one-up-the-tree by using
 * good old conbar().
 *
 * If the level is not zero then we construct at that stage using conbar()
 * but to obtain the Cs we recurse. 
 */

if (level != 0) {

    /* Get C's at this level by asking the next level down. */

    genericC = av_basis(wst, wstC, nlevels, level-1, 2*ix1, 2*ix1+1,
            H, LengthH, error);

    if (*error != 0)
        return(NULL); 

    /* Get D's straight from the wst matrix */

    genericD = getpacket(wst, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    conbar(genericC, LengthCin, 0, 
           genericD, LengthCin, 0, 
           H, LengthH,
           cl, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    free((void *)genericC);
    free((void *)genericD);

    /* Now do the RHS */
    
    genericC = av_basis(wst, wstC, nlevels, level-1, 2*ix2, 2*ix2+1,
        H, LengthH, error);

    if (*error != 0)
        return(NULL); 

    /* Get D's straight from the wst matrix */

    genericD = getpacket(wst, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    conbar(genericC, LengthCin, 0, 
           genericD, LengthCin, 0,
           H, LengthH,
           cr, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    /* Rotate the RHS back */

    rotateback(cr, LengthC);

    /* Can get rid of generics now */

    free((void *)genericC);
    free((void *)genericD);
    }

else    {
    /* Have to really do it! */

    genericC = getpacket(wstC, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    genericD = getpacket(wst, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    conbar(genericC, LengthCin, 0, 
           genericD, LengthCin, 0, 
           H, LengthH,
           cl, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    free((void *)genericC);
    free((void *)genericD);

    genericC = getpacket(wstC, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    genericD = getpacket(wst, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    conbar(genericC, LengthCin, 0,
           genericD, LengthCin, 0,
           H, LengthH,
           cr, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    /* Rotate the RHS back */

    rotateback(cr, LengthC);

    free((void *)genericC);
    free((void *)genericD);
    }

for(i=0; i<LengthC; ++i)
    *(cl+i) = ((double)0.5)*( *(cl+i) + *(cr+i) );

/*
 *  Return the answer in cl (which has to be freed later)
 *  Destroy pointer to cr, as it is not needed now
 *
 */

free((void *)cr);
return(cl);
}

/* Wrapper for av_basis */
void av_basisWRAP(double *wst, double *wstC, int *LengthData, int *level,
	double *H, int *LengthH, double *answer, int *error)
{
register int i;
int nlevels;
double *acopy;
double *av_basis(double *wst, double *wstC, int nlevels, int level,
	int ix1, int ix2, double *H, int LengthH, int *error);

nlevels = 2 + (int)*level;

acopy =  av_basis(wst, wstC, nlevels, (int)*level, 0, 1, H,
        (int)*LengthH, error);

for(i=0; i< (int)*LengthData; ++i)
    *(answer+i) = *(acopy+i);

free((void *)acopy);
}
