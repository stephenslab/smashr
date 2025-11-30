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
