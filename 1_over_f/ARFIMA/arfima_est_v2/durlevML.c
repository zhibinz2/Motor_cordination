#include "mex.h" 
#include "matrix.h"
#define square(p) ((p)*(p))
#define inv(q) (1/(q))
/* Durbin-Levinson algorithm for linear stationary AR(FI)MA(p,d,q) processes 
   Slightly altered for the maximum likelihood estimation  
   (C) György Inzelt 2011                                                    */
void levinson_recur1(double* v,double* L, int N, double* gammas,int step)
{
 int i,k;   
 
 if(step==0)
 {
   *(v + step) = *(gammas + step);
   *(L + step) = 1;
   for(k = step+1;k < N;++k)
   {  
     *(L + k) = 0;
   }
 }
 else if(step > 0 && step < N)
 {
    //phi_tt 
    *(L + step*N ) =  (-1.00)* *(gammas + step);
    if(step > 1)
    {
        for(i = 1;i < step ;++i)
        {
          *(L + step*N) -=  *(L + (step-1)*N + (step -1)  - i ) * *(gammas + step - i)   ;
        }
    }
     *(L +step*N)  *= inv( *(v + step-1) );
    //v_t
    *(v + step) = *(v + step-1)*(1- square( *(L + step*N) ));
    //phi_tj
    for(i =1; i < step; ++i)
    {
     *(L + step*N + step - i) =  *(L + (step-1)*N + (step -1) - i) + *(L  + step*N  ) * *(L + (step-1)*N + i -1 ) ;
    }
    //filling L with zeros and ones
    *(L + step*N + step ) = 1; 
    if(step != N-1)
    {
        for(k = step*N +step+1 ;k < step*N + N ;++k)
        {  
        *(L + k) =0;
        } 
    }
 }
   if(step < N-1)
     levinson_recur1(v,L,N,gammas,++step);
}
/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
int step=0;
int N;
double *gammas,*v,*L;
// getting the autocovariances 
gammas = mxGetPr(prhs[0]);
N = mxGetM(prhs[0]);
// v
plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
mxSetM(plhs[0],N); 
mxSetN(plhs[0],1);
mxSetData(plhs[0], mxMalloc(sizeof(double)*N*1));
// L
plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
mxSetM(plhs[1],square(N)); 
mxSetN(plhs[1],1);
mxSetData(plhs[1], mxMalloc(sizeof(double)*square(N)*1));
//
v = mxGetPr(plhs[0]);
L = mxGetPr(plhs[1]);
//
levinson_recur1(v, L, N,gammas,step);
//
return;
}