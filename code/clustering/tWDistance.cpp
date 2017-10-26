/*==========================================================
 * twDistance.cpp
 * The calling syntax is:
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2010 The MathWorks, Inc.
 *
 *========================================================*/
/* $Revision: 1.1.10.3 $ */

#include "mex.h"
#include <math.h>
#include <float.h>

int max(int a, int b) {return (a>b)?a:b;}
int min(int a, int b) {return (a<b)?a:b;}
double max(double a, double b) {return (a>b)?a:b;}
double min(double a, double b) {return (a<b)?a:b;}
//double fabs(double a) {return (a>0)?a:-a;}
int i2d(int a, int b, int outer) {return a + outer * b;}
int i3d(int a, int b, int c, int outer, int mid) {return a + outer * (b + mid * c);}

/* The computational routine */
void DPdist(const int nRows, double *A, double *B, const int lenA, const int lenB,
        double *dist, double *featureDists, double *matchIndex);

// inputs are two arrays/sequences normalized
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int nFeats;
    double *totalDist, *featureDists /* may be unneeded */;
    double *matchInterval;
    double *sampleA, *sampleB;
    size_t lenA, lenB;
    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two output required.");
    }
    
    /* check that number of rows in first two input arguments are equal */
    if(mxGetM(prhs[0])!=mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("twDistance:unequalRows","Inputs must have equal number of rows.");
    }
    
    /* create a pointer to the real data in the input matrix  */
    
    /* get dimensions of the input matrix */
    nFeats = (int) mxGetM(prhs[0]);
    
    /* assign samples according to size, A is smaller than B */
    
    if(mxGetN(prhs[0]) < mxGetN(prhs[1])) {
        sampleA = mxGetPr(prhs[0]);
        sampleB = mxGetPr(prhs[1]);
        lenA = mxGetN(prhs[0]);
        lenB = mxGetN(prhs[1]);
    } else {
        sampleA = mxGetPr(prhs[1]);
        sampleB = mxGetPr(prhs[0]);
        lenA = mxGetN(prhs[1]);
        lenB = mxGetN(prhs[0]);
    }
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,nFeats,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,lenB, mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    totalDist = mxGetPr(plhs[0]);
    featureDists = mxGetPr(plhs[1]);
    matchInterval = mxGetPr(plhs[2]);
    DPdist(nFeats, sampleA, sampleB, lenA, lenB, totalDist, featureDists, matchInterval);
}

void DPdist(const int nRows, double *A, double *B, const int lenA, const int lenB,
        double *dist, double *featureDists, double *matchIndex)
{
    const int dLen = abs(lenA - lenB);
    
    const int margin = 5;
    int i,j,k;
    enum Path {diagonal, onlyA, onlyB};
    double * staticCostMatrix = (double*) malloc(nRows * lenA * lenB * sizeof(double)); //[nRows][lenA][lenB];
    double * costOverFeatures = (double*) malloc(lenA * lenB * sizeof(double)); //[lenA][lenB];
    double * dynCostMatrix    = (double*) malloc(lenA * lenB * sizeof(double));  //[lenA][lenB];
    Path * bestCand = (Path*) malloc(lenA * lenB * sizeof(Path)); //[lenA][lenB];
    for(i = 0; i < lenA * lenB; i++) {
        costOverFeatures[i] = DBL_MAX;
        dynCostMatrix[i] = DBL_MAX;
    }
    
    const double warpingCost = 1.0;
    const bool L2 = true;
    bool breakOut = false;
    double cumul;
    double candA, candAB, candB, minCand;
    int revTravA, revTravB;
    
    for(i = 0; i < nRows; i++) featureDists[i] = 0;
    
    for(j = 0; j < lenA; j++) {
        for(k = max(0, j - margin); k < min(lenB, j + dLen + margin); k++) {
            cumul = 0;
            for(i = 0; i < nRows; i++)
            {
                if(L2)
                    staticCostMatrix[i3d(i,j,k, nRows, lenA)] =
                            (A[i2d(i, j, nRows)] - B[i2d(i, k, nRows)]) * (A[i2d(i, j, nRows)] - B[i2d(i, k, nRows)]);
                else
                    staticCostMatrix[i3d(i,j,k, nRows, lenA)] =
                            fabs(A[i2d(i, j, nRows)] - B[i2d(i, k, nRows)]);
                cumul += staticCostMatrix[i3d(i,j,k, nRows, lenA)];
            }
            costOverFeatures[i2d(j,k,lenA)] = cumul;
        }
    }

    dynCostMatrix[i2d(0,0,lenA)] = costOverFeatures[i2d(0,0,lenA)];
    bestCand[i2d(0,0,lenA)] = diagonal;
    
    // along the warping on short end
    if(margin >= 1) 
    {
        dynCostMatrix[i2d(1,0,lenA)] = dynCostMatrix[i2d(1,0,lenA)] +
                costOverFeatures[i2d(0,0,lenA)];
        bestCand[i2d(1,0,lenA)] = diagonal;
    }    

    for(i = 2; i < min(lenA, margin + 1); i++)  {
        dynCostMatrix[i2d(i,0,lenA)] = dynCostMatrix[i2d(i-1,0,lenA)] +
                costOverFeatures[i2d(i,0,lenA)];                
        bestCand[i2d(i,0,lenA)] = onlyA;
    }
    
    // along the warping on long end 
    for(j = 1; j < min(lenB, dLen + margin + 1); j++) {
        dynCostMatrix[i2d(0,j,lenA)] = dynCostMatrix[i2d(0,j-1,lenA)] +
                costOverFeatures[i2d(0,j,lenA)] + warpingCost;
        bestCand[i2d(0,j,lenA)] = onlyB;
    }
    // start the DP
    for(i = 1; i < lenA; i++) {
        for(j = max(1, i - (dLen + margin)); j < min(lenB, i + dLen + margin); j++)
        {
            candAB = dynCostMatrix[i2d(i-1,j-1,lenA)];
            candA  = dynCostMatrix[i2d(i-1,j  ,lenA)] + warpingCost;
            candB  = dynCostMatrix[i2d(i  ,j-1,lenA)] + warpingCost;
            
            if(candAB < candA && candAB < candB)
            {
                minCand = candAB;
                bestCand[i2d(i,j,lenA)] = diagonal;
            }
            else if(candA < candB)
            {
                minCand = candA;
                bestCand[i2d(i,j,lenA)] = onlyA;
            }
            else
            {
                minCand = candB;
                bestCand[i2d(i,j,lenA)] = onlyB;
            } //bestCand bad access traversal block
            dynCostMatrix[i2d(i,j,lenA)] = minCand + costOverFeatures[i2d(i,j,lenA)];
        }
    }
    if(L2)
        dist[0] = sqrt(dynCostMatrix[i2d(lenA-1,lenB-1,lenA)]);
    else
        dist[0] = dynCostMatrix[i2d(lenA-1,lenB-1,lenA)];
    
    // traverse backwards
    revTravA = lenA - 1; revTravB = lenB - 1;
    while(!(revTravA == 0 && revTravB == 0) && !breakOut)
    {
        for(i = 0; i < nRows; i++)
            featureDists[i] += staticCostMatrix[i3d(i,revTravA,revTravB,nRows, lenA)];
        switch(bestCand[i2d(revTravA,revTravB,lenA)]) {
            case diagonal: {
                matchIndex[revTravB] = revTravA;
                revTravA--; revTravB--;
                break;
            }
            case onlyA: {
                revTravA--;
                for(i = 0; i < nRows; i++) featureDists[i] += warpingCost/nRows;
                break;
            }
            case onlyB: {
                matchIndex[revTravB] = revTravA;
                revTravB--;
                for(i = 0; i < nRows; i++) featureDists[i] += warpingCost/nRows;
                break;
            }
            default: {mexPrintf("ERROR: out of bounds\n"); breakOut = true;}
        }
    }
    for(i = 0; i < revTravB; i++)
        
        if(L2) for(i = 0; i < nRows; i++) featureDists[i] = sqrt(featureDists[i]);
    
    free(staticCostMatrix);
    free(dynCostMatrix);
    free(costOverFeatures);
    free(bestCand);
    
}
