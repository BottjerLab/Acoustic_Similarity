/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier)
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2010 The MathWorks, Inc.
 *
 *========================================================*/
/* $Revision: 1.1.10.3 $ */

#include "mex.h"

/* The computational routine */
double DPdist(const int nRows, double *A, double *B, int & *matchIndex)
{
    // TODO: enforce that A is shorter than B
    const int lenA = sizeof(A)/sizeof(double)/nRows,
            lenB = sizeof(B)/sizeof(double)/nRows,
            dLen = abs(lenA - lenB);
    
    const int margin = 5;
    int i,j,k;
    double staticCostMatrix[nRows][lenA][lenB];
    double costOverFeatures[lenA][lenB];
    double dynCostMatrix[lenA][lenB];
    const double warpingCost = 1.0;
    bool L2 = true;
    double cum;
    int bestCand[lenA][lenB];
    double candA, candAB, candB, minCand;
    int revTravA = lenA, revTravB = lenB, step;
    
    for(j = 0; j < lenA; j++)
        for(k = max(0,j- margin); k < min(nRows, j+dLen + margin); k++) {
        cum = 0;
        for(i = 0; i < nRows; i++)
        {
            if(L2)
                staticCostMatrix[i][j][k] = (A[i][j] - B[i][k]) * (A[i][j] - B[i][k]);
            else
                staticCostMatrix[i][j][k] = abs(A[i][j] - B[i][k]);
            cum += staticCostMatrix[i][j][k];
        }
        costOverFeatures[j][k] = cum;
        }
    
    dynCostMatrix[0][0] = costOverFeatures[0][0]; bestCand[0][0] = 0;
    for(j = 1; j <= dLen; j++)  {
        dynCostMatrix[j][0] = costOverFeatures[j-1][0] + warpingCost;
        bestCand[j][0] = 0;
    }
    for(j = 1; j <= dLen + margin; j++) {
        dynCostMatrix[0][j] = costOverFeatures[0][j-1] + warpingCost;
        bestCand[0][j] = 0;
    }
    
    // start the DP
    for(i = 1; i < lenA; i++) {
        for(j = max(1, i - (dLen + margin)); j < min(lenB, i + dLen); j++)
        {
            candAB = dynCostMatrix[i-1][j-1];
            candA = dynCostMatrix[i-1][j] + warpingCost;
            candB = dynCostMatrix[i][j-1] + warpingCost;
            
            if(candAB < candA && candAB < candB)
            {
                minCand = candAB;
                bestCand[i][j] = 0;
            }
            else if(candA < candB)
            {
                minCand = candA;
                bestCand[i][j] = 1;
            }
            else
            {
                minCand = candB;
                bestCand[i][j] = 2;
            }
            dynCostMatrix[i][j] = minCand + costOverFeatures[i][j];
        }
        
    }
    dist = dynCostMatrix[lenA-1][lenB-1];
    if(L2)
        dist = sqrt(dist);
    
    matchIndex = malloc(sizeof(int) * lenA);
    while(revTravA > 0 && revTravB > 0)
    {
        step = bestCand[revTravA][revTravB];
        if (step < 0 || step > 2)
            break; // throw error
        
        if(step != 0)
        {
                matchIndex[revTravA] = revTravB;
                revTravA--;
        }
        if (step != 2)
            revTravB--;
    }
    for(i = 0; i < revTravA; i++) 
        matchIndex[i] = revTravB;
    
    return dist;
}

// inputs are two arrays/sequences normalized
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int i, nSamples, nFeats;
    double *totalDist, *featureDists /* may be unneeded */, *matchInterval;
    double *sampleA, *sampleB;
    
    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Three output required.");
    }
    
    /* check that number of rows in first two input arguments are equal */
    if(mxGetM(prhs[0])!=mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("twDistance:unequalRows","Inputs must have equal number of rows.");
    }
    
    /* create a pointer to the real data in the input matrix  */
    nSamples = (int) mxGetN(prhs[1]);
    
    /* get dimensions of the input matrix */
    nFeats = (int) mxGetM(prhs[0]);
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,nFeats,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,2,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    totalDist = mxGetPr(plhs[0]);
    featureDists = mxGetPr(plhs[1]);
    matchInterval = mxGetPr(plhs[2]);
    
    /* try a quick return */
    
    sampleA =mxGetPr(prhs[0]);
    sampleB =mxGetPr(prhs[1]);
    for(i = 0; i < nFeats; i++)
    {
        featureDists[i] = fabs(sampleA[i] - sampleB[i]);
        (*totalDist) += featureDists[i];
    }
    matchInterval[0] = 0;
    matchInterval[1] = 1;
}

/*
 * function [totalDist, featureDists, matchInterval] = timeWarpedDistance(spectrum1, spectrum2, params, varargin)
 * %TIMEWARPEDDISTANCE Dynamic time-warped distance of two spectral feature sets
 * %
 * %  [totalDist, featureDists] = timeWarpedDistance(SPECTRUM1, SPECTRUM2)
 * %  operates as standardDistance does, but takes into account some local time
 * %  warping measure as well, which allows two sequences do have different time courses.
 * %
 * %  A warping cost can be implemented, as well as a cost for the difference in duration.
 * %
 * % totalDist is the distance (L2) over all the features, and featureDists are the
 * % component (L1 difference) distance.
 *
 * % matchInterval returns the times (in sec) of the shorter interval (with relation to
 * % that interval) that best matches the other spectrum
 * %% argument handling
 * if nargin < 3
 * params = defaultParams;
 * end
 * params = processArgs(params,varargin{:});
 *
 * L2 = true && true;
 *
 * %% order two spectrograms by length
 * % convert lengths
 * len1 = length(spectrum1.times);
 * len2 = length(spectrum2.times);
 * if(len1 > len2)
 * temp = spectrum1;
 * spectrum1 = spectrum2;
 * spectrum2 = temp;
 * len1 = length(spectrum1.times);
 * len2 = length(spectrum2.times);
 * end
 * dLen = len2 - len1;
 * % len1 < len2 is true after this
 * % part of locality enforcement
 * margin = params.maxWarpAllowed;
 * %% name all features that should be compared
 * featCatalog = params.featureCatalog;
 *
 * featureSel = false(1,numel(featCatalog));
 * for ii = 1:numel(featCatalog)
 * featureSel(ii) = any(strcmp(fieldnames(spectrum1), featCatalog(ii).name)) && ...
 * any(strcmp(fieldnames(spectrum2), featCatalog(ii).name));
 * end
 * featCatalog = featCatalog(featureSel);
 * nF = numel(featCatalog);
 *
 * %% convert structure to array to speed indexing
 * specArray1 = zeros(nF,len1);
 * specArray2 = zeros(nF,len2); % rows are different features, columns are time points
 * for ii = 1:nF
 * if featCatalog(ii).doLog
 * specArray1(ii,:) = log(spectrum1.(featCatalog(ii).name));
 * specArray2(ii,:) = log(spectrum2.(featCatalog(ii).name));
 * else
 * specArray1(ii,:) = spectrum1.(featCatalog(ii).name);
 * specArray2(ii,:) = spectrum2.(featCatalog(ii).name);
 * end
 * end
 * featureMADs = [featCatalog.MAD]'; % column vector
 *
 * %% score differences between each time point
 * staticCostMatrix = zeros(nF, len1, len2);
 * % use linear indices for taking differences
 *
 * % generate acceptable indices according to warping boundaries:
 * [iPoss, jPoss] = meshgrid(1:len1, 1:len2);
 * dPoss = jPoss-iPoss;
 * inBounds = (dPoss >= -margin & dPoss <= margin + dLen);
 * iPoss = iPoss(inBounds); jPoss = jPoss(inBounds);
 * nPoss = numel(iPoss);
 *
 * % duplicate once for each feature, and get feature subscripts
 * fPoss = ones(nPoss,1) * (1:nF);
 * iPoss = iPoss(:,ones(nF,1));
 * jPoss = jPoss(:,ones(nF,1));
 *
 * % set up linear indices
 * sLin = sub2ind(size(staticCostMatrix), fPoss, iPoss, jPoss);
 * aLin = sub2ind(size(specArray1),       fPoss, iPoss       );
 * bLin = sub2ind(size(specArray2),       fPoss,        jPoss);
 *
 * % calculate cost in one pass
 * staticCostMatrix(sLin) = abs(specArray1(aLin) - specArray2(bLin)) ./ featureMADs(fPoss);
 * if L2
 * staticCostMatrix(sLin) = staticCostMatrix(sLin).^2;
 * end
 *
 * costOverFeatures = permute(sum(staticCostMatrix,1),[2 3 1]);
 *
 * % define warping cost
 * warpingCost = params.warpingCost;
 * if L2
 * warpingCost = warpingCost^2;
 * end
 *
 * %% dynamically search among paths using DP
 * dynCostMatrix = inf(len1, len2);
 * % score extreme first column: compare first spectrum against 1st sample of
 * % second
 * dynCostMatrix(1,1) = costOverFeatures(1,1);
 * dynCostMatrix(1:margin, 1) = cumsum(costOverFeatures(1:margin,1),2) + warpingCost;
 * % score extreme first row: compare second spectrum against 1st sample of
 * % first - we can allow the start to be different here.
 *
 * dynCostMatrix(1, 2:dLen + margin) = costOverFeatures(1,2:dLen + margin);
 *
 * % now fill out in DP style (L1 accumulation over features, L2 over time)
 * % warpingCost determines how much we don't like warping
 *
 * % allow a minimum timewarping window
 * %TODO: (5 samples corresponds to about how many ms?)
 *
 * % about bestcand:
 * %                 value = 1: step up seq1, step up len2;
 * %                 value = 2: const   seq1, step up len2;
 * %                 value = 3: step up seq1, const   len2;
 *
 * % keep track of trajectory
 * bestCand = NaN(len1,len2);
 * bestCand(1,:) = 2;
 *
 *
 * % pre-vectorize warp cost matrix
 * costMatrix = [0 warpingCost warpingCost];
 *
 * % some optimizations here for DP operations
 * for ii = 2:len1
 * % some cost-cutting: enforce some locality (can't warp past a certain
 * % distance off the diagonal)
 * for jj = max(2,ii-(dLen+margin)):min(len2,ii + dLen)
 * %if ii-jj > dLen || ii-jj < - dLen, continue; end;
 * %if ii - jj > dLen || ii - jj < -(dLen + margin), continue; end;
 *
 * % calculate the minimum distances
 * % can optimize further by precomputing the sums
 * %idxs = sub2ind([len1 len2], [ii-1 ii ii-1], [jj-1 jj-1 jj], but fast
 * idxs = (ii + (jj-2) * len1) + [-1, 0, len1];
 *
 * candidates = dynCostMatrix(idxs);
 * if warpingCost > 0
 * candidates = candidates + costMatrix;
 * end
 *
 * %[minCand,bestCand(ii,jj)]= min(candidates), but fast...
 * if candidates(1) < candidates(2) && candidates(1) < candidates(3),
 * minCand = candidates(1); bestCand(ii,jj) = 1;
 * elseif candidates(2) < candidates(3),
 * minCand = candidates(2); bestCand(ii,jj) = 2;
 * else
 * minCand = candidates(3); bestCand(ii,jj) = 3;
 * end
 *
 * % update the dynamic matrix
 * dynCostMatrix(ii,jj) = minCand + costOverFeatures(ii, jj);
 * end
 * end
 *
 * % duration cost doesn't make sense here because we are already counting it
 * % in the warping cost (and for each feature)
 * %durDev = 40; % this is an empirical number based on Lb189_4_25_5
 * %durationDist = (len2-len1) / (durDev / 1000);
 * totalDist = dynCostMatrix(end,end); % durationDist];
 * if L2
 * totalDist = sqrt(totalDist);
 * end
 *
 * % now also figure out what subset of the larger one best matches the shorter one (trimming):
 * % traverse back each path from the end point through bestCand
 * if nargout == 3
 * trav1 = len1; trav2 = len2;
 * matchIndex = zeros(1,len1);
 * while(trav1 > 1 && trav2 > 1)
 * step = bestCand(trav1,trav2);
 * if isnan(step) % something went wrong
 * break;
 * end
 * if step == 1 || step == 3
 * matchIndex(trav1) = trav2;
 * trav1 = trav1 - 1;
 * end
 * if step == 1 || step == 2
 * trav2 = trav2 - 1;
 * end
 * end
 * matchIndex(1:trav1) = trav2;
 * if all(matchIndex ~= 0)
 * matchTimes = spectrum2.times(matchIndex);
 * matchInterval = matchTimes([1,end]);
 * else
 * matchInterval = spectrum2.times(1,end);
 * end
 * end
 *
 * function compDists = getDistance(ind1, ind2)
 * M = numel(ind1); N = numel(ind2);
 * if M == N
 * compDists = abs(specArray1(:,ind1) - specArray2(:,ind2)) ./ featureMADs(:,ones(M,1));
 * elseif M == 1
 * compDists = abs(specArray1(:,ind1*ones(N,1)) - specArray2(:,ind2)) ./ featureMADs(:,ones(N,1));
 * elseif N == 1
 * compDists = abs(specArray1(:,ind1) - specArray2(:,ind2*ones(M,1))) ./ featureMADs(:,ones(M,1));
 * else
 * error('length of ind1 and ind2 must either be equal or at least one must be scalar');
 * end
 * end
 * end
 **/