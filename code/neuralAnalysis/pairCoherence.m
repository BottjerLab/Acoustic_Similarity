function ret = pairCoherence(clip1, clip2, specParams)

    %under heavy construction
    nWindow = specParams.windowSize / 1000 * specParams.fs;
    nWindow2 = 2.^ceil(log2(nWindow));
    nOverlap = floor(nWindow2 * specParams.nOverlap / specParams.windowSize); 
    ret = cpsd(clip1, clip2, hann(nWindow2), nOverlap, specParams.NfreqBands);  % problem: these are not the same length
end