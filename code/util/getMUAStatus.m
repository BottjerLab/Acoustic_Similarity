function isMUAUnit = getMUAStatus(spikeFiles, nNeuronsPerFile, MUAmeta)
% utility, only called by loadSpikeData
    cumIndex = [0 cumsum(nNeuronsPerFile)];
    isMUAUnit = false(cumIndex(end),1);
    for kk = 1:numel(spikeFiles)
        unitsFromFile = (cumIndex(kk)+1):cumIndex(kk+1);                
        % MUA/SUA-ness of each 'cluster' depends on the metadata
        [~,thisTrodeID] = fileparts(spikeFiles{kk});
        lookupMUA = strcmp(thisTrodeID, strcat({MUAmeta.trodeID}, '_times'));       
        if any(lookupMUA)
            idxInFile = MUAmeta(lookupMUA).MUAunits;
            idxInFile(idxInFile > nNeuronsPerFile(kk)) = [];       
            idxInSession = unitsFromFile(idxInFile);
            isMUAUnit(idxInSession) = true;
        end
    end
end