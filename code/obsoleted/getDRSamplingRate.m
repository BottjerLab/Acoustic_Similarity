function fs = getDRSamplingRate(region)
[st pa] = fiieparts(region.file)
metaFile = [pa filesep 'meta-' st];
fs = 1/getfield(getfield(load(metaFile, 'metaStruct'),'metaStruct'),'interval');