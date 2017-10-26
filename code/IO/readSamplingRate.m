function fs = readSamplingRate(manifest)
fs = 1/getfield(loadFromManifest(manifest, 'metaStruct'),'interval');
