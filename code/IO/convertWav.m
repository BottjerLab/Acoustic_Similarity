function sAudio = convertWav(birdID, sessionID)
fil     = [pwd filesep 'data' filesep birdID filesep sessionID '.mat'];
metafil = [pwd filesep 'data' filesep birdID filesep 'meta-' sessionID '.mat'];
wavfil  = [pwd filesep 'data' filesep birdID filesep sessionID '.wav'];

sAudio = [];
metaStruct = [];
load(fil);
load(metafil);
fprintf('Writing wav file (max, min sound = %0.2f,%0.2f) %s...\n', max(sAudio), min(sAudio), wavfil);

wavwrite(sAudio, 1/metaStruct.interval, wavfil);

