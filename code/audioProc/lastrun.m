% quickie temporary script
% supposed to get all the features for every syllable, not always
% guaranteed to be one to one

load('syllableBounds','qqq');
for ii = 1:numel(qqq)
    [newFeats,newSylls]=getFeatures(Lb277_3_27_4_Ch1,qqq(ii),noiseProfile,...
        'plot',false,'verbose',true,'playsample',false); 
    feats = [feats newFeats];
    sylls = [sylls newSylls];
end
save('newSyllableBounds','feats','sylls');

