function baseEvents = constructBaselineSilence(songStruct, bLength, syllables, params)
%Jenny 09/09/13

fs = params.fs;
baseEvents = [];


regions = stepSpectrogram(songStruct,params); %finds sounds
silentRegions = eventFromTimes ([regions(1:end-1).stop],[regions(2:end).start],fs); %need to switch starts and stops

% construct base events into segments of bLength length
for i = 1:length(silentRegions)
    howMany = floor((silentRegions(i).stop - silentRegions(i).start)/bLength);
    for j = 1:howMany
        baseEvent = eventFromTimes(silentRegions(i).start + ((bLength *j) - bLength), silentRegions(i).start + (bLength *j), fs);
        baseEvents = [baseEvents; baseEvent];
    end
end

%need to remove ones that are within 2 seconds of singing
syllAdd = addPrePost(syllables, params, 'preroll', 1000, 'postroll', 1000); % adding one second on each side
overlapBaseWithAny = findOverlapsBi(baseEvents,syllAdd');
baseEvents(overlapBaseWithAny(:,1)) = [];
fprintf('Deleted %d baseline events due to overlap...\n',length(overlapBaseWithAny));

%want to remove baseline regions that have some sound (even if not
%singing)
tempSoundEvents = {length(baseEvents)};
%progressbar('Checking baseline regions for complete silence');
for ii = 1:length(baseEvents)
    [~, tempSound] = findSilenceBaseline(songStruct,baseEvents(ii), params, 'plot', false, 'silenceThresh', 0.5);
    
    tempSoundEvents{ii} = tempSound;
 %   progressbar(ii/length(baseEvents));
end

isNoiseTime = false(length(baseEvents),1);
for j = 1:length(baseEvents)
    if ~isempty(tempSoundEvents{j})
        isNoiseTime(j,1) = true;
    end
end
baseEvents(isNoiseTime) = [];

end





