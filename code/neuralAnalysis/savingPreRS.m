function [ neuronData ] = savingPreRS( baselinePeriodsS, spikes, params, aS, fAS, lAS, preTime, neuronData )
%savingPreRS sends spikes off to calculate RS for and significance and saves it in
%neuronData. RS is for pre syllable onset and offsets
%for designated pre-onset/offset times
%   JMA 12/5/13
time = num2str(preTime);
preTime = preTime/1000; % convert to seconds
preAS = eventFromTimes([aS.start] - preTime, [aS.start],params.fs);
[preAS.type] = deal(['approvedSyllables_',time]);
postAS = eventFromTimes([aS.stop] - preTime, [aS.stop],params.fs);
[postAS.type] = deal(['postApprovedSyllables_',time]);
preFAS = eventFromTimes([fAS.start] - preTime, [fAS.start],params.fs);
[preFAS.type] = deal(['firstApprovedSyllables_',time]);
postLAS = eventFromTimes([lAS.start] - preTime, [lAS.start],params.fs);
[postLAS.type] = deal(['lastApprovedSyllables_',time]);

tmpNData = getRS([preAS; baselinePeriodsS], spikes, ['approvedSyllables_',time], 'silence', ['p_p_',time], ['preMeanRS_',time], ['preMeanRSI_',time], ['preStandRS_',time], ['preVarFR_',time], params, 'verbose',false); %onset reponse strength
neuronData = copyPartialStruct(tmpNData, neuronData, ...
    {['p_p_',time], ['preMeanRS_',time], ['preMeanRSI_',time], ['FR_approvedSyllables_',time], ['preStandRS_',time], ['preVarFR_',time]});

tmpNData = getRS([postAS; baselinePeriodsS],spikes, ['postApprovedSyllables_',time], 'silence', ['pp_p_',time], ['postMeanRS_',time], ['postMeanRSI_',time], ['postStandRS_',time],['postVarFR_',time], params, 'verbose',false); %offset response strength
neuronData = copyPartialStruct(tmpNData, neuronData, ...
    {['pp_p_',time], ['postMeanRS_',time], ['postMeanRSI_',time], ['FR_postApprovedSyllables_',time], ['postStandRS_',time], ['postVarFR_',time]});

tmpNData = getRS([preFAS; baselinePeriodsS],spikes, ['firstApprovedSyllables_',time], 'silence', ['fp_p_',time], ['preFMeanRS_',time], ['preFMeanRSI_',time], ['preFStandRS_',time] ,['preFVarFR_',time], params, 'verbose',false); %onset response strength for first syllable in motif
neuronData = copyPartialStruct(tmpNData, neuronData, ...
    {['fp_p_',time], ['preFMeanRS_',time], ['preFMeanRSI_',time], ['FR_firstApprovedSyllables_',time], ['preFStandRS_',time], ['preFVarFR_',time]});

tmpNData = getRS([postLAS; baselinePeriodsS],spikes, ['lastApprovedSyllables_',time], 'silence', ['lp_p_',time], ['postLMeanRS_',time], ['postLMeanRSI_',time], ['postLStandRS_',time], ['postLVarFR_',time], params, 'verbose',false); %offset response strength for last syllable in motif
neuronData = copyPartialStruct(tmpNData, neuronData, ...
    {['lp_p_',time], ['postLMeanRS_',time], ['postLMeanRSI_',time], ['FR_lastApprovedSyllables_',time], ['postLStandRS_',time], ['postLVarFR_',time]});
end

