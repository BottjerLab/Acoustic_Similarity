function [corrSig, sigLevel, corrSigP]=multLinearRegress(xStruct, yStruct, varargin)
    
% remove 
% convert xStruct to column data
[X, xNames] = structArrayToColumn(xStruct); 

% convert yStruct to column data
[allY, yNames] = structArrayToColumn(yStruct); 

% clean zeroed data (is this an exact match?)
Xclean = X;
Xclean(:,any(X==0,1)) = [];
xNames(any(X==0,1)) = [];
%keyboard

% significance levels
sig(1) = 0.05;
sig(2) = sig(1) / numel(xNames);
sig(3) = sig(2) / numel(yNames);

% run the linear regression
for ii = 1:numel(yNames)
    [b,dev,stats]=glmfit(Xclean,allY(:,ii));
    posFiring = (0 < allY(:,ii));
    [bP, devP, statsP] = glmfit(Xclean(posFiring,:),allY(posFiring,ii));
    corrSig.constant(ii) = stats.p(1);
    corrSigP.constant(ii) = statsP.p(1);
    for jj = 1:numel(xNames)
        corrSig.(xNames{jj})(ii) = stats.p(jj+1);
        corrSigP.(xNames{jj})(ii) = statsP.p(jj+1);

        sigLevel.(xNames{jj})(ii) = 0; 
        for kk = 1:numel(sig)
            sigLevel.(xNames{jj})(ii) = sigLevel.(xNames{jj})(ii) + [stats.p(jj+1) < sig(kk)];
        end
        if sigLevel.(xNames{jj})(ii) == numel(sig) % unquestionably significant after mult. t-tests
            figure; plot(Xclean(:,jj), allY(:,ii),'r.'); 
            xlabel(sprintf('Similarity to feature %s', nos(xNames{jj}))); 
            ylabel(sprintf('Firing rate within syllable (Hz), %s',nos(yNames{ii}))); 
            fprintf('With zero firing rate, p = %f, without, p = %f\n',...
                corrSig.(xNames{jj})(ii), corrSigP.(xNames{jj})(ii))
            title(sprintf('With zero firing rate, p = %f, without, p = %f',...
                corrSig.(xNames{jj})(ii), corrSigP.(xNames{jj})(ii)));
            drawnow;
            pause;
        end
    end
end


end

function str = nos(str)
    str = strrep(str,'_', ' ');
end
function [colArray, names] = structArrayToColumn(structArray)
names = fieldnames(structArray);
colArray = zeros(numel(structArray), numel(names));
for kk = 1:numel(names)
    colArray(:,kk) = [structArray.(names{kk})];
end
end