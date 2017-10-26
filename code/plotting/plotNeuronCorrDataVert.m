function plotNeuronCorrDataVert(params, varargin)
if nargin < 1 || isempty(params)
    params = defaultParams;
end
params = processArgs(params, varargin{:});

% plot excited difference between firing rates for near tutor/far from tutor
% and also p-values for correlations between neurons and firing rates
allNeuronCorrData = [];
load('data/allNeuronCorrelations.mat');

%% flags
isCore        = [allNeuronCorrData.isCore];
isMUA         = [allNeuronCorrData.isMUA];
isPlastic     = [allNeuronCorrData.isPlastic];
isSignificant = [allNeuronCorrData.sigResponse];
isExcited     = [allNeuronCorrData.isExcited];

isPresel = isSignificant & ~isMUA;
RSdiffBins = -10:0.2:10;

% subplot rows / columns
nR = 3; nC = 2;
% measures
distanceTypes = {'tutor', 'intra', 'inter', 'consensus', 'central','humanMatch'}';
distanceDescriptions = {'closest tutor', 'cluster center', 'normed center', ...
    'closest tutor to cluster consensus', 'closest tutor to cluster center', 'expert-designated tutor'};
dFields = [strcat(distanceTypes, '_nearMeanRS') strcat(distanceTypes, '_farMeanRS')];
eiTitle = {'Significantly inhibited single unit-syllable pairs', ...
           'Significantly excited single unit-syllable pairs', ...
           'All significant single unit-syllable pairs'};
xlabels = strcat({'RS for near - far to '}, distanceDescriptions);
filsuff = {'inh','exc','all'};
%%
ycommlims = [-10 10];
for hh = 1:3 % inhibited, excited, all
    for ii = 1:numel(distanceTypes) % six of them
        diffTutorMeanRS = [allNeuronCorrData.(dFields{ii,1})] - [allNeuronCorrData.(dFields{ii,2})];
        if hh < 3
             coreDiffRS = diffTutorMeanRS(isPresel & isExcited == hh-1 &  isCore);
            shellDiffRS = diffTutorMeanRS(isPresel & isExcited == hh-1 & ~isCore);

             coreSubDiffRS = diffTutorMeanRS(isPresel & isExcited == hh-1 &  isCore & ~isPlastic);
            shellSubDiffRS = diffTutorMeanRS(isPresel & isExcited == hh-1 & ~isCore & ~isPlastic);
            
             corePlastDiffRS = diffTutorMeanRS(isPresel & isExcited == hh-1 &  isCore & isPlastic);
            shellPlastDiffRS = diffTutorMeanRS(isPresel & isExcited == hh-1 & ~isCore & isPlastic);

        else
            coreDiffRS  = diffTutorMeanRS(isPresel &  isCore);
            shellDiffRS = diffTutorMeanRS(isPresel & ~isCore);

             coreSubDiffRS = diffTutorMeanRS(isPresel &  isCore & ~isPlastic);
            shellSubDiffRS = diffTutorMeanRS(isPresel & ~isCore & ~isPlastic);
            
             corePlastDiffRS = diffTutorMeanRS(isPresel &  isCore & isPlastic);
            shellPlastDiffRS = diffTutorMeanRS(isPresel & ~isCore & isPlastic);
        end
        pc = signrank( coreDiffRS);
        ps = signrank(shellDiffRS);
        pMannU = ranksum(coreDiffRS(~isnan(coreDiffRS)), shellDiffRS(~isnan(shellDiffRS)));
        fprintf(['%s-%s:\n\tsign-rank test p-value for core: %0.3f' ...
                 '      \n\tsign-rank test p-value for shell: %0.3f',...
                 '      \n\tMann-Whitney U test p-value for core v shell: %0.3f\n'],...
            eiTitle{hh},xlabels{ii},pc,ps,pMannU);

        figure
        % plot the core/shell contrast for single units for this 
        % RS contrast
        subplot(1,2,1);        
        plotInterlaceBarsVert(coreDiffRS, shellDiffRS, RSdiffBins, pMannU < 0.05);
        
        legend(sprintf('CORE:  n = %d', numel( coreDiffRS(~isnan( coreDiffRS)))),...
               sprintf('SHELL: n = %d', numel(shellDiffRS(~isnan(shellDiffRS)))))
        xlabel('Count');
        ylabel(xlabels{ii});
        ylim(ycommlims);
        set(gca,'Box','off');
        set(gca, 'FontSize', 14);
        set(get(gca,'XLabel'),'FontSize', 14);
        set(get(gca,'YLabel'),'FontSize', 14);
        set(get(gca,'Title' ),'FontSize', 14);
        
        title(eiTitle{hh});
        set(gcf,'Color',[1 1 1]);
        
        subplot(1,2,2)
        % plot the core/core subsong/core plastic//
        %          shell/shell subsong/shell plastic individual bars
        [      coreMean,      coreVar] = meanVar(     coreDiffRS);
        [   coreSubMean,   coreSubVar] = meanVar(  coreSubDiffRS);
        [ corePlastMean, corePlastVar] = meanVar(corePlastDiffRS);
        
        [     shellMean,      shellVar] = meanVar(     shellDiffRS);
        [  shellSubMean,   shellSubVar] = meanVar(  shellSubDiffRS);
        [shellPlastMean, shellPlastVar] = meanVar(shellPlastDiffRS);
        xlim([0.5 6.5])
        hold on;
        plotHorzErrorBarVert(1,      coreMean,      coreVar, [0.5 0.5 0.5]);
        plotHorzErrorBarVert(2,   coreSubMean,   coreSubVar, [0.5 0.5 0.5]);
        plotHorzErrorBarVert(3, corePlastMean, corePlastVar, [0.5 0.5 0.5]);
        
        plotHorzErrorBarVert(4,      shellMean,      shellVar, [0.5 0.5 0.5]);
        plotHorzErrorBarVert(5,   shellSubMean,   shellSubVar, [0.5 0.5 0.5]);
        plotHorzErrorBarVert(6, shellPlastMean, shellPlastVar, [0.5 0.5 0.5]);
        ylim(ycommlims);
        
        set(gca,'Box','off', 'XTick', 1:6,...
            'XTickLabel',{'Core', 'Core/Subsong','Core/Plastic',...
            'Shell','Shell/Subsong','Shell/Plastic'});
        xticklabel_rotate([],45, [], 'FontSize', 14);
        hold on;
        keyboard
    end
    if params.saveplot
        saveCurrFigure(sprintf('figures/distanceCorrelations/RSdiffs-SUA-%s.jpg', filsuff{hh}));
    end
end

%{
pFields = strcat(distanceTypes, 'Distance_p');
xlabels = strcat({'Linear trend p-values of FR to '}, distanceDescriptions);
   
for hh = 1:3 % inhibited, excited, all
    figure(hh+3);
    for ii = 1:numel(distanceTypes)
        subplot(nR,nC,ii)
        corrPVals = [allNeuronCorrData.(pFields{ii})];
        if hh < 3
             coreCPVs = corrPVals(isPresel & isExcited == hh-1 &  isCore);
            shellCPVs = corrPVals(isPresel & isExcited == hh-1 & ~isCore);
        else
             coreCPVs = corrPVals(isPresel &  isCore);
            shellCPVs = corrPVals(isPresel & ~isCore);
        end
        pMannU = ranksum(coreCPVs(~isnan(coreCPVs)), shellCPVs(~isnan(shellCPVs)));
        fprintf('%s - %s: Mann-Whitney U p-value for core v shell: %0.3f\n',...
            eiTitle{hh}, xlabels{ii},pMannU);
        pBins = logspace(-3,0,30);
        plotInterlaceBarsVert(coreCPVs, shellCPVs, pBins, pMannU < 0.05);
        
        legend(sprintf('CORE:  n = %d', numel( coreCPVs(~isnan( coreCPVs)))),...
               sprintf('SHELL: n = %d', numel(shellCPVs(~isnan(shellCPVs)))));
        xlabel(xlabels{ii});
        ylabel('Count');
        xlim([0 1])
        set(gca, 'XTick', [0.01 0.05 0.1 0.2 0.4 0.6 0.8]);
        set(gca, 'Box', 'off');
        set(gca, 'FontSize', 14);
        set(get(gca,'XLabel'),'FontSize', 14);
        
        set(get(gca,'YLabel'),'FontSize', 14);
        set(get(gca,'Title' ),'FontSize', 14);
    end
    subplot(nR,nC,1);
    title(eiTitle{hh});
    set(gcf,'Color',[1 1 1]);
    if params.saveplot
        saveCurrFigure(sprintf('figures/distanceCorrelations/neuronDistanceCorr-SUA-%s.jpg', filsuff{hh}));
    end
end
%}
end

function [m, v] = meanVar(set1)
m = nanmean(set1); 
v = nanstd(set1 )/sqrt(numel(set1)-1);
end

function plotInterlaceBarsVert(setCore, setShell, bins, sigLevel)
% core plotted in gray, shell plotted in red

% run mann-whitney u test
[p,h] = ranksum(setCore, setShell);
binTol = 1e-5;
hCore  = histc(setCore , bins);
hShell = histc(setShell, bins);
[m1 sem1] = meanVar(setCore ); 
[m2 sem2] = meanVar(setShell);

if all(diff(bins) - (bins(2) - bins(1)) < binTol)
    bw = bins(2) - bins(1);
    barh(bins, hCore, 0.5, 'FaceColor', [0.5 0.5 0.5]); 
    hold on; 
    barh(bins+bw/2, hShell, 0.5, 'r');
else
    bw = diff(bins); bw = [bins(1)/2 bw bw(end)];
    % make visible legend groups
    ghCore = hggroup; ghShell = hggroup;
    set(get(get(ghCore , 'Annotation'),'LegendInformation'), 'IconDisplayStyle','on');
    set(get(get(ghShell, 'Annotation'),'LegendInformation'), 'IconDisplayStyle','on');
    for ii = 1:numel(bins) % draw histogram bin by bin
        % core bin
        yd = bins(ii) - bw(ii)/2; yu = bins(ii);
        xl = 0; xr = hCore(ii);
        patch([xl xr; xr xr; xl xl], [yd yd; yd yu; yu yu], [0.5 0.5 0.5],...
            'EdgeColor','none', 'Parent', ghCore);
        hold on;
        % shell bins
        yd = bins(ii); yu = bins(ii) + bw(ii+1)/2;
        xl = 0; xr = hShell(ii);
        patch([xl xr; xr xr; xl xl], [yd yd; yd yu; yu yu], [1 0 0],...
            'EdgeColor','none', 'Parent', ghShell);                
    end
end
xlims = xlim; 
right = xlims(2) + 2;
xlim([0 right]);
plotHorzErrorBarVert(right - 0.8, m1, sem1, [0.5 0.5 0.5]);
hold on;
plotHorzErrorBarVert(right - 1.2, m2, sem2, [1 0 0]);
if sigLevel > 0
    plot(right-1.0, mean([m1 m2]), 'k*','MarkerSize',6);
end
hold off;
end

function plotHorzErrorBarVert(x, y, yerr, col)
    xw = 0.02*diff(xlim);
    % the connecting line, and the bottom and top end lines
    xx = [x      x      NaN x-xw   x+xw   NaN x-xw   x+xw  ]; 
    yy = [y-yerr y+yerr NaN y-yerr y-yerr NaN y+yerr y+yerr];
    plot(xx,yy, '-','Color', col, 'LineWidth', 1.5);        
end