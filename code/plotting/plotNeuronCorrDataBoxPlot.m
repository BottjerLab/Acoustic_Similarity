function plotNeuronCorrData(params, varargin)
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
%isSignificant = true(1,numel(allNeuronCorrData));
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
for hh = 1:3 % inhibited, excited, all
    for ii = 1:numel(distanceTypes) % six of them
        figure;
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

        % clunky way just to get the top histogram value
        plotInterlaceBars(coreDiffRS, shellDiffRS, RSdiffBins, pMannU < 0.05);
        ytop = ylim * [0 1]';
        % plot boxplots above - b/c of boxplot behavior, this eliminates
        % the histograms        
        
        grps = [ones(size(coreDiffRS)) 2*ones(size(coreSubDiffRS)) 3*ones(size(corePlastDiffRS)) ,...
            4*ones(size(shellDiffRS)), 5*ones(size(shellSubDiffRS)) 6*ones(size(shellPlastDiffRS))];
        dats = [coreDiffRS coreSubDiffRS corePlastDiffRS shellDiffRS shellSubDiffRS shellPlastDiffRS];
        boxplot(gca, dats, grps, 'orientation', 'horizontal', 'notch','on','positions', ytop:ytop+5, ...
            'colors', [0.5 * ones(3); [1 1 1]' * [1 0 0]]);
        
        % plot histograms underneath again
        hold on;
        plotInterlaceBars(coreDiffRS, shellDiffRS, RSdiffBins, pMannU < 0.05);
        ylim([0 ytop+6]);
        plot([0 0], ylim,'k--','HandleVisibility','off');
        legend(sprintf('CORE:  n = %d', numel( coreDiffRS(~isnan( coreDiffRS)))),...
               sprintf('SHELL: n = %d', numel(shellDiffRS(~isnan(shellDiffRS)))))
        
       
%{
        plotSEMBar(      coreDiffRS, ytop  , [0.5 0.5 0.5]); 
        plotSEMBar(   coreSubDiffRS, ytop+1, [0.5 0.5 0.5]);
        plotSEMBar( corePlastDiffRS, ytop+2, [0.5 0.5 0.5]);
        plotSEMBar(     shellDiffRS, ytop+3, [  1   0   0]); 
        plotSEMBar(  shellSubDiffRS, ytop+4, [  1   0   0]);
        plotSEMBar(shellPlastDiffRS, ytop+5, [  1   0   0]);
        %}
        
        % redo y axis labels
        
        yt = get(gca,'YTick');
        yt = [yt(yt < ytop) ytop:ytop+5];
        ytl = cellfun(@(x) sprintf('%d',x),num2cell(yt),'UniformOutput',false);
        ytl(end-5:end) = {'Core','Core/Subsong','Core/Plastic','Shell','Shell/Subsong','Shell/Plastic'};
        set(gca,'YTick',yt,'YTickLabel',ytl);
        
        
        % figure formatting
        xlabel(xlabels{ii});
        ylabel('Count');
        xlim([-10 10]);
        set(gca,'Box','off');
        set(gca, 'FontSize', 14);
        set(get(gca,'XLabel'),'FontSize', 14);
        set(get(gca,'YLabel'),'FontSize', 14);
        set(get(gca,'Title' ),'FontSize', 14);
        title(eiTitle{hh});
        set(gcf,'Color',[1 1 1]);
        hold off;
       
        if params.saveplot
            saveCurrFigure(sprintf('figures/distanceCorrelations/RSdiffs-SUA-box-%s-%s.jpg', distanceTypes{ii}, filsuff{hh}));
        end
    end

end

pFields = strcat(distanceTypes, 'Distance_p');
xlabels = strcat({'Linear trend p-values of FR to '}, distanceDescriptions);
   
for hh = 1:3 % inhibited, excited, all
    figure;
    for ii = 1:numel(distanceTypes)
        subplot(nR,nC,ii)
        corrPVals = [allNeuronCorrData.(pFields{ii})];
        if hh < 3
             coreCPVs = corrPVals(isPresel & isExcited == hh-1 &  isCore);
            shellCPVs = corrPVals(isPresel & isExcited == hh-1 & ~isCore);
        
             %coreSubCPVs = corrPVals(isPresel & isExcited == hh-1 &  isCore & ~isPlastic);
            %shellSubCPVs = corrPVals(isPresel & isExcited == hh-1 & ~isCore & ~isPlastic);
            
             %corePlastCPVs = corrPVals(isPresel & isExcited == hh-1 &  isCore & isPlastic);
            %shellPlastCPVs = corrPVals(isPresel & isExcited == hh-1 & ~isCore & isPlastic);
        else
             coreCPVs = corrPVals(isPresel &  isCore);
            shellCPVs = corrPVals(isPresel & ~isCore);

             %coreSubCPVs = corrPVals(isPresel &  isCore & ~isPlastic);
            %shellSubCPVs = corrPVals(isPresel & ~isCore & ~isPlastic);
            
             %corePlastCPVs = corrPVals(isPresel &  isCore & isPlastic);
            %shellPlastCPVs = corrPVals(isPresel & ~isCore & isPlastic);
        end
        pMannU = ranksum(coreCPVs(~isnan(coreCPVs)), shellCPVs(~isnan(shellCPVs)));
        fprintf('%s - %s: Mann-Whitney U p-value for core v shell: %0.3f\n',...
            eiTitle{hh}, xlabels{ii},pMannU);
        pBins = logspace(-3,0,30);
        plotInterlaceBars(coreCPVs, shellCPVs, pBins, pMannU < 0.05);
        
        legend(sprintf('CORE:  n = %d', numel( coreCPVs(~isnan( coreCPVs)))),...
               sprintf('SHELL: n = %d', numel(shellCPVs(~isnan(shellCPVs)))));
        xlabel(xlabels{ii});
        ylabel('Count');
        xlim([0 1])
        set(gca, 'XTick', [0.01 0.05 0.1 0.2 0.4 0.6 0.8]);
        set(gca, 'Box', 'off');
        set(gca, 'FontSize', 11);
        set(get(gca,'XLabel'),'FontSize', 14);        
        set(get(gca,'YLabel'),'FontSize', 14);
        set(get(gca,'Title' ),'FontSize', 14);

        ytop = ylim * [0 1]'; ylim([0 ytop+2]);
        plotSEMBar( coreCPVs, ytop, [0.5 0.5 0.5]); 
        plotSEMBar(shellCPVs, ytop, [1 0 0]); 
    end
    subplot(nR,nC,1);
    title(eiTitle{hh});
    set(gcf,'Color',[1 1 1]);
    if params.saveplot
        saveCurrFigure(sprintf('figures/distanceCorrelations/neuronDistanceCorr-SUA-%s.jpg', filsuff{hh}));
    end

end
end

function plotInterlaceBars(setCore, setShell, bins, sigLevel)
% core plotted in gray, shell plotted in red

% run mann-whitney u test
[p,h] = ranksum(setCore, setShell);
binTol = 1e-5;
hCore  = histc(setCore , bins);
hShell = histc(setShell, bins);
%[m1 sem1] = meanSEM( setCore);
%[m2 sem2] = meanSEM(setShell);

holdState = ishold;
if all(diff(bins) - (bins(2) - bins(1)) < binTol)
    bw = bins(2) - bins(1);
    bar(bins, hCore, 0.5, 'FaceColor', [0.5 0.5 0.5]); 
    hold on; 
    bar(bins+bw/2, hShell, 0.5, 'r');
else
    bw = diff(bins); bw = [bins(1)/2 bw bw(end)];
    % make visible legend groups
    ghCore = hggroup; ghShell = hggroup;
    set(get(get(ghCore , 'Annotation'),'LegendInformation'), 'IconDisplayStyle','on');
    set(get(get(ghShell, 'Annotation'),'LegendInformation'), 'IconDisplayStyle','on');
    for ii = 1:numel(bins) % draw histogram bin by bin
        % core bin
        xl = bins(ii) - bw(ii)/2; xr = bins(ii);
        yd = 0; yu = hCore(ii);
        patch([xl xl; xl xr; xr xr], [yd yu; yu yu; yd yd], [0.5 0.5 0.5],...
            'EdgeColor','none', 'Parent', ghCore);
        hold on;
        % shell bins
        xl = bins(ii); xr = bins(ii) + bw(ii+1)/2;
        yd = 0; yu = hShell(ii);
        patch([xl xl; xl xr; xr xr], [yd yu; yu yu; yd yd], [1 0 0],...
            'EdgeColor','none', 'Parent', ghShell);        
    end    
end

% plot the error bars w/ SEM on top
%{
ylims = ylim; 
top = ylims(2) + 2;
ylim([0 top]);
plotHorzErrorBar(m1, top - 0.8, sem1, [0.5 0.5 0.5]);
plotHorzErrorBar(m2, top - 1.2, sem2, [1 0 0]);
if sigLevel > 0
    plot(mean([m1 m2]), top-1, 'k*','MarkerSize',8);
end
%}
if ~ishold
    hold off
end
end

function [m, v] = meanSEM(set1)
m = nanmean(set1); 
v = nanstd(set1 )/sqrt(numel(set1)-1);
end

function plotSEMBar(set, y, col)
[m,v] = meanSEM(set);
plotHorzErrorBar(m,y,v,col);

end

function plotHorzErrorBar(x, y, xwidth, col)
    yh = 0.02*diff(ylim);
    xx = [x-xwidth x+xwidth NaN x-xwidth x-xwidth NaN x+xwidth x+xwidth];
    yy = [y y NaN y-yh y+yh NaN y-yh y+yh];
    plot(xx,yy, '-','Color', col, 'LineWidth', 1.5);    
    %hold on;
    %plot(x,y,'.','MarkerSize', 16, 'Color', col);
end