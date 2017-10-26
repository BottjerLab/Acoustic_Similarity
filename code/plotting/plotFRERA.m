function plotFRERA(subset1, subset2, field, bins)
% plots a firing rate comparison

if ~isempty(subset1)
    set1Hist     = sum(horzcat(subset1.(field)),2)/length(subset1);
    set1HistSem  = std(horzcat(subset1.(field)),[],2)/sqrt(length(subset1)); 
    plotLineTransparentSEM(bins(1:end-3), set1Hist, set1HistSem, [0.5 0.5 0.5]);%removing bins because removed last three points in eraFind
    hold on;
end
if ~isempty(subset2)
    set2Hist     = sum(horzcat(subset2.(field)),2)/length(subset2);
    set2HistSem  = std(horzcat(subset2.(field)),[],2)/sqrt(length(subset2)); 
    plotLineTransparentSEM(bins(1:end-3), set2Hist, set2HistSem, [1 0.2 0.2]); %removing bins because removed last three points in eraFind
end
if ~isempty(subset1)
    hold off;
end

% let's look for significance
N1 = length(subset1);
N2 = length(subset2);

t_stat = tinv(0.975, N1-1); %0.05 two-tailed critical value 
if any((set1Hist - set1HistSem * t_stat) > 0) 
    sig_set1IdxPos = find((set1Hist - set1HistSem * t_stat) > 0);
    sig_set1xxPos = insertNonConsecNans(bins(sig_set1IdxPos), sig_set1IdxPos);
    hold on;
    ybar = ylim * [0.08 0.92]'; % 92% of the way to the top
    plot(sig_set1xxPos, ybar * ones(1,numel(sig_set1xxPos)), '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); 
end

if any((set1Hist + set1HistSem * t_stat) < 0)    
    sig_set1IdxNeg = find((set1Hist + set1HistSem * t_stat) < 0);    
    sig_set1xxNeg = insertNonConsecNans(bins(sig_set1IdxNeg), sig_set1IdxNeg);
    hold on;
    ybar = ylim * [0.92 0.08]'; % 92% of the way to the bottom
    plot(sig_set1xxNeg, ybar * ones(1,numel(sig_set1xxNeg)), '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);     
end

t_stat = tinv(0.975, N2-1); %0.05 two-tailed critical value 
if any((set2Hist - set2HistSem * t_stat) > 0)    
    sig_set2IdxPos = find((set2Hist - set2HistSem * t_stat) > 0);        
    sig_set2xxPos = insertNonConsecNans(bins(sig_set2IdxPos), sig_set2IdxPos);
    hold on;
    ybar = ylim * [0.12 0.88]'; % 88% of the way to the top
    plot(sig_set2xxPos, ybar * ones(1,numel(sig_set2xxPos)), 'r-', 'LineWidth', 2);     
end

if any((set2Hist + set2HistSem * t_stat) < 0)
    sig_set2IdxNeg = find((set2Hist + set2HistSem * t_stat) < 0);
    sig_set2xxNeg = insertNonConsecNans(bins(sig_set2IdxNeg), sig_set2IdxNeg);   
    hold on;
    ybar = ylim * [0.88 0.12]'; % 88% of the way to the bottom
    plot(sig_set2xxNeg, ybar * ones(1,numel(sig_set2xxNeg)), 'r-', 'LineWidth', 2); 
end
end
function ret = insertNonConsecNans(lookupvec, idxvec)
    nonConsec = find(diff(idxvec)>1);
    ret = lookupvec;
    for jj = 1:numel(nonConsec)
        ret = [ret(1:nonConsec(jj)); NaN; ...
            ret(nonConsec(jj)+1:end)];
        nonConsec = nonConsec+1;
    end
    xlim([-150 150]);
end