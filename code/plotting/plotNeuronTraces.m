function plotNeuronTraces(neurons)

% plot rasters in a quick way, one row per column vector in the 
% cell array 'neurons'

nNeurons = numel(neurons);
neuronInds = cell(1,nNeurons);
for ii = 1:nNeurons
    neuronInds{ii} = ii * ones(numel(neurons{ii}),1);
end
xx = vertcat(neurons{:});
yy = vertcat(neuronInds{:});
plot(xx,yy,'s','MarkerFaceColor', 'k', 'MarkerSize', 2);
ylim([0.5 nNeurons+0.5]);
xlim([0 max(xx)+1]);
end