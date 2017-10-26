% Step1: reading Data from the file
file_data = load('Ionosphere.txt');
Data = file_data(:,1:end-1)';
Labels = file_data(:, end)';
Labels = Labels*2 - 1;

MaxIter = 100; % boosting iterations

% Step2: splitting data to training and control set
TrainData   = Data(:,1:2:end);
TrainLabels = Labels(1:2:end);

ControlData   = Data(:,2:2:end);
ControlLabels = Labels(2:2:end);

% Step3: constructing weak learner
weak_learner = tree_node_w(3); % pass the number of tree splits to the constructor

% Step4: training with Gentle AdaBoost
[RLearners RWeights] = RealAdaBoost(weak_learner, TrainData, TrainLabels, MaxIter);

% Step5: training with Modest AdaBoost
[MLearners MWeights] = ModestAdaBoost(weak_learner, TrainData, TrainLabels, MaxIter);

% Step6: evaluating on control set
ResultR = sign(Classify(RLearners, RWeights, ControlData));

ResultM = sign(Classify(MLearners, MWeights, ControlData));

% Step7: calculating error
ErrorR  = sum(ControlLabels ~= ResultR) / length(ControlLabels)

ErrorM  = sum(ControlLabels ~= ResultM) / length(ControlLabels)