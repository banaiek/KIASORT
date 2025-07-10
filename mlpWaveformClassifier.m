function net = mlpWaveformClassifier(X, labels, hp)

labels = categorical(labels);
uniqueLabels = categories(labels);
numClasses = numel(uniqueLabels);

% Split data into training and validation sets
% ---------------------------%
N = size(X, 1);  % Number of observations

valRatio = hp.validationRatio;
trainRatio = 1 - valRatio;
numTrain = floor(trainRatio * N);
numVal = N - numTrain;

% Shuffle the data and labels
rng(141);  % For reproducibility
shuffledIndices = randperm(N);
XShuffled = X(shuffledIndices, :);
YShuffled = labels(shuffledIndices);

% Split the Data
XTrain = XShuffled(1:numTrain, :);
YTrain = YShuffled(1:numTrain);
XVal = XShuffled(numTrain+1:end, :);
YVal = YShuffled(numTrain+1:end);

% Define MLP model architecture
% ---------------------------%
M = size(XTrain, 2);  % Number of input features

% Define input layer
layers = [
    featureInputLayer(M, 'Name', 'input')
    ];

% Add hidden layers based on numHidLayer
for i = 1:hp.numHidLayer
    layers = [
        layers
        fullyConnectedLayer(hp.hiddenLayerSize(i), 'Name', ['fc' num2str(i)])
        batchNormalizationLayer('Name', ['bn' num2str(i)], 'MeanDecay', 0.05, 'VarianceDecay', 0.05)
        geluLayer('Name', ['gelu' num2str(i)])
        dropoutLayer(hp.dropoutRate(i), 'Name', ['dropout' num2str(i)])
        ];
end

% Add the final output layers
layers = [
    layers
    fullyConnectedLayer(numClasses, 'Name', 'fc_output')  % Output layer with numClasses neurons
    softmaxLayer('Name','softmax')                        % Softmax layer for probabilities    
    ];

% Training Options
% ---------------------------%
options = trainingOptions("adam", ...
    "MaxEpochs", hp.maxEpochs, ...
    "MiniBatchSize", hp.miniBatchSize, ...
    "Shuffle", "every-epoch", ...
    "ValidationData", {XVal, YVal}, ...
    "ValidationFrequency", floor(numTrain / 64), ...
    "Verbose", true, ...
    "LearnRateSchedule", "piecewise", ...
    "LearnRateDropFactor", hp.learningRateDropFactor, ...
    "LearnRateDropPeriod", hp.learningRateDropPeriod, ...
    "InitialLearnRate", hp.initialLearningRate, ...
    "L2Regularization", hp.l2Reg, ...
    "GradientThreshold", hp.gradientThreshold, ...
    "ExecutionEnvironment", hp.executionEnvironment);

% Train the MLP model and the network using trainnet
% ---------------------------%
net = trainnet(XTrain, YTrain, layers,'crossentropy', options);

end