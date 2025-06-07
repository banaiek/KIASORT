function net = mlpDimReduction(X, Y, hp)


%% Split the data into training and validation sets

N = size(X,1);
% Define split Ratios and calculating split sizes
trainRatio = 1 - hp.validationRatio;
numTrain = floor(trainRatio * N);
numVal = N - numTrain;

% Shuffle the Data
rng(42);  % For reproducibility
shuffledIndices = randperm(N);
XShuffled = X(shuffledIndices, :);
YShuffled = Y(shuffledIndices, :);

% Split the Data
XTrain = XShuffled(1:numTrain, :);
YTrain = YShuffled(1:numTrain, :);

XVal = XShuffled(numTrain+1:end, :);
YVal = YShuffled(numTrain+1:end, :);


%% Setting MLP Model Architecture

M = size(X,2);
P = size(Y,2);
% Initialize the layers with the input layer
layers = [
    featureInputLayer(M, 'Name', 'input')
    ];

% Dynamically add hidden layers based on numHidLayer
for i = 1:hp.numHidLayer
    layers = [
        layers
        fullyConnectedLayer(hp.hiddenLayerSize(i), 'Name', ['fc' num2str(i)])
        batchNormalizationLayer('Name', ['bn' num2str(i)], 'MeanDecay', 0.05, 'VarianceDecay', 0.05)
        geluLayer('Name', ['gelu' num2str(i)])
        dropoutLayer(hp.dropoutRate(i), 'Name', ['dropout' num2str(i)])
        ];
end

% Append the final output layer
layers = [
    layers
    fullyConnectedLayer(P, 'Name', 'fc_output')
    ];


%% Specifying training options and training the MLP model

% Training Options
options = trainingOptions("adam", ...
    "MaxEpochs", hp.maxEpochs, ...
    "MiniBatchSize", hp.miniBatchSize, ...
    "Shuffle","every-epoch", ...
    "ValidationData",{XVal, YVal}, ...
    "ValidationFrequency", floor(numTrain/64), ...
    "Verbose",true, ...
    "Plots","training-progress", ...
    "LearnRateSchedule","piecewise", ...
    "LearnRateDropFactor", hp.learningRateDropFactor, ...
    "LearnRateDropPeriod", hp.learningRateDropPeriod, ...
    "InitialLearnRate", hp.initialLearningRate, ...
    "L2Regularization", hp.l2Reg, ...
    "GradientThreshold", hp.gradientThreshold, ...
    "ExecutionEnvironment", hp.executionEnvironment);

net = trainnet(XTrain, YTrain, layers, hp.mlpLossFcn, options);

end