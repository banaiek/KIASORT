function net = cnnWaveformClassifier(X, labels, hp)


labels = categorical(labels);
uniqueLabels = categories(labels);
numClasses = numel(uniqueLabels);

%  Split data
% ---------------------------%
[N, C, T] = size(X);
valRatio = hp.validationRatio;
trainRatio = 1 - valRatio;
numTrain = floor(trainRatio * N);
numVal = N - numTrain;

rng(141);  % For reproducibility
shuffledIndices = randperm(N);
XShuffled = X(shuffledIndices, :, :);   % 2D structure for CNN
YShuffled = labels(shuffledIndices);

XTrain = XShuffled(1:numTrain, :, :);
YTrain = YShuffled(1:numTrain);
XVal   = XShuffled(numTrain+1:end, :, :);
YVal   = YShuffled(numTrain+1:end);

% Reshape data for CNN
% ---------------------------%
% CNN input: (Height, Width, Channels, NumberOfObservations) (C, T, 1, N)
XTrainCNN = permute(XTrain, [2, 3, 1]);   % becomes (C, T, N)
XTrainCNN = reshape(XTrainCNN, [C, T, 1, size(XTrainCNN,3)]);

XValCNN   = permute(XVal, [2, 3, 1]);     % becomes (C, T, N_val)
XValCNN   = reshape(XValCNN, [C, T, 1, size(XValCNN,3)]);



% Define CNN model architecture
% ---------------------------%
layers = [
    imageInputLayer([C T 1], 'Name','input','Normalization','none')
];

% Add convolution + pooling blocks
for i = 1:hp.numConvBlocks
    convLayerName = sprintf('conv%d', i);
    bnLayerName   = sprintf('bn%d', i);
    reluLayerName = sprintf('relu%d', i);
    poolLayerName = sprintf('pool%d', i);
    dropoutName   = sprintf('dropout%d', i);

    layers = [
        layers
        convolution2dLayer(hp.filterSize, hp.numFilters(i), 'Padding','same','Name',convLayerName)
        batchNormalizationLayer('Name', bnLayerName)
        reluLayer('Name', reluLayerName)
        maxPooling2dLayer(hp.poolSize, 'Stride', 2, 'Name', poolLayerName,'Padding','same')
        dropoutLayer(hp.dropoutRate(i), 'Name', dropoutName)
    ];
end

% Add final classification layers
layers = [
    layers
    fullyConnectedLayer(numClasses, 'Name','fc_output')
    softmaxLayer('Name','softmax')
];


% Training Options
% ---------------------------%
options = trainingOptions("adam", ...
    "MaxEpochs", hp.maxEpochs, ...
    "MiniBatchSize", hp.miniBatchSize, ...
    "Shuffle", "every-epoch", ...
    "ValidationData", {XValCNN, YVal}, ...
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
net = trainnet(XTrainCNN, YTrain, layers, hp.cnnLossFcn, options);


end