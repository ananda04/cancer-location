%% load Data
%unzip("NT_Brain.zip")
imds = imageDatastore("NT_Brain", ...
    'LabelSource','foldernames');
label = [1,0,1,0,1,0,1,0]
label = num2cell(transpose(label))
a = [imds.Files, label]
numTrainFiles = 8;
c = cvpartition(a,"Holdout",0.25)

%[imdsTrain,imdsValidation] = splitEachLabel(a,numTrainFiles,'randomized');

%% Network Architecture 
inputSize = [1260 974 3];
numClasses = 10;

layers = [
    imageInputLayer(inputSize)
    convolution2dLayer(5,20)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];
%% Train Network
options = trainingOptions('sgdm', ...
    'MaxEpochs',4, ...
    'ValidationData',imdsValidation, ...
    'ValidationFrequency',30, ...
    'Verbose',false, ...
    'Plots','training-progress');

net = trainNetwork(imdsTrain,layers,options);