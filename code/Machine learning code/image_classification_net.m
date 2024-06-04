%% load Data
unzip("Brains.zip")
imds = imageDatastore("Brains","IncludeSubfolders",1,'LabelSource','foldernames');
[imdsTrain,imdsValidation] = splitEachLabel(imds,0.75);
%% Network Architecture 
inputSize = [974 1260 3];
numClasses = 2;

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
%%
YTest = classify(net,imdsValidation);