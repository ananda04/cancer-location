% Create TQC Image
% Created by Stefan Stryker

% run in folder with single .mat file contianing image data of interest
files = dir('*.mat');
names = files.name; 

%% User Selectable Parameters

minQ = 0.05;              % where to trim spectra for mean q calculation (in 1/A)
maxQ = 0.305;             % ||
minVal = 0.163;           % what to window TQC image to  % CHANGE THIS WINDOWING OF minVal and maxVal FOR SHIFTING COLOR RANGE
maxVal = 0.198;           % ||  
tBP = 10;                 % transmission bottom percentile for display
tTP = 100;               % transmission top percentile for display

%% Open Reconstructed Data, create qvalue range

load(names) 

qvals = qvals'; % variable in data file

% finding index of trim locations
[~,minQInd] = min(abs(qvals-minQ));
[~,maxQInd] = min(abs(qvals-maxQ));

%% Compute Weighted Mean Q

meanQMat = zeros(size(reconstructedData,1),size(reconstructedData,2));
for i = 1 : size(reconstructedData,1)
    for j = 1 : size(reconstructedData,2)

        spectra = squeeze(reconstructedData(i,j,minQInd:maxQInd));  % seeing if smooth improves stability
        meanQMat(i,j) = sum((spectra .* qvals(minQInd:maxQInd))) ./ sum(spectra);         
    end
end

%% Crop Transmission Image

startScat = scatSliceLocations(1);
endScat = scatSliceLocations(end);

% accounting for motor moving 0.74mm vs 1mm
startScat = ((158*startScat+1541)-1582)/117;  % floor()
endScat = ((158*endScat+1541)-1582)/117; % ceil()

cutToRight = round((startScat - transSliceLocations(1)) * 16); % 16 pixels per fan was +1 in the parenethesis
cutToLeft = round((transSliceLocations(end) - endScat) * 16) - 1;

transCropped = transScan;
transCropped(:,(end-cutToLeft):end) = [];
transCropped(:,1:cutToRight) = [];

transCropped(2861:end,:) = [];     % removing lower rows not used for recon

detPixelsMapping = round((((1:2860) - 1430.5) .* (312/632.6)) .* (0.1*3) + (72.5*3)); % for super-res XRD images
transCropped = mat2gray(transCropped,[prctile(transCropped(:),tBP),prctile(transCropped(:),tTP)]); % normalizing between 0 to 1
scatWidth = size(transCropped,2)/size(meanQMat,2); % used for image to map color on to trans

transMask = repmat(transCropped,1,1,3);

%% Project Mean Q into Detector Space

detMeanQMat = [];
for i = 1 : size(meanQMat,1)    
    numRows = length(find(detPixelsMapping == i)); % pulling rows to give color to 

    detMeanQTemp = [];
    for j = 1 : size(meanQMat,2)

        % calculate correct amount of pixels to pull for decimal amount
        startCol = round(scatWidth * j - scatWidth + 1);
        endCol = round(scatWidth * j);

        detMeanQTemp = [detMeanQTemp, ones(numRows,length(startCol:endCol)).*meanQMat(i,j)];     
    end

    detMeanQMat = [detMeanQMat;detMeanQTemp];
end

%% Create tXRD Images Plot

qGrid = minVal : ((maxVal-minVal)/(115)) : maxVal;

rgb1 = [0 0.9 0.85];
rgb2 = [1 1 0];
rgb3 = [1 0.2 0];
rgb = [interp1([1,size(qGrid,2)/2], [rgb1;rgb2], 1:size(qGrid,2)/2);...
    interp1([1,size(qGrid,2)/2], [rgb2;rgb3], 1:size(qGrid,2)/2)];

[~,index] = min(abs(qGrid-detMeanQMat(:)),[],2);  
rgbReconImg = rgb(index,:);
rgbReconImg = reshape(rgbReconImg,[size(detMeanQMat,1),size(detMeanQMat,2),3]);
rgbTransImg = rgbReconImg .* transMask; 

figure()
imshow(transMask)

figure();
imhandle = imshow(rgbTransImg);    % showing smoothed transmission image
colormap(rgb)
h = colorbar('Location','southoutside');
h.YTick = [0 0.5 1];
h.YTickLabel = {num2str(round(minVal,3)), num2str(round((minVal+maxVal)/2,3)), num2str(round(maxVal,3))};
set(gca,'fontweight','bold')
ylabel(h, 'Mean q  (1/Å)')
set(h,'fontweight','bold','FontSize',20)
brighten(0.3)
