%% Data for NT_176
load('ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat')
load("NT_176.mat")
redChannel = NT176(:, :, 1);
greenChannel = NT176(:, :, 2);
blueChannel = NT176(:, :, 3);
    
mask1 = blueChannel > 245, redChannel > 250, greenChannel > 248;
figure(2); imshow(mask1)
hold on;
[x2 y2] = find(mask1 ==0)
l = length(x2)
healallspec1 = []
for k1 = 1:l
    healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(x2(k1),y2(k1),:)))
end
magnitude1 = sqrt(sum(healallspec1.^2))
normhealth1 = healallspec1./magnitude1
label0 = zeros(length(normhealth1),1)
%% Data for GBMT
load("ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat")
load("cancerblock.mat")
[x y] = find(cancerblock ~= 0)
l = length(x)
canallspec = []
for k1 = 1:l
    canallspec = cat(2, canallspec, squeeze(reconstructedData(x(k1),y(k1),:)))
end
magnitude = sqrt(sum(canallspec.^2))
normcan1 = canallspec./magnitude
%normcan1 = interp1(normcan1, 0:90, "nearest")
%qvals1 = 0.04:0.0098:0.5
qvals2 = 0.04:0.0051:0.5
normcan1 = interp1(qvals, normcan1, qvals2)
label1 = ones(length(normcan1),1)

%%

response =[label0,
           label1];
response = transpose(response)
completenorms = [normhealth1,normcan1, 
                 response]


