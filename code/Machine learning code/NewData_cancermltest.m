% Description: only new data classificcation trainer
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
L = length(x2)
healallspec1 = []
for k1 = 1:L
    healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(x2(k1),y2(k1),:)))
end
magnitude1 = sqrt(sum(healallspec1.^2))
normhealth1 = healallspec1./magnitude1
label0 = zeros(length(normhealth1),1)
%% data for cancerous slice (t1187)
load("ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat")
load("T1187_data.mat")
l = length(r)
canallspec = []
for k = 1:l
    canallspec = cat(2, canallspec, squeeze(reconstructedData(r(k),c(k),:)))
end
magnitude = sqrt(sum(canallspec.^2))
normcan1 = canallspec./magnitude
label1 = ones(length(normcan1),1)
% labeled data
response =[label0,
           label1];
response = transpose(response)
completenorms = [normhealth1,normcan1, 
                 response]