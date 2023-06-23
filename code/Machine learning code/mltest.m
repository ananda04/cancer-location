%% Data for NT_187
load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
load("NT_187.mat")
redChannel = NT187(:, :, 1);
greenChannel = NT187(:, :, 2);
blueChannel = NT187(:, :, 3);
    
mask = blueChannel > 240, redChannel > 240, greenChannel > 240;
figure(2); imshow(mask)
hold on;

[r c] = find(mask == 0)
l = length(r)
healallspec = []
for k1 = 1:l
    healallspec = cat(2, healallspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(healallspec.^2))
normtest187 = healallspec./magnitude

%% Data for Normal
load('ReconResults_Model_3_BrainScan_Normal_FullFanExtent_400iter_M3.mat')
load('u1.mat')
redChannel = u1(:, :, 1);
greenChannel = u1(:, :, 2);
blueChannel = u1(:, :, 3);

mask2 = blueChannel > 250,  greenChannel > 250;
figure(4); imshow(mask2)
hold on;

[x y] = find(mask2 == 0)
l = length(x)
healallspec2 = []
for k1 = 1:l
    healallspec2 = cat(2, healallspec2, squeeze(reconstructedData(x(k1),y(k1),:)))
end
magnitude2 = sqrt(sum(healallspec2.^2))
normhealth2 = healallspec2./magnitude2
normtest = interp1(normhealth2,0:90)



%% load Data: t1_158
load('ReconResults_Brain_158_T1_20s_300iter_M3_Try1.mat')
load("t1_158.mat")
% Mask t2_177
redChannel = T1158(:, :, 1);
greenChannel = T1158(:, :, 2);
blueChannel = T1158(:, :, 3);

mask = blueChannel > 250, redChannel > 250, greenChannel > 250;
figure(3); imshow(mask)
[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
norm158can = allspec./magnitude

