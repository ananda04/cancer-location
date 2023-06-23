%% load Data: t2_177
load('ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat')
load("t2_177.mat")

% Mask t2_177
redChannel = T2177(:, :, 1);
greenChannel = T2177(:, :, 2);
blueChannel = T2177(:, :, 3);

mask = blueChannel > 251, redChannel > 255, greenChannel > 254;
figure(3); imshow(mask)
hold on;

[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude
anomalyfunc(normallspec, 0.7, mask,r,c)

%% load Data: t2_158 -
load('ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat')
load("t2_158.mat")
% Mask t2_187
redChannel = T2158(:, :, 1);
greenChannel = T2158(:, :, 2);
blueChannel = T2158(:, :, 3);

mask = blueChannel > 251, redChannel > 255, greenChannel > 254;
figure(3); imshow(mask)
hold on;

[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude
anomalyfunc(normallspec, 0.7,mask, r,c)

%% load Data: t2_176 -
load('ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat')
load("t2_176.mat")

% Mask t2_187
redChannel = t2176(:, :, 1);
greenChannel = t2176(:, :, 2);
blueChannel = t2176(:, :, 3);

mask = blueChannel > 250, redChannel > 252, greenChannel > 254;
figure(3); imshow(mask)
hold on;

[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude
anomalyfunc(normallspec, 0.7,mask, r,c)

%% load Data: t1_187 -
load('ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat')
load("t1_187.mat")
% Mask t1_187
redChannel = t1187(:, :, 1);
greenChannel = t1187(:, :, 2);
blueChannel = t1187(:, :, 3);

mask = blueChannel > 240, redChannel > 200, greenChannel > 200;
figure(3); imshow(mask)
hold on;

[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude
anomalyfunc(normallspec, 0.68, mask, r,c)


%% load Data: t2_187 -
load('ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat')
load("t2_187.mat")
% Mask t2_187
redChannel = t2187(:, :, 1);
greenChannel = t2187(:, :, 2);
blueChannel = t2187(:, :, 3);

mask = blueChannel > 253, redChannel > 253, greenChannel > 254;
figure(1); imshow(mask)
hold on;

[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude
anomalyfunc(normallspec, 0.7,mask, r,c)

%% load Data: t1_176
load('ReconResults_Brain_176_T1_20s_300iter_M3_Try1.mat')
load("t1_176.mat")
% Mask t1_176
redChannel = T1176(:, :, 1);
greenChannel = T1176(:, :, 2);
blueChannel = T1176(:, :, 3);
    
mask = blueChannel > 248, redChannel > 248, greenChannel > 248;
hold on;
[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude
anomalyfunc(normallspec, 0.7,mask, r,c)

%% load Data: t1_158
load('ReconResults_Brain_158_T1_20s_300iter_M3_Try1.mat')
load("t1_158.mat")
mask = imbinarize(T1158)
[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude
anomalyfunc(normallspec, 0.7,mask, r,c)
