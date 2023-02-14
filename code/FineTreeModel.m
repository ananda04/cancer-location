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
normT1158 = allspec./magnitude
T1158 = testcancer.predictFcn(normT1158);
%% load Data: t1_176
load('ReconResults_Brain_176_T1_20s_300iter_M3_Try1.mat')
load("t1_176.mat")
% Mask t2_177
redChannel = T1176(:, :, 1);
greenChannel = T1176(:, :, 2);
blueChannel = T1176(:, :, 3);

mask = blueChannel > 248, redChannel > 248, greenChannel > 248;
figure(3); imshow(mask)

hold on;[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normT1176 = allspec./magnitude
T1176 = testcancer.predictFcn(normT1176);
%% load Data: t1_187
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
normt1187 = allspec./magnitude
t1187 = testcancer.predictFcn(normt1187);
%% load Data: t2_177
load('ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat')
load("t2_177.mat")

% Mask t2_177
redChannel = T2177(:, :, 1);
greenChannel = T2177(:, :, 2);
blueChannel = T2177(:, :, 3);

mask = blueChannel > 254, redChannel > 254, greenChannel > 254;
figure(3); imshow(mask)
hold on;

[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normT2177 = allspec./magnitude
T2177 = testcancer.predictFcn(normT2177);

%% load Data: t2_187
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
normt2187 = allspec./magnitude
t2187 = testcancer.predictFcn(normt2187);
%% load Data: t2_158
load('ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat')
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
normT2158 = allspec./magnitude
T2158 = testcancer.predictFcn(normT2158);
%% load Data: t2_176
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
normt2176 = allspec./magnitude
t2176 = testcancer.predictFcn(normt2176);
%% load Data: NT_187
    load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
    load("NT_187.mat")
    % Mask NT_187
    redChannel = NT187(:, :, 1);
    greenChannel = NT187(:, :, 2);
    blueChannel = NT187(:, :, 3);
    
    mask4 = blueChannel > 240, redChannel > 240, greenChannel > 240;
    figure(2); imshow(mask4)
    hold on;

    [x y] = find(mask4 == 0)
    l = length(x)
    healallspec = []
    for k1 = 1:l
        healallspec = cat(2, healallspec, squeeze(reconstructedData(x(k1),y(k1),:)))
    end
    magnitude = sqrt(sum(healallspec.^2))
    normNT187 = healallspec./magnitude
    NT187 = testcancer.predictFcn(normNT187);
    %% load Data: NT_158
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')
    load("NT_158.mat")
        % Mask NT_158
    redChannel = NT158(:, :, 1);
    greenChannel = NT158(:, :, 2);
    blueChannel = NT158(:, :, 3);
    
    mask5 = blueChannel > 245, redChannel > 250, greenChannel > 248;
    figure(2); imshow(mask5)
    hold on;
    [x1 y1] = find(mask5 ==0)
    l = length(x1)
    healallspec1 = []
    for k1 = 1:l
        healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(x1(k1),y1(k1),:)))
    end
    magnitude1 = sqrt(sum(healallspec1.^2))
    normNT158 = healallspec1./magnitude1
    NT158 = testcancer.predictFcn(normNT158);
