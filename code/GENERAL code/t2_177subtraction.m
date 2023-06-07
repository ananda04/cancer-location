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