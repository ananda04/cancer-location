%% internal data labels for cancerous and noncancerous data
% created by Arnav Nanda, Rising Freshman at Duke University 
%% load in NT Slices and label as zero 
% NT_158
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')
    load("NT_158.mat")
    redChannel = NT158(:, :, 1);
    greenChannel = NT158(:, :, 2);
    blueChannel = NT158(:, :, 3);
        
    mask = blueChannel > 245, redChannel > 250, greenChannel > 248;
    figure(2); imshow(mask)
    hold on;
    [x y] = find(mask ==0)
    l = length(x)
    healallspec = []
    for i = 1:l
        healallspec = cat(2, healallspec, squeeze(reconstructedData(x(i),y(i),:)))
    end
    magnitude = sqrt(sum(healallspec.^2))
    normhealth = healallspec./magnitude
    label0_158 = zeros(length(normhealth),1)
    label0_158 = transpose(label0_158)
    NT158_resp = [normhealth,
                  label0_158]
% NT_176
    load('ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat')
    load("NT_176.mat")
    redChannel = NT176(:, :, 1);
    greenChannel = NT176(:, :, 2);
    blueChannel = NT176(:, :, 3);
        
    mask1 = blueChannel > 245, redChannel > 250, greenChannel > 248;
    figure(2); imshow(mask1)
    hold on;
    [x1 y1] = find(mask1 ==0)
    l = length(x1)
    healallspec1 = []
    for j = 1:l
        healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(x1(j),y1(j),:)))
    end
    magnitude1 = sqrt(sum(healallspec1.^2))
    normhealth1 = healallspec1./magnitude1
    label0_176 = zeros(length(normhealth1),1)
    label0_176 = transpose(label0_176)
    NT176_resp = [normhealth1,
                  label0_176]
% NT_187
    load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
    load("NT_187.mat")
    redChannel = NT187(:, :, 1);
    greenChannel = NT187(:, :, 2);
    blueChannel = NT187(:, :, 3);
        
    mask2 = blueChannel > 245, redChannel > 250, greenChannel > 248;
    figure(2); imshow(mask2)
    hold on;
    [x2 y2] = find(mask2 ==0)
    l = length(x2)
    healallspec2 = []
    for k = 1:l
        healallspec2 = cat(2, healallspec2, squeeze(reconstructedData(x2(k),y2(k),:)))
    end
    magnitude2 = sqrt(sum(healallspec2.^2))
    normhealth2 = healallspec2./magnitude2
    label0_187 = zeros(length(normhealth2),1)
    label0_187 = transpose(label0_187)
    NT187_resp = [normhealth2,
                  label0_187]
% T1_158
    load('ReconResults_Brain_158_T1_20s_300iter_M3_Try1.mat')
    T1158spec = [squeeze(reconstructedData(88,28,:)),squeeze(reconstructedData(90,36,:)),squeeze(reconstructedData(93,28,:)),squeeze(reconstructedData(91,25,:)),squeeze(reconstructedData(86,30,:)),squeeze(reconstructedData(87,28,:)),squeeze(reconstructedData(89,33,:)),squeeze(reconstructedData(88,24,:)),squeeze(reconstructedData(94,21,:)),squeeze(reconstructedData(89,20,:))]
    magT1158 = sqrt(sum(T1158spec.^2))
    normT1158 = T1158spec./magT1158
    label1_T1158 = ones(10,1)
    label1_T1158 = transpose(label1_T1158)
    T1158_resp = [normT1158,
                  label1_T1158]
% T1_176
    load('ReconResults_Brain_176_T1_20s_300iter_M3_Try1.mat') 
    T1176spec = [squeeze(reconstructedData(48,44,:)),squeeze(reconstructedData(47,48,:)),squeeze(reconstructedData(51,37,:)),squeeze(reconstructedData(54,27,:)),squeeze(reconstructedData(89,56,:)),squeeze(reconstructedData(89,45,:)),squeeze(reconstructedData(45,58,:)),squeeze(reconstructedData(47,55,:)),squeeze(reconstructedData(46,46,:)),squeeze(reconstructedData(48,42,:))]
    magT1176 = sqrt(sum(T1176spec.^2))
    normT1176 = T1176spec./magT1176
    label1_T1176 = ones(10,1)
    label1_T1176 = transpose(label1_T1176)
    T1176_resp = [normT1176,
                  label1_T1176]
% T1_187
    load('ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat')
    T1187spec = [squeeze(reconstructedData(82,26,:)),squeeze(reconstructedData(86,29,:)),squeeze(reconstructedData(82,29,:)),squeeze(reconstructedData(82,27,:)),squeeze(reconstructedData(75,27,:)),squeeze(reconstructedData(77,26,:)),squeeze(reconstructedData(70,26,:)),squeeze(reconstructedData(70,24,:)),squeeze(reconstructedData(71,22,:)),squeeze(reconstructedData(72,32,:))]
    magT1187 = sqrt(sum(T1187spec.^2))
    normT1187 = T1187spec./magT1187
    label1_T1187 = ones(10,1)
    label1_T1187 = transpose(label1_T1187)
    T1187_resp = [normT1187,
                  label1_T1187]
% T2_158
    load('ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat')
    T2158spec = [squeeze(reconstructedData(81,17,:)),squeeze(reconstructedData(74,16,:)),squeeze(reconstructedData(74,19,:)),squeeze(reconstructedData(83,14,:)),squeeze(reconstructedData(84,16,:)),squeeze(reconstructedData(82,14,:)),squeeze(reconstructedData(79,13,:)),squeeze(reconstructedData(79,15,:)),squeeze(reconstructedData(72,18,:)),squeeze(reconstructedData(76,22,:))]
    magT2158 = sqrt(sum(T2158spec.^2))
    normT2158 = T2158spec./magT2158
    label1_T2158 = ones(10,1)
    label1_T2158 = transpose(label1_T2158)
    T2158_resp = [normT2158,
                  label1_T2158]
% T2_176
    load('ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat')
    T2176spec = [squeeze(reconstructedData(70,48,:)),squeeze(reconstructedData(70,35,:)),squeeze(reconstructedData(70,52,:)),squeeze(reconstructedData(66,52,:)),squeeze(reconstructedData(72,51,:)),squeeze(reconstructedData(74,45,:)),squeeze(reconstructedData(64,48,:)),squeeze(reconstructedData(69,45,:)),squeeze(reconstructedData(68,44,:)),squeeze(reconstructedData(67,44,:))]
    magT2176 = sqrt(sum(T2176spec.^2))
    normT2176 = T2176spec./magT2176
    label1_T2176 = ones(10,1)
    label1_T2176 = transpose(label1_T2176)
    T2176_resp = [normT2176,
                  label1_T2176]
% T2_187
    load('ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat')
    T2187spec = [squeeze(reconstructedData(77,35,:)),squeeze(reconstructedData(73,20,:)),squeeze(reconstructedData(81,25,:)),squeeze(reconstructedData(74,26,:)),squeeze(reconstructedData(69,29,:)),squeeze(reconstructedData(77,43,:)),squeeze(reconstructedData(84,43,:)),squeeze(reconstructedData(84,32,:)),squeeze(reconstructedData(80,29,:)),squeeze(reconstructedData(84,26,:))]
    magT2187 = sqrt(sum(T2187spec.^2))
    normT2187 = T2187spec./magT2187
    label1_T2187 = ones(10,1)
    label1_T2187 = transpose(label1_T2187)
    T2187_resp = [normT2187,
                  label1_T2187]

%% different combos
% T1158
    NT_T1158 = [NT158_resp,T1158_resp]
%T2158
    NT_T2158 = [NT158_resp,T2158_resp]
% T1176
    NT_T1176 = [NT176_resp,T1176_resp]
% T2176
    NT_T2176 = [NT176_resp,T2176_resp]
% T1187
    NT_T1187 = [NT187_resp,T1187_resp]
% T2187
    NT_T2187 = [NT187_resp,T2187_resp]