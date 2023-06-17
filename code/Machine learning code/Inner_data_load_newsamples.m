%% internal data labels for cancerous and noncancerous data
% created by Arnav Nanda, Rising Freshman at Duke University 
%% load in NT Slices and label as zero 
% NT_158
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')
    [x y] = bmask("NT_158")
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
    [x1 y1] = bmask("NT_176")
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
    [x2 y2] = bmask("NT_187")
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
%NT_177
    load("ReconResults_Brain_177_NT_30s_300iter_M3_Try1.mat")
    [x3 y3] = bmask("NT_177")
    l = length(x3)
    healallspec3 = []
    for k = 1:l
        healallspec3 = cat(2, healallspec3, squeeze(reconstructedData(x3(k),y3(k),:)))
    end
    magnitude3 = sqrt(sum(healallspec3.^2))
    normhealth3 = healallspec3./magnitude3
    label0_177 = zeros(length(normhealth3),1)
    label0_177 = transpose(label0_177)
    NT177_resp = [normhealth3,
                  label0_177]
% T1_158
    load('ReconResults_Brain_158_T1_20s_300iter_M3_Try1.mat')
    T1158spec = [squeeze(reconstructedData(56,57,:)),squeeze(reconstructedData(53,59,:)),squeeze(reconstructedData(58,61,:)),squeeze(reconstructedData(61,57,:)),squeeze(reconstructedData(66,56,:)),squeeze(reconstructedData(69,50,:)),squeeze(reconstructedData(58,45,:)),squeeze(reconstructedData(58,47,:)),squeeze(reconstructedData(58,35,:)),squeeze(reconstructedData(61,42,:))]
    magT1158 = sqrt(sum(T1158spec.^2))
    normT1158 = T1158spec./magT1158
    label1_T1158 = ones(10,1)
    label1_T1158 = transpose(label1_T1158)
    T1158_resp = [normT1158,
                  label1_T1158]
% T1_176
    load('ReconResults_Brain_176_T1_20s_300iter_M3_Try1.mat') 
    T1176spec = [squeeze(reconstructedData(88,60,:)),squeeze(reconstructedData(84,57,:)),squeeze(reconstructedData(85,52,:)),squeeze(reconstructedData(84,55,:)),squeeze(reconstructedData(84,62,:)),squeeze(reconstructedData(83,64,:)),squeeze(reconstructedData(83,60,:)),squeeze(reconstructedData(85,53,:)),squeeze(reconstructedData(88,52,:)),squeeze(reconstructedData(82,52,:))]
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
%T1_177
    load("ReconResults_Brain_177_T1_30s_300iter_M3_Try1.mat")
    T1177spec = [squeeze(reconstructedData(90,29,:)),squeeze(reconstructedData(90,28,:)),squeeze(reconstructedData(92,27,:)),squeeze(reconstructedData(87,27,:)),squeeze(reconstructedData(87,24,:)),squeeze(reconstructedData(87,29,:)),squeeze(reconstructedData(90,27,:)),squeeze(reconstructedData(86,27,:)),squeeze(reconstructedData(90,71,:)),squeeze(reconstructedData(86,68,:))]
    magT1177 = sqrt(sum(T1177spec.^2))
    normT1177 = T1177spec./magT1177
    label1_T1177 = ones(10,1)
    label1_T1177 = transpose(label1_T1177)
    T1177_resp = [normT1177,
                  label1_T1177]
% T2_158
    load('ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat')
    T2158spec = [squeeze(reconstructedData(51,18,:)),squeeze(reconstructedData(57,20,:)),squeeze(reconstructedData(58,22,:)),squeeze(reconstructedData(62,24,:)),squeeze(reconstructedData(58,32,:)),squeeze(reconstructedData(61,31,:)),squeeze(reconstructedData(61,21,:)),squeeze(reconstructedData(59,20,:)),squeeze(reconstructedData(57,16,:)),squeeze(reconstructedData(59,14,:))]
    magT2158 = sqrt(sum(T2158spec.^2))
    normT2158 = T2158spec./magT2158
    label1_T2158 = ones(10,1)
    label1_T2158 = transpose(label1_T2158)
    T2158_resp = [normT2158,
                  label1_T2158]
% T2_176
    load('ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat')
    T2176spec = [squeeze(reconstructedData(87,42,:)),squeeze(reconstructedData(87,37,:)),squeeze(reconstructedData(89,37,:)),squeeze(reconstructedData(89,39,:)),squeeze(reconstructedData(89,33,:)),squeeze(reconstructedData(90,35,:)),squeeze(reconstructedData(55,35,:)),squeeze(reconstructedData(59,38,:)),squeeze(reconstructedData(56,39,:)),squeeze(reconstructedData(57,34,:))]
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
% T2_177
    load("ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat")
    T2177spec = [squeeze(reconstructedData(90,30,:)),squeeze(reconstructedData(88,26,:)),squeeze(reconstructedData(86,25,:)),squeeze(reconstructedData(88,25,:)),squeeze(reconstructedData(90,25,:)),squeeze(reconstructedData(89,23,:)),squeeze(reconstructedData(89,22,:)),squeeze(reconstructedData(87,24,:)),squeeze(reconstructedData(87,26,:)),squeeze(reconstructedData(89,25,:))]
    magT2177 = sqrt(sum(T2177spec.^2))
    normT2177 = T2177spec./magT2177
    label1_T2177 = ones(10,1)
    label1_T2177 = transpose(label1_T2177)
    T2177_resp = [normT2177,
                  label1_T2177]
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
% T1177
    NT_T1177 = [NT177_resp,T1177_resp]
% T2177
    NT_T2177 = [NT177_resp,T2177_resp]

%% Run ML Code
% load Models
    load("n_inner_T1158_FSVM.mat")
    load("n_inner_T1158_FT.mat")
    load("n_inner_T1158_WNN.mat")
    load("n_inner_T1176_FSVM.mat")
    load("n_inner_T1176_FT.mat")
    load("n_inner_T1176_WNN.mat")
    load("n_inner_T1177_FSVM.mat")
    load("n_inner_T1177_FT.mat")
    load("n_inner_T1177_WNN.mat")
    load("n_inner_T1187_FSVM.mat")
    load("n_inner_T1187_FT.mat")
    load("n_inner_T1187_WNN.mat")
    load("n_inner_T2158_FSVM.mat")
    load("n_inner_T2158_FT.mat")
    load("n_inner_T2158_WNN.mat")
    load("n_inner_T2176_FSVM.mat")
    load("n_inner_T2176_FT.mat")
    load("n_inner_T2176_WNN.mat")
    load("n_inner_T2177_FSVM.mat")
    load("n_inner_T2177_FT.mat")
    load("n_inner_T2177_WNN.mat")
    load("n_inner_T2187_FSVM.mat")
    load("n_inner_T2187_FT.mat")
    load("n_inner_T2187_WNN.mat")
%% load Data: t1_158
    load('ReconResults_Brain_158_T1_20s_300iter_M3_Try1.mat')
    % Mask t2_177
    [mask r c] = bmask("T1_158")
    L = length(r)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT1158 = allspec./magnitude
    T1158NN = n_inner_t1158_WNN.predictFcn(normT1158);
    T1158FT = n_inner_t1158_FT.predictFcn(normT1158);
    T1158SVM = n_inner_t1158_FSVM.predictFcn(normT1158);
    L = length(normT1158)
    tissue_assigner(T1158NN,mask,r,c, L,1,1, "WNN")
    tissue_assigner(T1158FT,mask,r,c, L,1,2, "FT")
    tissue_assigner(T1158SVM,mask,r,c,L,1,3, "FSVM")
    figure(1); sgtitle("Classificaion of T1158")
%% load Data: t1_176
    load('ReconResults_Brain_176_T1_20s_300iter_M3_Try1.mat')
    % Mask t1_176
    [mask r1 c1] = bmask("T1_176")
    L = length(r1)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r1(k1),c1(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT1176 = allspec./magnitude
    T1176NN = n_inner_T1176_WNN.predictFcn(normT1176);
    T1176FT = n_inner_T1176_FT.predictFcn(normT1176);
    T1176SVM = n_inner_T1176_FSVM.predictFcn(normT1176);
    L = length(normT1176)
    tissue_assigner(T1176NN,mask,r1,c1, L,2,1, "WNN")
    tissue_assigner(T1176FT,mask,r1,c1, L,2,2, "FT")
    tissue_assigner(T1176SVM,mask,r1,c1,L,2,3, "FSVM")
    figure(2); sgtitle("Classificaion of T1176")
%% load Data: t1_187
    load('ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat')
    % Mask t1_187
    [mask r2 c2] = bmask("T1_187")
    L = length(r2)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r2(k1),c2(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT1187 = allspec./magnitude
    T1187NN = n_inner_T1187_WNN.predictFcn(normT1187);
    T1187FT = n_inner_T1187_FT.predictFcn(normT1187);
    T1187SVM = n_inner_T1187_FSVM.predictFcn(normT1187);
    L = length(normT1187)
    tissue_assigner(T1187NN,mask,r2,c2, L,3,1, "WNN")
    tissue_assigner(T1187FT,mask,r2,c2, L,3,2, "FT")
    tissue_assigner(T1187SVM,mask,r2,c2,L,3,3, "FSVM")
    figure(3); sgtitle("Classificaion of t1187")
%% load Data: t1_177
    load("ReconResults_Brain_177_T1_30s_300iter_M3_Try1.mat")
    % Mask t1_177
    [mask r3 c3] = bmask("T1_177")
    L = length(r3)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r3(k1),c3(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT1177 = allspec./magnitude
    T1177NN = n_inner_T1177_WNN.predictFcn(normT1177);
    T1177FT = n_inner_T1177_FT.predictFcn(normT1177);
    T1177SVM = n_inner_T1177_FSVM.predictFcn(normT1177);
    L = length(normT1177)
    tissue_assigner(T1177NN,mask,r3,c3, L,4,1, "WNN")
    tissue_assigner(T1177FT,mask,r3,c3, L,4,2, "FT")
    tissue_assigner(T1177SVM,mask,r3,c3,L,4,3, "FSVM")
    figure(4); sgtitle("Classificaion of t1177")
%% load Data: t2_158
    load('ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat')
    [mask r4 c4] = bmask("T2_158")
    L = length(r4)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r4(k1),c4(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2158 = allspec./magnitude
    T2158NN = n_inner_T2158_WNN.predictFcn(normT2158);
    T2158FT = n_inner_T2158_FT.predictFcn(normT2158);
    T2158SVM = n_inner_T2158_FSVM.predictFcn(normT2158);
    L = length(normT2158)
    tissue_assigner(T2158NN,mask,r4,c4, L,5,1, "WNN")
    tissue_assigner(T2158FT,mask,r4,c4, L,5,2, "FT")
    tissue_assigner(T2158SVM,mask,r4,c4,L,5,3, "FSVM")
    figure(5); sgtitle("Classificaion of t2158")
%% load Data: t2_176
    load('ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat')
    [mask r5 c5] = bmask("T2_176")
    L = length(r5)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r5(k1),c5(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2176 = allspec./magnitude
    T2176NN = n_inner_T2176_WNN.predictFcn(normT2176);
    T2176FT = n_inner_T2176_FT.predictFcn(normT2176);
    T2176SVM = n_inner_T2176_FSVM.predictFcn(normT2176);
    L = length(normT2176)
    tissue_assigner(T2176NN,mask,r5,c5, L,6,1, "WNN")
    tissue_assigner(T2176FT,mask,r5,c5, L,6,2, "FT")
    tissue_assigner(T2176SVM,mask,r5,c5,L,6,3, "FSVM")
    figure(6); sgtitle("Classificaion of t2176")
%% load Data: t2_187
    load('ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat')
    [mask r6 c6] = bmask("T2_187")
    L = length(r6)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r6(k1),c6(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2187 = allspec./magnitude
    T2187NN = n_inner_T2187_WNN.predictFcn(normT2187);
    T2187FT = n_inner_T2187_FT.predictFcn(normT2187);
    T2187SVM = n_inner_T2187_FSVM.predictFcn(normT2187);
    L = length(normT2187)
    tissue_assigner(T2187NN,mask,r6,c6, L,7,1, "NN")
    tissue_assigner(T2187FT,mask,r6,c6, L,7,2, "FT")
    tissue_assigner(T2187SVM,mask,r6,c6,L,7,3, "SVM")
    figure(7); sgtitle("Classificaion of t2187")
%% load Data: t2_177
    load("ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat")
    % Mask t1_177
    [mask r7 c7] = bmask("T2_177")
    L = length(r7)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r7(k1),c7(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2177 = allspec./magnitude
    T2177NN = n_inner_T2177_WNN.predictFcn(normT2177);
    T2177FT = n_inner_T2177_FT.predictFcn(normT2177);
    T2177SVM = n_inner_T2177_FSVM.predictFcn(normT2177);
    L = length(normT2177)
    tissue_assigner(T2177NN,mask,r7,c7, L,8,1, "WNN")
    tissue_assigner(T2177FT,mask,r7,c7, L,8,2, "FT")
    tissue_assigner(T2177SVM,mask,r7,c7,L,8,3, "FSVM")
    figure(8); sgtitle("Classificaion of t2177")