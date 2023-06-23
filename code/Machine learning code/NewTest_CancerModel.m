    load('new_FSVM_model.mat')
    load('new_WNN_model.mat')
    load('new_FT_model.mat')
%% load Data: t1_158
    load('ReconResults_Brain_158_T1_20s_300iter_M3_Try1.mat')
    load("t1_158.mat")
    % Mask t2_177
    redChannel = T1158(:, :, 1);
    greenChannel = T1158(:, :, 2);
    blueChannel = T1158(:, :, 3);
    
    mask = blueChannel > 250, redChannel > 250, greenChannel > 250;
    [r c] = find(mask == 0)
    L = length(r)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT1158 = allspec./magnitude
    T1158NN = new_cancermodel_WNN.predictFcn(normT1158);
    T1158FT = new_cancermodel_FT.predictFcn(normT1158);
    T1158SVM = new_cancermodel_FSVM.predictFcn(normT1158);
    L = length(normT1158)
    tissue_assigner(T1158NN,mask,r,c, L,1,1, "NN")
    tissue_assigner(T1158FT,mask,r,c, L,1,2, "FT")
    tissue_assigner(T1158SVM,mask,r,c,L,1,3, "SVM")
    figure(1); sgtitle("Classificaion of t1158")
%% load Data: t1_176
    load('ReconResults_Brain_176_T1_20s_300iter_M3_Try1.mat')
    load("t1_176.mat")
    % Mask t1_176
    redChannel = T1176(:, :, 1);
    greenChannel = T1176(:, :, 2);
    blueChannel = T1176(:, :, 3);
    
    mask1 = blueChannel > 248, redChannel > 248, greenChannel > 248;
    
    hold on;
    [r1 c1] = find(mask1 == 0)
    L1 = length(r1)
    allspec = []
    for k1 = 1:L1
        allspec = cat(2, allspec, squeeze(reconstructedData(r1(k1),c1(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT1176 = allspec./magnitude
    T1176NN = new_cancermodel_WNN.predictFcn(normT1176);
    T1176FT = new_cancermodel_FT.predictFcn(normT1176);
    T1176SVM = new_cancermodel_FSVM.predictFcn(normT1176);
    hold on;
    L1 = length(normT1176)
    tissue_assigner(T1176NN,mask1,r1,c1, L1,2,1, "NN")
    tissue_assigner(T1176FT,mask1,r1,c1, L1,2,2, "FT")
    tissue_assigner(T1176SVM,mask1,r1,c1,L1,2,3, "SVM")
    figure(2); sgtitle("Classificaion of t1176")
%% load Data: t2_177
    load('ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat')
    load("t2_177.mat")
    
    % Mask t2_177
    redChannel = T2177(:, :, 1);
    greenChannel = T2177(:, :, 2);
    blueChannel = T2177(:, :, 3);
    
    mask3 = blueChannel > 254, redChannel > 250, greenChannel > 254;
    
    % deleting extraneous variables  (black lines that shouldnt be there)
    [r3 c3] = find(mask3 == 0)
    for k1 = 1:length(r3)
        mask3(r3(k1),80)=1
    end
    % reconfiguring slice
    
    hold on;
    [r3 c3] = find(mask3 == 0)
    
    L = length(r3)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r3(k1),c3(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2177 = allspec./magnitude
    T2177NN = new_cancermodel_WNN.predictFcn(normT2177);
    T2177FT = new_cancermodel_FT.predictFcn(normT2177);
    T2177SVM = new_cancermodel_FSVM.predictFcn(normT2177);
    hold on;
    L3 = length(normT2177)
    tissue_assigner(T2177NN,mask3,r3,c3, L3,4,1, "NN")
    tissue_assigner(T2177FT,mask3,r3,c3, L3,4,2, "FT")
    tissue_assigner(T2177SVM,mask3,r3,c3,L3,4,3, "SVM")
    figure(4); sgtitle("Classificaion of T2177")
%% load Data: t2_187
    load('ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat')
    load("t2_187.mat")
    
    % Mask t2_187
    redChannel = t2187(:, :, 1);
    greenChannel = t2187(:, :, 2);
    blueChannel = t2187(:, :, 3);
    
    mask4 = blueChannel > 253, redChannel > 253, greenChannel > 254;
    hold on;
    
    [r4 c4] = find(mask4 == 0)
    L = length(r4)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r4(k1),c4(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normt2187 = allspec./magnitude
    t2187NN = new_cancermodel_WNN.predictFcn(normt2187);
    t2187FT = new_cancermodel_FT.predictFcn(normt2187);
    t2187SVM = new_cancermodel_FSVM.predictFcn(normt2187);
    hold on;
    L4 = length(normt2187)
    tissue_assigner(t2187NN,mask4,r4,c4, L4,5,1, "NN")
    tissue_assigner(t2187FT,mask4,r4,c4, L4,5,2, "FT")
    tissue_assigner(t2187SVM,mask4,r4,c4,L4,5,3, "SVM")
    figure(5); sgtitle("Classificaion of t2187")
% load Data: t2_158
    load('ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat')
    load("t2_158.mat")
    
    % Mask t2_187
    redChannel = T2158(:, :, 1);
    greenChannel = T2158(:, :, 2);
    blueChannel = T2158(:, :, 3);
    
    mask5 = blueChannel > 251, redChannel > 255, greenChannel > 254;
    hold on;
    
    [r5 c5] = find(mask5 == 0)
    L = length(r5)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r5(k1),c5(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2158 = allspec./magnitude
    T2158NN = new_cancermodel_WNN.predictFcn(normT2158);
    T2158FT = new_cancermodel_FT.predictFcn(normT2158);
    T2158SVM = new_cancermodel_FSVM.predictFcn(normT2158);
    hold on;
    L5 = length(normT2158)
    tissue_assigner(T2158NN,mask5,r5,c5, L5,6,1, "NN")
    tissue_assigner(T2158FT,mask5,r5,c5, L5,6,2, "FT")
    tissue_assigner(T2158SVM,mask5,r5,c5,L5,6,3, "SVM")
    figure(6); sgtitle("Classificaion of T2158")
% load Data: t2_176
    load('ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat')
    load("t2_176.mat")
    
    % Mask t2_187
    redChannel = t2176(:, :, 1);
    greenChannel = t2176(:, :, 2);
    blueChannel = t2176(:, :, 3);
    
    mask6 = blueChannel > 250, redChannel > 252, greenChannel > 254;
    hold on;
    
    [r6 c6] = find(mask6 == 0)
    L = length(r6)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r6(k1),c6(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normt2176 = allspec./magnitude
    t2176NN = new_cancermodel_WNN.predictFcn(normt2176);
    t2176FT = new_cancermodel_FT.predictFcn(normt2176);
    t2176SVM = new_cancermodel_FSVM.predictFcn(normt2176);
    hold on;
    L6 = length(normt2176)
    tissue_assigner(t2176NN,mask6,r6,c6, L6,7,1, "NN")
    tissue_assigner(t2176FT,mask6,r6,c6, L6,7,2, "FT")
    tissue_assigner(t2176SVM,mask6,r6,c6,L6,7,3, "SVM")
    figure(7); sgtitle("Classificaion of t2176")
% load Data: NT_187
    load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
    load("NT_187.mat")
    % Mask NT_187
    redChannel = NT187(:, :, 1);
    greenChannel = NT187(:, :, 2);
    blueChannel = NT187(:, :, 3);
    
    mask7 = blueChannel > 240, redChannel > 240, greenChannel > 240;
    hold on;

    [r7 c7] = find(mask7 == 0)
    l = length(r7)
    healallspec = []
    for k1 = 1:l
        healallspec = cat(2, healallspec, squeeze(reconstructedData(r7(k1),c7(k1),:)))
    end
    magnitude = sqrt(sum(healallspec.^2))
    normNT187 = healallspec./magnitude
    NT187NN = new_cancermodel_WNN.predictFcn(normNT187);
    NT187FT = new_cancermodel_FT.predictFcn(normNT187);
    NT187SVM = new_cancermodel_FSVM.predictFcn(normNT187);
    hold on;
    L7 = length(normNT187)
    tissue_assigner(NT187NN,mask7,r7,c7, L7,8,1, "NN")
    tissue_assigner(NT187FT,mask7,r7,c7, L7,8,2, "FT")
    tissue_assigner(NT187SVM,mask7,r7,c7,L7,8,3, "SVM")
    figure(8); sgtitle("Classificaion of NT187")
% load Data: NT_158
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')
    load("NT_158.mat")
        % Mask NT_158
    redChannel = NT158(:, :, 1);
    greenChannel = NT158(:, :, 2);
    blueChannel = NT158(:, :, 3);
    
    mask8 = blueChannel > 245, redChannel > 250, greenChannel > 248;
    hold on;
    [r8 c8] = find(mask8 ==0)
    l = length(r8)
    healallspec1 = []
    for k1 = 1:l
        healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(r8(k1),c8(k1),:)))
    end
    magnitude1 = sqrt(sum(healallspec1.^2))
    normNT158 = healallspec1./magnitude1
    NT158NN = new_cancermodel_WNN.predictFcn(normNT158);
    NT158FT = new_cancermodel_FT.predictFcn(normNT158);
    NT158SVM = new_cancermodel_FSVM.predictFcn(normNT158);
    hold on;
    L8 = length(normNT158)
    tissue_assigner(NT158NN,mask8,r8,c8, L8,9,1, "NN")
    tissue_assigner(NT158FT,mask8,r8,c8, L8,9,2, "FT")
    tissue_assigner(NT158SVM,mask8,r8,c8,L8,9,3, "SVM")
    figure(9); sgtitle("Classificaion of NT158")