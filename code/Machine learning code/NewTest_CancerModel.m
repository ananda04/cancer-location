% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

    %load('new_FSVM_model.mat')
    %load('new_WNN_model.mat')
    %load('new_FT_model.mat')
    load('reviewed_cancerModel_FSVM.mat')
    load('reviewed_cancerModel_WNN.mat')
    load('reviewed_cancerModel_FT.mat')
    reviewed_cancerModel_WNN
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
    T1158NN = reviewed_cancerModel_WNN.predictFcn(normT1158);
    T1158FT = reviewed_cancerModel_FT.predictFcn(normT1158);
    T1158SVM = reviewed_cancerModel_FSVM.predictFcn(normT1158);
    L = length(normT1158)
    tissue_assigner(T1158NN,mask,r,c, L,1,1, "WNN")
    tissue_assigner(T1158FT,mask,r,c, L,1,2, "FT")
    tissue_assigner(T1158SVM,mask,r,c,L,1,3, "FSVM")
    figure(1); sgtitle("Classificaion of T1158")
% load Data: t1_176
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
    T1176NN = reviewed_cancerModel_WNN.predictFcn(normT1176);
    T1176FT = reviewed_cancerModel_FT.predictFcn(normT1176);
    T1176SVM = reviewed_cancerModel_FSVM.predictFcn(normT1176);
    L = length(normT1176)
    tissue_assigner(T1176NN,mask,r1,c1, L,2,1, "WNN")
    tissue_assigner(T1176FT,mask,r1,c1, L,2,2, "FT")
    tissue_assigner(T1176SVM,mask,r1,c1,L,2,3, "FSVM")
    figure(2); sgtitle("Classificaion of T1176")
% load Data: t1_187
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
    T1187NN = reviewed_cancerModel_WNN.predictFcn(normT1187);
    T1187FT = reviewed_cancerModel_FT.predictFcn(normT1187);
    T1187SVM = reviewed_cancerModel_FSVM.predictFcn(normT1187);
    L = length(normT1187)
    tissue_assigner(T1187NN,mask,r2,c2, L,3,1, "WNN")
    tissue_assigner(T1187FT,mask,r2,c2, L,3,2, "FT")
    tissue_assigner(T1187SVM,mask,r2,c2,L,3,3, "FSVM")
    figure(3); sgtitle("Classificaion of t1187")
% load Data: t1_177
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
    T1177NN = reviewed_cancerModel_WNN.predictFcn(normT1177);
    T1177FT = reviewed_cancerModel_FT.predictFcn(normT1177);
    T1177SVM = reviewed_cancerModel_FSVM.predictFcn(normT1177);
    L = length(normT1177)
    tissue_assigner(T1177NN,mask,r3,c3, L,4,1, "WNN")
    tissue_assigner(T1177FT,mask,r3,c3, L,4,2, "FT")
    tissue_assigner(T1177SVM,mask,r3,c3,L,4,3, "FSVM")
    figure(4); sgtitle("Classificaion of t1177")
% load Data: t2_158
    load('ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat')
    [mask r4 c4] = bmask("T2_158")
    L = length(r4)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r4(k1),c4(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2158 = allspec./magnitude
    T2158NN = reviewed_cancerModel_WNN.predictFcn(normT2158);
    T2158FT = reviewed_cancerModel_FT.predictFcn(normT2158);
    T2158SVM = reviewed_cancerModel_FSVM.predictFcn(normT2158);
    L = length(normT2158)
    tissue_assigner(T2158NN,mask,r4,c4, L,5,1, "WNN")
    tissue_assigner(T2158FT,mask,r4,c4, L,5,2, "FT")
    tissue_assigner(T2158SVM,mask,r4,c4,L,5,3, "FSVM")
    figure(5); sgtitle("Classificaion of t2158")

% load Data: t2_187
    load('ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat')
    [mask r6 c6] = bmask("T2_187")
    L = length(r6)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r6(k1),c6(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2187 = allspec./magnitude
    T2187NN = reviewed_cancerModel_WNN.predictFcn(normT2187);
    T2187FT = reviewed_cancerModel_FT.predictFcn(normT2187);
    T2187SVM = reviewed_cancerModel_FSVM.predictFcn(normT2187);
    L = length(normT2187)
    tissue_assigner(T2187NN,mask,r6,c6, L,7,1, "NN")
    tissue_assigner(T2187FT,mask,r6,c6, L,7,2, "FT")
    tissue_assigner(T2187SVM,mask,r6,c6,L,7,3, "SVM")
    figure(7); sgtitle("Classificaion of t2187")
% load Data: t2_177
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
    T2177NN = reviewed_cancerModel_WNN.predictFcn(normT2177);
    T2177FT = reviewed_cancerModel_FT.predictFcn(normT2177);
    T2177SVM = reviewed_cancerModel_FSVM.predictFcn(normT2177);
    L = length(normT2177)
    tissue_assigner(T2177NN,mask,r7,c7, L,8,1, "WNN")
    tissue_assigner(T2177FT,mask,r7,c7, L,8,2, "FT")
    tissue_assigner(T2177SVM,mask,r7,c7,L,8,3, "FSVM")
    figure(8); sgtitle("Classificaion of t2177")
% load Data: t2_176
    load('ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat')
    [mask r5 c5] = bmask("T2_176")
    L = length(r5)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r5(k1),c5(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2176 = allspec./magnitude
    T2176NN = reviewed_cancerModel_WNN.predictFcn(normT2176);
    T2176FT = reviewed_cancerModel_FT.predictFcn(normT2176);
    T2176SVM = reviewed_cancerModel_FSVM.predictFcn(normT2176);
    L = length(normT2176)
    tissue_assigner(T2176NN,mask,r5,c5, L,6,1, "WNN")
    tissue_assigner(T2176FT,mask,r5,c5, L,6,2, "FT")
    tissue_assigner(T2176SVM,mask,r5,c5,L,6,3, "FSVM")
    figure(6); sgtitle("Classificaion of t2176")
