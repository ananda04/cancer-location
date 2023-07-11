%% load data
load('BrainData_andErrorBars.mat',"Cancer_Spec")
load("BrainData_andErrorBars.mat","GreyMatter_Spec")
load("BrainData_andErrorBars.mat","WhiteMatter_Spec")
cancerData = [Cancer_Spec, 
              0]
whiteData = [WhiteMatter_Spec,
             1]
grayData = [GreyMatter_errorbars,
            2]
allData = [cancerData,whiteData,grayData]
%% Test model 
load("Cancer_gray_white_FSVM.mat")

% NT_187
    load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
    [mask x2 y2] = bmask("NT_187")
    l = length(x2)
    healallspec = []
    for k = 1:l
        healallspec = cat(2, healallspec, squeeze(reconstructedData(x2(k),y2(k),:)))
    end
    magnitude2 = sqrt(sum(healallspec.^2))
    normNT187 = healallspec./magnitude2
    %%
    NT187SVM = Cancer_gray_white_FSVM.predictFcn(normNT187);
    l = length(normNT187)
    tissue_assigner(NT187SVM,mask,x2,y2,l,1,3, "FSVM")
    figure(1); sgtitle("Classificaion of NT187")
% NT_176
    load('ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat')
    [mask x1 y1] = bmask("NT_176")
    L = length(x1)
    healallspec1 = []
    for j = 1:L
        healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(x1(j),y1(j),:)))
    end
    magnitude1 = sqrt(sum(healallspec1.^2))
    normhealth1 = healallspec1./magnitude1
    NT176FSVM = Cancer_gray_white_FSVM.predictFcn(normhealth1)
    l = length(NT176FSVM)
    tissue_assigner(NT176FSVM,mask,x2,y2,l,1,3, "FSVM")
    figure(1); sgtitle("Classificaion of NT176")
% NT_177
    load("ReconResults_Brain_177_NT_30s_300iter_M3_Try1.mat")
    [mask2 x3 y3] = bmask("NT_177")
    L = length(x3)
    healallspec3 = []
    for k = 1:L
        healallspec3 = cat(2, healallspec3, squeeze(reconstructedData(x3(k),y3(k),:)))
    end
    magnitude3 = sqrt(sum(healallspec3.^2))
    normhealth3 = healallspec3./magnitude3
    NT177FSVM = Cancer_gray_white_FSVM.predictFcn(normhealth3)
    l = length(NT177FSVM)
    tissue_assigner(NT177FSVM,mask,x2,y2,l,1,3, "FSVM")
    figure(1); sgtitle("Classificaion of NT177")
% NT_158
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')
    [mask3 x y] = bmask("NT_158")
    L = length(x)
    healallspec = []
    for i = 1:L
        healallspec = cat(2, healallspec, squeeze(reconstructedData(x(i),y(i),:)))
    end
    magnitude = sqrt(sum(healallspec.^2))
    normhealth = healallspec./magnitude
    NT158FSVM = Cancer_gray_white_FSVM.predictFcn(normhealth)
    l = length(NT158FSVM)
    tissue_assigner(NT158FSVM,mask,x2,y2,l,1,3, "FSVM")
    figure(1); sgtitle("Classificaion of NT158")
% T1_187
    load("ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat")
    [mask x1 y1] = bmask("T1_187")
    l = length(x1)
    healallspec1 = []
    for k = 1:l
        healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(r(k),c(k),:)))
    end
    magnitude2 = sqrt(sum(healallspec1.^2))
    normT1187 = healallspec1./magnitude2
    T1187SVM = Cancer_gray_white_FSVM.predictFcn(normT1187);
    L = length(normT1187)
    tissue_assigner(T1187SVM,mask,r,c,L,1,3, "FSVM")
    figure(1); sgtitle("Classificaion of T1187")

% t2_187
    load("ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat")
    load("T2187_cancer.mat")
    [mask x y] = bmask("T2_187")
    healallspec2 = []
    for k = 1:l
        healallspec2 = cat(2, healallspec2, squeeze(reconstructedData(r(k),c(k),:)))
    end
    magnitude2 = sqrt(sum(healallspec2.^2))
    normT2187 = healallspec2./magnitude2
    T2187SVM = Cancer_gray_white_FSVM.predictFcn(normT2187);
    L = length(normT2187)
    tissue_assigner(T2187SVM,mask,r,c,L,1,3, "FSVM")
    figure(1); sgtitle("Classificaion of T2187")

 % t1_158
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
    T1158SVM = Cancer_gray_white_FSVM.predictFcn(normT1158);
    L = length(normT1158)
    tissue_assigner(T1158SVM,mask,r,c,L,1,3, "FSVM")
    figure(1); sgtitle("Classificaion of T1158")

% t1_176
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
    T1176SVM = Cancer_gray_white_FSVM.predictFcn(normT1176);
    L = length(normT1176)
    tissue_assigner(T1176SVM,mask,r1,c1,L,2,3, "FSVM")
    figure(2); sgtitle("Classificaion of T1176")

% t1_187
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
    T1187SVM = Cancer_gray_white_FSVM.predictFcn(normT1187);
    L = length(normT1187)
    tissue_assigner(T1187SVM,mask,r2,c2,L,3,3, "FSVM")
    figure(3); sgtitle("Classificaion of t1187")

% t1_177
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
    T1177SVM = Cancer_gray_white_FSVM.predictFcn(normT1177);
    L = length(normT1177)
    tissue_assigner(T1177SVM,mask,r3,c3,L,4,3, "FSVM")
    figure(4); sgtitle("Classificaion of t1177")

% t2_158
    load('ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat')
    [mask r4 c4] = bmask("T2_158")
    L = length(r4)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r4(k1),c4(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2158 = allspec./magnitude
    T2158SVM = Cancer_gray_white_FSVM.predictFcn(normT2158);
    L = length(normT2158)
    tissue_assigner(T2158SVM,mask,r4,c4,L,5,3, "FSVM")
    figure(5); sgtitle("Classificaion of t2158")

% t2_187
    load('ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat')
    [mask r6 c6] = bmask("T2_187")
    L = length(r6)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r6(k1),c6(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2187 = allspec./magnitude
    T2187SVM = Cancer_gray_white_FSVM.predictFcn(normT2187);
    L = length(normT2187)
    tissue_assigner(T2187SVM,mask,r6,c6,L,7,3, "SVM")
    figure(7); sgtitle("Classificaion of t2187")

% t2_177
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
    T2177SVM = Cancer_gray_white_FSVM.predictFcn(normT2177);
    L = length(normT2177)
    tissue_assigner(T2177SVM,mask,r7,c7,L,8,3, "FSVM")
    figure(8); sgtitle("Classificaion of t2177")

% t2_176
    load('ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat')
    [mask r5 c5] = bmask("T2_176")
    L = length(r5)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r5(k1),c5(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2176 = allspec./magnitude
    T2176SVM = Cancer_gray_white_FSVM.predictFcn(normT2176);
    L = length(normT2176)
    tissue_assigner(T2176SVM,mask,r5,c5,L,6,3, "FSVM")
    figure(6); sgtitle("Classificaion of t2176")