% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

%% load data
load('BrainData_andErrorBars.mat',"Cancer_Spec")
load("BrainData_andErrorBars.mat","GreyMatter_Spec")
load("BrainData_andErrorBars.mat","WhiteMatter_Spec")
Cancer_Spec = Cancer_Spec - 0.164382933871042
magnitudeB = sqrt(sum(Cancer_Spec.^2))
Cancer_Spec = Cancer_Spec/magnitudeB
cancerData = [Cancer_Spec, 
              0]
WhiteMatter_Spec = WhiteMatter_Spec - 0.190099192902471
magnitudeB = sqrt(sum(WhiteMatter_Spec.^2))
WhiteMatter_Spec = WhiteMatter_Spec/magnitudeB
whiteData = [WhiteMatter_Spec,
             1]
GreyMatter_Spec = GreyMatter_Spec - 0.172140096190729
magnitudeB = sqrt(sum(GreyMatter_Spec.^2))
GreyMatter_Spec = GreyMatter_Spec/magnitudeB
grayData = [GreyMatter_Spec,
            2]
allData = [cancerData,whiteData,grayData]
%% Test model 
load("bruker_CSVM.mat")
load("bruker_FSVM.mat")
load("bruker_MSVM.mat")
%% NT_187
    load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
    [mask r c] = bmask("NT_187")
    l = length(r)
    healallspec = []
    for k = 1:l
        healallspec = cat(2, healallspec, squeeze(reconstructedData(r(k),c(k),:)))
    end
    magnitude2 = sqrt(sum(healallspec.^2))
    normNT187 = healallspec./magnitude2
    qvals2 = 0.05:0.45/460:0.5
    normNT187 =interp1(qvals,normNT187,qvals2)
    NT187CSVM = bruker_CSVM.predictFcn(normNT187);
    NT187FSVM = bruker_FSVM.predictFcn(normNT187);
    NT187MSVM = bruker_CSVM.predictFcn(normNT187);
    l = length(normNT187)
    tissue_assigner(NT187CSVM,mask,r,c,L,1,1, "CSVM")
    tissue_assigner(NT187FSVM,mask,r,c,L,1,2, "FSVM")
    tissue_assigner(NT187MSVM,mask,r,c,L,1,3, "MSVM")
    %hold on;
    figure(1); sgtitle("Classificaion of NT187")
%% NT_176
    load('ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat')
    [mask r c] = bmask("NT_176")
    L = length(r)
    healallspec1 = []
    for j = 1:L
        healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(r(j),c(j),:)))
    end
    magnitude1 = sqrt(sum(healallspec1.^2))
    normhealth1 = healallspec1./magnitude1
     qvals2 = 0.05:0.45/460:0.5
    normhealth1 =interp1(qvals,normhealth1,qvals2)
    NT176CSVM = bruker_CSVM.predictFcn(normhealth1)
    NT176FSVM = bruker_FSVM.predictFcn(normhealth1)
    NT176MSVM = bruker_MSVM.predictFcn(normhealth1)

    l = length(NT176FSVM)
    tissue_assigner(NT176CSVM,mask,r,c,L,2,1, "CSVM")
    tissue_assigner(NT176FSVM,mask,r,c,L,2,2, "FSVM")
    tissue_assigner(NT176MSVM,mask,r,c,L,2,3, "MSVM")
    %hold on;
    figure(2); sgtitle("Classificaion of NT176")
%% NT_177
    load("ReconResults_Brain_177_NT_30s_300iter_M3_Try1.mat")
    [mask r c] = bmask("NT_177")
    L = length(r)
    healallspec3 = []
    for k = 1:L
        healallspec3 = cat(2, healallspec3, squeeze(reconstructedData(r(k),c(k),:)))
    end
    magnitude3 = sqrt(sum(healallspec3.^2))
    normhealth3 = healallspec3./magnitude3
     qvals2 = 0.05:0.45/460:0.5
    normhealth3 =interp1(qvals,normhealth3,qvals2)
    NT177CSVM = bruker_CSVM.predictFcn(normhealth3)
    NT177FSVM = bruker_FSVM.predictFcn(normhealth3)
    NT177MSVM = bruker_MSVM.predictFcn(normhealth3)
    l = length(NT177FSVM)
    tissue_assigner(NT177CSVM,mask,r,c,L,3,1, "CSVM")
    tissue_assigner(NT177FSVM,mask,r,c,L,3,2, "FSVM")
    tissue_assigner(NT177MSVM,mask,r,c,L,3,3, "MSVM")
    %hold on;
    figure(3); sgtitle("Classificaion of NT177")
%% NT_158
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')
    [mask r c] = bmask("NT_158")
    L = length(r)
    healallspec = []
    for i = 1:L
        healallspec = cat(2, healallspec, squeeze(reconstructedData(r(i),c(i),:)))
    end
    magnitude = sqrt(sum(healallspec.^2))
    normhealth = healallspec./magnitude
     qvals2 = 0.05:0.45/460:0.5
    normhealth =interp1(qvals,normhealth,qvals2)
    NT158CSVM = bruker_CSVM.predictFcn(normhealth)
    NT158FSVM = bruker_FSVM.predictFcn(normhealth)
    NT158MSVM = bruker_MSVM.predictFcn(normhealth)
    l = length(NT158FSVM)
    tissue_assigner(NT158CSVM,mask,r,c,L,4,1, "CSVM")
    tissue_assigner(NT158FSVM,mask,r,c,L,4,2, "FSVM")
    tissue_assigner(NT158MSVM,mask,r,c,L,4,3, "MSVM")
    hold on;
    figure(4);sgtitle("Classificaion of NT158")

 %% t1_158
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
     qvals2 = 0.05:0.45/460:0.5
    normT1158 =interp1(qvals,normT1158,qvals2)
    T1158CSVM = bruker_CSVM.predictFcn(normT1158);
    T1158FSVM = bruker_FSVM.predictFcn(normT1158);
    T1158MSVM = bruker_MSVM.predictFcn(normT1158);
    L = length(normT1158)
    tissue_assigner(T1158CSVM,mask,r,c,L,5,1, "CSVM")
    tissue_assigner(T1158FSVM,mask,r,c,L,5,2, "FSVM")
    tissue_assigner(T1158MSVM,mask,r,c,L,5,3, "MSVM")
    hold on;
    figure(5); sgtitle("Classificaion of T1158")

%% t1_176
    load('ReconResults_Brain_176_T1_20s_300iter_M3_Try1.mat')
    % Mask t1_176
    [mask r c] = bmask("T1_176")
    L = length(r)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT1176 = allspec./magnitude
     qvals2 = 0.05:0.45/460:0.5
    normT1176 =interp1(qvals,normT1176,qvals2)
    T1176CSVM = bruker_CSVM.predictFcn(normT1176);
    T1176FSVM = bruker_FSVM.predictFcn(normT1176);
    T1176MSVM = bruker_MSVM.predictFcn(normT1176);
    L = length(normT1176)
    tissue_assigner(T1176CSVM,mask,r,c,L,6,1, "CSVM")
    tissue_assigner(T1176FSVM,mask,r,c,L,6,2, "FSVM")
    tissue_assigner(T1176MSVM,mask,r,c,L,6,3, "MSVM")
    hold on;
    figure(8); sgtitle("Classificaion of T1176")

% t1_187
    load('ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat')
    % Mask t1_187
    [mask r c] = bmask("T1_187")
    L = length(r)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT1187 = allspec./magnitude
    qvals2 = 0.05:0.45/460:0.5
    normT1187 =interp1(qvals,normT1187,qvals2)
    T1187CSVM = bruker_CSVM.predictFcn(normT1187);
    T1187FSVM = bruker_FSVM.predictFcn(normT1187);
    T1187MSVM = bruker_MSVM.predictFcn(normT1187);
    L = length(normT1187)
    tissue_assigner(T1187CSVM,mask,r,c,L,7,1, "CSVM")
    tissue_assigner(T1187FSVM,mask,r,c,L,7,2, "FSVM")
    tissue_assigner(T1187MSVM,mask,r,c,L,7,3, "MSVM")
    hold on;
    figure(9); sgtitle("Classificaion of t1187")

% t1_177
    load("ReconResults_Brain_177_T1_30s_300iter_M3_Try1.mat")
    % Mask t1_177
    [mask r c] = bmask("T1_177")
    L = length(r)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT1177 = allspec./magnitude
    qvals2 = 0.05:0.45/460:0.5
    normT1177 =interp1(qvals,normT1177,qvals2)
    T1177CSVM = bruker_CSVM.predictFcn(normT1177);
    T1177FSVM = bruker_FSVM.predictFcn(normT1177);
    T1177MSVM = bruker_MSVM.predictFcn(normT1177);
    L = length(normT1177)
    tissue_assigner(T1177CSVM,mask,r,c,L,8,1, "CSVM")
    tissue_assigner(T1177FSVM,mask,r,c,L,8,2, "FSVM")
    tissue_assigner(T1177MSVM,mask,r,c,L,8,3, "MSVM")
    hold on;
    figure(10); sgtitle("Classificaion of t1177")

 %t2_158
    load('ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat')
    [mask r c] = bmask("T2_158")
    L = length(r)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2158 = allspec./magnitude
    qvals2 = 0.05:0.45/460:0.5
    normT2158 =interp1(qvals,normT2158,qvals2)
    T2158CSVM = bruker_CSVM.predictFcn(normT2158);
    T2158FSVM = bruker_FSVM.predictFcn(normT2158);
    T2158MSVM = bruker_MSVM.predictFcn(normT2158);
    L = length(normT2158)
    tissue_assigner(T2158CSVM,mask,r,c,L,9,1, "CSVM")
    tissue_assigner(T2158FSVM,mask,r,c,L,9,2, "FSVM")
    tissue_assigner(T2158MSVM,mask,r,c,L,9,3, "MSVM")
    hold on;
    figure(11); sgtitle("Classificaion of t2158")

% t2_187
    load('ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat')
    [mask r c] = bmask("T2_187")
    L = length(r)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2187 = allspec./magnitude
    qvals2 = 0.05:0.45/460:0.5
    normT2187 =interp1(qvals,normT2187,qvals2)
    T2187CSVM = bruker_CSVM.predictFcn(normT2187);
    T2187FSVM = bruker_FSVM.predictFcn(normT2187);
    T2187MSVM = bruker_MSVM.predictFcn(normT2187);
    L = length(normT2187)
    tissue_assigner(T2187CSVM,mask,r,c,L,10,1, "CSVM")
    tissue_assigner(T2187FSVM,mask,r,c,L,10,2, "FSVM")
    tissue_assigner(T2187MSVM,mask,r,c,L,10,3, "MSVM")
    hold on;
    figure(7); sgtitle("Classificaion of t2187")

% t2_177
    load("ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat")
    % Mask t1_177
    [mask r c] = bmask("T2_177")
    L = length(r7)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2177 = allspec./magnitude
    qvals2 = 0.05:0.45/460:0.5
    normT2177 =interp1(qvals,normT2177,qvals2)
    T2177CSVM = bruker_CSVM.predictFcn(normT2177);
    T2177FSVM = bruker_FSVM.predictFcn(normT2177);
    T2177MSVM = bruker_MSVM.predictFcn(normT2177);
    L = length(normT2177)
    tissue_assigner(T2177CSVM,mask,r,c,L,11,1, "CSVM")
    tissue_assigner(T2177FSVM,mask,r,c,L,11,2, "FSVM")
    tissue_assigner(T2177MSVM,mask,r,c,L,11,3, "MSVM")
    hold on;
    figure(13); sgtitle("Classificaion of t2177")

% t2_176
    load('ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat')
    [mask r c] = bmask("T2_176")
    L = length(r)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normT2176 = allspec./magnitude
    qvals2 = 0.05:0.45/460:0.5
    normT2176 =interp1(qvals,normT2176,qvals2)
    T2176CSVM = bruker_CSVM.predictFcn(normT2176);
    T2176FSVM = bruker_FSVM.predictFcn(normT2176);
    T2176MSVM = bruker_MSVM.predictFcn(normT2176);
    L = length(normT2176)
    tissue_assigner(T2176CSVM,mask,r,c,L,12,1, "CSVM")
    tissue_assigner(T2176FSVM,mask,r,c,L,12,2, "FSVM")
    tissue_assigner(T2176MSVM,mask,r,c,L,12,3, "MSVM")
    hold on;
    figure(14); sgtitle("Classificaion of t2176")
