load("T1187_data.mat")
load("ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat")
normT1187 = ReconResults(reconstructedData,r,c)
T1187_samp = datasample(normT1187,10,2)

load("ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat")
[mask r c] = bmask("NT_187")
normNT187 = ReconResults(reconstructedData,r,c)
%NT187_samp = datasample(normNT187,20,2)
T1187_label = [T1187_samp, 
                transpose(ones(10,1))]
NT187_label = [normNT187, 
                transpose(zeros(length(normNT187),1))]
allData = [NT187_label,T1187_label]
%%
%load("WNN_10p.mat")
%load("MSVM_10p.mat")
%load("FKNN_10p.mat")
%load("FKNN_2010.mat")
%load("MSVM_2010.mat")
%load("WNN_2010.mat")
%load("SD_2010.mat")

load("WNN_cp_10.mat")
load("MSVM_cp_10.mat")
load("FKNN_cp_10.mat")
load("SD_cp_10.mat")
    load("ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat")
    [mask r c] = bmask("NT_187")
    normNT187 = ReconResults(reconstructedData,r,c)
    %NT187FKNN = FKNN_10p.predictFcn(normNT187)
    %NT187MSVM = MSVM_10p.predictFcn(normNT187)
    %NT187WNN = WNN_10p.predictFcn(normNT187)
    NT187FKNN = FKNN_cp_10.predictFcn(normNT187)
    NT187MSVM = MSVM_cp_10.predictFcn(normNT187)
    NT187WNN = WNN_cp_10.predictFcn(normNT187)
    NT187SD = SD_cp_10.predictFcn(normNT187)
    L = length(r)
    tissue_assigner(NT187FKNN,mask,r,c,L,1,1, "FKNN")
    tissue_assigner(NT187MSVM,mask,r,c,L,1,2, "MSVM")
    tissue_assigner(NT187WNN,mask,r,c,L,1,3, "WNN")
    tissue_assigner(NT187SD,mask,r,c,L,1,4, "SD")
    hold on; figure(1); sgtitle("Classification for NT187")

    load("ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat")
    [mask r c] = bmask("T1_187")
    normT1187 = ReconResults(reconstructedData,r,c)
    %T1187FKNN = FKNN_10p.predictFcn(normT1187)
    %T1187MSVM = MSVM_10p.predictFcn(normT1187)
    %T1187WNN = WNN_10p.predictFcn(normT1187)
    T1187FKNN = FKNN_cp_10.predictFcn(normT1187)
    T1187MSVM = MSVM_cp_10.predictFcn(normT1187)
    T1187WNN = WNN_cp_10.predictFcn(normT1187)
    T1187SD = SD_cp_10.predictFcn(normT1187)
    L = length(r)
    tissue_assigner(T1187FKNN,mask,r,c,L,2,1, "FKNN")
    tissue_assigner(T1187MSVM,mask,r,c,L,2,2, "MSVM")
    tissue_assigner(T1187WNN,mask,r,c,L,2,3, "WNN")
    tissue_assigner(T1187SD,mask,r,c,L,2,4, "SD")
    hold on; figure(2); sgtitle("Classification for T1187")

    load("ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat")
    [mask r c] = bmask("T2_187")
    normT2187 = ReconResults(reconstructedData,r,c)
    %T2187FKNN = FKNN_10p.predictFcn(normT2187)
    %T2187MSVM = MSVM_10p.predictFcn(normT2187)
    %T2187WNN = WNN_10p.predictFcn(normT2187)
    T2187FKNN = FKNN_cp_10.predictFcn(normT2187)
    T2187MSVM = MSVM_cp_10.predictFcn(normT2187)
    T2187WNN = WNN_cp_10.predictFcn(normT2187)
    T2187SD = SD_cp_10.predictFcn(normT2187)
    L = length(r)
    tissue_assigner(T2187FKNN,mask,r,c,L,3,1, "FKNN")
    tissue_assigner(T2187MSVM,mask,r,c,L,3,2, "MSVM")
    tissue_assigner(T2187WNN,mask,r,c,L,3,3, "WNN")
    tissue_assigner(T2187SD,mask,r,c,L,3,4, "SD")
    hold on; figure(3); sgtitle("Classification for T2187")

    load("ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat")
    [mask r c] = bmask("T2_177")
    normT2177 = ReconResults(reconstructedData,r,c)
    %T2177FKNN = FKNN_10p.predictFcn(normT2177)
    %T2177MSVM = MSVM_10p.predictFcn(normT2177)
    %T2177WNN = WNN_10p.predictFcn(normT2177)
    T2177FKNN = FKNN_cp_10.predictFcn(normT2177)
    T2177MSVM = MSVM_cp_10.predictFcn(normT2177)
    T2177WNN = WNN_cp_10.predictFcn(normT2177)
    T2177SD = SD_cp_10.predictFcn(normT2177)
    L = length(r)
    tissue_assigner(T2177FKNN,mask,r,c,L,4,1, "FKNN")
    tissue_assigner(T2177MSVM,mask,r,c,L,4,2, "MSVM")
    tissue_assigner(T2177WNN,mask,r,c,L,4,3, "WNN")
    tissue_assigner(T2177SD,mask,r,c,L,4,4, "SD")
    hold on; figure(4); sgtitle("Classification for T2177")

    load("ReconResults_Brain_177_T1_30s_300iter_M3_Try1.mat")
    [mask r c] = bmask("T1_177")
    normT1177 = ReconResults(reconstructedData,r,c)
    %T1177FKNN = FKNN_10p.predictFcn(normT1177)
    %T1177MSVM = MSVM_10p.predictFcn(normT1177)
    %T1177WNN = WNN_10p.predictFcn(normT1177)
    T1177FKNN = FKNN_cp_10.predictFcn(normT1177)
    T1177MSVM = MSVM_cp_10.predictFcn(normT1177)
    T1177WNN = WNN_cp_10.predictFcn(normT1177)
    T1177SD = SD_cp_10.predictFcn(normT1177)
    L = length(r)
    tissue_assigner(T1177FKNN,mask,r,c,L,5,1, "FKNN")
    tissue_assigner(T1177MSVM,mask,r,c,L,5,2, "MSVM")
    tissue_assigner(T1177WNN,mask,r,c,L,5,3, "WNN")
    tissue_assigner(T1177SD,mask,r,c,L,5,4, "SD")
    hold on; figure(5); sgtitle("Classification for T1177")
    
    load("ReconResults_Brain_177_NT_30s_300iter_M3_Try1.mat")
    [mask r c] = bmask("NT_177")
    normNT177 = ReconResults(reconstructedData,r,c)
    %NT177FKNN = FKNN_10p.predictFcn(normNT177)
    %NT177MSVM = MSVM_10p.predictFcn(normNT177)
    %NT177WNN = WNN_10p.predictFcn(normNT177)
    NT177FKNN = FKNN_cp_10.predictFcn(normNT177)
    NT177MSVM = MSVM_cp_10.predictFcn(normNT177)
    NT177WNN = WNN_cp_10.predictFcn(normNT177)
    NT177SD = SD_cp_10.predictFcn(normNT177)
    L = length(r)
    tissue_assigner(NT177FKNN,mask,r,c,L,6,1, "FKNN")
    tissue_assigner(NT177MSVM,mask,r,c,L,6,2, "MSVM")
    tissue_assigner(NT177WNN,mask,r,c,L,6,3, "WNN")
    tissue_assigner(NT177SD,mask,r,c,L,6,4, "SD")
    hold on; figure(6); sgtitle("Classification for NT177")

    load("ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat")
    [mask r c] = bmask("T2_176")
    normT2176 = ReconResults(reconstructedData,r,c)
    %T2176FKNN = FKNN_10p.predictFcn(normT2176)
    %T2176MSVM = MSVM_10p.predictFcn(normT2176)
    %T2176WNN = WNN_10p.predictFcn(normT2176)
    T2176FKNN = FKNN_cp_10.predictFcn(normT2176)
    T2176MSVM = MSVM_cp_10.predictFcn(normT2176)
    T2176WNN = WNN_cp_10.predictFcn(normT2176)
    T2176SD = SD_cp_10.predictFcn(normT2176)
    L = length(r)
    tissue_assigner(T2176FKNN,mask,r,c,L,7,1, "FKNN")
    tissue_assigner(T2176MSVM,mask,r,c,L,7,2, "MSVM")
    tissue_assigner(T2176WNN,mask,r,c,L,7,3, "WNN")
    tissue_assigner(T2176SD,mask,r,c,L,7,4, "SD")
    hold on; figure(7); sgtitle("Classification for T2176")
    
    load("ReconResults_Brain_176_T1_20s_300iter_M3_Try1.mat")
    [mask r c] = bmask("T1_176")
    normT1176 = ReconResults(reconstructedData,r,c)
    %T1176FKNN = FKNN_10p.predictFcn(normT1176)
    %T1176MSVM = MSVM_10p.predictFcn(normT1176)
    %T1176WNN = WNN_10p.predictFcn(normT1176)
    T1176FKNN = FKNN_cp_10.predictFcn(normT1176)
    T1176MSVM = MSVM_cp_10.predictFcn(normT1176)
    T1176WNN = WNN_cp_10.predictFcn(normT1176)
    T1176SD = SD_cp_10.predictFcn(normT1176)
    L = length(r)
    tissue_assigner(T1176FKNN,mask,r,c,L,8,1, "FKNN")
    tissue_assigner(T1176MSVM,mask,r,c,L,8,2, "MSVM")
    tissue_assigner(T1176WNN,mask,r,c,L,8,3, "WNN")
    tissue_assigner(T1176SD,mask,r,c,L,8,4, "SD")
    hold on; figure(8); sgtitle("Classification for T1176")
    
    load("ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat")
    [mask r c] = bmask("NT_176")
    normNT176 = ReconResults(reconstructedData,r,c)
    %NT176FKNN = FKNN_10p.predictFcn(normNT176)
    %NT176MSVM = MSVM_10p.predictFcn(normNT176)
    %NT176WNN = WNN_10p.predictFcn(normNT176)
    NT176FKNN = FKNN_cp_10.predictFcn(normNT176)
    NT176MSVM = MSVM_cp_10.predictFcn(normNT176)
    NT176WNN = WNN_cp_10.predictFcn(normNT176)
    NT176SD = SD_cp_10.predictFcn(normNT176)
    L = length(r)
    tissue_assigner(NT176FKNN,mask,r,c,L,9,1, "FKNN")
    tissue_assigner(NT176MSVM,mask,r,c,L,9,2, "MSVM")
    tissue_assigner(NT176WNN,mask,r,c,L,9,3, "WNN")
    tissue_assigner(NT176SD,mask,r,c,L,9,4, "SD")
    hold on; figure(9); sgtitle("Classification for NT176")
    
    load("ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat")
    [mask r c] = bmask("T2_158")
    normT2158 = ReconResults(reconstructedData,r,c)
    %T2158FKNN = FKNN_10p.predictFcn(normT2158)
    %T2158MSVM = MSVM_10p.predictFcn(normT2158)
    %T2158WNN = WNN_10p.predictFcn(normT2158)
    T2158FKNN = FKNN_cp_10.predictFcn(normT2158)
    T2158MSVM = MSVM_cp_10.predictFcn(normT2158)
    T2158WNN = WNN_cp_10.predictFcn(normT2158)
    T2158SD = SD_cp_10.predictFcn(normT2158)
    L = length(r)
    tissue_assigner(T2158FKNN,mask,r,c,L,10,1, "FKNN")
    tissue_assigner(T2158MSVM,mask,r,c,L,10,2, "MSVM")
    tissue_assigner(T2158WNN,mask,r,c,L,10,3, "WNN")
    tissue_assigner(T2158SD,mask,r,c,L,10,4, "SD")
    hold on; figure(10); sgtitle("Classification for T2158")
    
    load("ReconResults_Brain_158_T1_20s_300iter_M3_Try1.mat")
    [mask r c] = bmask("T1_158")
    normT1158 = ReconResults(reconstructedData,r,c)
    %T1158FKNN = FKNN_10p.predictFcn(normT1158)
    %T1158MSVM = MSVM_10p.predictFcn(normT1158)
    %T1158WNN = WNN_10p.predictFcn(normT1158)
    T1158FKNN = FKNN_cp_10.predictFcn(normT1158)
    T1158MSVM = MSVM_cp_10.predictFcn(normT1158)
    T1158WNN = WNN_cp_10.predictFcn(normT1158)
    T1158SD = SD_cp_10.predictFcn(normT1158)
    L = length(r)
    tissue_assigner(T1158FKNN,mask,r,c,L,11,1, "FKNN")
    tissue_assigner(T1158MSVM,mask,r,c,L,11,2, "MSVM")
    tissue_assigner(T1158WNN,mask,r,c,L,11,3, "WNN")
    tissue_assigner(T1158SD,mask,r,c,L,11,4, "SD")
    hold on; figure(11); sgtitle("Classification for T1158")
    %%
    load("ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat")
    [mask r c] = bmask("NT_158")
    normNT158 = ReconResults(reconstructedData,r,c)
    %NT158FKNN = FKNN_10p.predictFcn(normNT158)
    %NT158MSVM = MSVM_10p.predictFcn(normNT158)
    %NT158WNN = WNN_10p.predictFcn(normNT158)
    NT158FKNN = FKNN_cp_10.predictFcn(normNT158)
    NT158MSVM = MSVM_cp_10.predictFcn(normNT158)
    NT158WNN = WNN_cp_10.predictFcn(normNT158)
    NT158SD = SD_cp_10.predictFcn(normNT158)
    L = length(r)
    tissue_assigner(NT158FKNN,mask,r,c,L,12,1, "FKNN")
    tissue_assigner(NT158MSVM,mask,r,c,L,12,2, "MSVM")
    tissue_assigner(NT158WNN,mask,r,c,L,12,3, "WNN")
    tissue_assigner(NT158SD,mask,r,c,L,12,4, "SD")
    hold on; figure(12); sgtitle("Classification for NT158")