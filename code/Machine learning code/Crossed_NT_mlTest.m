% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

%% NT_158 -  Gray  
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')    
    graySpec = [squeeze(reconstructedData(82,13,:)),squeeze(reconstructedData(79,13,:)),squeeze(reconstructedData(76,12,:)),squeeze(reconstructedData(70,11,:)),squeeze(reconstructedData(66,13,:)),squeeze(reconstructedData(65,17,:)),squeeze(reconstructedData(69,20,:)),squeeze(reconstructedData(72,22,:)),squeeze(reconstructedData(63,20,:)),squeeze(reconstructedData(58,22,:)),squeeze(reconstructedData(53,30,:)),squeeze(reconstructedData(48,35,:)),squeeze(reconstructedData(51,39,:)),squeeze(reconstructedData(54,40,:)),squeeze(reconstructedData(53,43,:)),squeeze(reconstructedData(48,43,:)),squeeze(reconstructedData(48,46,:)),squeeze(reconstructedData(50,48,:)),squeeze(reconstructedData(58,57,:)),squeeze(reconstructedData(71,54,:))]
    magGray = sqrt(sum(graySpec.^2))
    normGray = graySpec./magGray
    label0_Gray = zeros(20,1)
    label0_Gray = transpose(label0_Gray)
    Gray_resp = [normGray,
                  label0_Gray]
% NT_158 - White
    whiteSpec = [squeeze(reconstructedData(60,47,:)),squeeze(reconstructedData(62,45,:)),squeeze(reconstructedData(64,43,:)),squeeze(reconstructedData(67,41,:)),squeeze(reconstructedData(70,40,:)),squeeze(reconstructedData(73,40,:)),squeeze(reconstructedData(73,36,:)),squeeze(reconstructedData(68,37,:)),squeeze(reconstructedData(76,42,:)),squeeze(reconstructedData(79,41,:)), squeeze(reconstructedData(78,37,:)),squeeze(reconstructedData(74,33,:)),squeeze(reconstructedData(82,35,:)),squeeze(reconstructedData(84,39,:)),squeeze(reconstructedData(86,38,:)),squeeze(reconstructedData(87,35,:)),squeeze(reconstructedData(86,26,:)),squeeze(reconstructedData(87,29,:)),squeeze(reconstructedData(81,23,:)),squeeze(reconstructedData(79,21,:))]
    magWhite = sqrt(sum(whiteSpec.^2))
    normWhite = whiteSpec./magWhite
    label1_White = ones(20,1)
    label1_White = transpose(label1_White)
    White_resp = [normWhite,
                  label1_White]
    Totalresp = [Gray_resp,White_resp]
%% load data
    load("crossed_nt_FKNN.mat")
    load("crossed_nt_FSVM.mat")
    load("crossed_nt_FT.mat")
    load("crossed_nt_MNN.mat")
    load("crossed_nt_SKNN.mat")
    load("crossed_nt_WNN.mat")
%% NT_176
    load('ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat')
    [mask x1 y1] = bmask("NT_176")
    L = length(x1)
    healallspec1 = []
    for j = 1:L
        healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(x1(j),y1(j),:)))
    end
    magnitude1 = sqrt(sum(healallspec1.^2))
    normhealth1 = healallspec1./magnitude1
    NT176WNN = crossed_nt_WNN.predictFcn(normhealth1)
    NT176SKNN = crossed_nt_SKNN.predictFcn(normhealth1)
    NT176MNN = crossed_nt_MNN.predictFcn(normhealth1)
    NT176FT = crossed_nt_FT.predictFcn(normhealth1)
    NT176FSVM = crossed_nt_FSVM.predictFcn(normhealth1)
    NT176FKNN = crossed_nt_FKNN.predictFcn(normhealth1)
    tissue_assigner(NT176WNN,mask,x1,y1, L,1,1, "WNN")
    tissue_assigner(NT176FT,mask,x1,y1, L,1,2, "FT")
    tissue_assigner(NT176FSVM,mask,x1,y1,L,1,3, "FSVM")
    tissue_assigner(NT176SKNN,mask,x1,y1, L,1,4, "SKNN")
    tissue_assigner(NT176FKNN,mask,x1,y1, L,1,5, "FKNN")
    tissue_assigner(NT176MNN,mask,x1,y1,L,1,6, "MNN")
    figure(1); sgtitle("Classificaion of NT176")
    
% NT_187
    load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
    [mask1 x2 y2] = bmask("NT_187")
    L = length(x2)
    healallspec2 = []
    for k = 1:L
        healallspec2 = cat(2, healallspec2, squeeze(reconstructedData(x2(k),y2(k),:)))
    end
    magnitude2 = sqrt(sum(healallspec2.^2))
    normhealth2 = healallspec2./magnitude2
    NT187WNN = crossed_nt_WNN.predictFcn(normhealth2)
    NT187SKNN = crossed_nt_SKNN.predictFcn(normhealth2)
    NT187MNN = crossed_nt_MNN.predictFcn(normhealth2)
    NT187FT = crossed_nt_FT.predictFcn(normhealth2)
    NT187FSVM = crossed_nt_FSVM.predictFcn(normhealth2)
    NT187FKNN = crossed_nt_FKNN.predictFcn(normhealth2)
    tissue_assigner(NT187WNN,mask1,x2,y2, L,2,1, "WNN")
    tissue_assigner(NT187FT,mask1,x2,y2, L,2,2, "FT")
    tissue_assigner(NT187FSVM,mask1,x2,y2,L,2,3, "FSVM")
    tissue_assigner(NT187SKNN,mask1,x2,y2, L,2,4, "SKNN")
    tissue_assigner(NT187FKNN,mask1,x2,y2, L,2,5, "FKNN")
    tissue_assigner(NT187MNN,mask1,x2,y2,L,2,6, "MNN")
    figure(2); sgtitle("Classificaion of NT187")
    
%NT_177
    load("ReconResults_Brain_177_NT_30s_300iter_M3_Try1.mat")
    [mask2 x3 y3] = bmask("NT_177")
    L = length(x3)
    healallspec3 = []
    for k = 1:L
        healallspec3 = cat(2, healallspec3, squeeze(reconstructedData(x3(k),y3(k),:)))
    end
    magnitude3 = sqrt(sum(healallspec3.^2))
    normhealth3 = healallspec3./magnitude3
    NT177WNN = crossed_nt_WNN.predictFcn(normhealth3)
    NT177SKNN = crossed_nt_SKNN.predictFcn(normhealth3)
    NT177MNN = crossed_nt_MNN.predictFcn(normhealth3)
    NT177FT = crossed_nt_FT.predictFcn(normhealth3)
    NT177FSVM = crossed_nt_FSVM.predictFcn(normhealth3)
    NT177FKNN = crossed_nt_FKNN.predictFcn(normhealth3)
    tissue_assigner(NT177WNN,mask2,x3,y3, L,3,1, "WNN")
    tissue_assigner(NT177FT,mask2,x3,y3, L,3,2, "FT")
    tissue_assigner(NT177FSVM,mask2,x3,y3,L,3,3, "FSVM")
    tissue_assigner(NT177SKNN,mask2,x3,y3, L,3,4, "SKNN")
    tissue_assigner(NT177FKNN,mask2,x3,y3, L,3,5, "FKNN")
    tissue_assigner(NT177MNN,mask2,x3,y3,L,3,6, "MNN")
    figure(3); sgtitle("Classificaion of NT177")
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
    NT158WNN = crossed_nt_WNN.predictFcn(normhealth)
    NT158SKNN = crossed_nt_SKNN.predictFcn(normhealth)
    NT158MNN = crossed_nt_MNN.predictFcn(normhealth)
    NT158FT = crossed_nt_FT.predictFcn(normhealth)
    NT158FSVM = crossed_nt_FSVM.predictFcn(normhealth)
    NT158FKNN = crossed_nt_FKNN.predictFcn(normhealth)
    tissue_assigner(NT158WNN,mask3,x,y, L,4,1, "WNN")
    tissue_assigner(NT158FT,mask3,x,y, L,4,2, "FT")
    tissue_assigner(NT158FSVM,mask3,x,y,L,4,3, "FSVM")
    tissue_assigner(NT158SKNN,mask3,x,y, L,4,4, "SKNN")
    tissue_assigner(NT158FKNN,mask3,x,y, L,4,5, "FKNN")
    tissue_assigner(NT158MNN,mask3,x,y,L,4,6, "MNN")
    figure(4); sgtitle("Classificaion of NT158")
