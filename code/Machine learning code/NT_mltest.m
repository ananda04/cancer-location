% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

%% NT_158 -  Gray  
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')    
    graySpec = [squeeze(reconstructedData(66,13,:)),squeeze(reconstructedData(70,12,:)),squeeze(reconstructedData(73,12,:)),squeeze(reconstructedData(77,12,:)),squeeze(reconstructedData(83,11,:)),squeeze(reconstructedData(58,22,:)),squeeze(reconstructedData(49,36,:)),squeeze(reconstructedData(48,44,:)),squeeze(reconstructedData(49,52,:)),squeeze(reconstructedData(49,56,:)), squeeze(reconstructedData(83,11,:)),squeeze(reconstructedData(86,12,:)),squeeze(reconstructedData(79,11,:)),squeeze(reconstructedData(64,15,:)),squeeze(reconstructedData(54,28,:)),squeeze(reconstructedData(50,34,:)),squeeze(reconstructedData(48,45,:)),squeeze(reconstructedData(54,58,:)),squeeze(reconstructedData(59,58,:)),squeeze(reconstructedData(64,57,:))]
    magGray = sqrt(sum(graySpec.^2))
    normGray = graySpec./magGray
    label0_Gray = zeros(20,1)
    label0_Gray = transpose(label0_Gray)
    Gray_resp = [normGray,
                  label0_Gray]
% NT_158 - White
    whiteSpec = [squeeze(reconstructedData(73,35,:)),squeeze(reconstructedData(86,26,:)),squeeze(reconstructedData(81,23,:)),squeeze(reconstructedData(77,21,:)),squeeze(reconstructedData(86,37,:)),squeeze(reconstructedData(80,36,:)),squeeze(reconstructedData(79,28,:)),squeeze(reconstructedData(62,41,:)),squeeze(reconstructedData(64,47,:)),squeeze(reconstructedData(79,45,:)),squeeze(reconstructedData(83,31,:)),squeeze(reconstructedData(81,40,:)),squeeze(reconstructedData(84,43,:)),squeeze(reconstructedData(70,33,:)),squeeze(reconstructedData(72,41,:)),squeeze(reconstructedData(65,37,:)),squeeze(reconstructedData(77,41,:)),squeeze(reconstructedData(74,30,:)),squeeze(reconstructedData(84,32,:)),squeeze(reconstructedData(68,39,:))]
    magWhite = sqrt(sum(whiteSpec.^2))
    normWhite = whiteSpec./magWhite
    label1_White = ones(20,1)
    label1_White = transpose(label1_White)
    White_resp = [normWhite,
                  label1_White]
    Totalresp = [Gray_resp,White_resp]
%% load data
    load("nt_FKNN.mat")
    load("nt_FSVM.mat")
    load("nt_FT.mat")
    load("nt_MNN.mat")
    load("nt_SKNN.mat")
    load("nt_WNN.mat")
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
    NT176WNN = ntModel_WNN.predictFcn(normhealth1)
    NT176SKNN = ntModel_SKNN.predictFcn(normhealth1)
    NT176MNN = ntModel_MNN.predictFcn(normhealth1)
    NT176FT = ntModel_FT.predictFcn(normhealth1)
    NT176FSVM = ntModel_FSVM.predictFcn(normhealth1)
    NT176FKNN = ntModel_FKNN.predictFcn(normhealth1)
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
    NT187WNN = ntModel_WNN.predictFcn(normhealth2)
    NT187SKNN = ntModel_SKNN.predictFcn(normhealth2)
    NT187MNN = ntModel_MNN.predictFcn(normhealth2)
    NT187FT = ntModel_FT.predictFcn(normhealth2)
    NT187FSVM = ntModel_FSVM.predictFcn(normhealth2)
    NT187FKNN = ntModel_FKNN.predictFcn(normhealth2)
    %%
    tissue_assigner(NT187WNN,mask1,x2,y2, L,2,1, "WNN")
    tissue_assigner(NT187FT,mask1,x2,y2, L,2,2, "FT")
    tissue_assigner(NT187FSVM,mask1,x2,y2,L,2,3, "FSVM")
    tissue_assigner(NT187SKNN,mask1,x2,y2, L,2,4, "SKNN")
    tissue_assigner(NT187FKNN,mask1,x2,y2, L,2,5, "FKNN")
    tissue_assigner(NT187MNN,mask1,x2,y2,L,2,6, "MNN")
    figure(2); sgtitle("Classificaion of NT_187")
    
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
    NT177WNN = ntModel_WNN.predictFcn(normhealth3)
    NT177SKNN = ntModel_SKNN.predictFcn(normhealth3)
    NT177MNN = ntModel_MNN.predictFcn(normhealth3)
    NT177FT = ntModel_FT.predictFcn(normhealth3)
    NT177FSVM = ntModel_FSVM.predictFcn(normhealth3)
    NT177FKNN = ntModel_FKNN.predictFcn(normhealth3)
    %%
    tissue_assigner(NT177WNN,mask2,x3,y3, L,3,1, "WNN")
    tissue_assigner(NT177FT,mask2,x3,y3, L,3,2, "FT")
    tissue_assigner(NT177FSVM,mask2,x3,y3,L,3,3, "FSVM")
    tissue_assigner(NT177SKNN,mask2,x3,y3, L,3,4, "SKNN")
    tissue_assigner(NT177FKNN,mask2,x3,y3, L,3,5, "FKNN")
    tissue_assigner(NT177MNN,mask2,x3,y3,L,3,6, "MNN")
    figure(3); sgtitle("Classificaion of NT_177")
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
    NT158WNN = ntModel_WNN.predictFcn(normhealth)
    NT158SKNN = ntModel_SKNN.predictFcn(normhealth)
    NT158MNN = ntModel_MNN.predictFcn(normhealth)
    NT158FT = ntModel_FT.predictFcn(normhealth)
    NT158FSVM = ntModel_FSVM.predictFcn(normhealth)
    NT158FKNN = ntModel_FKNN.predictFcn(normhealth)
    %%
    tissue_assigner(NT158WNN,mask3,x,y, L,4,1, "WNN")
    tissue_assigner(NT158FT,mask3,x,y, L,4,2, "FT")
    tissue_assigner(NT158FSVM,mask3,x,y,L,4,3, "FSVM")
    tissue_assigner(NT158SKNN,mask3,x,y, L,4,4, "SKNN")
    tissue_assigner(NT158FKNN,mask3,x,y, L,4,5, "FKNN")
    tissue_assigner(NT158MNN,mask3,x,y,L,4,6, "MNN")
    figure(3); sgtitle("Classificaion of NT158")
