%% NT_158 -  Gray  
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')    
    graySpec = [squeeze(reconstructedData(66,13,:)),squeeze(reconstructedData(70,12,:)),squeeze(reconstructedData(73,12,:)),squeeze(reconstructedData(77,12,:)),squeeze(reconstructedData(83,11,:)),squeeze(reconstructedData(58,22,:)),squeeze(reconstructedData(49,36,:)),squeeze(reconstructedData(48,44,:)),squeeze(reconstructedData(49,52,:)),squeeze(reconstructedData(49,56,:))]
    magGray = sqrt(sum(graySpec.^2))
    normGray = graySpec./magGray
    label0_Gray = zeros(10,1)
    label0_Gray = transpose(label0_Gray)
    Gray_resp = [normGray,
                  label0_Gray]
% NT_158 - White
    whiteSpec = [squeeze(reconstructedData(73,35,:)),squeeze(reconstructedData(86,26,:)),squeeze(reconstructedData(73,12,:)),squeeze(reconstructedData(77,12,:)),squeeze(reconstructedData(86,37,:)),squeeze(reconstructedData(80,36,:)),squeeze(reconstructedData(79,28,:)),squeeze(reconstructedData(62,41,:)),squeeze(reconstructedData(64,47,:)),squeeze(reconstructedData(79,45,:))]
    magWhite = sqrt(sum(whiteSpec.^2))
    normWhite = whiteSpec./magWhite
    label1_White = ones(10,1)
    label1_White = transpose(label1_White)
    White_resp = [normWhite,
                  label1_White]
    Totalresp = [Gray_resp,White_resp]
%% NT_176
    load('ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat')
    [x1 y1] = bmask("NT_176")
    l = length(x1)
    healallspec1 = []
    for j = 1:l
        healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(x1(j),y1(j),:)))
    end
    magnitude1 = sqrt(sum(healallspec1.^2))
    normhealth1 = healallspec1./magnitude1
    
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
    