% NT_158
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')
    [mask x y] = bmask("NT_158")
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
%% NT_176
    load('ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat')
    [mask x1 y1] = bmask("NT_176")
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
%% NT_187
    load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
    [mask x2 y2] = bmask("NT_187")
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
    [mask x3 y3] = bmask("NT_177")
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
    T1158spec = [squeeze(reconstructedData(61,51,:)),squeeze(reconstructedData(64,41,:)),squeeze(reconstructedData(58,43,:)),squeeze(reconstructedData(59,36,:)),squeeze(reconstructedData(69,35,:)),squeeze(reconstructedData(72,41,:)),squeeze(reconstructedData(65,22,:)),squeeze(reconstructedData(69,18,:)),squeeze(reconstructedData(76,20,:)),squeeze(reconstructedData(61,12,:))]
    magT1158 = sqrt(sum(T1158spec.^2))
    normT1158 = T1158spec./magT1158
    label1_T1158 = ones(10,1)
    label1_T1158 = transpose(label1_T1158)
    T1158_resp = [normT1158,
                  label1_T1158]
    [mask] = bmask("T1_158")
    figure(1);imshow(mask)
    hold on;
    figure(1);plot(51,61,"rX")
    hold on;
    figure(1); plot(41,64,"rX")
    hold on;
    figure(1); plot(43,58,"rX")
    hold on;
    figure(1); plot(36,59,"rX")
    hold on;
    figure(1); plot(35,69,"rX")
    hold on;
    figure(1); plot(41,72,"rX")
    hold on;
    figure(1); plot(18,69,"rX")
    hold on;
    figure(1); plot(20,76,"rX")
    hold on;
    figure(1); plot(12,61,"rX")
    hold on;
    figure(1); plot(22,65,"rX")
    hold on;
    figure(1); title("t1158 data collected from approximated regions")
% T1_176
    load('ReconResults_Brain_176_T1_20s_300iter_M3_Try1.mat') 
    T1176spec = [squeeze(reconstructedData(84,48,:)),squeeze(reconstructedData(83,51,:)),squeeze(reconstructedData(84,54,:)),squeeze(reconstructedData(85,58,:)),squeeze(reconstructedData(88,49,:)),squeeze(reconstructedData(88,51,:)),squeeze(reconstructedData(88,54,:)),squeeze(reconstructedData(88,57,:)),squeeze(reconstructedData(87,59,:)),squeeze(reconstructedData(85,61,:))]
    magT1176 = sqrt(sum(T1176spec.^2))
    normT1176 = T1176spec./magT1176
    label1_T1176 = ones(10,1)
    label1_T1176 = transpose(label1_T1176)
    T1176_resp = [normT1176,
                  label1_T1176]
    [mask] = bmask("T1_176")
    figure(1);imshow(mask)
    hold on;
    figure(1);plot(48,84,"rX")
    hold on;
    figure(1); plot(51,83,"rX")
    hold on;
    figure(1); plot(54,84,"rX")
    hold on;
    figure(1); plot(58,85,"rX")
    hold on;
    figure(1); plot(49,88,"rX")
    hold on;
    figure(1); plot(51,88,"rX")
    hold on;
    figure(1); plot(54,88,"rX")
    hold on;
    figure(1); plot(57,88,"rX")
    hold on;
    figure(1); plot(59,87,"rX")
    hold on;
    figure(1); plot(61,85,"rX")
    hold on;
    figure(1); title("t1176 data collected from approximated regions")
% T1_187
    load('ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat')
    T1187spec = [squeeze(reconstructedData(75,22,:)),squeeze(reconstructedData(73,28,:)),squeeze(reconstructedData(72,32,:)),squeeze(reconstructedData(75,36,:)),squeeze(reconstructedData(77,31,:)),squeeze(reconstructedData(79,26,:)),squeeze(reconstructedData(83,29,:)),squeeze(reconstructedData(82,33,:)),squeeze(reconstructedData(80,38,:)),squeeze(reconstructedData(86,38,:))]
    magT1187 = sqrt(sum(T1187spec.^2))
    normT1187 = T1187spec./magT1187
    label1_T1187 = ones(10,1)
    label1_T1187 = transpose(label1_T1187)
    T1187_resp = [normT1187,
                  label1_T1187]
    [mask] = bmask("T1_187")
    figure(1);imshow(mask)
    hold on;
    figure(1);plot(22,75,"rX")
    hold on;
    figure(1); plot(28,73,"rX")
    hold on;
    figure(1); plot(32,72,"rX")
    hold on;
    figure(1); plot(36,75,"rX")
    hold on;
    figure(1); plot(31,77,"rX")
    hold on;
    figure(1); plot(26,79,"rX")
    hold on;
    figure(1); plot(29,83,"rX")
    hold on;
    figure(1); plot(33,82,"rX")
    hold on;
    figure(1); plot(38,80,"rX")
    hold on;
    figure(1); plot(38,86,"rX")
    hold on;
    figure(1); title("t1187 data collected from approximated regions")
%T1_177
    load("ReconResults_Brain_177_T1_30s_300iter_M3_Try1.mat")
    T1177spec = [squeeze(reconstructedData(91,19,:)),squeeze(reconstructedData(89,20,:)),squeeze(reconstructedData(86,21,:)),squeeze(reconstructedData(84,23,:)),squeeze(reconstructedData(83,25,:)),squeeze(reconstructedData(92,21,:)),squeeze(reconstructedData(90,21,:)),squeeze(reconstructedData(89,21,:)),squeeze(reconstructedData(87,22,:)),squeeze(reconstructedData(91,22,:))]
    magT1177 = sqrt(sum(T1177spec.^2))
    normT1177 = T1177spec./magT1177
    label1_T1177 = ones(10,1)
    label1_T1177 = transpose(label1_T1177)
    T1177_resp = [normT1177,
                  label1_T1177]
    [mask] = bmask("T1_177")
    figure(1);imshow(mask)
    hold on;
    figure(1);plot(19,91,"rX")
    hold on;
    figure(1); plot(20,89,"rX")
    hold on;
    figure(1); plot(21,86,"rX")
    hold on;
    figure(1); plot(23,84,"rX")
    hold on;
    figure(1); plot(25,83,"rX")
    hold on;
    figure(1); plot(21,92,"rX")
    hold on;
    figure(1); plot(21,90,"rX")
    hold on;
    figure(1); plot(21,89,"rX")
    hold on;
    figure(1); plot(22,87,"rX")
    hold on;
    figure(1); plot(22,91,"rX")
    hold on;
    figure(1); title("t1177 data collected from approximated regions")
% T2_158
    load('ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat')
    T2158spec = [squeeze(reconstructedData(62,65,:)),squeeze(reconstructedData(62,63,:)),squeeze(reconstructedData(62,61,:)),squeeze(reconstructedData(61,59,:)),squeeze(reconstructedData(63,57,:)),squeeze(reconstructedData(64,59,:)),squeeze(reconstructedData(65,62,:)),squeeze(reconstructedData(64,65,:)),squeeze(reconstructedData(67,65,:)),squeeze(reconstructedData(75,62,:))]
    magT2158 = sqrt(sum(T2158spec.^2))
    normT2158 = T2158spec./magT2158
    label1_T2158 = ones(10,1)
    label1_T2158 = transpose(label1_T2158)
    T2158_resp = [normT2158,
                  label1_T2158]
    [mask] = bmask("T2_158")
    figure(1);imshow(mask)
    hold on;
    figure(1);plot(65,62,"rX")
    hold on;
    figure(1); plot(63,62,"rX")
    hold on;
    figure(1); plot(61,62,"rX")
    hold on;
    figure(1); plot(59,61,"rX")
    hold on;
    figure(1); plot(57,63,"rX")
    hold on;
    figure(1); plot(59,64,"rX")
    hold on;
    figure(1); plot(62,65,"rX")
    hold on;
    figure(1); plot(65,64,"rX")
    hold on;
    figure(1); plot(65,67,"rX")
    hold on;
    figure(1); plot(62,75,"rX")
    hold on;
    figure(1); title("t2158 data collected from approximated regions")

%% T2_176
    load('ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat')
    T2176spec = [squeeze(reconstructedData(68,10,:)),squeeze(reconstructedData(66,11,:)),squeeze(reconstructedData(65,11,:)),squeeze(reconstructedData(64,10,:)),squeeze(reconstructedData(63,10,:)),squeeze(reconstructedData(66,9,:)),squeeze(reconstructedData(68,8,:)),squeeze(reconstructedData(67,7,:)),squeeze(reconstructedData(66,7,:)),squeeze(reconstructedData(65,8,:))]
    magT2176 = sqrt(sum(T2176spec.^2))
    normT2176 = T2176spec./magT2176
    label1_T2176 = ones(10,1)
    label1_T2176 = transpose(label1_T2176)
    T2176_resp = [normT2176,
                  label1_T2176]
    [mask] = bmask("T2_176")
    figure(1);imshow(mask)
    hold on;
    figure(1);plot(10,68,"rX")
    hold on;
    figure(1); plot(11,66,"rX")
    hold on;
    figure(1); plot(11,65,"rX")
    hold on;
    figure(1); plot(10,64,"rX")
    hold on;
    figure(1); plot(10,63,"rX")
    hold on;
    figure(1); plot(9,66,"rX")
    hold on;
    figure(1); plot(8,68,"rX")
    hold on;
    figure(1); plot(7,67,"rX")
    hold on;
    figure(1); plot(7,66,"rX")
    hold on;
    figure(1); plot(8,65,"rX")
    hold on;
    figure(1); title("t2176 data collected from approximated regions")
%% T2_187
    load('ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat')
    T2187spec = [squeeze(reconstructedData(83,27,:)),squeeze(reconstructedData(79,27,:)),squeeze(reconstructedData(76,27,:)),squeeze(reconstructedData(76,32,:)),squeeze(reconstructedData(76,36,:)),squeeze(reconstructedData(78,37,:)),squeeze(reconstructedData(81,37,:)),squeeze(reconstructedData(84,40,:)),squeeze(reconstructedData(83,33,:)),squeeze(reconstructedData(79,32,:))]
    magT2187 = sqrt(sum(T2187spec.^2))
    normT2187 = T2187spec./magT2187
    label1_T2187 = ones(10,1)
    label1_T2187 = transpose(label1_T2187)
    T2187_resp = [normT2187,
                  label1_T2187]
    [mask] = bmask("T2_187")
    figure(1);imshow(mask)
    hold on;
    figure(1);plot(27,83,"rX")
    hold on;
    figure(1); plot(27,79,"rX")
    hold on;
    figure(1); plot(27,76,"rX")
    hold on;
    figure(1); plot(32,76,"rX")
    hold on;
    figure(1); plot(36,76,"rX")
    hold on;
    figure(1); plot(37,78,"rX")
    hold on;
    figure(1); plot(37,81,"rX")
    hold on;
    figure(1); plot(40,84,"rX")
    hold on;
    figure(1); plot(33,83,"rX")
    hold on;
    figure(1); plot(32,79,"rX")
    hold on;
    figure(1); title("t2187 data collected from approximated regions")
% T2_177
    load("ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat")
    T2177spec = [squeeze(reconstructedData(89,17,:)),squeeze(reconstructedData(88,18,:)),squeeze(reconstructedData(87,21,:)),squeeze(reconstructedData(85,24,:)),squeeze(reconstructedData(86,27,:)),squeeze(reconstructedData(88,28,:)),squeeze(reconstructedData(89,29,:)),squeeze(reconstructedData(89,26,:)),squeeze(reconstructedData(89,21,:)),squeeze(reconstructedData(88,24,:))]
    magT2177 = sqrt(sum(T2177spec.^2))
    normT2177 = T2177spec./magT2177
    label1_T2177 = ones(10,1)
    label1_T2177 = transpose(label1_T2177)
    T2177_resp = [normT2177,
                  label1_T2177]
    [mask] = bmask("T2_177")
    figure(1);imshow(mask)
    hold on;
    figure(1);plot(17,89,"rX")
    hold on;
    figure(1); plot(18,88,"rX")
    hold on;
    figure(1); plot(21,87,"rX")
    hold on;
    figure(1); plot(24,85,"rX")
    hold on;
    figure(1); plot(27,86,"rX")
    hold on;
    figure(1); plot(28,88,"rX")
    hold on;
    figure(1); plot(29,89,"rX")
    hold on;
    figure(1); plot(26,89,"rX")
    hold on;
    figure(1); plot(21,89,"rX")
    hold on;
    figure(1); plot(24,88,"rX")
    hold on;
    figure(1); title("t2177 data collected from approximated regions")
    %% different combos
% T1158
    NT_T1158 = [NT158_resp,T1158_resp]
%T2158
    NT_T2158 = [NT158_resp,T2158_resp]
% T1176
    NT_T1176 = [NT176_resp,T1176_resp]
%% T2176
    NT_T2176 = [NT176_resp,T2176_resp]
%% T1187
    NT_T1187 = [NT187_resp,T1187_resp]
% T2187
    NT_T2187 = [NT187_resp,T2187_resp]
% T1177
    NT_T1177 = [NT177_resp,T1177_resp]
% T2177
    NT_T2177 = [NT177_resp,T2177_resp]
%% Run ML Code
% load Models
    load("cmi_T1158_FSVM.mat")
    load("cmi_T1158_FT.mat")
    load("cmi_T1158_WNN.mat")
    load("cmi_T1176_FSVM.mat")
    load("cmi_T1176_FT.mat")
    load("cmi_T1176_WNN.mat")
    load("cmi_T1177_FSVM.mat")
    load("cmi_T1177_FT.mat")
    load("cmi_T1177_WNN.mat")
    load("cmi_T1187_FSVM.mat")
    load("cmi_T1187_FT.mat")
    load("cmi_T1187_WNN.mat")
    load("cmi_T2158_FSVM.mat")
    load("cmi_T2158_FT.mat")
    load("cmi_T2158_WNN.mat")
    load("cmi_T2177_FSVM.mat")
    load("cmi_T2177_FT.mat")
    load("cmi_T2177_WNN.mat")
    load("cmi_T2187_FSVM.mat")
    load("cmi_T2187_FT.mat")
    load("cmi_T2187_WNN.mat")
    load("cmi_T2176_FSVM.mat")
    load("cmi_T2176_WNN.mat")
    load("cmi_T2176_FT")
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
    T1158NN = cmi_T1158_WNN.predictFcn(normT1158);
    T1158FT = cmi_T1158_FT.predictFcn(normT1158);
    T1158SVM = cmi_T1158_FSVM.predictFcn(normT1158);
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
    T1176NN = cmi_T1176_WNN.predictFcn(normT1176);
    T1176FT = cmi_T1176_FT.predictFcn(normT1176);
    T1176SVM = cmi_T1176_FSVM.predictFcn(normT1176);
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
    T1187NN = cmi_T1187_WNN.predictFcn(normT1187);
    T1187FT = cmi_T1187_FT.predictFcn(normT1187);
    T1187SVM = cmi_T1187_FSVM.predictFcn(normT1187);
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
    T1177NN = cmi_T1177_WNN.predictFcn(normT1177);
    T1177FT = cmi_T1177_FT.predictFcn(normT1177);
    T1177SVM = cmi_T1177_FSVM.predictFcn(normT1177);
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
    T2158NN = cmi_T2158_WNN.predictFcn(normT2158);
    T2158FT = cmi_T2158_FT.predictFcn(normT2158);
    T2158SVM = cmi_T2158_FSVM.predictFcn(normT2158);
    L = length(normT2158)
    tissue_assigner(T2158NN,mask,r4,c4, L,5,1, "WNN")
    tissue_assigner(T2158FT,mask,r4,c4, L,5,2, "FT")
    tissue_assigner(T2158SVM,mask,r4,c4,L,5,3, "FSVM")
    figure(5); sgtitle("Classificaion of t2158")

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
    T2187NN = cmi_T2187_WNN.predictFcn(normT2187);
    T2187FT = cmi_T2187_FT.predictFcn(normT2187);
    T2187SVM = cmi_T2187_FSVM.predictFcn(normT2187);
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
    T2177NN = cmi_T2177_WNN.predictFcn(normT2177);
    T2177FT = cmi_T2177_FT.predictFcn(normT2177);
    T2177SVM = cmi_T2177_FSVM.predictFcn(normT2177);
    L = length(normT2177)
    tissue_assigner(T2177NN,mask,r7,c7, L,8,1, "WNN")
    tissue_assigner(T2177FT,mask,r7,c7, L,8,2, "FT")
    tissue_assigner(T2177SVM,mask,r7,c7,L,8,3, "FSVM")
    figure(8); sgtitle("Classificaion of t2177")
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
    T2176NN = cmi_T2176_WNN.predictFcn(normT2176);
    T2176FT = cmi_T2176_FT.predictFcn(normT2176);
    T2176SVM = cmi_T2176_FSVM.predictFcn(normT2176);
    L = length(normT2176)
    tissue_assigner(T2176NN,mask,r5,c5, L,6,1, "WNN")
    tissue_assigner(T2176FT,mask,r5,c5, L,6,2, "FT")
    tissue_assigner(T2176SVM,mask,r5,c5,L,6,3, "FSVM")
    figure(6); sgtitle("Classificaion of t2176")