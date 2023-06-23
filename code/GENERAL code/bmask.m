%% more efficient and accurate binary mask function 
% created by Arnav Nanda
function [mask r c] = bmask(String)
if String == "GBMT_Normal" 
    load("ReconResults_Model_3_BrainScan_Normal_FullFanExtent_400iter_M3.mat")
    load("normal_mask.mat")
    redChannel = normal(:, :, 1);
    greenChannel = normal(:, :, 2);
    blueChannel = normal(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "GBMT"
    load("ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat")
    load("gbmt_mask.mat")
    redChannel = gbmt(:, :, 1);
    greenChannel = gbmt(:, :, 2);
    blueChannel = gbmt(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "T2_187"
    load("ReconResults_Brain_187_T2_20s_300iter_M3_Try1")
    load("t2187_mask.mat")
    redChannel = t2187(:, :, 1);
    greenChannel = t2187(:, :, 2);
    blueChannel = t2187(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "T1_187"
    load("ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat")
    load("t1187_mask.mat")
    redChannel = t1187(:, :, 1);
    greenChannel = t1187(:, :, 2);
    blueChannel = t1187(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "NT_187"
    load("ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat")
    load("nt187_mask.mat")
    redChannel = nt187(:, :, 1);
    greenChannel = nt187(:, :, 2);
    blueChannel = nt187(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "T2_177"
    load("ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat")
    load("t2177_mask.mat")
    redChannel = t2177(:, :, 1);
    greenChannel = t2177(:, :, 2);
    blueChannel = t2177(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "T1_177"
    load("ReconResults_Brain_177_T1_30s_300iter_M3_Try1.mat")
    load("t1177_mask.mat")
    redChannel = t1177(:, :, 1);
    greenChannel = t1177(:, :, 2);
    blueChannel = t1177(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "NT_177"
    load("ReconResults_Brain_177_NT_30s_300iter_M3_Try1.mat")
    load("nt177_mask.mat")
    redChannel = nt177(:, :, 1);
    greenChannel = nt177(:, :, 2);
    blueChannel = nt177(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "T2_176"
    load("ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat")
    load("t2176_mask.mat")
    redChannel = t2176(:, :, 1);
    greenChannel = t2176(:, :, 2);
    blueChannel = t2176(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "T1_176"
    load("ReconResults_Brain_176_T1_20s_300iter_M3_Try1.mat")
    load("t1176_mask.mat")
    redChannel = t1176(:, :, 1);
    greenChannel = t1176(:, :, 2);
    blueChannel = t1176(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "NT_176"
    load("ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat")
    load("nt176_mask.mat")
    redChannel = nt176(:, :, 1);
    greenChannel = nt176(:, :, 2);
    blueChannel = nt176(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "T2_158"
    load("ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat")
    load("t2158_mask.mat")
    redChannel = t2158(:, :, 1);
    greenChannel = t2158(:, :, 2);
    blueChannel = t2158(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "T1_158"
    load("ReconResults_Brain_158_T1_20s_300iter_M3_Try1.mat")
    load("t1158_mask.mat")
    redChannel = t1158(:, :, 1);
    greenChannel = t1158(:, :, 2);
    blueChannel = t1158(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 
if String == "NT_158"
    load("ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat")
    load("nt158_mask.mat")
    redChannel = nt158(:, :, 1);
    greenChannel = nt158(:, :, 2);
    blueChannel = nt158(:, :, 3);
    mask = blueChannel > 1, redChannel > 1, greenChannel > 1;
    [r c] = find(mask == 1)
end 