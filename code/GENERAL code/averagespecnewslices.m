[mask r c] = bmask("NT_176")
load("ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat")
normNT176 = ReconResults(reconstructedData,r,c)
stdNT176 = std(normNT176)
normNT176 = mean(normNT176,2)

[mask r c] = bmask("NT_187")
load("ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat")
normNT187 = ReconResults(reconstructedData,r,c)
stdNT187 = std(normNT187)
normNT187 = mean(normNT187,2)

[mask r c] = bmask("NT_177")
load("ReconResults_Brain_177_NT_30s_300iter_M3_Try1.mat")
normNT177 = ReconResults(reconstructedData,r,c)
stdNT177 = std(normNT177)
normNT177 = mean(normNT177,2)

[mask r c] = bmask("NT_158")
load("ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat")
normNT158 = ReconResults(reconstructedData,r,c)
stdNT158 = std(normNT158)
normNT158 = mean(normNT158,2)

[mask r c] = bmask("T1_158")
load("ReconResults_Brain_158_T1_20s_300iter_M3_Try1.mat")
normT1158 = ReconResults(reconstructedData,r,c)
stdT1158 = std(normT1158)
normT1158 = mean(normT1158,2)

[mask r c] = bmask("T1_187")
load("ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat")
normT1187 = ReconResults(reconstructedData,r,c)
stdT1187 = std(normT1187)
normT1187 = mean(normT1187,2)

[mask r c] = bmask("T1_177")
load("ReconResults_Brain_177_T1_30s_300iter_M3_Try1.mat")
normT1177 = ReconResults(reconstructedData,r,c)
stdT1177 = std(normT1177)
normT1177 = mean(normT1177,2)

[mask r c] = bmask("T1_176")
load("ReconResults_Brain_176_T1_20s_300iter_M3_Try1.mat")
normT1176 = ReconResults(reconstructedData,r,c)
stdT1176 = std(normT1176)
normT1176 = mean(normT1176,2)

[mask r c] = bmask("T2_158")
load("ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat")
normT2158 = ReconResults(reconstructedData,r,c)
stdT1158 = std(normT2158)
normT2158 = mean(normT2158,2)

[mask r c] = bmask("T2_187")
load("ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat")
normT2187 = ReconResults(reconstructedData,r,c)
stdT2187 = std(normT2187)
normT2187 = mean(normT2187,2)

[mask r c] = bmask("T2_176")
load("ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat")
normT2176 = ReconResults(reconstructedData,r,c)
stdT2176 = std(normT2176)
normT2176 = mean(normT2176,2)

[mask r c] = bmask("T2_177")
load("ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat")
normT2177 = ReconResults(reconstructedData,r,c)
stdT2177 = std(normT2177,2)
normT2177 = mean(normT2177,2)

% Plot graphs with error bars
figure(1); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
e1 = errorbar(qvals, normNT158,stdNT158); 
e1.Color = 'b'
hold on;
figure(1); e2 = errorbar(qvals, normT1158,stdT1158); 
e2.Color = 'g'
hold on;
figure(1); e3 = errorbar(qvals, normT2158,stdT2158); 
e3.Color = 'r'
legend([e1(1),e2(1),e3(1)],'No-Cancer', 'T1 slice','T2 Slice') 

figure(2); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
e1 = errorbar(qvals, normNT187,stdNT187); 
e1.Color = 'b'
hold on;
figure(2); e2 = errorbar(qvals, normT1187,stdT1187); 
e2.Color = 'g'
hold on;
figure(2); e3 = errorbar(qvals, normT2187,stdT2187); 
e3.Color = 'r'
legend([e1(1),e2(1),e3(1)],'No-Cancer', 'T1 slice','T2 Slice')

figure(3); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
e1 = errorbar(qvals, normNT176,stdNT176); 
e1.Color = 'b'
hold on;
figure(3); e2 = errorbar(qvals, normT1176,stdT1176); 
e2.Color = 'g'
hold on;
figure(3); e3 = errorbar(qvals, normT2176,stdT2176); 
e3.Color = 'r'
legend([e1(1),e2(1),e3(1)],'No-Cancer', 'T1 slice','T2 Slice')

figure(4); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
e1 = errorbar(qvals, normNT177,stdNT177); 
e1.Color = 'b'
hold on;
figure(4); e2 = errorbar(qvals, normT1177,stdT1177); 
e2.Color = 'g'
hold on;
figure(4); e3 = errorbar(qvals, normT2177,stdT2177); 
e3.Color = 'r'
legend([e1(1),e2(1),e3(1)],'No-Cancer', 'T1 slice','T2 Slice') 
