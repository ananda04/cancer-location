% Description: This script draws the spectra line graphs by plotting mean and error from GBMT reconstructedData.
%              Creating graphs for Cancer, White and Gray Matter spectra points
%              Cancer, White and Gray Matter regions were determined through Dr.Greenberg & coordinating Oncologist.
% Author:      Arnav Nanda 
%              Apex High School 
%              Fitzpatrick Center for Interdisciplinary Engineering, Medicine and Applied Science (CIEMAS), Duke University  
              
load('ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat') 

figure(5); imagesc(squeeze(reconstructedData(:, :, 17)))

% Drawing cancerous spectra lines 
figure(1); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Cancer Spectra")
hold on; plot(qvals, squeeze(reconstructedData(79,23,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(82,23,:))) 
hold on;figure(1); plot(qvals, squeeze(reconstructedData(81,21,:)))
hold on; figure(1); plot(qvals, squeeze(reconstructedData(79,19,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(79,22,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(78,25,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(78,27,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(76,26,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(77,26,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(80,33,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(75,32,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(74,30,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(73,23,:)))

% Drawing White Matter spectra lines 
figure(2); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("White Matter Spectra"); 
hold on;figure(2); plot(qvals, squeeze(reconstructedData(88,11,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(86,10,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(87,13,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(85,12,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(57,25,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(59,31,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(63,37,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(86,61,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(85,59,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(81,56,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(77,15,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(76,13,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(75,14,:)))

% Drawing Gray Matter spectra lines 
figure(3); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Gray Matter Spectra")
hold on;figure(3); plot(qvals, squeeze(reconstructedData(79,10,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(76,9,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(78,10,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(77,9,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(60,13,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(60,14,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(59,15,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(58,19,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(56,21,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(63,44,:))); 
hold on;figure(3); plot(qvals, squeeze(reconstructedData(56,62,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(74,65,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(77,64,:)));

% calculating averge spectraline for 
blueline = squeeze(reconstructedData(79,23,:))+squeeze(reconstructedData(73,23,:))+squeeze(reconstructedData(82,23,:))+squeeze(reconstructedData(81,21,:))+squeeze(reconstructedData(79,19,:))+squeeze(reconstructedData(79,22,:))+squeeze(reconstructedData(78,25,:))+squeeze(reconstructedData(78,27,:))+squeeze(reconstructedData(76,26,:))+squeeze(reconstructedData(77,26,:))+squeeze(reconstructedData(80,33,:))+squeeze(reconstructedData(75,32,:))+squeeze(reconstructedData(74,30,:))
greenline = squeeze(reconstructedData(88,11,:))+squeeze(reconstructedData(86,10,:))+squeeze(reconstructedData(87,13,:))+squeeze(reconstructedData(85,12,:))+squeeze(reconstructedData(57,25,:))+squeeze(reconstructedData(59,31,:))+squeeze(reconstructedData(63,37,:))+squeeze(reconstructedData(86,61,:))+squeeze(reconstructedData(85,59,:))+squeeze(reconstructedData(81,56,:))+squeeze(reconstructedData(77,15,:))+squeeze(reconstructedData(76,13,:))+squeeze(reconstructedData(75,14,:))
redline = squeeze(reconstructedData(79,10,:))+squeeze(reconstructedData(76,9,:))+squeeze(reconstructedData(78,10,:))+squeeze(reconstructedData(77,9,:))+squeeze(reconstructedData(60,13,:))+squeeze(reconstructedData(60,14,:))+squeeze(reconstructedData(59,15,:))+squeeze(reconstructedData(58,19,:))+squeeze(reconstructedData(56,21,:))+squeeze(reconstructedData(63,44,:))+squeeze(reconstructedData(56,62,:))+squeeze(reconstructedData(74,65,:))+squeeze(reconstructedData(77,64,:))

avgblueline = blueline/13
avggreenline = greenline/13
avgredline = redline/13

errblue = erf(avgblueline)
errgreen = erf(avggreenline)
errred = erf(avgredline)

% Plot graphs with error bars
figure(5); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
e1 = errorbar(qvals, avgblueline,errblue); 
e1.Color = 'b'
hold on;
figure(5); e2 = errorbar(qvals, avggreenline,errgreen); 
e2.Color = 'g'
hold on;
figure(5); e3 = errorbar(qvals, avgredline, errred); 
e3.Color = 'r'
legend([e1(1),e2(1),e3(1)],'Cancer', 'White Matter','Gray Matter' )