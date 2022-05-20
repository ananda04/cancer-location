% Description: This script draws the spectra line graphs by plotting mean and error from GBMT reconstructedData.
%              Creating graphs for Cancer, White and Gray Matter spectra points
%              Cancer, White and Gray Matter regions were determined through Dr.Greenberg & coordinating Oncologist.
% Author:      Arnav Nanda 
%              Apex High School 
%              Fitzpatrick Center for Interdisciplinary Engineering, Medicine and Applied Science (CIEMAS), Duke University  
              
load('ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat') 

figure(4); imagesc(squeeze(reconstructedData(:, :, 17)))

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
hold on;figure(1); plot(qvals, squeeze(reconstructedData(83,34,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(79,33,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(77,33,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(74,31,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(77,30,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(79,28,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(86,27,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(87,35,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(84,37,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(81,36,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(80,34,:)))
hold on;figure(1); plot(qvals, squeeze(reconstructedData(86,31,:)))

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
hold on;figure(2); plot(qvals, squeeze(reconstructedData(74,49,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(75,56,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(79,51,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(80,50,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(64,40,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(64,30,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(64,24,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(86,56,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(64,33,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(64,34,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(65,25,:)))
hold on;figure(2); plot(qvals, squeeze(reconstructedData(65,34,:)))

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
hold on;figure(3); plot(qvals, squeeze(reconstructedData(83,63,:))); 
hold on;figure(3); plot(qvals, squeeze(reconstructedData(55,61,:))); 
hold on;figure(3); plot(qvals, squeeze(reconstructedData(57,61,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(55,64,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(55,65,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(62,63,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(65,64,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(68,63,:))); 
hold on;figure(3); plot(qvals, squeeze(reconstructedData(70,65,:)));
hold on;figure(3); plot(qvals, squeeze(reconstructedData(58,15,:))); 
hold on;figure(3); plot(qvals, squeeze(reconstructedData(59,18,:))); 
hold on; figure(3); plot(qvals, squeeze(reconstructedData(55,23,:))); 

% calculating averge spectraline  
blueline = squeeze(reconstructedData(79,23,:))+squeeze(reconstructedData(73,23,:))+squeeze(reconstructedData(82,23,:))+squeeze(reconstructedData(81,21,:))+squeeze(reconstructedData(79,19,:))+squeeze(reconstructedData(79,22,:))+squeeze(reconstructedData(78,25,:))+squeeze(reconstructedData(78,27,:))+squeeze(reconstructedData(76,26,:))+squeeze(reconstructedData(77,26,:))+squeeze(reconstructedData(80,33,:))+squeeze(reconstructedData(75,32,:))+squeeze(reconstructedData(74,30,:))+squeeze(reconstructedData(83,34,:))+squeeze(reconstructedData(79,33,:))+squeeze(reconstructedData(77,33,:))+squeeze(reconstructedData(74,31,:))+squeeze(reconstructedData(77,30,:))+squeeze(reconstructedData(79,28,:))+squeeze(reconstructedData(86,27,:))+squeeze(reconstructedData(87,35,:))+squeeze(reconstructedData(84,37,:))+squeeze(reconstructedData(81,36,:))+squeeze(reconstructedData(80,34,:))+squeeze(reconstructedData(86,31,:))
greenline = squeeze(reconstructedData(88,11,:))+squeeze(reconstructedData(86,10,:))+squeeze(reconstructedData(87,13,:))+squeeze(reconstructedData(85,12,:))+squeeze(reconstructedData(57,25,:))+squeeze(reconstructedData(59,31,:))+squeeze(reconstructedData(63,37,:))+squeeze(reconstructedData(86,61,:))+squeeze(reconstructedData(85,59,:))+squeeze(reconstructedData(81,56,:))+squeeze(reconstructedData(77,15,:))+squeeze(reconstructedData(76,13,:))+squeeze(reconstructedData(75,14,:))+squeeze(reconstructedData(74,49,:))+squeeze(reconstructedData(75,56,:))+squeeze(reconstructedData(79,51,:))+squeeze(reconstructedData(80,50,:))+squeeze(reconstructedData(64,40,:))+squeeze(reconstructedData(64,30,:))+squeeze(reconstructedData(64,24,:))+squeeze(reconstructedData(86,56,:))+squeeze(reconstructedData(64,33,:))+squeeze(reconstructedData(64,34,:))+squeeze(reconstructedData(65,25,:))+squeeze(reconstructedData(65,34,:))
redline = squeeze(reconstructedData(79,10,:))+squeeze(reconstructedData(76,9,:))+squeeze(reconstructedData(78,10,:))+squeeze(reconstructedData(77,9,:))+squeeze(reconstructedData(60,13,:))+squeeze(reconstructedData(60,14,:))+squeeze(reconstructedData(59,15,:))+squeeze(reconstructedData(58,19,:))+squeeze(reconstructedData(56,21,:))+squeeze(reconstructedData(63,44,:))+squeeze(reconstructedData(56,62,:))+squeeze(reconstructedData(74,65,:))+squeeze(reconstructedData(77,64,:))+squeeze(reconstructedData(83,63,:))+squeeze(reconstructedData(55,61,:))+squeeze(reconstructedData(57,61,:))+squeeze(reconstructedData(55,64,:))+squeeze(reconstructedData(55,65,:))+squeeze(reconstructedData(62,63,:))+squeeze(reconstructedData(65,64,:))+squeeze(reconstructedData(68,63,:))+squeeze(reconstructedData(70,65,:))+squeeze(reconstructedData(58,15,:))+squeeze(reconstructedData(59,18,:))+squeeze(reconstructedData(55,23,:))

avgblueline = blueline/25
avggreenline = greenline/25
avgredline = redline/25

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
xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Average Spectra")