% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

load('ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat') 

cancerarray = [squeeze(reconstructedData(79,23,:)),
     squeeze(reconstructedData(82,23,:)),
     squeeze(reconstructedData(81,21,:)),
     squeeze(reconstructedData(79,19,:)),
     squeeze(reconstructedData(79,22,:)),
     squeeze(reconstructedData(78,25,:)),
     squeeze(reconstructedData(78,27,:)),
     squeeze(reconstructedData(76,26,:)),
     squeeze(reconstructedData(77,26,:)),
     squeeze(reconstructedData(80,33,:)),
     squeeze(reconstructedData(75,32,:)),
     squeeze(reconstructedData(74,30,:)),
     squeeze(reconstructedData(73,23,:))]

% Drawing White Matter spectra lines 
whitematterarray = [squeeze(reconstructedData(58,15,:)), 
 squeeze(reconstructedData(59,18,:)),
 squeeze(reconstructedData(55,65,:)),
 squeeze(reconstructedData(56,62,:)),
 squeeze(reconstructedData(65,64,:)),
 squeeze(reconstructedData(57,25,:)),
 squeeze(reconstructedData(57,65,:)),
 squeeze(reconstructedData(55,61,:)),
 squeeze(reconstructedData(68,63,:)),
 squeeze(reconstructedData(62,63,:)),
 squeeze(reconstructedData(57,61,:)),
 squeeze(reconstructedData(57,64,:)),
 squeeze(reconstructedData(70,65,:))]
%hold on;figure(2); plot(qvals, squeeze(reconstructedData(63,37,:)))

% Drawing Gray Matter spectra lines 
 graymatterarray = [squeeze(reconstructedData(64,40,:)), 
 squeeze(reconstructedData(75,50,:)),
 squeeze(reconstructedData(79,51,:)),
 squeeze(reconstructedData(81,49,:)),
 squeeze(reconstructedData(86,49,:)),
 squeeze(reconstructedData(78,48,:)),
 squeeze(reconstructedData(63,44,:)), 
 squeeze(reconstructedData(74,49,:)), 
 squeeze(reconstructedData(78,44,:)),
 squeeze(reconstructedData(82,50,:)),
 squeeze(reconstructedData(81,48,:)),
 squeeze(reconstructedData(83,63,:)), 
 squeeze(reconstructedData(80,50,:))] 

% Making matrices
%largearray = [cancerarray;whitematterarray;graymatterarray]

% RUN UMAP
import umapFileExchange.*;
[reduction] = run_umap(cancerarray)
[reduction1] = run_umap(whitematterarray)
[reduction2] = run_umap(graymatterarray)
figure(4); plot(reduction(:,1), reduction(:,2),'r.')
hold on;
figure(4); plot(reduction1(:,1), reduction1(:,2),'g.')
hold on;
figure(4); plot(reduction2(:,1), reduction2(:,2),'k.')
legend('Cancer', 'White Matter','Gray Matter' )
