% Description: This script uses UMAP to create graphs of cancer, white and gray matter spectra points from GBMT reconstructedData.
%              Creating graphs for Cancer, White and Gray Matter spectra points
%              Cancer, White and Gray Matter regions were determined through Dr.Greenberg & coordinating Oncologist.
% Author:      Arnav Nanda 
%              Apex High School 
%              Fitzpatrick Center for Interdisciplinary Engineering, Medicine and Applied Science (CIEMAS), Duke University

load('ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat')

% Creating Cancer Array
A = [squeeze(reconstructedData(79,23,:)),squeeze(reconstructedData(73,23,:)),squeeze(reconstructedData(82,23,:)),squeeze(reconstructedData(81,21,:)),squeeze(reconstructedData(79,19,:)),squeeze(reconstructedData(79,22,:)),squeeze(reconstructedData(78,25,:)),squeeze(reconstructedData(78,27,:)),squeeze(reconstructedData(76,26,:)),squeeze(reconstructedData(77,26,:)),squeeze(reconstructedData(80,33,:)),squeeze(reconstructedData(75,32,:)),squeeze(reconstructedData(74,30,:)),squeeze(reconstructedData(83,34,:)),squeeze(reconstructedData(79,33,:)),squeeze(reconstructedData(77,33,:)),squeeze(reconstructedData(74,31,:)),squeeze(reconstructedData(77,30,:)),squeeze(reconstructedData(79,28,:)),squeeze(reconstructedData(86,27,:)),squeeze(reconstructedData(87,35,:)),squeeze(reconstructedData(84,37,:)),squeeze(reconstructedData(81,36,:)),squeeze(reconstructedData(80,34,:)),squeeze(reconstructedData(86,31,:))]
cancerarray = transpose(A)

% Creating White Matter array
B = [squeeze(reconstructedData(88,11,:)),squeeze(reconstructedData(86,10,:)),squeeze(reconstructedData(87,13,:)),squeeze(reconstructedData(85,12,:)),squeeze(reconstructedData(57,25,:)),squeeze(reconstructedData(59,31,:)),squeeze(reconstructedData(63,37,:)),squeeze(reconstructedData(86,61,:)),squeeze(reconstructedData(85,59,:)),squeeze(reconstructedData(81,56,:)),squeeze(reconstructedData(77,15,:)),squeeze(reconstructedData(76,13,:)),squeeze(reconstructedData(75,14,:)),squeeze(reconstructedData(74,49,:)),squeeze(reconstructedData(75,56,:)),squeeze(reconstructedData(79,51,:)),squeeze(reconstructedData(80,50,:)),squeeze(reconstructedData(64,40,:)),squeeze(reconstructedData(64,30,:)),squeeze(reconstructedData(64,24,:)),squeeze(reconstructedData(86,56,:)),squeeze(reconstructedData(64,33,:)),squeeze(reconstructedData(64,34,:)),squeeze(reconstructedData(65,25,:)),squeeze(reconstructedData(65,34,:))]

whitematterarray = transpose(B) 

% Creating Gray Matter array
C = [squeeze(reconstructedData(79,10,:)),squeeze(reconstructedData(76,9,:)),squeeze(reconstructedData(78,10,:)),squeeze(reconstructedData(77,9,:)),squeeze(reconstructedData(60,13,:)),squeeze(reconstructedData(60,14,:)),squeeze(reconstructedData(59,15,:)),squeeze(reconstructedData(58,19,:)),squeeze(reconstructedData(56,21,:)),squeeze(reconstructedData(63,44,:)),squeeze(reconstructedData(56,62,:)),squeeze(reconstructedData(74,65,:)),squeeze(reconstructedData(77,64,:)),squeeze(reconstructedData(83,63,:)),squeeze(reconstructedData(55,61,:)),squeeze(reconstructedData(57,61,:)),squeeze(reconstructedData(55,64,:)),squeeze(reconstructedData(55,65,:)),squeeze(reconstructedData(62,63,:)),squeeze(reconstructedData(65,64,:)),squeeze(reconstructedData(68,63,:)),squeeze(reconstructedData(70,65,:)),squeeze(reconstructedData(58,15,:)),squeeze(reconstructedData(59,18,:)),squeeze(reconstructedData(55,23,:))]
graymatterarray = transpose(C)

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
legend('Cancer', 'White Matter','Gray Matter' );
title('Transpose')