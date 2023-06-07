% Description: This script uses UMAP to create graphs of cancer, white and gray matter spectra points from GBMT reconstructedData.
%              Creating graphs for Cancer, White and Gray Matter spectra points
%              Cancer, White and Gray Matter regions were determined through Dr.Greenberg & coordinating Oncologist.
% Author:      Arnav Nanda 
%              Apex High School 
%              Fitzpatrick Center for Interdisciplinary Engineering, Medicine and Applied Science (CIEMAS), Duke University

load('ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat')

% Creating Cancer Array
A = [squeeze(reconstructedData(79,23,:)),squeeze(reconstructedData(73,23,:)),squeeze(reconstructedData(82,23,:)),squeeze(reconstructedData(81,21,:)),squeeze(reconstructedData(79,19,:)),squeeze(reconstructedData(79,22,:)),squeeze(reconstructedData(78,25,:)),squeeze(reconstructedData(78,27,:)),squeeze(reconstructedData(76,26,:)),squeeze(reconstructedData(77,26,:)),squeeze(reconstructedData(80,33,:)),squeeze(reconstructedData(75,32,:)),squeeze(reconstructedData(74,30,:))]
cancerarray = transpose(A)

% Creating White Matter array
B = [squeeze(reconstructedData(58,15,:)),squeeze(reconstructedData(75,56,:)),squeeze(reconstructedData(63,25,:)),squeeze(reconstructedData(56,62,:)),squeeze(reconstructedData(64,40,:)),squeeze(reconstructedData(57,25,:)),squeeze(reconstructedData(64,33,:)),squeeze(reconstructedData(86,56,:)),squeeze(reconstructedData(80,50,:)),squeeze(reconstructedData(64,30,:)),squeeze(reconstructedData(64,24,:)),squeeze(reconstructedData(64,34,:)),squeeze(reconstructedData(79,51,:))]

whitematterarray = transpose(B) 

% Creating Gray Matter array
C = [squeeze(reconstructedData(55,23,:)),squeeze(reconstructedData(57,61,:)),squeeze(reconstructedData(62,63,:)),squeeze(reconstructedData(70,65,:)),squeeze(reconstructedData(59,18,:)),squeeze(reconstructedData(55,65,:)),squeeze(reconstructedData(63,44,:)),squeeze(reconstructedData(61,55,:)),squeeze(reconstructedData(51,64,:)),squeeze(reconstructedData(58,15,:)),squeeze(reconstructedData(68,63,:)),squeeze(reconstructedData(83,63,:)),squeeze(reconstructedData(65,64,:))]
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