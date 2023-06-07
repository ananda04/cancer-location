% Description: This script plot spatial points onto reconstructed image from GBMT reconstructedData.
%              Cancer, White and Gray Matter regions were determined through Dr.Greenberg & coordinating Oncologist.
% Author:      Arnav Nanda 
%              Apex High School 
%              Fitzpatrick Center for Interdisciplinary Engineering, Medicine and Applied Science (CIEMAS), Duke University  

load('ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat')

% Cancer data
x = [23,23,21,19,22,25,27,26,26,33,32,30,23]
y = [79,82,81,79,79,78,78,76,77,80,75,74,73]

% White Matter data
x1 = [11,10,13,12,15,13,14,25,31,37,61,59,56]
y1 = [88,86,87,85,77,76,75,57,59,63,86,85,81]

% Gray Matter data
x2 = [10, 9,10, 9,13,14,15,19,21,44,61,65,64]
y2 = [79,76,78,77,60,60,59,58,56,63,56,74,77]

% Picture
figure(4); imagesc(squeeze(reconstructedData(:, :, 17)))
hold on; figure(4); plot(x,y,'rx')
hold on; figure(4); plot(x1,y1, 'gx')
hold on; figure(4); plot(x2,y2, 'kx')
legend('Cancer', 'White Matter','Gray Matter' );
title("Spatial Locations")