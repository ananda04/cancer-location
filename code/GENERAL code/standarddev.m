% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

A = [squeeze(reconstructedData(87,11,:)),squeeze(reconstructedData(89,15,:)),squeeze(reconstructedData(81,20,:)),squeeze(reconstructedData(85,20,:)),squeeze(reconstructedData(88,19,:)),squeeze(reconstructedData(89,14,:)),squeeze(reconstructedData(72,28,:)),squeeze(reconstructedData(84,25,:)),squeeze(reconstructedData(84,21,:)),squeeze(reconstructedData(89,16,:))]
B = transpose(A)
C = mean(B)
V = std(C)
