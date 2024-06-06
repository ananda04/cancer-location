% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

load("T1187_data.mat")
load("ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat")
[mask x2 y2] = bmask("NT_187")
    l = length(x2)
    healallspec = []
    for k = 1:l
        healallspec = cat(2, healallspec, squeeze(reconstructedData(x2(k),y2(k),:)))
    end
    magnitude2 = sqrt(sum(healallspec.^2))
    normT1187 = healallspec./magnitude2

    a = dot(Cancer_spec,normT1187)
