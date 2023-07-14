load("BrainData_andErrorBars.mat")
magnitudeB = sqrt(sum(Cancer_Spec.^2))
normB = Cancer_Spec/magnitudeB
load("ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat")
load("T1187_data.mat")
    l = length(r)
    healallspec1 = []
    for k = 1:l
        healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(r(k),c(k),:)))
    end
    magnitude1 = sqrt(sum(healallspec1.^2))
    normT1187 = healallspec1./magnitude1
    qvals2 = 0.05:0.45/460:0.5
    normT1187 =interp1(qvals,normT1187,qvals2)
%%
    a = []
    for j = 1:l
        a = cat(2,a, sum(normB.*normT1187(:,j)))
    end