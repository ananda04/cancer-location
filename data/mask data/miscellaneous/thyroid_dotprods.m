function [t] = thyroid_dotprods(pixel)
    load('ReconResults_Thyroid_300iter_M3_Try1.mat')
    load("Thyroid_noncancer_regions.mat")
    load('Thyroid_Cancer_regions.mat')

   %gt known noncancerous data
healthspec = []
l = length(noncancerthy)
for s1 = 1:l
    healthspec = cat(2, healthspec, squeeze(reconstructedData(noncancerthy(s1,2),noncancerthy(s1,1),:)))
end 
healthspec = sum(healthspec)
L = length(Cancer)
cancerspec = []
for s1 = 1:L
    cancerspec = cat(2, cancerspec, squeeze(reconstructedData(Cancer(s1,2),Cancer(s1,1),:)))
end 
cancerspec = sum(cancerspec)

    
    %average spectra
    avgcancer = mean(cancerspec)
    avghealth = mean(healthspec)
    
    %magnitude
    magavgcancer = sqrt(sum(avgcancer.^2))
    magavghealth = sqrt(sum(avghealth.^2))
    
    %normalized
    normavgcancer = avgcancer./magavgcancer
    normavghealth = avghealth./magavghealth

    % denoising 
    %denoise = wdenoise(normallspec)
    %denoise(denoise<0)=0
    %figure(6);plot(qvals, denoise,'r')

    % abundace

    i1 = []
    j1 = []
    for i = 0:0.1:1
        for j = 0:0.1:1
            if i+j == 1 
                i1 = [i1, i]
                j1 = [j1, j]
            end 
        end 
    end

    z = []
    l = length(i1)
    for m = 1:l
        z  = cat(2, z, (i1(:,m).*normavgcancer+j1(:,m).*normavghealth))
    end 

     t = []
    for b = 1:l
        t = cat(2,t, sum(pixel.*z(:,b)))
    end 
end