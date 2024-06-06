% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

load('ReconResults_Thyroid_300iter_M3_Try1.mat')
load('thyroid.mat')
%load('gtdatanoncancer.mat')
load("Thyroid_noncancer_regions.mat")
%Mask
figure(1); imshow(AA)
mask = im2bw(AA,0)
figure(1); imshow(mask)
% collect all data 
[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude

%gt known noncancerous data
%healthspec = []
% l = length(noncancerthy)
%for s1 = 1:l
 %   healthspec = cat(2, healthspec, squeeze(reconstructedData(noncancerthy(s1,2),noncancerthy(s1,1),:)))
%end

[rn cn] = find(noncancer == 1)
l = length(rn)
healthspec = []
for s1 = 1:l
    healthspec = cat(2, healthspec, squeeze(reconstructedData(rn(s1),cn(s1),:)))
end 
magavghealth = sqrt(sum(healthspec.^2))
normhealth = healthspec./magavghealth
%%
L = length(normallspec)
similarity = []
for i = 1:L
    for j = 1:80
        dotprod = sum(normallspec(:,i).*normhealth(:,j))
    end 
    if dotprod > 0.92
       similarity = [similarity, 1]
    else
        similarity = [similarity, 2]
    end 
end 
%%
similar1 = transpose(similarity)
figure(1); imshow(mask)
hold on;
for f1 = 1:L
    if similar1(f1) == 1
        figure(1); plot(c(f1),r(f1),'bx')
        hold on;
    elseif similar1(f1) == 2
        figure(1); plot(c(f1),r(f1),'rx')
        hold on;
    end
end


