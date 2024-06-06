% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

%% load Data: t1_158
load('ReconResults_Brain_158_T1_20s_300iter_M3_Try1.mat')
load("t1_158.mat")
% Mask t2_177
redChannel = T1158(:, :, 1);
greenChannel = T1158(:, :, 2);
blueChannel = T1158(:, :, 3);

mask = blueChannel > 250, redChannel > 250, greenChannel > 250;
figure(3); imshow(mask)
[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude

load('ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat')
load("NT_176.mat")
% Mask NT_176
redChannel = NT176(:, :, 1);
greenChannel = NT176(:, :, 2);
blueChannel = NT176(:, :, 3);

mask6 = blueChannel > 245, redChannel > 250, greenChannel > 248;
figure(2); imshow(mask6)
hold on;
[x2 y2] = find(mask6 ==0)
l = length(x2)
healallspec2 = []
for k1 = 1:l
    healallspec2 = cat(2, healallspec2, squeeze(reconstructedData(x2(k1),y2(k1),:)))
end
magnitude2 = sqrt(sum(healallspec2.^2))
normhealth2 = healallspec2./magnitude2

%% Similarity
similarity = []
for i = 1:L
    for j = 1:l
        dotprod = sum(normhealth2(:,j).*normallspec(:,i))
    end 
    if dotprod > 0.95
       similarity = [similarity, 1]
    else
       similarity = [similarity, 2]
    end 
end 

similar1 = transpose(similarity)
figure(3); imshow(mask)
hold on;
for f1 = 1:L
    if similar1(f1) == 1
        figure(3); plot(c(f1),r(f1),'bx')
        hold on;
    elseif similar1(f1) == 2
        figure(3); plot(c(f1),r(f1),'rx')
        hold on;
    end
end
