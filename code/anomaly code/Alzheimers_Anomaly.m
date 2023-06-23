load('ReconResults_AD_Brain_300iter_M3_BrainScan_7mm_V2_Try2.mat')
load("reconNormal.mat")
load('Alzheimers.mat')
load("u1.mat")
figure(1); imshow(AA)
figure(2);imshow(u1)

% Create a mask of the background only for Alzheimers slies
redChannel = AA(:, :, 1);
greenChannel = AA(:, :, 2);
blueChannel = AA(:, :, 3);

mask = blueChannel > 252, redChannel > 210, greenChannel > 210;
figure(3); imshow(mask)
hold on;

% collect all data 
[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude

%create a mask for non-cancerous slice 
redChannel = u1(:, :, 1);
greenChannel = u1(:, :, 2);
blueChannel = u1(:, :, 3);

mask1 = blueChannel > 250,  greenChannel > 250;
figure(4); imshow(mask1)
hold on;

[x y] = find(mask1 == 0)
l = length(x)
healallspec = []
for k1 = 1:l
    healallspec = cat(2, healallspec, squeeze(reconstructedData3(x(k1),y(k1),:)))
end
magnitude = sqrt(sum(healallspec.^2))
normhealth = healallspec./magnitude

%% Similarity
similarity = []
for i = 1:l
    for j = 1:L
        dotprod = sum(normhealth(:,i).*normallspec(:,j))
    end 
    if dotprod > 0.7
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
