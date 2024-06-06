% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

%% load Data: t2_158
load('ReconResults_Brain_158_T2_30s_300iter_M3_Try1.mat')
load("t2_158.mat")

% Mask t2_187
redChannel = T2158(:, :, 1);
greenChannel = T2158(:, :, 2);
blueChannel = T2158(:, :, 3);

mask = blueChannel > 251, redChannel > 255, greenChannel > 254;
figure(3); imshow(mask)
hold on;

[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude

%% load Data: NT_158
load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')
load("NT_158.mat")
% Mask NT_187
redChannel = NT158(:, :, 1);
greenChannel = NT158(:, :, 2);
blueChannel = NT158(:, :, 3);

mask3 = blueChannel > 245, redChannel > 250, greenChannel > 248;
figure(2); imshow(mask3)
hold on;

[x y] = find(mask3 == 0)
l = length(x)
healallspec = []
for k1 = 1:l
    healallspec = cat(2, healallspec, squeeze(reconstructedData(x(k1),y(k1),:)))
end
magnitude = sqrt(sum(healallspec.^2))
normhealth = healallspec./magnitude
%% Similarity
similarity = []
for i = 1:L
    for j = 1:l
        dotprod = sum(normhealth(:,j).*normallspec(:,i))
    end 
    if dotprod > 0.7
       similarity = [similarity, 1]
    else
       similarity = [similarity, 2]
    end 
end 
%%
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

