%% load Data: t2_177
load('ReconResults_Brain_177_T2_30s_300iter_M3_Try1.mat')
load("t2_177.mat")

% Mask t2_177
redChannel = T2177(:, :, 1);
greenChannel = T2177(:, :, 2);
blueChannel = T2177(:, :, 3);

mask = blueChannel > 254, redChannel > 254, greenChannel > 254;
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

mask5 = blueChannel > 245, redChannel > 250, greenChannel > 248;
figure(2); imshow(mask5)
hold on;
[x1 y1] = find(mask5 ==0)
l = length(x1)
healallspec1 = []
for k1 = 1:l
    healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(x1(k1),y1(k1),:)))
end
magnitude1 = sqrt(sum(healallspec1.^2))
normhealth1 = healallspec1./magnitude1
%% Similarity
similarity = []
for i = 1:L
    for j = 1:l
        dotprod = sum(normhealth1(:,j).*normallspec(:,i))
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

