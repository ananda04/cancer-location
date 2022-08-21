load('ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat')
load("reconNormal.mat")
load('e2.mat')
load("u1.mat")
figure(1);imshow(e2)
figure(2);imshow(u1)


% Create a mask of the background only for cancerous slice (e2)
redChannel = e2(:, :, 1);
greenChannel = e2(:, :, 2);
blueChannel = e2(:, :, 3);

mask = blueChannel > 200, redChannel > 200, greenChannel > 200;
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


similarity = []
for i = 1:L
    for j = 1:l
        dotprod = sum(normhealth(:,j).*normallspec(:,i))
    end 
    if dotprod > 0.62
       similarity = [similarity, 1]
    elseif dotprod>0.6 && dotprod <0.62
        similarity = [similarity, 2]
    else 
        similarity = [similarity, 3]
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
        figure(3); plot(c(f1),r(f1),'gx')
        hold on;
    else
        figure(3); plot(c(f1),r(f1),'rx')
        hold on;
    end
end
%%
cancercheck = []
for f1 = 1:l
    if similar1(f1) ==3
        cancercheck = cat(2, cancercheck, squeeze(reconstructedData3(r(f1),c(f1),:)))
    end
end
figure(11); plot(qvals, cancercheck)
