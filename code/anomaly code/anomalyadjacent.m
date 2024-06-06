% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

load('ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat')
load("reconNormal.mat")
load('e2.mat')
load("allspec.mat")
load("healthyspectra.mat")

redChannel = e2(:, :, 1);
greenChannel = e2(:, :, 2);
blueChannel = e2(:, :, 3);

mask = blueChannel > 200, redChannel > 200, greenChannel > 200;
figure(3); imshow(mask)
hold on;
L = length(normallspec)
l = length(normhealth)

similarity = []
for i = 1:L
    for j = 1:l
        dotprod = sum(normhealth(:,j).*normallspec(:,i))
    end 
    if dotprod > 0.63
       similarity = [similarity, 1]
    elseif dotprod>0.6 && dotprod <0.62
        similarity = [similarity, 2]
    else 
        similarity = [similarity, 3]
    end 
end 

similar1 = transpose(similarity)
figure(3); imshow(mask)
[r c] = find(mask == 0)
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
