load('ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat')
load('e2.mat')
figure(1);imshow(e2)
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
pixel = allspec./magnitude

%% mapping
w= []
for s1 = 1:2256
    w = [w,abundancefunc(pixel(:,s1))]
end 
%%
figure(1); imshow(mask)
hold on;
for s1 =1:2256
    if w(1,s1) == 1
        figure(1);plot(c(s1),r(s1),'rx')
        hold on;
    elseif  w(2,s1) == 1
        figure(1);plot(c(s1),r(s1), 'bx')
        hold on;
    else  w(3,s1) == 1
        figure(1);plot(c(s1),r(s1), 'gx')
        hold on;
    end 
end
%%
mask = []
for s1 = 1:2256
    mask(r(s1),c(s1),:) = [w(1,s1),w(2,s1),w(3,s1)]
end
figure(2);imshow(mask)




