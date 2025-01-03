% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

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
%%
%g = sgolayfilt(0,3)
%allspec1 = []
%for s1 = 1:2256
%    allspec1 = [allspec1, conv(allspec(:,s1),g(:,1).','same')]
%end 
magnitude = sqrt(sum(allspec.^2))
pixel = allspec./magnitude

%% mapping
w= []
w = [w,dotprods(pixel)]
%%
t =  reshape(w, [217 2256])
z = abundvals(t)
%%
mask = []
mask1 = []
mask2 = []
for s1 = 1:2256
    mask(r(s1),c(s1),:) = [z(1,s1),0,0]
    mask1(r(s1),c(s1),:) = [0,z(2,s1),0]
    mask2(r(s1),c(s1),:) = [0,0,z(3,s1)]
end
figure(2);imshow(mask)
figure(4);imshow(mask1)
figure(5);imshow(mask2)
figure(11); imshow(mask+mask1+mask2)
