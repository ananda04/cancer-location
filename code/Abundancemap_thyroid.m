load('ReconResults_Thyroid_300iter_M3_Try1.mat')
load('thyroid.mat')
figure(1); imshow(AA)
mask = im2bw(AA,0)
figure(1); imshow(mask)

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
w = [w,thyroid_dotprods(pixel)]
%%
t =  reshape(w, [217 2256])
z = thyroid_abundvals(t)
%%
mask = []
mask1 = []
for s1 = 1:2256
    mask(r(s1),c(s1),:) = [z(1,s1),0,0]
    mask1(r(s1),c(s1),:) = [0,z(2,s1),0]
end
figure(2);imshow(mask)
figure(4);imshow(mask1)
figure(11); imshow(mask+mask1)