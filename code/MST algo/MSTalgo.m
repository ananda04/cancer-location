% Cancer Image 
A = imread('HnE.png');
lab_hne = rgb2lab(A);
ab = lab_hne(:,:,2:3);
ab = im2single(ab);

numColors = 3;
L = imsegkmeans(ab,numColors)
B = labeloverlay(A,L);
imshow(B)
% filted_I = CoherenceFilter(A,struct('T',5,'rho',2,'Scheme','I', 'sigma', 1));
% L = graphSeg(filted_I, 0.5, 50, 2, 0);
% Lnn = graphSeg(filted_I, 0.5/sqrt(3), 50, 10, 1);
% subplot(3, 1, 1), imshow(I_gray, []), title('original image');
% subplot(3, 1, 2), imshow(label2rgb(L)), title('adjacent neighborhood based segmentation');
% subplot(3, 1, 3), imshow(label2rgb(Lnn)), title('k nearest neighborhood based segmentation');
