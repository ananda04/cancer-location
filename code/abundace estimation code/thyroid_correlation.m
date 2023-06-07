load('ReconResults_Thyroid_300iter_M3_Try1.mat')
load('thyroid.mat')

%Mask
figure(1); imshow(AA)
mask = im2bw(AA,0)
figure(1); imshow(mask)
% collect all data 
[r c] = find(mask = 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude


%library values 
cancer = squeeze(reconstructedData(38,37,:))+squeeze(reconstructedData(39,37,:))+squeeze(reconstructedData(41,37,:))+squeeze(reconstructedData(44,15,:))+squeeze(reconstructedData(44,14,:))+squeeze(reconstructedData(40,37,:))+squeeze(reconstructedData(84,10,:))+squeeze(reconstructedData(85,11,:))+squeeze(reconstructedData(76,86,:))+squeeze(reconstructedData(54,15,:))+squeeze(reconstructedData(53,9,:))+squeeze(reconstructedData(51,9,:))+squeeze(reconstructedData(55,14,:))+squeeze(reconstructedData(48,84,:))+squeeze(reconstructedData(44,67,:))+squeeze(reconstructedData(45,66,:))+squeeze(reconstructedData(49,86,:))+squeeze(reconstructedData(47,87,:))+squeeze(reconstructedData(76,85,:))+squeeze(reconstructedData(86,27,:))+squeeze(reconstructedData(44,14,:))+squeeze(reconstructedData(52,9,:))+squeeze(reconstructedData(50,90,:))+squeeze(reconstructedData(48,90,:))+squeeze(reconstructedData(51,89,:))
nocancer = squeeze(reconstructedData(46,51,:))+squeeze(reconstructedData(48,77,:))+squeeze(reconstructedData(67,74,:))+squeeze(reconstructedData(86,79,:))+squeeze(reconstructedData(83,54,:))+squeeze(reconstructedData(80,29,:))+squeeze(reconstructedData(66,28,:))+squeeze(reconstructedData(69,52,:))+squeeze(reconstructedData(61,70,:))+squeeze(reconstructedData(65,75,:))+squeeze(reconstructedData(46,52,:))+squeeze(reconstructedData(65,49,:))+squeeze(reconstructedData(72,48,:))+squeeze(reconstructedData(68,44,:))+squeeze(reconstructedData(75,25,:))+squeeze(reconstructedData(63,70,:))+squeeze(reconstructedData(75,71,:))+squeeze(reconstructedData(88,49,:))+squeeze(reconstructedData(85,44,:))+squeeze(reconstructedData(66,43,:))+squeeze(reconstructedData(69,27,:))+squeeze(reconstructedData(96,53,:))+squeeze(reconstructedData(71,74,:))+squeeze(reconstructedData(88,84,:))+squeeze(reconstructedData(61,38,:))
l1 = length(cancer)
l2 = length(nocancer)

% Average
avgcancer = cancer/l1
avgnocancer = nocancer/l2

%magnitude
magavgcancer = sqrt(sum(avgcancer.^2))
magavgnocancer = sqrt(sum(avgnocancer.^2))

%normalized
normavgcancer = avgcancer./magavgcancer
normavgnocancer = avgnocancer./magavgnocancer

similar = []
for d1 = 1:L
    simcancer = dot(normallspec(:,d1), normavgcancer)
    simnocancer = dot(normallspec(:,d1),normavgnocancer)

    if simcancer > simnocancer
       similar = [similar, 1]
    else
    end
    if simnocancer > simcancer
       similar = [similar, 2]
    else
    end
end 

similar1 = transpose(similar)
%%
xcan = []
ycan = []
xnocan = []
ynocan = []
figure(3); imshow(AA)
hold on;
for f1 = 1:L
    if similar1(f1) == 1
        figure(3); plot(c(f1),r(f1),'rx')
        xcan = [xcan, c(f1)]
        ycan = [ycan, r(f1)]
        hold on;
    elseif similar1(f1) == 2
        figure(3); plot(c(f1),r(f1),'bx')
        hold on;
        xnocan = [xnocan, c(f1)]
        ynocan = [ynocan, r(f1)]
    end
end
l = length(xcan)
canspec = []
for k1 = 1:l
    canspec = cat(2, canspec, squeeze(reconstructedData(ycan(k1),xcan(k1),:)))
end
l = length(xnocan)
nocanspec = []
for k1 = 1:l
    nocanspec = cat(2, nocanspec, squeeze(reconstructedData(ynocan(k1),xnocan(k1),:)))
end
%%
canspec = transpose(canspec)
nocanspec = transpose(nocanspec)

avgcanspec = mean(canspec)
avgnocanspec = mean(nocanspec)

errblue = std(canspec)
errgreen = std(nocanspec)

% Plot graphs with error bars
figure(5); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
hold on;
e1 = errorbar(qvals, avgcanspec,errblue); 
e1.Color = 'r'
hold on;
figure(5); e2 = errorbar(qvals, avgnocanspec,errgreen); 
e2.Color = 'b'
hold on;
legend([e1(1),e2(1)],'Cancer', 'No cancer')

% UMAP reduction 

[reduction] = run_umap(canspec)
[reduction1] = run_umap(nocanspec)
figure(4); plot(reduction(:,1), reduction(:,2),'r.')
hold on;
figure(4); plot(reduction1(:,1), reduction1(:,2),'b.')
hold on;
legend('Cancer','No cancer')



