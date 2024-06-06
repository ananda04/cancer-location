% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

load('ReconResults_Prostate_V2_300iter_M3_Try1.mat')
import umapFileExchange.*;

% Create a mask of the background only for prostate slies
mask = imbinarize(squeeze(reconstructedData(:, :, 22)))

%load('prostate.mat')
%figure(1);imshow(pro)
%redChannel = pro(:, :, 1);
%greenChannel = pro(:, :, 2);
%blueChannel = pro(:, :, 3);
%mask = blueChannel > 254, redChannel > 254, greenChannel > 254;
%figure(3); imshow(mask)
%hold on;

% collect all data 
[r c] = find(mask == 1)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end
magnitude = sqrt(sum(allspec.^2))
normallspec = allspec./magnitude

%library values 
cancer = squeeze(reconstructedData(73,87,:))+squeeze(reconstructedData(76,88,:))+squeeze(reconstructedData(74,86,:))+squeeze(reconstructedData(76,87,:))+squeeze(reconstructedData(74,49,:))+squeeze(reconstructedData(71,47,:))+squeeze(reconstructedData(76,53,:))+squeeze(reconstructedData(77,51,:))+squeeze(reconstructedData(76,86,:))+squeeze(reconstructedData(78,88,:))+squeeze(reconstructedData(72,49,:))+squeeze(reconstructedData(76,84,:))+squeeze(reconstructedData(75,84,:))+squeeze(reconstructedData(77,50,:))+squeeze(reconstructedData(76,49,:))+squeeze(reconstructedData(75,48,:))+squeeze(reconstructedData(73,84,:))+squeeze(reconstructedData(73,85,:))+squeeze(reconstructedData(76,85,:))+squeeze(reconstructedData(86,27,:))+squeeze(reconstructedData(73,51,:))+squeeze(reconstructedData(76,51,:))+squeeze(reconstructedData(76,48,:))+squeeze(reconstructedData(73,52,:))+squeeze(reconstructedData(72,84,:))
nocancer = squeeze(reconstructedData(51,65,:))+squeeze(reconstructedData(53,36,:))+squeeze(reconstructedData(56,104,:))+squeeze(reconstructedData(39,130,:))+squeeze(reconstructedData(108,113,:))+squeeze(reconstructedData(74,26,:))+squeeze(reconstructedData(69,62,:))+squeeze(reconstructedData(63,103,:))+squeeze(reconstructedData(57,126,:))+squeeze(reconstructedData(111,117,:))+squeeze(reconstructedData(73,22,:))+squeeze(reconstructedData(68,55,:))+squeeze(reconstructedData(66,100,:))+squeeze(reconstructedData(40,132,:))+squeeze(reconstructedData(75,25,:))+squeeze(reconstructedData(109,115,:))+squeeze(reconstructedData(50,66,:))+squeeze(reconstructedData(47,97,:))+squeeze(reconstructedData(39,126,:))+squeeze(reconstructedData(90,109,:))+squeeze(reconstructedData(60,15,:))+squeeze(reconstructedData(49,59,:))+squeeze(reconstructedData(48,91,:))+squeeze(reconstructedData(51,102,:))+squeeze(reconstructedData(106,109,:))
calcification = squeeze(reconstructedData(54,92,:))+squeeze(reconstructedData(53,92,:))+squeeze(reconstructedData(54,92,:))+squeeze(reconstructedData(52,93,:))+squeeze(reconstructedData(53,93,:))+squeeze(reconstructedData(54,93,:))+squeeze(reconstructedData(54,92,:))+squeeze(reconstructedData(51,93,:))+squeeze(reconstructedData(50,94,:))+squeeze(reconstructedData(49,95,:))
l = length(calcification)
l1 = length(cancer)
l2 = length(nocancer)

% Average
avgcancer = cancer/l1
avgnocancer = nocancer/l2
avggcalcification = calcification/l

%magnitude
magavgcancer = sqrt(sum(avgcancer.^2))
magavgnocancer = sqrt(sum(avgnocancer.^2))
magavgcalc = sqrt(sum(avggcalcification.^2))

%normalized
normavgcancer = avgcancer./magavgcancer
normavgnocancer = avgnocancer./magavgnocancer
normavgcalc = avggcalcification./magavgcalc

similar = []
for d1 = 1:L
    simcancer = dot(normallspec(:,d1), normavgcancer)
    simnocancer = dot(normallspec(:,d1),normavgnocancer)
    simcalc = dot(normallspec(:,d1),normavgcalc)

    if simcancer > simnocancer
        if simcancer > simcalc
           similar = [similar, 1]
        else
        end
    end
    
    if simnocancer > simcancer
        if simnocancer > simcalc
           similar = [similar, 2]
        else
        end
    end
    
    if simcalc > simcancer
        if simcalc > simnocancer
           similar = [similar, 3]
        else
        end
    end
end
xcan = []
ycan = []
xnocan= []
ynocan = []
xcalc = []
ycalc = []
similar1 = transpose(similar)
figure(3); imshow(mask)
hold on;
for f1 = 1:L
    if similar1(f1) == 1
        figure(3); plot(c(f1),r(f1),'rx')
        xcan = [xcan,c(f1)]
        ycan = [ycan, r(f1)]
        hold on;
    elseif similar1(f1) == 2
        figure(3); plot(c(f1),r(f1),'bx')
        xnocan = [xnocan,c(f1)]
        ynocan = [ynocan, r(f1)]
        hold on;
    elseif similar1(f1) == 3
        figure(3); plot(c(f1),r(f1),'yx')
        xcalc = [xcalc,c(f1)]
        ycalc = [ycalc,r(f1)]
        hold on;
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
l = length(xcalc)
calcspec = []
for k1 = 1:l
    calcspec = cat(2, calcspec, squeeze(reconstructedData(ycalc(k1),xcalc(k1),:)))
end

%% UMAP reduction 
import umapFileExchange.*;
[reduction] = run_umap(canspec)
[reduction1] = run_umap(nocanspec)
[reduction2] = run_umap(calcspec)
figure(4); plot(reduction(:,1), reduction(:,2),'r.')
hold on;
figure(4); plot(reduction1(:,1), reduction1(:,2),'b.')
hold on;
figure(4); plot(reduction2(:,1), reduction2(:,2),'y.')
legend('Cancer', 'No cancer','Calcification' )

%%
canspec = transpose(canspec)
nocanspec = transpose(nocanspec)
calcspec = transpose(calcspec)

avgcanspec = mean(canspec)
avgnocanspec = mean(nocanspec)
avgcalcspec = mean(calcspec)

errblue = std(canspec)
errgreen = std(nocanspec)
errred = std(calcspec)

% Plot graphs with error bars
figure(5); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
e1 = errorbar(qvals, avgcanspec,errblue); 
e1.Color = 'r'
hold on;
figure(5); e2 = errorbar(qvals, avgnocanspec,errgreen); 
e2.Color = 'b'
hold on;
figure(5); e3 = errorbar(qvals, avgcalcspec, errred); 
e3.Color = 'y'
legend([e1(1),e2(1),e3(1)],'Cancer', 'No cancer','Calcification' )



