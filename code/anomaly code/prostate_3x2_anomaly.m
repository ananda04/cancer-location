% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

load('ReconResults_Prostate_V3_300iter_M3_Try1.mat')

% Create a mask of the background only for prostate slies
mask = imbinarize(squeeze(reconstructedData(:, :, 28)))


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

similar1 = transpose(similar)
figure(3); imshow(mask)
hold on;
for f1 = 1:L
    if similar1(f1) == 1
        figure(3); plot(c(f1),r(f1),'rx')
        hold on;
    elseif similar1(f1) == 2
        figure(3); plot(c(f1),r(f1),'bx')
        hold on;
    elseif similar1(f1) == 3
        figure(3); plot(c(f1),r(f1),'yx')
        hold on;
    end
end
legend(['rx','bx','yx'],'Cancer', 'White Matter','Gray Matter' )

figure(4); imagesc(squeeze(reconstructedData(:, :, 22)))
hold on;
figure(4); plot(c,r,'gx')
