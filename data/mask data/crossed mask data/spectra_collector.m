load("ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat")
load("gbmt_annotated_mask.mat")
load("noncan.mat")
[r c] = find(mask == 1)
noncanSpec = ReconResults(reconstructedData, r,c)
import umapFileExchange.*;


%25 points
normSpec = ReconResults(reconstructedData, r,c)
sample = datasample(normSpec,25,2)
tsample = transpose(sample)
stand = std(tsample)
norm = mean(tsample,1)

Nsample = datasample(noncanSpec,25,2)
tNsample = transpose(Nsample)
Nstd = std(tNsample)
Nnorm = mean(tNsample,1)

[reduction] = run_umap(tsample)
[reduction1] = run_umap(tNsample)
figure(1); plot(reduction(:,1), reduction(:,2),'r.')
hold on;
figure(1); plot(reduction1(:,1), reduction1(:,2),'b.')
hold on;
legend('Cancer', 'Non-cacner' );
title('25 Points')

% 50 points
sample1 = datasample(normSpec,50,2)
tsample1 = transpose(sample1)
std1 = std(tsample1)
norm1 = mean(tsample1,1)

Nsample1 = datasample(noncanSpec,50,2)
tNsample1 = transpose(Nsample1)
Nstd1 = std(tNsample1)
Nnorm1 = mean(tNsample1,1)

[reduction] = run_umap(tsample1)
[reduction1] = run_umap(tNsample1)
figure(2); plot(reduction(:,1), reduction(:,2),'r.')
hold on;
figure(2); plot(reduction1(:,1), reduction1(:,2),'b.')
hold on;
legend('Cancer', 'Non-cacner' );
title('50 Points')

% 75 points
sample2 = datasample(normSpec,75,2)
tsample2 = transpose(sample2)
std2 = std(tsample2)
norm2 = mean(tsample2,1)

Nsample2 = datasample(noncanSpec,75,2)
tNsample2 = transpose(Nsample2)
Nstd2 = std(tNsample2)
Nnorm2 = mean(tNsample2,1)

[reduction] = run_umap(tsample2)
[reduction1] = run_umap(tNsample2)
figure(3); plot(reduction(:,1), reduction(:,2),'r.')
hold on;
figure(3); plot(reduction1(:,1), reduction1(:,2),'b.')
hold on;
legend('Cancer', 'Non-cacner' );
title('75 Points')

% 100 points
sample3 = datasample(normSpec,100,2)
tsample3 = transpose(sample3)
std3 = std(tsample3)
norm3 = mean(tsample3,1)

Nsample3 = datasample(noncanSpec,100,2)
tNsample3 = transpose(Nsample3)
Nstd3 = std(tNsample3)
Nnorm3 = mean(tNsample3,1)

[reduction] = run_umap(tsample3)
[reduction1] = run_umap(tNsample3)
figure(4); plot(reduction(:,1), reduction(:,2),'r.')
hold on;
figure(4); plot(reduction1(:,1), reduction1(:,2),'b.')
hold on;
legend('Cancer', 'Non-cacner' );
title('100 Points')
%%
figure(5); e2 = errorbar(qvals, norm, stand); 
e2.Color = 'r'
hold on;
figure(5); e2 = errorbar(qvals, Nnorm, Nstd); 
e2.Color = 'k'

legend("Cancer", "Non-Cancer");
figure(5); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars with 25 points")

figure(6); e2 = errorbar(qvals, norm1, std1); 
e2.Color = 'r'
hold on;
figure(6); e2 = errorbar(qvals, Nnorm1, Nstd1); 
e2.Color = 'k'

legend("Cancer", "Non-Cancer");
figure(6); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars with 50 points")

figure(7); e2 = errorbar(qvals, norm2, std2); 
e2.Color = 'r'
hold on;
figure(7); e2 = errorbar(qvals, Nnorm2, Nstd2); 
e2.Color = 'k'

legend("Cancer", "Non-Cancer");
figure(7); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars with 75 points")

figure(8); e2 = errorbar(qvals, norm3, std3); 
e2.Color = 'r'
hold on;
figure(8); e2 = errorbar(qvals, Nnorm3, Nstd3); 
e2.Color = 'k'

legend("Cancer", "Non-Cancer");
figure(8); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars with 100 points")
%%
[sr,p] = ttest(stand, Nstd)
[sr1,p] = ttest(std1, Nstd1)
[sr2,p] = ttest(std2, Nstd2)
[sr3,p] = ttest(std3, Nstd3)