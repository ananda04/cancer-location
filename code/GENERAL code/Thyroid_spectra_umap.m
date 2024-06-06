% Author:      Arnav Nanda 
%              Apex High School 
%              Fitzpatrick Center for Interdisciplinary Engineering, Medicine and Applied Science (CIEMAS), Duke University

load('ReconResults_Thyroid_300iter_M3_Try1.mat')
load('Thyroid_Cancer_regions.mat')
load("Thyroid_noncancer_regions.mat")

[r c] = find(Cancer == 1)
L = length(r)
cancerspec = []
for k1 = 1:L
    cancerspec = cat(2, cancerspec, squeeze(reconstructedData(r(k1),c(k1),:)))
end

%gt known noncancerous data
[rn cn] = find(noncancer == 1)
l = length(rn)
healthspec = []
for s1 = 1:l
    healthspec = cat(2, healthspec, squeeze(reconstructedData(rn(s1),cn(s1),:)))
end 
%%
avgcancerspec = mean(cancerspec,2)
avghealthspec = mean(healthspec,2)

errcancer = std(cancerspec,[],2)
errnoncancer = std(healthspec, [],2)


% Plot graphs with error bars
figure(5); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
hold on;
e1 = errorbar(qvals, avgcancerspec,errcancer); 
e1.Color = 'r'
hold on;
figure(5); e2 = errorbar(qvals, avghealthspec,errnoncancer); 
e2.Color = 'b'
hold on;
legend([e1(1),e2(1)],'Cancer', 'No cancer')

cancerspec1 = transpose(cancerspec)
healthspec1 = transpose(healthspec)
% UMAP reduction 
[reduction] = run_umap(cancerspec1)
[reduction1] = run_umap(healthspec1)
figure(4); plot(reduction(:,1), reduction(:,2),'r.')
hold on;
figure(4); plot(reduction1(:,1), reduction1(:,2),'b.')
hold on;
legend('Cancer','No cancer')



