% Author:      Arnav Nanda 
%              Apex High School 
%              Fitzpatrick Center for Interdisciplinary Engineering, Medicine and Applied Science (CIEMAS), Duke University

%% load data
load('ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat')

%% Training data section 
M1 = []
for t=1:25
    %  validation for cancer / training cancer data
    canarray = [squeeze(reconstructedData(79,23,:)), squeeze(reconstructedData(73,23,:))  ,squeeze(reconstructedData(82,23,:))  ,squeeze(reconstructedData(81,21,:))  ,squeeze(reconstructedData(79,19,:))  ,squeeze(reconstructedData(79,22,:))  ,squeeze(reconstructedData(78,25,:))  ,squeeze(reconstructedData(78,27,:))  ,squeeze(reconstructedData(76,26,:))  ,squeeze(reconstructedData(77,26,:))  ,squeeze(reconstructedData(80,33,:))  ,squeeze(reconstructedData(75,32,:))  ,squeeze(reconstructedData(74,30,:))  ,squeeze(reconstructedData(83,34,:))  ,squeeze(reconstructedData(79,33,:))  ,squeeze(reconstructedData(77,33,:))  ,squeeze(reconstructedData(74,31,:))  ,squeeze(reconstructedData(77,30,:))  ,squeeze(reconstructedData(79,28,:))  ,squeeze(reconstructedData(86,27,:))  ,squeeze(reconstructedData(87,35,:))  ,squeeze(reconstructedData(84,37,:))  ,squeeze(reconstructedData(81,36,:))  ,squeeze(reconstructedData(80,34,:))  ,squeeze(reconstructedData(86,31,:))]
    valcancer = [canarray(:,t)]
    canarray(:,t) = []
    cancer = sum(canarray, 2)
    
    %  validation for whitem / training whitem data
    
    wmarray = [squeeze(reconstructedData(88,11,:)),squeeze(reconstructedData(86,10,:)) ,squeeze(reconstructedData(87,13,:)) ,squeeze(reconstructedData(85,12,:)) ,squeeze(reconstructedData(57,25,:)) ,squeeze(reconstructedData(59,31,:)) ,squeeze(reconstructedData(63,37,:)) ,squeeze(reconstructedData(86,61,:)) ,squeeze(reconstructedData(85,59,:)) ,squeeze(reconstructedData(81,56,:)) ,squeeze(reconstructedData(77,15,:)) ,squeeze(reconstructedData(76,13,:)) ,squeeze(reconstructedData(75,14,:)) ,squeeze(reconstructedData(74,49,:)) ,squeeze(reconstructedData(75,56,:)) ,squeeze(reconstructedData(79,51,:)) ,squeeze(reconstructedData(80,50,:)) ,squeeze(reconstructedData(64,40,:)) ,squeeze(reconstructedData(64,30,:)) ,squeeze(reconstructedData(64,24,:)) ,squeeze(reconstructedData(86,56,:)) ,squeeze(reconstructedData(64,33,:)) ,squeeze(reconstructedData(64,34,:)) ,squeeze(reconstructedData(65,25,:)) ,squeeze(reconstructedData(65,34,:))]
    valwhitem = [wmarray(:,t)]
    wmarray(:,t) = []
    whitem = sum(wmarray, 2)
    
    gmarray = [squeeze(reconstructedData(79,10,:)),squeeze(reconstructedData(76,9,:)) ,squeeze(reconstructedData(78,10,:)) ,squeeze(reconstructedData(77,9,:)) ,squeeze(reconstructedData(60,13,:)) ,squeeze(reconstructedData(60,14,:)) ,squeeze(reconstructedData(59,15,:)) ,squeeze(reconstructedData(58,19,:)) ,squeeze(reconstructedData(56,21,:)) ,squeeze(reconstructedData(63,44,:)) ,squeeze(reconstructedData(56,62,:)) ,squeeze(reconstructedData(74,65,:)) ,squeeze(reconstructedData(77,64,:)) ,squeeze(reconstructedData(83,63,:)) ,squeeze(reconstructedData(55,61,:)) ,squeeze(reconstructedData(57,61,:)) ,squeeze(reconstructedData(55,64,:)) ,squeeze(reconstructedData(55,65,:)) ,squeeze(reconstructedData(62,63,:)) ,squeeze(reconstructedData(65,64,:)) ,squeeze(reconstructedData(68,63,:)) ,squeeze(reconstructedData(70,65,:)) ,squeeze(reconstructedData(58,15,:)) ,squeeze(reconstructedData(59,18,:)) ,squeeze(reconstructedData(55,23,:))]
    valgraym = [gmarray(:,t)]
    gmarray(:,t) = []
    graym = sum(gmarray, 2)
    
    %average spectra
    avgcancer = cancer/24
    avgwhitem = whitem/24
    avggraym = graym/24
    
    %magnitude
    magavgcancer = sqrt(sum(avgcancer.^2))
    magavgwhitem = sqrt(sum(avgwhitem.^2))
    magavggraym = sqrt(sum(avggraym.^2))
    
    %normalized
    normavgcancer = avgcancer./magavgcancer
    normavgwhitem = avgwhitem./magavgwhitem
    normavggraym = avggraym./magavggraym
    
    % given weights
    weightc = 0.9.*normavgcancer
    weightwm = 0.*normavgwhitem
    weightgm = 0.1*normavggraym
    weightedspec = [weightc,weightwm,weightgm]
    
    gauscan = max(weightc+normrnd(0, 0.1),0)
    gauswm = max(weightwm+normrnd(0, 0.1),0)
    gausgm = max(weightgm+normrnd(0, 0.1),0)

    
    i1 = []
    j1 = []
    k1 = []
    for i = 0:0.1:1
        for j = 0:0.1:1
            for k = 0:0.1:1
                if i+j+k == 1 
                    i1 = [i1, i]
                    j1 = [j1, j]
                    k1 = [k1, k]
                end 
            end 
        end 
    end
    ijk = [i1,
        j1,
        k1]
    r = []
    l = length(i1)
    for m = 1:l
        r  = cat(2, r, (i1(:,m).*normavgcancer+j1(:,m).*normavgwhitem+k1(:,m).*normavggraym))
    end 
    
    
    t=[]
    t1=[]
    t2=[]
    for a = 1:l
        t = [t,sum(gauscan.*r(:,a))]
    end
    for b = 1:l
        t1 = [t1,sum(gauswm.*r(:,b))]
    end 
    for c = 1:l
        t2 = [t2, sum(gausgm.*r(:,c))]
    end 
       
    u = [t,
        t1,
        t2]
    M = max(u,[],2)
    M1 = [M1, M]

end 

%% Validation 
    canarray = [squeeze(reconstructedData(79,23,:)), squeeze(reconstructedData(73,23,:))  ,squeeze(reconstructedData(82,23,:))  ,squeeze(reconstructedData(81,21,:))  ,squeeze(reconstructedData(79,19,:))  ,squeeze(reconstructedData(79,22,:))  ,squeeze(reconstructedData(78,25,:))  ,squeeze(reconstructedData(78,27,:))  ,squeeze(reconstructedData(76,26,:))  ,squeeze(reconstructedData(77,26,:))  ,squeeze(reconstructedData(80,33,:))  ,squeeze(reconstructedData(75,32,:))  ,squeeze(reconstructedData(74,30,:))  ,squeeze(reconstructedData(83,34,:))  ,squeeze(reconstructedData(79,33,:))  ,squeeze(reconstructedData(77,33,:))  ,squeeze(reconstructedData(74,31,:))  ,squeeze(reconstructedData(77,30,:))  ,squeeze(reconstructedData(79,28,:))  ,squeeze(reconstructedData(86,27,:))  ,squeeze(reconstructedData(87,35,:))  ,squeeze(reconstructedData(84,37,:))  ,squeeze(reconstructedData(81,36,:))  ,squeeze(reconstructedData(80,34,:))  ,squeeze(reconstructedData(86,31,:))]
    wmarray = [squeeze(reconstructedData(88,11,:)),squeeze(reconstructedData(86,10,:)) ,squeeze(reconstructedData(87,13,:)) ,squeeze(reconstructedData(85,12,:)) ,squeeze(reconstructedData(57,25,:)) ,squeeze(reconstructedData(59,31,:)) ,squeeze(reconstructedData(63,37,:)) ,squeeze(reconstructedData(86,61,:)) ,squeeze(reconstructedData(85,59,:)) ,squeeze(reconstructedData(81,56,:)) ,squeeze(reconstructedData(77,15,:)) ,squeeze(reconstructedData(76,13,:)) ,squeeze(reconstructedData(75,14,:)) ,squeeze(reconstructedData(74,49,:)) ,squeeze(reconstructedData(75,56,:)) ,squeeze(reconstructedData(79,51,:)) ,squeeze(reconstructedData(80,50,:)) ,squeeze(reconstructedData(64,40,:)) ,squeeze(reconstructedData(64,30,:)) ,squeeze(reconstructedData(64,24,:)) ,squeeze(reconstructedData(86,56,:)) ,squeeze(reconstructedData(64,33,:)) ,squeeze(reconstructedData(64,34,:)) ,squeeze(reconstructedData(65,25,:)) ,squeeze(reconstructedData(65,34,:))]
    gmarray = [squeeze(reconstructedData(79,10,:)),squeeze(reconstructedData(76,9,:)) ,squeeze(reconstructedData(78,10,:)) ,squeeze(reconstructedData(77,9,:)) ,squeeze(reconstructedData(60,13,:)) ,squeeze(reconstructedData(60,14,:)) ,squeeze(reconstructedData(59,15,:)) ,squeeze(reconstructedData(58,19,:)) ,squeeze(reconstructedData(56,21,:)) ,squeeze(reconstructedData(63,44,:)) ,squeeze(reconstructedData(56,62,:)) ,squeeze(reconstructedData(74,65,:)) ,squeeze(reconstructedData(77,64,:)) ,squeeze(reconstructedData(83,63,:)) ,squeeze(reconstructedData(55,61,:)) ,squeeze(reconstructedData(57,61,:)) ,squeeze(reconstructedData(55,64,:)) ,squeeze(reconstructedData(55,65,:)) ,squeeze(reconstructedData(62,63,:)) ,squeeze(reconstructedData(65,64,:)) ,squeeze(reconstructedData(68,63,:)) ,squeeze(reconstructedData(70,65,:)) ,squeeze(reconstructedData(58,15,:)) ,squeeze(reconstructedData(59,18,:)) ,squeeze(reconstructedData(55,23,:))]
    cancer = sum(canarray, 2)
    whitem = sum(wmarray, 2)
    graym = sum(gmarray, 2)
%average spectra
    avgcancer = cancer/25
    avgwhitem = whitem/25
    avggraym = graym/25
    
    %magnitude
    magavgcancer = sqrt(sum(avgcancer.^2))
    magavgwhitem = sqrt(sum(avgwhitem.^2))
    magavggraym = sqrt(sum(avggraym.^2))
    
    %normalized
    normavgcancer = avgcancer./magavgcancer
    normavgwhitem = avgwhitem./magavgwhitem
    normavggraym = avggraym./magavggraym

w = []
for s = 1:25
    w = cat(2,w,M1(1,s).*normavgcancer+M1(2,s).*normavgwhitem+M1(3,s).*normavggraym)
end 

  z=[]
    z1=[]
    z2=[]
    for a = 1:25
        z = [z,sum(gauscan.*w(:,a))]
    end
    for b = 1:25
        z1 = [z1,sum(gauswm.*w(:,b))]
    end 
    for c = 1:25
        z2 = [z2, sum(gausgm.*w(:,c))]
    end 

    y = [z,
        z1,
        z2]
    [f,I] = max(y,[],2)
    M1(:,I(1,:))


