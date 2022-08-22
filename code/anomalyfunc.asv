function anomalyfunc(normallspec, s) 
    %% load Data: NT_187
    load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
    load("NT_187.mat")
    mask = imbinarize(NT187)


    [x y] = find(mask =~ 1)
    l = length(x)
    healallspec = []
    for k1 = 1:l
        healallspec = cat(2, healallspec, squeeze(reconstructedData3(x(k1),y(k1),:)))
    end
    magnitude = sqrt(sum(healallspec.^2))
    normhealth = healallspec./magnitude
    %% load Data: NT_158
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')
    load("NT_158.mat")
    mask1 = imbinarize(NT158)

    [x1 y1] = find(mask1 =~ 1)
    l = length(x1)
    healallspec1 = []
    for k1 = 1:l
        healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(x(k1),y(k1),:)))
    end
    magnitude1 = sqrt(sum(healallspec1.^2))
    normhealth1 = healallspec1./magnitude1
    %%  
    load('ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat')
    load("NT_176.mat")
    mask1 = imbinarize(NT176)

    [x2 y2] = find(mask1 =~ 1)
    l = length(x2)
    healallspec2 = []
    for k1 = 1:l
        healallspec2 = cat(2, healallspec2, squeeze(reconstructedData(x(k1),y(k1),:)))
    end
    magnitude2 = sqrt(sum(healallspec2.^2))
    normhealth2 = healallspec2./magnitude2

    %% large NT dataset
    normallhealth = [normhealth,normhealth1, normhealth2]
    l=length(normallhealth)

        %% Similarity
    similarity = []
    for i = 1:l
        for j = 1:L
            dotprod = sum(normallhealth(:,i).*normallspec(:,j))
        end 
        if dotprod > s
           similarity = [similarity, 1]
        else
           similarity = [similarity, 2]
        end 
    end 
    
    similar1 = transpose(similarity)
    figure(3); imshow(mask)
    hold on;
    for f1 = 1:L
        if similar1(f1) == 1
            figure(3); plot(c(f1),r(f1),'bx')
            hold on;
        elseif similar1(f1) == 2
            figure(3); plot(c(f1),r(f1),'rx')
            hold on;
        end
    end
end 
