% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

function anomalyfunc(normallspec, s, mask,r,c) 
    %% load Data: NT_187
    load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
    load("NT_187.mat")
    % Mask NT_187
    redChannel = NT187(:, :, 1);
    greenChannel = NT187(:, :, 2);
    blueChannel = NT187(:, :, 3);
    
    mask4 = blueChannel > 240, redChannel > 240, greenChannel > 240;
    figure(2); imshow(mask4)
    hold on;

    [x y] = find(mask4 == 0)
    l = length(x)
    healallspec = []
    for k1 = 1:l
        healallspec = cat(2, healallspec, squeeze(reconstructedData(x(k1),y(k1),:)))
    end
    magnitude = sqrt(sum(healallspec.^2))
    normhealth = healallspec./magnitude
    %% load Data: NT_158
    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')
    load("NT_158.mat")
        % Mask NT_158
    redChannel = NT158(:, :, 1);
    greenChannel = NT158(:, :, 2);
    blueChannel = NT158(:, :, 3);
    
    mask5 = blueChannel > 245, redChannel > 250, greenChannel > 248;
    figure(2); imshow(mask5)
    hold on;
    [x1 y1] = find(mask5 ==0)
    l = length(x1)
    healallspec1 = []
    for k1 = 1:l
        healallspec1 = cat(2, healallspec1, squeeze(reconstructedData(x1(k1),y1(k1),:)))
    end
    magnitude1 = sqrt(sum(healallspec1.^2))
    normhealth1 = healallspec1./magnitude1
    %%  
    load('ReconResults_Brain_176_NT_20s_300iter_M3_Try1.mat')
    load("NT_176.mat")
    % Mask NT_176
    redChannel = NT176(:, :, 1);
    greenChannel = NT176(:, :, 2);
    blueChannel = NT176(:, :, 3);
    
    mask6 = blueChannel > 245, redChannel > 250, greenChannel > 248;
    figure(2); imshow(mask6)
    hold on;
    [x2 y2] = find(mask6 ==0)
    l = length(x2)
    healallspec2 = []
    for k1 = 1:l
        healallspec2 = cat(2, healallspec2, squeeze(reconstructedData(x2(k1),y2(k1),:)))
    end
    magnitude2 = sqrt(sum(healallspec2.^2))
    normhealth2 = healallspec2./magnitude2

    %% large NT dataset
    %normallhealth = [normhealth,normhealth1, normhealth2]
    %l=length(normallhealth)
    L = length(normallspec)
    l = length(normhealth)
    l1 = length(normhealth1)
    l2 = length(normhealth2)
    %% Similarity
    %%NT_187
    similarity = []
    for i = 1:L
        for j = 1:l
            dotprod = sum(normhealth(:,j).*normallspec(:,i))
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
    %%NT_158
    similarity1 = []
    for i = 1:L
        for j = 1:l1
            dotprod = sum(normhealth1(:,j).*normallspec(:,i))
        end 
        if dotprod > s
           similarity1 = [similarity1, 1]
        else
           similarity1 = [similarity1, 2]
        end 
    end 
    
    similar2 = transpose(similarity1)
    figure(4); imshow(mask)
    hold on;
    for f1 = 1:L
        if similar1(f1) == 1
            figure(4); plot(c(f1),r(f1),'bx')
            hold on;
        elseif similar1(f1) == 2
            figure(4); plot(c(f1),r(f1),'rx')
            hold on;
        end
    end

    %%NT_176
    similarity2 = []
    for i = 1:L
        for j = 1:l2
            dotprod = sum(normhealth2(:,j).*normallspec(:,i))
        end 
        if dotprod > s
           similarity2 = [similarity2, 1]
        else
           similarity2 = [similarity2, 2]
        end 
    end 
    
    similar3 = transpose(similarity2)
    figure(5); imshow(mask)
    hold on;
    for f1 = 1:L
        if similar1(f1) == 1
            figure(5); plot(c(f1),r(f1),'bx')
            hold on;
        elseif similar1(f1) == 2
            figure(5); plot(c(f1),r(f1),'rx')
            hold on;
        end
    end
end 
