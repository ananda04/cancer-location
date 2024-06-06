% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

function anomaly_detection(r,r1) 
    [r c] = find(mask == 0)
    L = length(r)
    allspec = []
    for k1 = 1:L
        allspec = cat(2, allspec, squeeze(reconstructedData(r(k1),c(k1),:)))
    end
    magnitude = sqrt(sum(allspec.^2))
    normallspec = allspec./magnitude
    %create a mask for non-cancerous slice 
    redChannel = u1(:, :, 1);
    greenChannel = u1(:, :, 2);
    blueChannel = u1(:, :, 3);
    
    mask1 = blueChannel > 250,  greenChannel > 250;
    figure(4); imshow(mask1)
    hold on;
    
    [x y] = find(mask1 == 0)
    l = length(x)
    healallspec = []
    for k1 = 1:l
        healallspec = cat(2, healallspec, squeeze(reconstructedData3(x(k1),y(k1),:)))
    end
    magnitude = sqrt(sum(healallspec.^2))
    normhealth = healallspec./magnitude
        %% Similarity
    similarity = []
    for i = 1:l
        for j = 1:L
            dotprod = sum(normhealth(:,i).*normallspec(:,j))
        end 
        if dotprod > 0.7
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
