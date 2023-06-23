
load('e2.mat')
imshow(e2)

redChannel = e2(:, :, 1);
greenChannel = e2(:, :, 2);
blueChannel = e2(:, :, 3);

% Create a mask of the background only.
mask = blueChannel > 200, redChannel > 200, greenChannel > 200;
figure(3); imshow(mask)
hold on;

[r c] = find(mask == 0)
L = length(r)
allspec = []
for k1 = 1:L
    allspec = cat(2, allspec, e2(r(k1),c(k1),:))
end
% creating array
allspec1 = [allspec(:,:,1)
            allspec(:,:,2)
            allspec(:,:,3)]
allspec2 = im2double(allspec1)

magnitude = sqrt(sum(allspec2.^2))
normallspec = allspec2./magnitude


% 25 choosen points 
ctable = e2(84,31,:)
canarray = [ctable(:,:,1), 
          ctable(:,:,2), 
          ctable(:,:,3)]
cancer = im2double(canarray)

wmtable = e2(66,33,:)
wmarray = [wmtable(:,:,1), 
          wmtable(:,:,2), 
          wmtable(:,:,3)]
whitem = im2double(wmarray)

gmtable = e2(64,10,:)
gmarray = [gmtable(:,:,1), 
         gmtable(:,:,2), 
         gmtable(:,:,3)]
graym = im2double(gmarray)

%magnitude
magavgcancer = sqrt(sum(cancer.^2))
magavgwhitem = sqrt(sum(whitem.^2))
magavggraym = sqrt(sum(graym.^2))

%normalized
normavgcancer = cancer./magavgcancer
normavgwhitem = whitem./magavgwhitem
normavggraym = graym./magavggraym


similar = []
for d1 = 1:L
    simcancer = sum(normallspec(:,d1).*normavgcancer)
    simwhite = sum(normallspec(:,d1).*normavgwhitem)
    simgray = sum(normallspec(:,d1).*normavggraym)

    if simcancer > simwhite
        if simcancer > simgray
           similar = [similar, 1]
        else
        end
    end
    
    if simwhite > simcancer
        if simwhite > simgray
           similar = [similar, 2]
        else
        end
    end
    
    if simgray > simcancer
        if simgray > simwhite
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
        figure(3); plot(c(f1),r(f1),'gx')
        hold on;
    elseif similar1(f1) == 3
        figure(3); plot(c(f1),r(f1),'yx')
        hold on;
    end
end



% Running UMAP
canarray1 = []
gmarray1 = []
wmarray1 = []
L = length(r)

for s1 = 1:L
    if similar1(s1) == 1
        canarray1 = cat(2, canarray1, e2(r(s1),c(s1),:))
    end 
    if similar1(s1) == 2
        gmarray1 = cat(2, gmarray1, e2(r(s1),c(s1), :))
    end
    if similar1(s1) == 3
        wmarray1 = cat(2, wmarray1, e2(r(s1),c(s1), :))
    end 
end
%cancer
cancer = [canarray1(:,:,1),
    canarray1(:,:,2),
     canarray1(:,:,3)]
cancer1 = im2double(cancer)
cancer2 = transpose(cancer1)
[creduction] = run_umap(cancer2)

%gray
graym = [gmarray1(:,:,1),
    gmarray1(:,:,2),
    gmarray1(:,:,3)]
graym1 = im2double(graym)
graym2 = transpose(graym1)
[gmreduction] = run_umap(graym2)

%white
whitem = [wmarray1(:,:,1),
    wmarray1(:,:,2),
    wmarray1(:,:,3)]
whitem1 = im2double(whitem)
whitem2 = transpose(whitem1)
[wmreduction] = run_umap(whitem2)



 figure(4); plot(creduction(:,1), creduction(:,2),'r.')
 hold on;
 figure(4); plot(gmreduction(:,1), gmreduction(:,2),'g.')
 hold on;
 figure(4); plot(wmreduction(:,1), wmreduction(:,2),'y.') 




