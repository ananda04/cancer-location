load('e2.mat')
imshow(e2)
import umapFileExchange.*;

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
[creduction] = run_umap(cancer1)

%gray
graym = [gmarray1(:,:,1),
    gmarray1(:,:,2),
    gmarray1(:,:,3)]
graym1 = im2double(graym)
[gmreduction] = run_umap(graym1)

%white
whitem = [wmarray1(:,:,1),
    wmarray1(:,:,2),
    wmarray1(:,:,3)]
whitem1 = im2double(whitem)
[wmreduction] = run_umap(whitem1)



     figure(4); plot(creduction(:,1), creduction(:,2),'r.')
     hold on;
     figure(4); plot(gmreduction(:,1), gmreduction(:,2),'g.')
     hold on;
     figure(4); plot(wmreduction(:,1), wmreduction(:,2),'y.') 