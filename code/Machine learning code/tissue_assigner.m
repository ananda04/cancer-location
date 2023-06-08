function tissue_assigner(type,mask,r,c,length,figure_number,subplot_position,String)

figure(figure_number); subplot(1,3,subplot_position)
figure(figure_number); imshow(mask)
hold all;
for f1 = 1:length
    if type(f1) == 0
        figure(figure_number); subplot(1,3,subplot_position)
        plot(c(f1),r(f1),'bx')
        hold on;
    elseif type(f1) == 1
        figure(figure_number);subplot(1,3,subplot_position)
        hold on;
        plot(c(f1),r(f1),'rx')
        hold on;
    end 
end 
subplot(1,3,subplot_position)
title(String)