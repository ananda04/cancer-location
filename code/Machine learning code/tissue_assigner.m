% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

function tissue_assigner(type,mask,r,c,length,figure_number,subplot_position,String)

%figure(figure_number); subplot(1,4,subplot_position)
figure(figure_number); imshow(mask)
hold on;
for f1 = 1:length
    if type(f1) == 0
        hold on;
        figure(figure_number);plot(c(f1),r(f1),'bx')

        %subplot(1,4,subplot_position)
        hold on;
    elseif type(f1) == 1
        hold on;
        figure(figure_number);plot(c(f1),r(f1),'rx')

        %subplot(1,4,subplot_position)
        %hold on;
        hold on;
    end
end 
%subplot(1,4,subplot_position)
%hold on;
%title(String)
