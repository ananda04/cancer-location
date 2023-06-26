% This Script file is written to extract data from a figure.
% First open your figure and then run this m-fle. your extracted data is
% saved to Data.mat file. you can change this name from line 11 of this
% code.
% -------------------------------------
% Good luck
% Davood Shaghaghi
% K.N.Toosi University of Technology
% davood.shaghagh@gmail.com
%--------------------------------------

h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
xdata = cell2mat(get(dataObjs, 'XData'));  %data from low-level grahics objects
ydata = cell2mat(get(dataObjs, 'YData'));
zdata = cell2mat(get(dataObjs, 'ZData'));
Data = [xdata' ydata' zdata'];
