function v = mse(x,y)
if any(size(x) ~= size(y))==1
    x=x';
end
v=mean((x-y).^2);