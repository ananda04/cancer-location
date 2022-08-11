function z = cdfcalc(x,disn,ind)

z=zeros(1,length(x));

for i=1:length(x)
    temp=max(disn(ind <= x(1,i)));
    if(isempty(temp))
        temp=0;
    end
    z(1,i)=temp;
end

end