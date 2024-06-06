% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

function z = thyroid_abundvals(t)

    i1 = []
    j1 = []
    for i = 0:0.1:1
        for j = 0:0.1:1
            if i+j == 1
                i1 = [i1, i]
                j1 = [j1, j]
            end 
        end 
    end
    ij = [i1,
           j1]
    [M, g] = max(t, [], 1)
    u = rem(g,217)
    y = find(u == 0)
    u(y) = u(y)+217
    
    l = length(u)
    z= []
    for d = 1:2256
        z = [z,ij(:,u(d))]
    end
end 




