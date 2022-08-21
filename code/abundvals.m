function z = abundvals(t)

    i1 = []
    j1 = []
    k1 = []
    for i = 0:0.1:1
        for j = 0:0.1:1
            for k = 0:0.1:1
                if i+j+k == 1
                    i1 = [i1, i]
                    j1 = [j1, j]
                    k1 = [k1, k]
                end 
            end 
        end 
    end
    ijk = [i1,
           j1,
           k1]
    [M, g] = max(t, [], 1)
    u = rem(g,217)
    y = find(u == 0)
    u(y) = u(y)+217
    
    l = length(u)
    z= []
    for d = 1:2256
        z = [z,ijk(:,u(d))]
    end
end 




