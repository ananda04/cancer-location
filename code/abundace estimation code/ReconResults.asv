% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

function normSpec = ReconResults(reconstructedData,r,c)
l = length(r)
spec = []
    for k = 1:l
        spec = cat(2, spec, squeeze(reconstructedData(r(k),c(k),:)))
    end
magnitude = sqrt(sum(spec.^2))
normSpec = spec./magnitude
end
