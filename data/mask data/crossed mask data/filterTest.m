%% load Data: NT_187
load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
%% Mask NT_187
[mask x2 y2] = bmask("NT_187")

% Determine base set 
ncBset = ReconResults(reconstructedData, x2,y2)


%%
index = []
for z = 1:5
    index = cat(2, index, randi(numel(x2)))
end
index1 = x2(index)
index2 = y2(index)
trainingSet = ReconResults(reconstructedData, index1, index2)
label0_187 = transpose(zeros(5,1))
trainingSet = [trainingSet,
                label0_187]
%%
load("T1187_data.mat")
load("ReconResults_Brain_187_T1_20s_300iter_M3_Try1.mat")
cindex = []
for z = 1:5
    cindex = cat(2, cindex, randi(numel(r)))
end
index1 = r(cindex)
index2 = c(cindex)
trainingSet1 = ReconResults(reconstructedData, index1, index2)
label1_187 = transpose(ones(5,1))
trainingSet1 = [trainingSet1,
                label1_187]
%%
Train = [trainingSet, trainingSet1]

%%
can  = ReconResults(reconstructedData, r, c)
yfit = filtertest_SVM_model.predictFcn(can)
%%
[mask x2 y2] = bmask("T1_187")
tissue_assigner(yfit,mask,r,c,564,1,"SVM")
%%
load("ReconResults_Brain_187_T2_20s_300iter_M3_Try1.mat")
[mask x2 y2] = bmask("T2_187")
can  = ReconResults(reconstructedData, x2, y2)
t2187 = filtertest_SVM_model.predictFcn(can)
%%
tissue_assigner(t2187,mask,x2,y2,length(x2),1,"SVM")