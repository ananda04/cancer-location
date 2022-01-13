classdef UMAP_extra_results < handle
    properties
        supervisorMatchedLabels; 
                        %vector of labels matched from a supervised template 
                        % 1 value per row of input data matrix
        fig;            % main figure with UMAP output plot
        qft;            % instance of qf tree object for input data. 
                        % The fig property contains the view of the tree.
        qftSupervisors; % instance of qf tree object for supervisor data. 
                        % The fig property
                        % contains the view of the tree.
        qfd;            % instance of qf dissimilarity object. Relevant 
                        % figure properties are fig for table view, 
                        % qHistFig for view of dismilarity score histogram
                        % fHistFig for view of f measure histogram
        
    end
end