function [reduction, umap, clusterIdentifiers, extras]=run_umap(varargin)
%%RUN_UMAP reduces data matrices with 3+ parameters down to fewer
%   parameters using the algorithm UMAP (Uniform Manifold Approximation and
%   Projection).
%
%   [reduction,umap,clusterIdentifiers,extras]=RUN_UMAP(csv_file_or_data,...
%   'NAME1',VALUE1, 'NAMEN',VALUEN) 
%   
%
%   OUTPUT ARGUMENTS
%   Invoking run_umap returns these values:
%   1)  reduction, the actual data that UMAP reduces from the data 
%       specified by the input argument csv_file_or_data; 
%   2)  umap, an instance of the UMAP class made ready for the invoker 
%       to save in a MATLAB file for further use as a template.
%   3)  clusterIdentifiers, identifiers of clusters found by dbscan 
%       or DBM methods when run on the reduction of umap.
%   4)  extras, an instance of the class UMAP_extra_results.
%       See properties comments in UMAP_extra_results.m.
%
%
%   REQUIRED INPUT ARGUMENT
%   The argument csv_file_or_data is either 
%   A) a char array identifying a CSV text file containing the data 
%      to be reduced. 
%   B) the actual data to be reduced; a numeric matrix.
%
%   If A) then the CSV file must have data column names in the first line.
%   These annotate the parameters which UMAP reduces.  If B) then
%   parameter names are needed by the input argument 'parameter_names'
%   when creating or running a template.
%
%   Invoke run_umap with no arguments to download CSV files that 
%   our examples below rely upon.
%
%
%   OPTIONAL INPUT ARGUMENTS
%   Some of these are identical to those in the original Python
%   implementation documented by the inventors in their document "Basic
%   UMAP parameters" which can be retrieved at
%   https://umap-learn.readthedocs.io/en/latest/parameters.html.
%   The optional argument name/value pairs are:
%
%   Name                    Value
%
%   'n_neighbors'           Controls local and global structure as does the
%                           same input argument for the original
%                           implementation.
%                           Default is 30. 
%   
%   'min_dist'              Controls how tightly UMAP is allowed to pack
%                           points together as does the same input argument
%                           for the original implementation. Modifying this
%                           value requires the Curve Fitting Toolbox.
%                           Default is 0.3.
%
%   'metric'                Controls how distance is computed in the
%                           ambient space as does the same input argument
%                           for the original implementation. Accepted
%                           values for metric include 'euclidean',
%                           'cosine', 'cityblock', 'seuclidean',
%                           'correlation', 'jaccard', 'spearman',
%                           'hamming'. These metrics are described in
%                           MATLAB's documentation for knnsearch.
%                           Default is 'euclidean'.
%
%   'randomize'             true/false.  If false run_umap invokes
%                           MATLAB's "rng default" command to ensure the
%                           same random sequence of numbers between
%                           invocations.
%                           Default is false.
%
%   'template_file'         This identifies a .mat file with a saved
%                           instance of the UMAP class that run_umap
%                           previously produced. The instance must be be a
%                           suitable "training set" for the current "test
%                           set" of data supplied by the argument
%                           csv_file_or_data. Template processing
%                           accelerates the UMAP reduction and augments
%                           reproducibility. run_umap prechecks the
%                           suitability of the template's training set for
%                           the test set by checking the name and standard
%                           deviation distance from the mean for each
%                           parameter (AKA data column).
%                           Default is empty ([]...no template).
%
%   'see_training'          true/false to see/hide plots of both the
%                           supervising data and the supervised data with
%                           label coloring and legend. This takes effect
%                           when applying a UMAP template of a supervised
%                           reduction and when the input argument
%                           verbose='graphic'. Examples 5, 10, 11, 12 and
%                           16 apply a supervised template.  Example 16
%                           illustrates this.
%                           Default is false.
%
%   'parameter_names'       Cell of char arrays to annotate each column
%                           of the data matrix specified by
%                           csv_file_or_data. This is only needed if a
%                           template is being used or saved.
%                           Default is {}.
%                           
%   'verbose'               Accepted values are 'graphic', 'text', or
%                           'none'. If verbose='graphic' then the data
%                           displays with probability coloring and contours
%                           as is conventional in flow cytometry analysis.
%                           If method='Java' or method='MEX', then the
%                           display refreshes as optimize_layout progresses
%                           and a progress bar is shown along with a handy
%                           cancel button. If verbose='text', the progress
%                           is displayed in the MATLAB console as textual
%                           statements.
%                           Default is 'graphic'.
%                           
%   'method'                Selects 1 of 7 implementations for UMAP's
%                           optimize_layout processing phase which does
%                           stochastic gradient descent.  Accepted values
%                           are 'MEX', 'C++', 'Java', 'C', 'C vectorized',
%                           'MATLAB' or 'MATLAB Vectorized'. 'MEX' is our
%                           fastest & most recent implementation. The
%                           source
%                           umap/sgdCpp_files/mexStochasticGradientDescent.cpp
%							provides an illustration of the simplicity and
%							power of MATLAB's C++ MEX programming
%							framework.
%                           The other methods are provided for educational
%                           value.  They represent our iterative history of
%                           speeding up the slowest area of our translation
%                           from Python. We found stochastic gradient
%                           descent to be the least vectorizable. 'C' and
%                           'C vectorized', produced by MATLAB's "C coder"
%                           app, were our first attempts to accelerate our
%                           MATLAB programming. We were surprised to find
%                           our next attempt with 'Java' was faster than
%                           the code produced by C Coder.  Thus we
%                           proceeded to speed up with the 'C++' and then
%                           'MEX' implementations. 'C++', our 2nd fastest,
%                           is a separate spawned executable.  The build
%                           script and cpp source file are found in
%                           umap/sgdCpp_files.
%							Note that MathWorks open source license
%                           prevented the 'C' and 'C vectorized' modules to
%							be distributed.  You can download them too from
%                           http://cgworkspace.cytogenie.org/GetDown2/demo/umapDistribution.zip
%                           MEX, Java and C++ support the progress plots
%                           and cancellation options given by argument
%                           verbose='graphic'.
%                           Default is 'MEX'.
%
%  'progress_callback'      A MATLAB function handle that run_umap
%                           invokes when method is 'Java', 'C++', or 'MEX'
%                           and verbose='graphic'. The input/output
%                           expected of this function is
%                           keepComputing=progress_report(objectOrString).
%                           The function returns true/false to tell the
%                           reduction to keep computing or stop computing.
%                           The objectOrString argument is either a 
%                           status description before stochastic
%                           gradient descent starts or an object
%                           with properties (getEmbedding, getEpochsDone 
%                           and getEpochsToDo) which convey the state of
%                           progress. The function function 
%                           progress_report here in run_umap.m exemplifies
%                           how to write a callback. .
%                           Default is the function progress_report
%                           in run_umap.m.
%
%   'ask_to_save_template'  true/false instructs run_umap to ask/not ask
%                           to save a template PROVIDING method='Java',
%                           verbose='graphic', and template_file is empty.
%                           Default is false.
%
%   'label_column'          number identifying the column in the input data
%                           matrix which contains numeric identifiers to
%                           label the data for UMAP supervision mode.                    
%   `                       Default is 0, which indicates no label column.
%
%   'label_file'            the name of a properties file that contains
%                           the label names.  The property name/value
%                           format is identifier=false.
%                           Default is [].
%
%   'n_components'          The dimension of the space into which to embed
%                           the data.
%                           Default is 2.
%
%   'epsilon'               The epsilon input argument used by MATLAB's
%                           dbscan algorithm. 
%                           Default is 0.6.
%
%   'minpts'                The minpts input argument used by MATLAB's 
%                           dbscan.
%                           Default is 5.
%
%   'dbscan_distance'       The distance input argument used by MATLAB's
%                           dbscan.
%                           Default is 'euclidean'.
%
%
%   'cluster_output'        Allowed values: 'none', 'numeric', 'graphic'.  
%                           When the value~='none' && nargout>2 
%                           cluster results are returned in the 3rd output
%                           argument clusterIdentifiers.
%                           Default is 'numeric'.
%
%   'cluster_2D_method'     Clustering method when n_components==2.
%                           Allowed values are 'dbscan' or 'dbm'.
%                           Default is our own method 'dbm'.
%
%   'cluster_detail'        Used when (nargout>2 and the input 
%                           argument 'cluster_output'~='none') OR 
%                           'cluster_output'=='graphic'.  
%                           Allowed values are 'very low', 'low', 'medium',
%                           'high', 'very high', 'most high', 'adaptive' or
%                           'nearest neighbor' or 'dbscan arguments'
%                           if 'dbscan arguments' then run_umap uses the 
%                           input arguments 'epsilon' and 'minpts' to
%                           determine cluster detail IF the dbscan method 
%                           is needed.  If needed and 'cluster_detail'
%                           value is 'adaptive' or 'nearest neighbor' 
%                           then run_umap replaces with 'dbscan arguments'.
%                           Default is 'very high'.
%
%   'save_template_file'    Fully qualified path of the file to save the
%                           resulting UMAP object as a template.  One
%                           can also save run_umap's 2nd output argument.
%                           Default is [].
%
%   'match_supervisors'     A number indicating how to relabel data points 
%                           in the embedding data if the UMAP reduction 
%                           is guided by a template that in turn is guided 
%                           by supervisory labels.
%                           0 matches supervised and supervising data
%                             groupings by distance of medians. Supervising
%                             groupings are data points in the template's
%                             embedding that have the same supervisory
%                             label. Supervised groupings are DBM clusters
%                             in the final template-guided embedding. The
%                             publication that introduces DBM clustering is
%                             http://cgworkspace.cytogenie.org/GetDown2/demo/dbm.pdf.
%                           1 (default) matches groupings by quadratic 
%                             form dissimilarity.  The publication 
%                             that introduces QF dissimilarity is
%                             https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5818510/.
%                           2 matches supervised DBM clusters by assigning
%                             the label of the closest supervising data 
%                             point to each supervised data point and then 
%                             choosing the most frequent label in the 
%                             cluster.  Closeness is based on euclidean 
%                             distance in the supervised and supervising
%                             embedding data spaces.
%                           3 is similar to 2 except it only uses closeness
%                             to the supervising data point to relabel the
%                             supervised data points without the aid of DBM
%                             clustering.  Thus supervised groupings in the
%                             embedding space may have small fragments
%                             of data points that occur in distant 
%                             islands/clusters of the embedding.
%
%   'match_3D_limit'        The lower limit for the # of data rows before
%                           3D progress plotting avoids supervisor label
%                           matching.  This applies only when reducing with
%                           a supervised template and the n_components>2
%                           and verbose=graphic. If > limit then supervisor
%                           matching ONLY occurs in the final plot
%                           ...otherwise supervisor label matching occurs
%                           during progress plotting before epochs finish.
%                           Default is 20000. 
%                           
%   'qf_dissimilarity'      Show QF dissimilarity scores between data
%                           groupings in the supervised and supervising 
%                           embeddings. The showing uses a sortable data 
%                           table as well as a histogram.
%                           Default is false.
%                           run_umap only consults this argument when it
%                           guides a reduction with a supervised
%                           template.
%                           
%   'qf_tree'               Show a dendrogram plot that represents the
%                           relatedness of data groupings in the
%                           supervising and supervised embeddings. The
%                           above documentation for the match_supervisors
%                           argument defines "data groupings". The
%                           publication that introduces the QF tree is
%                           https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/.
%                           This uses phytree from MATLAB's Bioinformatics
%                           Toolbox, hence changing this value to true
%                           requires the Bioinformatics Toolbox.
%                           Default is false.
%                           run_umap only consults this argument when it
%                           guides a reduction with a supervised
%                           template.
%                           
%   'joined_transform'      true/false for a new transform method to avoid
%                           false positives when applying a template whose
%                           training data differs too much from test set
%                           data. This feature is not part of UMAP's
%                           original Python implementation. 
%                           Currently this not supported when
%                           method=C++.  Support for C++ is coming soon.
%                           Default is false.
%
%   'python'                true/false to use UMAP's original
%                           implementation written in Python instead of
%                           this one written in MATLAB, C++ and Java.  The
%                           Python implementation is from Leland McInnes,
%                           John Healy, and James Melville.
%                           If true then certain arguments are ignored:
%                           joined_transform, method, verbose, 
%                           and progress_callback.
%                           Default is false.
%
%   EXAMPLES 
%   Note these examples assume your current MATLAB folder is where
%   run_umap.m is stored.
%
%   1.  Download the example CSV files and run sample10k.csv.
%
%       run_umap
%
%   2.  Reduce parameters for sample30k.csv and save as UMAP template (ut).
%
%       run_umap('sample30k.csv', 'save_template_file', 'utBalbc2D.mat');
%
%   3.  Reduce parameters for sample130k.csv using prior template.
%
%       run_umap('sample130k.csv', 'template_file', 'utBalbc2D.mat');
%
%   4.  Reduce parameters for sampleBalbcLabeled55k.csv supervised by
%       labels produced by EPP and save as a UMAP supervised template 
%       (ust), EPP is a conservative clustering technique described at
%       https://www.nature.com/articles/s42003-019-0467-6. EPP stands for
%       "Exhaustive Projection Pursuit".  By clustering exhaustively in 2
%       dimension pairs, this technique steers more carefully away from the
%       curse of dimensionality than does UMAP or t-SNE.
%
%       To use EPP you can download AutoGate from cytogenie.org which
%       contains tutorials on using EPP.
%
%       run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'save_template_file', 'ustBalbc2D.mat');
%
%   5.  Reduce parameters for sampleRag148k.csv using template that is
%       supervised by EPP.  This takes the clusters created by EPP on the
%       lymphocytes of a normal mouse strain (BALB/c) and applies them via
%       a template to a mouse strain (RAG) that has neither T cells nor B
%       cells.
%
%       run_umap('sampleRag148k.csv', 'template_file', 'ustBalbc2D.mat');
%
%   6.  Reduce parameters for sample30k.csv and return & plot cluster 
%       identifiers using density-based merging described at 
%       http://cgworkspace.cytogenie.org/GetDown2/demo/dbm.pdf.
%
%       [~,~, clusterIds]=run_umap('sample30k.csv', 'cluster_output', 'graphic');
%
%   7.  Repeat sample 2 but for 3D output and return cluster identifiers
%       and save the result as 3D template.
%
%       [~, ~, clusterIds]=run_umap('sample30k.csv', 'n_components', 3, 'save_template_file', 'utBalbc3D.mat');
%
%   8.  Repeat example 3 in 3D.
%
%       run_umap('sample130k.csv', 'template_file', 'utBalbc3D.mat');
%
%   9.  Reduce parameters and save template for sampleRagLabeled60k.csv
%       using labels produced by an expert biologist drawing manual gate 
%       sequences on lymphocyte data taken from a RAG mouse strain which 
%       has no T cells or B cells.
%
%       run_umap('sampleRagLabeled60k.csv', 'label_column', 11, 'label_file', 'ragLabels.properties', 'save_template_file', 'ustRag2D.mat');
%
%   10. Reduce parameters for lymphocyte data taken from a BALB/c mouse
%       strain using template created in example 9.  This takes the
%       clusters created on the lymphocyte data of a knockout mouse strain
%       (RAG) with no B cells or T cells and applies them to a normal mouse
%       strain (BALB/c) which has both cell types.  This illustrates logic
%       to prevent false positives for data not seen when training/creating
%       supervised templates.  Choose to re-supervise to see effect.
%
%       run_umap('sample30k.csv', 'template_file', 'ustRag2D.mat');
%
%   11. Repeat example 10 but use joined_transform.  Currently 'method'==
%       'Java' is the only support for this.
%
%        run_umap('sample30k.csv', 'template_file', 'ustRag2D.mat', 'method', 'Java', 'joined_transform', true);
%
%   12. Run example 5 again showing training/test set plot pair, QF tree
%       and QF dissimilarity plots.
%
%       run_umap('sampleRag148k.csv', 'template_file', 'ustBalbc2D.mat', 'qf_tree', true, 'qf_dissimilarity', true, 'see_training', true);
%
%   13. Compare our implementation to the original Python implementation by
%       repeating example 2 as follows.
%       
%       run_umap('sample30k.csv');
%       run_umap('sample30k.csv', 'python', true);
%
%   14. Compare our implementation with C++ method to the original Python 
%       implementation by repeating example 4 as follows.
%
%       NOTE: Since February 2020, there is an incompatibility between the
%       Python packages "umap-learn" and "numba". You may find that the
%       numba package will throw an error while doing supervised UMAP. One
%       temporary workaround is to add a # symbol at the start of the line
%       "@numba.njit()" before the function "fast_intersection" (line 629
%       in version 0.4.3 of umap-learn) in the file (Python
%       folder)\Lib\site-packages\umap\umap_.py. This will disable numba
%       for one function; note that it will slow down the Python code,
%       making our time comparison less competitive.
%
%       run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'save_template_file', 'ustBalbc2D.mat');
%       run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'python', true, 'save_template_file', 'pyUstBalbc2D.mat');
%
%   15. Compare our implementation to the original Python 
%       implementation by repeating example 5 as follows.
%
%       run_umap('sampleRag148k.csv', 'template_file', 'ustBalbc2D.mat');
%       run_umap('sampleRag148k.csv', 'template_file', 'pyUstBalbc2D.mat');
%
%   16. Combining aspects of previous examples, this one creates
%       a UMAP supervised template for 3D output, then applies this
%       template to a different example. The final run_umap returns all
%       possible outputs, including the extras argument that contains
%       supervisor matching labels (1 per row of input data matrix) and
%       qf_tree and qf_dissimilarity arguments. The main plot shows
%       training/test set plot pair.
%
%       run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'qf_tree', true, 'n_components', 3, 'save_template_file', 'ustBalbc3D.mat');
%       [reduction, umap, clusterIdentifiers,extras]=run_umap('sample10k.csv', 'template_file', 'ustBalbc3D.mat', 'qf_tree', true, 'qf_dissimilarity', true, 'see_training', true, 'cluster_output', 'graphic');
%
%   NOTE that you can do supervised UMAP and templates with n_components
%       ...but we have not had time to update the 3D GUI to show where the 
%       supervised regions fall.
%
%
%   REQUIRED PATHS
%   This distribution has 2 folders:  umap and util.  
%   You must set paths to these folders plus the java inside of umap.jar.
%   Assume you have put these 2 folders under /Users/Stephen.
%   The commands that MATLAB requires would be:
%
%   addpath /Users/Stephen/umap
%   addpath /Users/Stephen/util
%   javaaddpath('/Users/Stephen/umap/umap.jar');
%
%
%   ALGORITHMS
%   UMAP is the invention of Leland McInnes, John Healy and James Melville
%   at Canada's Tutte Institute for Mathematics and Computing.  See
%   https://umap-learn.readthedocs.io/en/latest/.
%
%   AUTHORSHIP
%   Primary Developer+math lead: Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer:  Stephen Meehan <swmeehan@stanford.edu> 
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University. 
%   License: BSD 3 clause
%
%   IMPLEMENTATION NOTES
%   This is a total rewrite of the original Python implementation from
%   Leland McInnes, John Healy and James Melville. This implementation is
%   written in MATLAB, C++, and Java. The source is distributed openly on
%   MathWorks File Exchange. This implementation follows a very similar
%   structure to the Python implementation, and many of the function
%   descriptions are nearly identical. Leland McInnes has looked over it
%   and considered it "a fairly faithful direct translation of the original
%   Python code (except for the nearest neighbor search)". If you have
%   UMAP's Python implementation you can check how faithful and fast this
%   re-implementation is by using the argument python=true. When python is
%   false and method is MEX we observe superior performance on our Mac and
%   Windows laptops in most cases. For the cases of template-guided and
%   supervised parameter reduction the performance is significantly faster
%   than the Python implementation regardless of data size.
%
%   If you wish to have a simple user GUI to run these UMAP features, 
%   download AutoGate at CytoGenie.org.
%   
%
clusterIdentifiers=[];
extras=UMAP_extra_results;
initJava;
pth=fileparts(mfilename('fullpath'));
pPth=fileparts(pth);
utilPath=fullfile(pPth, 'util');
addpath(utilPath);
globals=BasicMap.Global;
try
    props=fullfile(File.Home, 'run_umap.mat');
    globals.load(props);
catch ex
end
reduction=[];
umap=[];
p=parseArguments();
parse(p,varargin{:});
args=p.Results;
if islogical(args.verbose)
    if args.verbose
        args.verbose = 'graphic';
    else
        args.verbose = 'none';
    end
end
plotting=strcmpi(args.verbose, 'graphic');
csv_file_or_data=args.csv_file_or_data;
save_template_file = args.save_template_file;
curAxes=[];
if plotting
    xLabel=[];
    yLabel=[];
    zLabel=[];
    fig=figure('name', 'Running UMAP ...');
    extras.fig=fig;
    if args.qf_tree 
        movegui(fig, 'south')
    elseif args.qf_dissimilarity
        movegui(fig, 'center')
    else
        movegui(fig, 'onscreen');
    end
    curAxes=gca;
    if isempty(csv_file_or_data)
        if askYesOrNo(Html.WrapHr(...
                ['Should run_umap.m download example csv files<br>',...
                'from the Herzenberg Lab @ Stanford University<br><br>', ...
                '.. and then run one of them?']))
            csv_file_or_data=downloadCsv;
        end
        if isempty(csv_file_or_data)
            if plotting
                delete(fig);
            end
            globals.save;
            return;
        end
        if ~askYesOrNo(Html.Wrap([...
                'Test csv files have been downloaded:<ol>'...
                ' <li>sample10k<li>sample30k<li>sampleBalbcLabeled55k'...
                '<li>sample130k<li>sampleRag148k.csv<li>sampleRag55k.csv'...
                '</ol><br><center>Run UMAP on <b>sample10k</b> now?<hr></center>']))
            if plotting
                delete(fig);
            end
            globals.save;
            return;
        end
    end
end
if nargout>=3 && args.n_components>2
    % check for presence of DBSCAN
    if ~Density.HasDbScan(false)
        if plotting
            if ~askYesOrNo(Html.WrapHr(['DBSCAN for clustering in 3+D is '...
                    '<br>not downloaded ...Continue?']))
                delete(fig);
                globals.save;
                return;
            end
        end
        dispNoDbScan;
    end
end
if ischar(csv_file_or_data)
    if ~exist(csv_file_or_data, 'file')
        showMsg(Html.WrapHr(['The text file "<b>' csv_file_or_data ...
            '</b>"<br><font color="red"><i>can not be found !!</i></font>']));
        if plotting
            delete(fig);
        end
        globals.save;
        return;
    end
    warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
    t=readtable(csv_file_or_data, 'ReadVariableNames', true);
    inData=table2array(t);
    parameter_names=File.CsvNames(t);
else
    inData=csv_file_or_data;
    parameter_names=args.parameter_names;
end
newSubsetIdxs=[];
template_file=args.template_file;
firstPlot=true;
if ~isempty(template_file)
    if ischar(template_file)
        if ~exist(template_file, 'file')
            showMsg(Html.WrapHr(['The template file "<b>' template_file ...
                '</b>"<br><font color="red"><i>can not be found !!</i></font>']));
            if plotting
                delete(fig);
            end
            globals.save;
            return;
        end
        if length(parameter_names)~=size(inData, 2)
            showMsg(Html.WrapHr(sprintf(['<b>Can not create '...
                'or use template</b> ...<br>'...
                '%d parameter_names... but data has %d parameters?'], ...
                length(parameter_names), size(inData,2))));
            if plotting
                delete(fig)
            end
            globals.save;
            return;
        end
        [umap, ~, canLoad, reOrgData, paramIdxs]=...
            Template.Get(inData, parameter_names, ...
            template_file, 3);
        if ~isempty(reOrgData)
            % column label order differed
            inData=reOrgData;
            if ~isempty(parameter_names) && ~isempty(paramIdxs)
                parameter_names=parameter_names(paramIdxs);
            end
        end
    elseif isa(template_file, 'UMAP')
        umap=template_file;
        canLoad = true;
    end
    if isempty(umap)
        if ~canLoad
            if plotting
                showMsg(Html.WrapHr(['No template data found in <br>"<b>', ...
                    template_file '</b>"']));
            else
                disp(['No template data found in ' template_file]);
            end
        end
        if plotting
            delete(fig);
        end
        globals.save;
        return;
    else
        args.n_components=umap.n_components;
        if ~isempty(umap.supervisors)
            %Connor's NEW joined_transform immunizes reduction from
            %false positives if items in the test set are TOO different
            %from the training set
            
            if ~args.joined_transform
                [percNewSubsets, unknownIdxs]=...
                    Template.CheckForUntrainedFalsePositives(umap, inData);
                if percNewSubsets>13 && plotting
                    [choice, cancelled]=Template.Ask(percNewSubsets);
                    if cancelled
                        if plotting
                            delete(fig);
                        end
                        globals.save;
                        return;
                    end
                    if choice==2
                        umap.clearLimits;
                        newSubsetIdxs=unknownIdxs;
                        template_file=[];
                    end
                end
            end
            umap.supervisors.initClustering(args.cluster_detail, ...
                args.cluster_method_2D, args.minpts, ...
                args.epsilon, args.dbscan_distance);
            umap.supervisors.initPlots(args.contour_percent);
        end
    end
else
    umap = UMAP;
    umap.dimNames=parameter_names;
end
[nRows, nCols]=size(inData);

isSupervising=isprop(umap, 'supervisors') && ~isempty(umap.supervisors);
if isSupervising
    supervisors=umap.supervisors;
    if args.n_components==2
        progressMatchType=0;
    else
        limit=args.match_3D_limit;
        if limit<nRows
            progressMatchType=-1;
        else
            progressMatchType=3;
        end
    end
    
end
umap.metric=args.metric;
umap.n_epochs=args.n_epochs;
umap.n_neighbors=args.n_neighbors;
umap.min_dist=args.min_dist;
umap.n_components=args.n_components;
if strcmpi('Java', args.method)
    if ~initJava
        args.method='MEX';
        showMsg(Html.WrapHr('Could not load umap.jar for Java method'), ...
            'Problem with JAVA...', 'south west', false, false);
    end
end
        
method=umap.setMethod(args.method);
umap.verbose=~strcmpi(args.verbose, 'none');
umap.random_state=~args.randomize;
%umap.negative_sample_rate=30;
tick=tic;

labelMap=[];
sCols=num2str(size(inData,2));
nParams=length(parameter_names);
good=nParams==0||nParams==nCols || (args.label_column>0 &&...
    (nParams==nCols-1 || nParams==nCols));
if ~good
    if args.label_column>0
        preAmble=sprintf(['# data columns=%d, # parameter_names=%d '...
            'since label_column=%d <br># parameter_names must be '...
            '%d or %d '],  nCols, nParams, args.label_column, ...
            nCols, nCols-1);
    else
        preAmble=sprintf(['# of data columns(%s) must equal '...
            '# of parameter_names(%d)'], sCols, nParams);
    end
    msg(Html.WrapHr(preAmble));
    assert(nParams==0||nParams==nCols || (args.label_column>0 &&...
        (nParams==nCols-1 || nParams==nCols)), preAmble);    
end
if ~isempty(newSubsetIdxs)
    hasLabels=true;
    labelCols=0;
    [labels, labelMap]=resupervise(umap, inData, newSubsetIdxs);
    nLabels=length(unique(labels));
elseif args.label_column>0
    if ~isempty(template_file)
        showMsg(Html.WrapHr(['Can not do supervised mode <br>'...
            'AND use prior template at<br>the same time!']));
        globals.save;
        return;
    end
    hasLabels=true;
    labelCols=1;
    good=args.label_column>0 && args.label_column<=nCols;
    if ~good
        msg(Html.WrapHr(['The input data has ' sCols ' columns ...<br>'...
            'THUS the label_column must be >=1 and <= ' sCols]));
        assert(args.label_column>0 && args.label_column<=nCols, [...
            'label_column must be >=1 and <= ' sCols]);
        
    end
    labels=inData(:, args.label_column);
    nLabels=length(unique(labels));
    if nLabels > .5*nRows
        preAmble=sprintf(...
            'WARNING:  %d is a LOT of distinct labels for a %dx%d matrix!',...
            nLabels, nRows, nCols);
        msg(preAmble);
        warning(preAmble);
    end
    inData(:, args.label_column)=[];
    if args.label_column<=nParams
        parameter_names(args.label_column)=[];
    end
    umap.dimNames=parameter_names;
    nLabels=length(unique(labels));
    if exist(args.label_file, 'file')
       map=java.util.Properties;
       try
           map.load(java.io.FileInputStream(args.label_file));
       catch ex
           showMsg(['Can not load ' args.label_file]);
           delete(fig);
           globals.save;
           return;
       end
       labelMap=map;
    elseif ~isempty(args.label_file)
        if ~askYesOrNo(['<html>Can not find label file<br><br><b>' ...
                BasicMap.Global.smallStart args.label_file ...
                BasicMap.Global.smallEnd '</b><br><br><center>'...
                'Continue WITHOUT supervision ??</center><hr></html>'], ...
                'No supervising labels...', 'center', false)
            msg(['Can not load ' args.label_file]);
            delete(fig);
            globals.save;
            return;
        end
   end
else
    hasLabels=false;
    labelCols=0;
    nLabels=0;
end
if any(isnan(inData(:)))
    if plotting
        if isequal('Yes', questdlg({...
                'Data matrix has NAN values which',...
                'which cause odd effects on UMAP!','', ...
                'Try to remove nan values?'}))
            allNanColumns=all(isnan(inData));
            if any(allNanColumns)
                inData=inData(:, ~allNanColumns);
            end
            allNanRows=all(isnan(inData'));
            if any(allNanRows)
                inData=inData(~allNanRows,:);
            end
            
        end
        if any(isnan(inData(:)))
            showMsg(Html.WrapHr(['Sorry...<br>can not proceed<br>'...
                '<br>NAN values exist... SIGH!']));
            globals.save;
            return;
        end
    else
        error('Can not proceed with NAN values');
    end
end

info=[String.encodeInteger(nRows) 'x' String.encodeInteger(nCols-labelCols)];
if ischar(csv_file_or_data)
    [~, fileName]=fileparts(csv_file_or_data);
    info=['UMAP on ' fileName ', ' info];
else
    info=['[UMAP on ' info];
end
if args.python
    info=[info ', Python'];
else
    info=[info ', ' method];
end
if plotting
    set(fig, 'NumberTitle', 'off', 'name', info );
    drawnow;
end
pause(.01);
info2=['(optimize\_layout method=' method ')'];
if strcmpi(method, 'C++')
    if ~StochasticGradientDescent.IsAvailable
        if ~askYesOrNo(Html.Wrap(...
                ['This C++ executable is missing or corrupt:'...
            '<br>"<b>' StochasticGradientDescent.GetCmd '</b>"'...
            '<br><br>Maybe try rebuilding by changing clang++ '...
            'to g++ in the build scripts in the same folder...<br>'...
            '<br><center>Try <b>method=Java</b> instead?</center><hr>']))
            return;
        end
        method='MEX';
    end
end
if strcmpi(method, 'Java') || strcmpi(method, 'C++') || strcmpi(method, 'MEX')
    if plotting
        umap.progress_callback=args.progress_callback;
        set(fig, 'NumberTitle', 'off', 'name', info);
        try
            nTh=edu.stanford.facs.swing.StochasticGradientDescent.EPOCH_REPORTS+3;
            figure(fig);
            if args.qf_tree
                puLocation='north++';
            else
                if args.see_training
                    puLocation='south east+';
                else
                    puLocation='south++';
                end
            end
            path=BasicMap.Path;
            pu=PopUp(Html.WrapHr(sprintf(['Using UMAP to reduce '...
                ' <b>%d</b> parameters down to ' ...
                num2str(args.n_components) '...'], nCols-labelCols)), ...
                puLocation, 'Reducing parameters...', false, true, ...
                fullfile(path, 'genieSearch.png'));
            pu.initProgress(nTh);
            pu.pb.setStringPainted(true);
            pu.setTimeSpentTic;
            drawnow;
        catch ex
            args.method='MEX';
            method=umap.setMethod(args.method);
            showMsg(Html.WrapHr(['Could not load umap.jar for Java method'...
                '<br><br>Switching optimize_layout method to "MEX" ']), ...
                'Problem with JAVA...', 'south west', false, false);
        end
    end
end
tc=tic;
if plotting
    if ispc
        left=.21;
        width=.61;
        height=.145;
        lbl=annotation(fig, 'textbox','String', {['\color{blue}Running '...
            info], ['\fontsize{9}' info2]}, 'units', 'normalized', ...
            'position', [left .4 width height], 'fontSize', 11, ...
            'HorizontalAlignment', 'center');
    else
        left=.21;
        width=.61;
        height=.131;
        lbl=annotation(fig, 'textbox','String', {['\color{blue}Running '...
            info], ['\fontsize{11}' info2]}, 'units', 'normalized', ...
            'position', [left .4 width height], 'fontSize', 13, ...
            'HorizontalAlignment', 'center');
    end
    updatePlot;
end
paramAnnotation=[];
if umap.verbose
    txt=sprintf(['n\\_neighbors=\\color{blue}%d\\color{black}, '...
        'min\\_dist=\\color{blue}%s\\color{black}, '...
        'metric=\\color{blue}%s\\color{black},'...
        'randomize=\\color{blue}%d\\color{black}, '...
        'labels=\\color{blue}%d'], ...
        umap.n_neighbors, num2str(umap.min_dist), umap.metric,...
        ~umap.random_state, nLabels); 
    %disp(txt);
    if plotting
        paramAnnotation=annotation(fig, 'textbox','String', txt,...
            'units', 'normalized', 'position', [.03 .94 .92 .05],...
            'fontSize', 9, 'HorizontalAlignment', 'center');
        drawnow;
    end
end
if ~isempty(template_file)
    if ~isempty(umap.pythonTemplate)
        args.python=true;
    end
    if ~args.python 
        if ~args.joined_transform
            reduction=umap.transform(inData);
        else
            reduction=umap.transform2(inData);
        end
    else
        if isempty(umap.pythonTemplate) || ~exist(umap.pythonTemplate, 'file')
            reduction=[];
            msg(Html.WrapHr('Python template file not found'),  8,...
                'south', 'Error...', 'error.png');
        else
            inFile=[tempname '.csv'];
            reduction=UmapPython.Go(inData,inFile, [], ...
                lower(umap.metric), umap.n_neighbors, ...
                umap.min_dist, umap.n_components, [], ...
                umap.pythonTemplate);
        end
    end
else
    if ~args.python
        if ~hasLabels
            reduction = umap.fit_transform(inData);
        else
            reduction = umap.fit_transform(inData, labels);
            if ~isempty(reduction)
                if ~isempty(labelMap)
                    umap.setSupervisors(labels, labelMap, curAxes);
                end
            end
        end
    else
        inFile=[tempname '.csv'];
        if ~hasLabels
            labels=[];
        end
        reduction=UmapPython.Go(inData,inFile, [], ...
            lower(umap.metric), umap.n_neighbors, ...
            umap.min_dist, umap.n_components, labels);
        pythonTemplate=fullfile(fileparts(inFile), ...
            [UmapPython.PYTHON_TEMPLATE '.umap']);
        if ~isempty(reduction)
            umap.embedding=reduction;
            umap.raw_data=inData;
            if ~isempty(labelMap)
                umap.setSupervisors(labels, labelMap, curAxes);
            end
        end
    end
end
if ~isempty(paramAnnotation)
    set(paramAnnotation, 'visible', 'on');
end
if ~isempty(reduction)
    if plotting
        figure(fig);
        if ~exist('pu', 'var')
            pu=PopUp('Updating plot', 'north east', [], false);
        end
        delete(lbl);
        if strcmpi(method, 'Java') || strcmpi(method, 'C++') || strcmpi(method, 'MEX')
            pu.pb.setString('All done');
        end
        updatePlot(reduction, true)
        annotation(fig, 'textbox', 'String', ['Compute time=\color{blue}' ...
            String.MinutesSeconds(toc(tick))],'units', 'normalized', ...
            'position', [.65 .01 .33 .05], 'fontSize', 9)
        if isempty(template_file) && args.ask_to_save_template
            if isequal('Yes', questdlg({'Save this UMAP reduction', ...
                    'as template to accelerate reduction', ...
                    'for compatible other data sets?'}))
                if length(parameter_names)~=size(inData, 2)
                    showMsg(Html.WrapHr(sprintf(['<b>Can not create '...
                        'template</b> ...<br>'...
                        '%d parameter_names ...but data has %d parameters?'], ...
                        length(parameter_names), size(inData,2))));
                else
                    umap.prepareForTemplate(curAxes);
                    if ischar(csv_file_or_data)
                        Template.Save(umap, csv_file_or_data);
                    else
                        Template.Save(umap, fullfile(pwd, 'template.csv'));
                    end
                end
            end
        end
    end
    if ~strcmpi(args.verbose, 'none')
        fprintf('UMAP reduction finished (cost %s)\n', ...
            String.MinutesSeconds(toc(tick)));
    end
else
    msg('Parameter reduction was cancelled or not done');
end
if plotting
    if exist('pu', 'var') && isa(pu, 'PopUp')
        pu.stop;
        pu.dlg.dispose;
    end
else
    if isSupervising
        pu=[];
        if nargout>3
            if ~strcmpi(args.verbose, 'none')
                disp('Setting supervisor labels');
            end
            nSupervisors=umap.supervisors.computeAndMatchClusters( ...
                reduction, args.match_supervisors, []);
            if nSupervisors>0
                extras.supervisorMatchedLabels=umap.supervisors.supervise(...
                    reduction, false, args.match_supervisors);
            end
        end
        doQfs(reduction);
    end
end
if (nargout>1 || ~isempty(save_template_file)) && ~isempty(reduction)
    if plotting
        umap.prepareForTemplate(curAxes);
    else
        umap.prepareForTemplate;
    end
    if args.python
        if  ~isempty(save_template_file)
            [f1, f2]=fileparts(save_template_file);
            if isempty(f1)
                f1=pwd;
            end
        elseif ischar(csv_file_or_data)
            [f1, f2]=fileparts(csv_file_or_data);
            if isempty(f1)
                f1=pwd;
            end
            f2=[f2 '.umap'];
        else
            f1=[];
        end
        if ~isempty(f1)
            umap.pythonTemplate=fullfile(f1, [f2 '.python']);
            movefile(pythonTemplate, umap.pythonTemplate, 'f');
        end
    end
    if  ~isempty(save_template_file)
        save(save_template_file, 'umap');
    end
end    

if nargout>2 || ~strcmpi(args.cluster_output, 'none')
    if nargout<3 && ~strcmpi(args.cluster_output, 'graphic')
        warning('No clusterIdentifiers output argument');
    else
        clusterIdentifiers=doClusters(reduction);
        if isempty(clusterIdentifiers)
            dispNoDbScan;
        end
    end
end
globals.save;

    function clues=doClusters(data)
        if isempty(data) || strcmpi('none', args.cluster_detail)
            clues=[];
        else
            [mins, maxs]=Supervisors.GetMinsMaxs(data);
            [nClues, clues]=Density.FindClusters(reduction, ...
                args.cluster_detail, args.cluster_method_2D, true,...
                args.epsilon, args.minpts, args.dbscan_distance, mins, maxs);
            if strcmpi(args.cluster_output, 'graphic')
                if ~isempty(clues)
                    if ~exist('xLabel', 'var')
                        dimInfo=sprintf('  %dD\\rightarrow%dD', nCols-labelCols, ...
                            args.n_components);
                        xLabel=['UMAP-X' dimInfo];
                        yLabel=['UMAP-Y' dimInfo];
                        zLabel=['UMAP-Z' dimInfo];
                    end
                    cp=ClusterPlots.Go([], data, clues, [], xLabel, ...
                        yLabel, zLabel, true, [], false, true, false);
                    if nCols==2
                        if isequal('dbscan', args.cluster_method_2D)
                            clue='dbscan';
                        else
                            clue='dbm';
                        end
                    else
                        clue='dbscan';
                    end
                    a=annotateClues(get(cp.ax, 'Parent'), args.cluster_detail, clue, args.epsilon, args.minpts, args.dbscan_distance);
                end
            end
        end
    end

    function lbl=annotateClues(fig, detail, clue, epsilon, minpts, dist)
        X=.005;
        Y=.872;
        W=.58;
        H=.115;
        [epsilon, minpts]=Density.GetDbscanParameters(detail, ...
            epsilon, minpts);
        info=['clue method="', clue '", detail="' detail '"'];
        info2=['epsilon=' num2str(epsilon) ', minpts=' num2str(minpts) ...
            ', dbscan distance="' dist '"'];
        lbl=annotation(fig, 'textbox','String', ...
            {['\color{blue} ' info], ['\fontsize{10} ' info2]}, ...
            'units', 'normalized', 'position', [X Y W H],...
            'fontSize', 11, 'HorizontalAlignment', 'center');
    end

    function updatePlot(data, lastCall)
        labelsDone=true;
        if nargin<2
            lastCall=false;
        end
        if nargin>0
            if isempty(xLabel)
                dimInfo=sprintf('  %dD\\rightarrow%dD', nCols-labelCols, ...
                    args.n_components);
                xLabel=['UMAP-X' dimInfo];
                yLabel=['UMAP-Y' dimInfo];
                if args.n_components>2
                    zLabel=['UMAP-Z' dimInfo];
                end
            end
            if args.n_components>2
                nD=size(data, 2);
                assert(nD==args.n_components);
                if ~plotLabels(data, lastCall)
                    labelsDone=false;
                    if args.frequencyDensity3D
                        Gui.PlotDensity3D(curAxes, data, 64, 'iso',...
                            xLabel, yLabel, zLabel);
                    else
                        Gui.PlotNeighDist3D(curAxes, data, ...
                            args.n_neighbors);
                    end
                end
                if args.n_components>3
                    title(curAxes, ['NOTE:  Only 3 of \color{red}' ...
                        num2str(args.n_components) ...
                        ' dimensions being shown...']);
                end
            else
                if plotLabels(data, lastCall)
                    if ~isempty(paramAnnotation)
                        set(paramAnnotation, 'visible', 'off');
                    end
                else
                    labelsDone=false;
                    if lastCall
                        ProbabilityDensity2.Draw(curAxes, data, ...
                            true, true, true, .05);
                    else
                        ProbabilityDensity2.Draw(curAxes, data);
                    end
                end
            end
        end
        if ~labelsDone
            if nargin>0
                umap.adjustLims(curAxes, data );
            else
                umap.adjustLims(curAxes);
            end
            xlabel(curAxes, xLabel);
            ylabel(curAxes, yLabel);
            if args.n_components>2
                zlabel(curAxes, zLabel);
            end
            grid(curAxes, 'on')
            set(curAxes, 'plotboxaspectratio', [1 1 1])
            if lastCall
                Gui.StretchLims(curAxes, data, .04);
            end
        end
        if lastCall
            if args.n_components==2
                if isSupervising
                    if ~hasLabels
                        supervisors.drawClusterBorders(curAxes);
                    end
                end
            end
        end
        drawnow;
    end

    function ok=plotLabels(data, lastCall)
        if hasLabels
            ok=true;
            if lastCall
                Supervisors.Plot(data, labels, labelMap, nCols-labelCols, ...
                    umap, curAxes, true, false, args.contour_percent);
                Gui.StretchLims(curAxes, data, .04);
                if ~isempty(umap.supervisors)
                    umap.supervisors.prepareForTemplate;
                    if args.qf_tree
                        umap.supervisors.inputData=umap.raw_data;
                        [~,qft]=umap.supervisors.qfTreeSupervisors(true, ...
                            [], 'UMAP training set');
                        if ~isempty(qft) && ~isempty(qft.fig)
                            extras.qft=qft;
                        end
                    end
                end
            else
                Supervisors.Plot(data, labels, labelMap, nCols-labelCols, ...
                    umap, curAxes, false, false, args.contour_percent);
            end
        elseif isSupervising
            ok=true;
            if lastCall
                if isempty(supervisors.embedding)
                    supervisors.embedding=umap.embedding;
                end
                if firstPlot && args.see_training
                    [curAxes, ~,~,extras.supervisorMatchedLabels]...
                        =supervisors.plotTrainingAndTestSets(...
                        data, curAxes, umap, pu, args.match_supervisors, ...
                        true);
                    
                else
                    [~,extras.supervisorMatchedLabels]=...
                        supervisors.plotTestSet(umap, curAxes, data, pu, ...
                        args.match_supervisors, true, true);
                end
                doQfs(data);
                drawnow;
            else
                if firstPlot && args.see_training
                    curAxes=supervisors.plotTrainingAndTestSets(...
                        data, curAxes, umap, pu, progressMatchType, false);
                else
                    supervisors.plotTestSet(umap, curAxes, ...
                        data, pu, progressMatchType, false);
                end
                firstPlot=false;
            end
        else
            ok=false;
        end
    end

    function doQfs(data)
        if args.qf_tree || args.qf_dissimilarity
            umap.supervisors.inputData=umap.raw_data;
            cascading={};
            hasFig=exist('fig', 'var');
            if hasFig
                scrFig=fig;
            else
                scrFig=[];
            end
            if args.qf_tree
                if ~strcmpi(args.verbose, 'none')
                    disp('Computing QF tree(s)');
                end
                [~,qft]=umap.supervisors.qfTreeSupervisors(false, pu);
                if ~isempty(qft) && ~isempty(qft.fig)
                    extras.qftSupervisors=qft;
                    cascading{end+1}=qft.fig;
                    Gui.CascadeFigs(cascading, false, true, 15, 2, ...
                        true, false, scrFig);
                    if hasFig
                        figure(scrFig);
                    end
                    [~,qft]=umap.supervisors.qfTreeSupervisees(data, ...
                        inData, false, pu);
                    if ~isempty(qft)
                        cascading{end+1}=qft.fig;
                        Gui.CascadeFigs(cascading, false, true, 40, 2, ...
                            true, false, scrFig);
                        extras.qft=qft;
                    end
                end
            end
            if args.qf_dissimilarity
                if ~strcmpi(args.verbose, 'none')
                    disp('Computing QF dissimilarity scores');
                end
                if hasFig
                    figure(scrFig);
                end
                [~,qfd]=umap.supervisors.qfDissimilarity(data, ...
                    inData, false, pu);
                if ~isempty(qfd)
                    cascading{end+1}=qfd.fig;
                    Gui.CascadeFigs(cascading, false, true, 40, 2, ...
                        true, false, scrFig);
                    Gui.Locate(qfd.qHistFig, qfd.fig, 'south+');
                    figure(qfd.qHistFig);
                    if ismac && 2==size(get(0, 'MonitorPositions'),1)
                        drawnow;
                        Gui.Locate(qfd.qHistFig, qfd.fig,...
                            'south+');
                    end
                    extras.qfd=qfd;
                end
                if hasFig
                    figure(scrFig);
                end
            end
        end
    end

    
    function keepComputing=progress_report(objectOrString)
        keepComputing=~pu.cancelled;
        if ischar(objectOrString)
            if ~isequal(objectOrString, ...
                    StochasticGradientDescent.FINDING_ISLANDS) 
                pu.pb.setValue(pu.pb.getValue+1);
            end
            pu.pb.setString(objectOrString);
            pu.pack;
            pu.showTimeSpent;
            return;
        end 
        done=objectOrString.getEpochsDone-1;
        toDo=objectOrString.getEpochsToDo;
        pu.pb.setValue(3+(pu.pb.getMaximum*(done/toDo)));
        pu.pb.setMaximum(toDo);
        pu.pb.setString(sprintf('%d/%d epochs done', done, toDo));
        if isvalid(lbl)
            delete(lbl);
        end
        updatePlot(objectOrString.getEmbedding);
        pu.showTimeSpent;
    end

    
    function file=downloadCsv
        if ispc
            prompt='Specify name & folder for saving zip file download';
        else
            prompt=Html.WrapHr(...
                ['Please specify the name and folder for the'...
                '<br>zip file being downloaded'...
                '<br>(which will be unzipped after)']);
        end
        [fldr, file]=FileBasics.UiPut(pwd, 'samplesFromHerzenbergLab.zip', ...
            prompt);
        if isnumeric(fldr)
            file=[];
            return;
        end
        pu=PopUp('Downloading & unzipping samples');        
        zipFile=fullfile(fldr, file);
        websave(zipFile, ...
            'http://cgworkspace.cytogenie.org/GetDown2/demo/samples.zip');
        unzip(zipFile);
        pu.close;
        file=fullfile(fldr, 'sample10k.csv');
    end

    function ok=validateCallback(x)
        ok=isequal('function_handle', class(x));
    end

    function ok=validateParameterNames(x)
        ok=false;
        if iscell(x)
            N=length(x);
            if N>0
                for i=1:N
                    if ~ischar(x{i})
                        ok=false;
                        return;
                    end
                end
                ok=true;
            end
        end
    end

    function ok=validateClusterDetail(x)
        ok=false;
        if ischar(x) && ~isempty(x)
            x=lower(x);
            idx=StringArray.IndexOf(Density.DETAILS, x);
            ok=idx>0;
        end
    end



    function p=parseArguments()
        p = inputParser;
        defaultMetric = 'euclidean';
        expectedMetric = {'precomputed', 'euclidean', 'l2', 'manhattan', 'l1',...
        'taxicab', 'cityblock', 'seuclidean', 'standardised_euclidean',...
        'chebychev', 'linfinity', 'linfty', 'linf', 'minkowski',...
        'mahalanobis', 'cosine', 'correlation', 'hamming', 'jaccard',...
        'spearman'};
        defaultVerbose= 'graphic';
        expectedVerbose = {'graphic','text','none'};
        defaultMethod='MEX';
        expectedMethod={'Java', 'C', 'C vectorized', 'MATLAB', 'MATLAB vectorized',...
            'MATLAB experimental', 'MATLAB experimental 2', 'C++', 'MEX'};
        addOptional(p,'csv_file_or_data',[],@(x) ischar(x) || isnumeric(x));
        addParameter(p,'save_template_file',[], @ischar);
        addParameter(p,'ask_to_save_template', false, @islogical);
        addParameter(p,'randomize', false, @islogical);
        addParameter(p,'template_file',[], @(x) ischar(x) || isa(x, 'UMAP'));
        addParameter(p,'n_neighbors', 30, @(x) isnumeric(x) && x>2 && x<200);
        addParameter(p,'min_dist', .3, @(x) isnumeric(x) && x>.05 && x<.8);
        addParameter(p,'metric', defaultMetric, ...
            @(x) any(validatestring(x,expectedMetric)));
        addParameter(p,'n_epochs',[], @(x) isnumeric(x) && x>4);
        addParameter(p,'verbose',defaultVerbose,...
            @(x)islogical(x) || any(validatestring(x,expectedVerbose)));
        addParameter(p,'method',defaultMethod,...
            @(x) any(validatestring(x,expectedMethod)));
        addParameter(p, 'parameter_names', {}, @validateParameterNames);
        addParameter(p, 'progress_callback', ...
            @(javaObject)progress_report(javaObject), @validateCallback);
        addParameter(p,'label_column',0,@(x) isnumeric(x) && x>0);
        addParameter(p,'label_file',[], @ischar);
        addParameter(p,'n_components', 2, @(x) isnumeric(x) && x>=2 && x<101);
        addParameter(p,'frequencyDensity3D', true, @islogical);
        addParameter(p,'match_supervisors', 1, @(x)isnumeric(x) && x>=0 && x<=3);
        addParameter(p,'match_3D_limit', 20000, @(x)isnumeric(x)&&x>=0);
        addParameter(p,'qf_dissimilarity', false, @islogical);
        addParameter(p,'qf_tree', false, @islogical);
        addParameter(p,'joined_transform', false, @islogical);
        addParameter(p,'python', false, @islogical);
        addParameter(p,'see_training', false, @islogical);
        addParameter(p,'cluster_detail', 'very high', @validateClusterDetail);
        expectedMethod={'dbm', 'dbscan'};
        addParameter(p, 'cluster_method_2D', 'dbm', ...
            @(x)any(validatestring(x,expectedMethod)));
        addParameter(p,'minpts', 5, @(x) isnumeric(x) && x>=3 && x<1501);
        addParameter(p,'epsilon', .6, @(x) isnumeric(x) && x>.1 && x<100);
        expectedClusterOutput={'graphic', 'numeric', 'none'};
        addParameter(p, 'cluster_output', 'none', ...
            @(x)any(validatestring(x,expectedClusterOutput)));
        addParameter(p,'dbscan_distance', 'euclidean', ...
            @(x) any(validatestring(x,expectedMetric)));
        addParameter(p,'contour_percent', 10, ...
            @(x) isnumeric(x) && x>=0 && x<=25);

    end
end

function dispNoDbScan
        warning(['dbscan for clustering in 3+D is not available ... '...
            '\nDownload from MathWorks File Exchange: '...
            'https://www.mathworks.com/matlabcentral/fileexchange/52905-dbscan-clustering-algorithm']);
end
