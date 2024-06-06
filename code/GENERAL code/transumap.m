% Author:      Arnav Nanda 
%              Duke University, Pratt School of Engineering 

    x = [251,350,464,355,444,563,355,375,315,399] 
    y = [1572,1631,1601,1527,1388,1433,1448,1374,1448,1314]
    A = [x;y]
    cancerarray = transpose(A)
    x1 = [1345,1340,1349,1369,1330,1320,1310,1270,1265,419]
    y1 = [1769,1695,1715,1730,1586,1542,1403,1369,1097,1106]
    B = [x1;y1]
    noncancerarray = transpose(B)

    % RUN UMAP
    import umapFileExchange.*;
    [reduction] = run_umap(cancerarray)
    [reduction1] = run_umap(noncancerarray)
    
    figure(5); plot(reduction(:,1), reduction(:,2),'r.')
    hold on;
    figure(5); plot(reduction1(:,1), reduction1(:,2),'b.')
    hold on;
    legend('Cancer', 'Non-cacner' );
    title('Transpose')
