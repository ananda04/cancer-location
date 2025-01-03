% Author:      Arnav Nanda 
%              Apex High School 
%              Fitzpatrick Center for Interdisciplinary Engineering, Medicine and Applied Science (CIEMAS), Duke University

    load('ReconResults_Brain_158_NT_30s_300iter_M3_Try1.mat')
    graySpec = [squeeze(reconstructedData(82,13,:)),squeeze(reconstructedData(79,13,:)),squeeze(reconstructedData(76,12,:)),squeeze(reconstructedData(70,11,:)),squeeze(reconstructedData(66,13,:)),squeeze(reconstructedData(65,17,:)),squeeze(reconstructedData(69,20,:)),squeeze(reconstructedData(72,22,:)),squeeze(reconstructedData(63,20,:)),squeeze(reconstructedData(58,22,:)),squeeze(reconstructedData(53,30,:)),squeeze(reconstructedData(48,35,:)),squeeze(reconstructedData(51,39,:)),squeeze(reconstructedData(54,40,:)),squeeze(reconstructedData(53,43,:)),squeeze(reconstructedData(48,43,:)),squeeze(reconstructedData(48,46,:)),squeeze(reconstructedData(50,48,:)),squeeze(reconstructedData(58,57,:)),squeeze(reconstructedData(71,54,:))]
    whiteSpec = [squeeze(reconstructedData(60,47,:)),squeeze(reconstructedData(62,45,:)),squeeze(reconstructedData(64,43,:)),squeeze(reconstructedData(67,41,:)),squeeze(reconstructedData(70,40,:)),squeeze(reconstructedData(73,40,:)),squeeze(reconstructedData(73,36,:)),squeeze(reconstructedData(68,37,:)),squeeze(reconstructedData(76,42,:)),squeeze(reconstructedData(79,41,:)), squeeze(reconstructedData(78,37,:)),squeeze(reconstructedData(74,33,:)),squeeze(reconstructedData(82,35,:)),squeeze(reconstructedData(84,39,:)),squeeze(reconstructedData(86,38,:)),squeeze(reconstructedData(87,35,:)),squeeze(reconstructedData(86,26,:)),squeeze(reconstructedData(87,29,:)),squeeze(reconstructedData(81,23,:)),squeeze(reconstructedData(79,21,:))]
    graySpec = transpose(graySpec)
    whiteSpec = transpose(whiteSpec)
    avg_graySpec = mean(graySpec,1)
    avg_whiteSpec = mean(whiteSpec,1)
    errGray = std(graySpec)
    errWhite = std(whiteSpec)
    % Plot graphs with error bars
    figure(1); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
    e1 = errorbar(qvals, avg_graySpec,errGray); 
    e1.Color = 'b'
    hold on;
    figure(1); e2 = errorbar(qvals, avg_whiteSpec, errWhite); 
    e2.Color = 'r'
    legend([e1(1),e2(1)], 'Gray Matter','White Matter' )
    xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Average Spectra NT158")
%%
    load('ReconResults_Brain_187_NT_20s_300iter_M3_Try1.mat')
    graySpec = [squeeze(reconstructedData(86,8,:)),squeeze(reconstructedData(86,9,:)),squeeze(reconstructedData(83,9,:)),squeeze(reconstructedData(80,9,:)),squeeze(reconstructedData(77,10,:)),squeeze(reconstructedData(77,14,:)),squeeze(reconstructedData(76,12,:)),squeeze(reconstructedData(73,10,:)),squeeze(reconstructedData(69,11,:)),squeeze(reconstructedData(66,11,:)),squeeze(reconstructedData(66,14,:)),squeeze(reconstructedData(68,16,:)),squeeze(reconstructedData(62,13,:)),squeeze(reconstructedData(59,14,:)),squeeze(reconstructedData(59,16,:)),squeeze(reconstructedData(51,17,:)),squeeze(reconstructedData(55,15,:)),squeeze(reconstructedData(49,23,:)),squeeze(reconstructedData(52,28,:)),squeeze(reconstructedData(55,30,:))]
    whiteSpec = [squeeze(reconstructedData(84,24,:)),squeeze(reconstructedData(80,25,:)),squeeze(reconstructedData(79,26,:)),squeeze(reconstructedData(75,26,:)),squeeze(reconstructedData(72,27,:)),squeeze(reconstructedData(69,27,:)),squeeze(reconstructedData(65,26,:)),squeeze(reconstructedData(60,24,:)),squeeze(reconstructedData(73,29,:)),squeeze(reconstructedData(77,30,:)),squeeze(reconstructedData(82,28,:)),squeeze(reconstructedData(85,30,:)),squeeze(reconstructedData(80,31,:)),squeeze(reconstructedData(82,32,:)),squeeze(reconstructedData(79,35,:)),squeeze(reconstructedData(78,38,:)),squeeze(reconstructedData(80,29,:)),squeeze(reconstructedData(86,26,:)),squeeze(reconstructedData(85,22,:)),squeeze(reconstructedData(85,34,:))] 
    graySpec = transpose(graySpec)
    whiteSpec = transpose(whiteSpec)
    avg_graySpec = mean(graySpec,1)
    avg_whiteSpec = mean(whiteSpec,1)
    errGray = std(graySpec)
    errWhite = std(whiteSpec)
    % Plot graphs with error bars
    figure(1); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
    e1 = errorbar(qvals, avg_graySpec,errGray); 
    e1.Color = 'b'
    hold on;
    figure(1); e2 = errorbar(qvals, avg_whiteSpec, errWhite); 
    e2.Color = 'r'
    legend([e1(1),e2(1)], 'Gray Matter','White Matter' )
    xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Average Spectra NT187")
%%
    load('ReconResults_Brain_176_T2_20s_300iter_M3_Try1.mat')
    graySpec = [squeeze(reconstructedData(73,47,:)),squeeze(reconstructedData(76,49,:)),squeeze(reconstructedData(75,46,:)),squeeze(reconstructedData(70,48,:)),squeeze(reconstructedData(68,49,:)),squeeze(reconstructedData(84,5,:)),squeeze(reconstructedData(81,8,:)),squeeze(reconstructedData(76,10,:)),squeeze(reconstructedData(73,11,:)),squeeze(reconstructedData(70,12,:)),squeeze(reconstructedData(69,13,:)),squeeze(reconstructedData(67,16,:)),squeeze(reconstructedData(82,49,:)),squeeze(reconstructedData(85,48,:)),squeeze(reconstructedData(84,47,:)),squeeze(reconstructedData(64,18,:)),squeeze(reconstructedData(60,21,:)),squeeze(reconstructedData(57,23,:)),squeeze(reconstructedData(55,25,:)),squeeze(reconstructedData(84,49,:))]
    whiteSpec =[squeeze(reconstructedData(82,23,:)),squeeze(reconstructedData(81,22,:)),squeeze(reconstructedData(84,23,:)),squeeze(reconstructedData(82,20,:)),squeeze(reconstructedData(85,25,:)),squeeze(reconstructedData(81,26,:)),squeeze(reconstructedData(84,28,:)),squeeze(reconstructedData(80,29,:)),squeeze(reconstructedData(83,31,:)),squeeze(reconstructedData(80,31,:)),squeeze(reconstructedData(80,33,:)),squeeze(reconstructedData(78,33,:)),squeeze(reconstructedData(75,35,:)),squeeze(reconstructedData(73,37,:)),squeeze(reconstructedData(70,38,:)),squeeze(reconstructedData(68,39,:)),squeeze(reconstructedData(65,41,:)),squeeze(reconstructedData(62,42,:)),squeeze(reconstructedData(77,25,:)),squeeze(reconstructedData(79,26,:))]
    graySpec = transpose(graySpec)
    whiteSpec = transpose(whiteSpec)
    avg_graySpec = mean(graySpec,1)
    avg_whiteSpec = mean(whiteSpec,1)
    errGray = std(graySpec)
    errWhite = std(whiteSpec)
    % Plot graphs with error bars
    figure(1); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
    e1 = errorbar(qvals, avg_graySpec,errGray); 
    e1.Color = 'b'
    hold on;
    figure(1); e2 = errorbar(qvals, avg_whiteSpec, errWhite); 
    e2.Color = 'r'
    legend([e1(1),e2(1)], 'Gray Matter','White Matter' )
    xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Average Spectra NT176")
%%
    load('ReconResults_Brain_177_NT_30s_300iter_M3_Try1.mat')
    graySpec =      [squeeze(reconstructedData(87,8,:)),squeeze(reconstructedData(82,10,:)),squeeze(reconstructedData(82,8,:)),squeeze(reconstructedData(78,5,:)),squeeze(reconstructedData(76,8,:)),squeeze(reconstructedData(79,10,:)),squeeze(reconstructedData(84,11,:)),squeeze(reconstructedData(87,10,:)),squeeze(reconstructedData(66,7,:)),squeeze(reconstructedData(62,11,:)),squeeze(reconstructedData(60,13,:)),squeeze(reconstructedData(57,17,:)),squeeze(reconstructedData(55,20,:)),squeeze(reconstructedData(54,25,:)),squeeze(reconstructedData(55,27,:)),squeeze(reconstructedData(59,27,:)),squeeze(reconstructedData(67,39,:)),squeeze(reconstructedData(68,40,:)),squeeze(reconstructedData(69,44,:)),squeeze(reconstructedData(72,49,:))]
    whiteSpec = [squeeze(reconstructedData(73,12,:)),squeeze(reconstructedData(72,14,:)),squeeze(reconstructedData(74,14,:)),squeeze(reconstructedData(75,15,:)),squeeze(reconstructedData(76,17,:)),squeeze(reconstructedData(76,19,:)),squeeze(reconstructedData(76,20,:)),squeeze(reconstructedData(76,22,:)),squeeze(reconstructedData(75,24,:)),squeeze(reconstructedData(75,26,:)),squeeze(reconstructedData(74,27,:)),squeeze(reconstructedData(74,29,:)),squeeze(reconstructedData(74,30,:)),squeeze(reconstructedData(74,33,:)),squeeze(reconstructedData(74,34,:)),squeeze(reconstructedData(74,36,:)),squeeze(reconstructedData(75,38,:)),squeeze(reconstructedData(85,30,:)),squeeze(reconstructedData(85,29,:)),squeeze(reconstructedData(86,25,:))]
    graySpec = transpose(graySpec)
    whiteSpec = transpose(whiteSpec)
    avg_graySpec = mean(graySpec,1)
    avg_whiteSpec = mean(whiteSpec,1)
    errGray = std(graySpec)
    errWhite = std(whiteSpec)
    % Plot graphs with error bars
    figure(1); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
    e1 = errorbar(qvals, avg_graySpec,errGray); 
    e1.Color = 'b'
    hold on;
    figure(1); e2 = errorbar(qvals, avg_whiteSpec, errWhite); 
    e2.Color = 'r'
    legend([e1(1),e2(1)], 'Gray Matter','White Matter' )
    xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Average Spectra NT187")
