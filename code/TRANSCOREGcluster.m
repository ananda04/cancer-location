% Description: This script creates the Co-registration image(MRI+Reconstructed image) 
%              and pulls cancerous and non-cancerous data from image to create cluster, 
%              spectra line graphs, and plot spatial locations on reconstructed image 
% Author:      Arnav Nanda 
%              Apex High School 
%              Fitzpatrick Center for Interdisciplinary Engineering, Medicine and Applied Science (CIEMAS), Duke University  



% Creating Co-registration images 

    % loading images
    a = imread('TRANSCOREG.jpg')
    
    % reshaping coreg
    a1 = imresize(a, [146 72])
    
    % Creating image
    figure(6); imshow(a1)

% Spectra line graphs     
    load('ReconResults_Model_3_BrainScan_GBMT_FullFanExtent_400iter_M3.mat') 
    
    figure(4); imagesc(squeeze(reconstructedData(:, :, 17)))
    
    % Drawing cancerous spectra lines 
    figure(1); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Cancer Spectra")
    hold on; plot(qvals, squeeze(reconstructedData(87,11,:)))
    hold on;figure(1); plot(qvals, squeeze(reconstructedData(89,15,:))) 
    hold on;figure(1); plot(qvals, squeeze(reconstructedData(81,20,:)))
    hold on; figure(1); plot(qvals, squeeze(reconstructedData(85,20,:)))
    hold on;figure(1); plot(qvals, squeeze(reconstructedData(88,19,:)))
    hold on;figure(1); plot(qvals, squeeze(reconstructedData(89,14,:)))
    hold on;figure(1); plot(qvals, squeeze(reconstructedData(72,28,:)))
    hold on;figure(1); plot(qvals, squeeze(reconstructedData(84,25,:)))
    hold on;figure(1); plot(qvals, squeeze(reconstructedData(84,21,:)))
    hold on;figure(1); plot(qvals, squeeze(reconstructedData(89,16,:)))

    
    % Drawing White Matter spectra lines 
    figure(2); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("No Cancer Spectra"); 
    hold on;figure(2); plot(qvals, squeeze(reconstructedData(78,57,:)))
    hold on;figure(2); plot(qvals, squeeze(reconstructedData(60,57,:)))
    hold on;figure(2); plot(qvals, squeeze(reconstructedData(65,48,:)))
    hold on;figure(2); plot(qvals, squeeze(reconstructedData(69,51,:)))
    hold on;figure(2); plot(qvals, squeeze(reconstructedData(76,48,:)))
    hold on;figure(2); plot(qvals, squeeze(reconstructedData(78,51,:)))
    hold on;figure(2); plot(qvals, squeeze(reconstructedData(59,21,:)))
    hold on;figure(2); plot(qvals, squeeze(reconstructedData(64,44,:)))
    hold on;figure(2); plot(qvals, squeeze(reconstructedData(71,44,:)))
    hold on;figure(2); plot(qvals, squeeze(reconstructedData(75,52,:)))
    
    %adding values
    redline = squeeze(reconstructedData(87,11,:))+squeeze(reconstructedData(89,15,:))+squeeze(reconstructedData(81,20,:))+squeeze(reconstructedData(85,20,:))+squeeze(reconstructedData(88,19,:))+squeeze(reconstructedData(89,14,:))+squeeze(reconstructedData(72,28,:))+squeeze(reconstructedData(84,25,:))+squeeze(reconstructedData(84,21,:))+squeeze(reconstructedData(89,16,:))
    blueline = squeeze(reconstructedData(78,57,:))+squeeze(reconstructedData(60,57,:))+squeeze(reconstructedData(65,48,:))+squeeze(reconstructedData(69,51,:))+squeeze(reconstructedData(76,48,:))+squeeze(reconstructedData(78,51,:))+squeeze(reconstructedData(59,21,:))+squeeze(reconstructedData(64,44,:))+squeeze(reconstructedData(71,44,:))+squeeze(reconstructedData(75,52,:))
    
    % averaging lines
    avgredline = redline/10
    avgblueline = blueline/10

    %plotting with error 
    errred = erf(avgredline)
    errblue = erf(avgblueline)

    % Plot graphs with error bars
    figure(3); xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Spectra Line Graph with Error Bars"); 
    figure(3); e1 = errorbar(qvals, avgredline, errred); 
    e1.Color = 'r'
    hold on;
    figure(3); e2 = errorbar(qvals, avgblueline,errblue); 
    e2.Color = 'b'
    hold on;
    legend([e1(1),e2(1)],'Cancer', 'Non-cancer' )
    xlabel('q [1/A]'); ylabel('XRD amplitude [arb]'); title("Average Spectra")

%Checking spatial locations

    % Cancer data
    x = [11,15,20,20,19,14,28,25,21,16]
    y = [87,89,81,85,88,89,72,84,84,89]
    
    % Non-cancerous data
    x1 = [57,57,48,51,48,51,21,44,44,52]
    y1 = [78,60,65,69,76,78,59,64,71,75]
    
    % Picture
    figure(4); imagesc(squeeze(reconstructedData(:, :, 17)))
    hold on; figure(4); plot(x,y,'rx')
    hold on; figure(4); plot(x1,y1, 'kx')
    legend('Cancer', 'Non-cacner');
    title("Spatial Locations")

% Creating Clusters

    % Creating Cancer Array
    A = [squeeze(reconstructedData(87,11,:)),squeeze(reconstructedData(89,15,:)),squeeze(reconstructedData(81,20,:)),squeeze(reconstructedData(85,20,:)),squeeze(reconstructedData(88,19,:)),squeeze(reconstructedData(89,14,:)),squeeze(reconstructedData(72,28,:)),squeeze(reconstructedData(84,25,:)),squeeze(reconstructedData(84,21,:)),squeeze(reconstructedData(89,16,:))]
    cancerarray = transpose(A)
    
    % Creating Non-Cancer Array
    B = [squeeze(reconstructedData(78,57,:)),squeeze(reconstructedData(60,57,:)),squeeze(reconstructedData(65,48,:)),squeeze(reconstructedData(69,51,:)),squeeze(reconstructedData(76,48,:)),squeeze(reconstructedData(78,51,:)),squeeze(reconstructedData(59,21,:)),squeeze(reconstructedData(64,44,:)),squeeze(reconstructedData(71,44,:)),squeeze(reconstructedData(75,52,:))]
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

% QVAL coloring - helps identify cancerous regions
%    figure(4); imagesc(qvals)
%
%    figure(6); imagesc(reconstructedData(:,:,15))
%    figure(13);imagesc(reconstructedData(:,:,16))
%    figure(14);imagesc(reconstructedData(:,:,17))
%    figure(9);imagesc(reconstructedData(:,:,18))
%    figure(10);imagesc(reconstructedData(:,:,19))
%    figure(11);imagesc(reconstructedData(:,:,20))
%    figure(12);imagesc(reconstructedData(:,:,21))
     