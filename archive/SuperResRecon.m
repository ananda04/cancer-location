%% Script for Super-Resolutioning Reconstructions

% set to folder containing all cases
files = dir('*.mat');
files = struct2cell(files);
files = files(1,:)';

for iter = 1 : length(files)

    load(files{iter})
    r = reconstructedData;
    dim = size(r);

    % making super resolution
    rS = zeros(dim(1)*3,dim(2)*3,dim(3));
    dimS = size(rS);
    for i = 1 : size(rS,1)
        for j = 1 : size(rS,2)

            % dealing with boundary pixels
            if i == 1
                rS(i,j,:) = r(ceil(i/3),ceil(j/3),:);
                continue
            elseif i == size(rS,1)
                rS(i,j,:) = r(ceil(i/3),ceil(j/3),:);
                continue
            elseif j == 1
                rS(i,j,:) = r(ceil(i/3),ceil(j/3),:);
                continue
            elseif j == size(rS,2)
                rS(i,j,:) = r(ceil(i/3),ceil(j/3),:);
                continue
            end

            % interpolating all other pixels spectra
            % UL
            if ismember(i,1:3:(size(rS,1)-2)) && ismember(j,1:3:(size(rS,2)-2))
                rS(i,j,:) = (2 .* r(ceil(i/3),ceil(j/3),:)) + r(ceil(i/3)-1,ceil(j/3)-1,:);
                continue  
            % U
            elseif ismember(i,1:3:(size(rS,1)-2)) && ismember(j,2:3:(size(rS,2)-1))
                rS(i,j,:) = (2 .* r(ceil(i/3),ceil(j/3),:)) + r(ceil(i/3)-1,ceil(j/3),:);
                continue  
            % UR
            elseif ismember(i,1:3:(size(rS,1)-2)) && ismember(j,3:3:size(rS,2))
                rS(i,j,:) = (2 .* r(ceil(i/3),ceil(j/3),:)) + r(ceil(i/3)-1,ceil(j/3)+1,:);
                continue 
            % L
            elseif ismember(i,2:3:(size(rS,1)-1)) && ismember(j,1:3:(size(rS,2)-2))
                rS(i,j,:) = (2 .* r(ceil(i/3),ceil(j/3),:)) + r(ceil(i/3),ceil(j/3)-1,:);
                continue 
            % C
            elseif ismember(i,2:3:(size(rS,1)-1)) && ismember(j,2:3:(size(rS,2)-1)) % central pixels set to old pixel spectra
                rS(i,j,:) = 3 .* r(ceil(i/3),ceil(j/3),:); % 3 just to make similar intensity to others if not normalizing
                continue
            % R
            elseif ismember(i,2:3:(size(rS,1)-1)) && ismember(j,3:3:size(rS,2))
                rS(i,j,:) = (2 .* r(ceil(i/3),ceil(j/3),:)) + r(ceil(i/3),ceil(j/3)+1,:);
                continue 
            % DL
            elseif ismember(i,3:3:size(rS,1)) && ismember(j,1:3:(size(rS,2)-2))
                rS(i,j,:) = (2 .* r(ceil(i/3),ceil(j/3),:)) + r(ceil(i/3)+1,ceil(j/3)-1,:);
                continue 
            % D
            elseif ismember(i,3:3:size(rS,1)) && ismember(j,2:3:(size(rS,2)-1))
                rS(i,j,:) = (2 .* r(ceil(i/3),ceil(j/3),:)) + r(ceil(i/3)+1,ceil(j/3),:);
                continue 
            % DR
            elseif ismember(i,3:3:size(rS,1)) && ismember(j,3:3:size(rS,2))
                rS(i,j,:) = (2 .* r(ceil(i/3),ceil(j/3),:)) + r(ceil(i/3)+1,ceil(j/3)+1,:);
                continue 
            end    
        end
    end

    % saving super-res
    reconstructedData = rS;
    save(files{iter},'qvals','reconstructedData','scatSliceLocations','transScan','transSliceLocations');

    disp(['Case ',num2str(iter),' finished.'])   
end
