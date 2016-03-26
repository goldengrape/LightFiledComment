%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Given a bunch of parameters, compute the sparse propagation matrix for
%   N 1D attenuation layers in log-space given a desired 4D light field
%
%   input:  drawMode -  0 draw nothing
%                       1 show progress bar
%                       2 show matrix as it is updated
%
%           basisFunctionType - 0 discrete light field ray positions and angles
%                               1 discrete light field ray positions, box
%                               area integration in angle
%                               2 box area integration around position and angle
%                               3 discrete light field ray positions, linear area integration in angle
%
%   Gordon Wetzstein [wetzste1@cs.ubc.ca]
%   PSM Lab | University of British Columbia
%   February 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/
% [row, col, val] = find(accumarray([row, col], val, [], @sum, [], true));
% OR
% [r, c, v] = find(accumarray({row, col}, val))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO:
%   - despite using a single global instance of T, the memory management
%       should be enhanced even more!
%   - linear interpolation for rays (basisFunctionType 1)
%   - linear interpolation for area integration (basisFunctionType 3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function precomputeSparsePropagationMatrixLayers3D( lightFieldAnglesY, lightFieldAnglesX, lightFieldSize, lightFieldResolution, lightFieldOrigin,...
                                                    layerResolution, layerSize, layerOrigin, layerDistance, ...
                                                    basisFunctionType, drawMode )

    % matrix is a global variable
    global T;

    % if large scale mode, momory will be saved by cycling through entire
    % loop and just computing the number of non-zero elements in the
    % matrix, then it'll allocate the memory and fill the matrix
    bLargeScale = true;
                                                    
    if (basisFunctionType<0) || (basisFunctionType>3)
        error(['Basis function type ' num2str(basisFunctionType) ' currently not supported!']);
    end
                                                                  
    % absolute x coordinates for light field pixel centers
    lightFieldPixelSize     = lightFieldSize ./ [lightFieldResolution(3) lightFieldResolution(4)];
    lightFieldPixelCentersX = lightFieldOrigin(2)+lightFieldPixelSize(2)/2:lightFieldPixelSize(2):lightFieldOrigin(2)+lightFieldSize(2)-lightFieldPixelSize(2)/2;
    lightFieldPixelCentersY = lightFieldOrigin(1)+lightFieldPixelSize(1)/2:lightFieldPixelSize(1):lightFieldOrigin(1)+lightFieldSize(1)-lightFieldPixelSize(1)/2;         
        
    % assuming that the angles are sampled at equal distances, get that
    % distance
    lightFieldAngleStep = [0 0];
    if lightFieldResolution(1)>1
        lightFieldAngleStep(1) = lightFieldAnglesY(2)-lightFieldAnglesY(1);
    end
    if lightFieldResolution(2)>1
        lightFieldAngleStep(2) = lightFieldAnglesX(2)-lightFieldAnglesX(1);
    end    
    
    if drawMode == 1
        h               = waitbar(0,'Generating propagation matrix ...');
        numWaitbarCalls = lightFieldResolution(1)*lightFieldResolution(2);
    end 
    
    % scale to be applied to the float-positions to get a basis function index
    lookupScaleY = (layerResolution(1))/layerSize(1);
    lookupScaleX = (layerResolution(2))/layerSize(2);
    
    % if basisFunctionType == 2 (area integration), the integrated area for
    % each spatial sample will have this size)
    areaIntegrationSize = lightFieldPixelSize;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cycle through the entire proceduce once to see how many non-zero
    % elements there are in the matrix        
    
    numNonzeroElements = 0;
    if bLargeScale
        
        % need to run through this twice
        numWaitbarCalls = 2*lightFieldResolution(1)*lightFieldResolution(2);
        
        % for all angles in the light field
        for vyIdx=1:lightFieldResolution(1)
            for vxIdx=1:lightFieldResolution(2)     

                % update waitbar
                if drawMode == 1
                    waitbar( (vxIdx+(vyIdx-1)*lightFieldResolution(2)) / numWaitbarCalls );
                end

                % actual angle in v units
                vy = lightFieldAnglesY(vyIdx);
                vx = lightFieldAnglesX(vxIdx);

                % if we use box-area integration filtering, use a range of angles
                if (basisFunctionType == 1) || (basisFunctionType == 2)
                    
                    % lower boundary of angular box
                    vy = lightFieldAnglesY(vyIdx) - lightFieldAngleStep(1)/2;
                    vx = lightFieldAnglesX(vxIdx) - lightFieldAngleStep(2)/2;

                    % upper boundary of angular box
                    vyUpper = lightFieldAnglesY(vyIdx) + lightFieldAngleStep(1)/2;     
                    vxUpper = lightFieldAnglesX(vxIdx) + lightFieldAngleStep(2)/2;
                    
                % linear angular sampling - full angular step!
                elseif basisFunctionType == 3
                    
                     % lower boundary of angular box
                    vy = lightFieldAnglesY(vyIdx) - lightFieldAngleStep(1);
                    vx = lightFieldAnglesX(vxIdx) - lightFieldAngleStep(2);

                    % upper boundary of angular box
                    vyUpper = lightFieldAnglesY(vyIdx) + lightFieldAngleStep(1);     
                    vxUpper = lightFieldAnglesX(vxIdx) + lightFieldAngleStep(2);
                    
                end

                % shift light field pixels to pixel centers of 1st volume slice
                rayPositionsX = lightFieldPixelCentersX - lightFieldOrigin(3)*vx + layerOrigin(3)*vx;
                rayPositionsY = lightFieldPixelCentersY - lightFieldOrigin(3)*vy + layerOrigin(3)*vy;  

                % if we use box-area integration filtering, use a range of angles
                if basisFunctionType > 0
                    % shift upper light field pixels to pixel centers of 1st layer
                    rayPositionsXUpper = lightFieldPixelCentersX - lightFieldOrigin(3)*vxUpper + layerOrigin(3)*vxUpper;
                    rayPositionsYUpper = lightFieldPixelCentersY - lightFieldOrigin(3)*vyUpper + layerOrigin(3)*vyUpper;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % iterate over all volume slices
                for layer = 1:layerResolution(3)

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % point sampling 
                    
                    % set values in propagation matrix            
                    if basisFunctionType == 0       % nearest interpolation
                                               
                        % get currently effected layer pixels - matrix column indices                       
                        layerPixelIndicesForRaysX = ceil(lookupScaleX .* (rayPositionsX-layerOrigin(2)));                        
                        layerPixelIndicesForRaysY = ceil(lookupScaleY .* (rayPositionsY-layerOrigin(1)));
                                                               
                        % kick out stuff that's outside
                        layerPixelIndicesForRaysX(layerPixelIndicesForRaysX>layerResolution(2)) = 0;
                        layerPixelIndicesForRaysX(layerPixelIndicesForRaysX<1) = 0;
                        layerPixelIndicesForRaysY(layerPixelIndicesForRaysY>layerResolution(1)) = 0;
                        layerPixelIndicesForRaysY(layerPixelIndicesForRaysY<1) = 0;       
                        
                        % update number of non-zero matrix elements
                        numNonzeroElements = numNonzeroElements + numel(find(layerPixelIndicesForRaysX))*numel(find(layerPixelIndicesForRaysY));

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % area integration in angle
                    elseif basisFunctionType == 1   
                    
                         % get currently effected layer pixels - matrix column indices for each ligh ray                      
                        layerPixelIndicesForRaysXLeft   = ceil(lookupScaleX .* (rayPositionsXUpper-layerOrigin(2)));                        
                        layerPixelIndicesForRaysYUp     = ceil(lookupScaleY .* (rayPositionsYUpper-layerOrigin(1)));                    
                        layerPixelIndicesForRaysXRight  = ceil(lookupScaleX .* (rayPositionsX-layerOrigin(2)));                        
                        layerPixelIndicesForRaysYDown   = ceil(lookupScaleY .* (rayPositionsY-layerOrigin(1)));                                                            

                        % this only works if the number of effected pixels is
                        % the same everywhere
                        numLayerPixelsPerRayPixelX = 1 + layerPixelIndicesForRaysXRight(1)-layerPixelIndicesForRaysXLeft(1);
                        numLayerPixelsPerRayPixelY = 1 + layerPixelIndicesForRaysYDown(1)-layerPixelIndicesForRaysYUp(1);

                        % check if it really is the same everywhere
                        if  (sum( 1 - ((layerPixelIndicesForRaysXRight-layerPixelIndicesForRaysXLeft)==(numLayerPixelsPerRayPixelX-1)) ) ~= 0) || ...
                            (sum( 1 - ((layerPixelIndicesForRaysYDown-layerPixelIndicesForRaysYUp)==(numLayerPixelsPerRayPixelY-1)) ) ~= 0)
                            error('Dimensions are wrong!');
                        end                                        
                        
                        % construct matrix columns
                        currentIndicesColumnsX  = zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                        currentIndicesColumnsY  = zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                        
                        for ky=1:numLayerPixelsPerRayPixelY
                            currentIndicesColumnsY(:,:,ky,:) = repmat( (layerPixelIndicesForRaysYUp')+ky-1, [1 lightFieldResolution(4) 1 numLayerPixelsPerRayPixelX]);
                        end

                        for kx=1:numLayerPixelsPerRayPixelX
                            currentIndicesColumnsX(:,:,:,kx) = repmat( layerPixelIndicesForRaysXLeft+kx-1, [lightFieldResolution(3) 1 numLayerPixelsPerRayPixelY 1]);
                        end

                        % remember all valid indices
                        validIndices = ones([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                        validIndices(currentIndicesColumnsY<1)                  = 0;
                        validIndices(currentIndicesColumnsY>layerResolution(1)) = 0;
                        validIndices(currentIndicesColumnsX<1)                  = 0;
                        validIndices(currentIndicesColumnsX>layerResolution(2)) = 0;
                        
                        % update number of non-zero matrix elements
                        numNonzeroElements = numNonzeroElements + sum(validIndices(:));
                    
                        
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % roectangular region around ray in space and angle

                    elseif basisFunctionType == 2  
                        
                        % get currently effected layer pixels - matrix column indices for each ligh ray                      
                        layerPixelIndicesForRaysXLeft   = ceil(lookupScaleX .* (rayPositionsXUpper-areaIntegrationSize(2)/2-layerOrigin(2)));                        
                        layerPixelIndicesForRaysYUp     = ceil(lookupScaleY .* (rayPositionsYUpper-areaIntegrationSize(1)/2-layerOrigin(1)));                    
                        layerPixelIndicesForRaysXRight  = ceil(lookupScaleX .* (rayPositionsX+areaIntegrationSize(2)/2-layerOrigin(2)));                        
                        layerPixelIndicesForRaysYDown   = ceil(lookupScaleY .* (rayPositionsY+areaIntegrationSize(1)/2-layerOrigin(1)));                                                            

                        % this only works if the number of effected pixels is
                        % the same everywhere
                        numLayerPixelsPerRayPixelX = 1 + layerPixelIndicesForRaysXRight(1)-layerPixelIndicesForRaysXLeft(1);
                        numLayerPixelsPerRayPixelY = 1 + layerPixelIndicesForRaysYDown(1)-layerPixelIndicesForRaysYUp(1);

                        % check if it really is the same everywhere
                        if  (sum( 1 - ((layerPixelIndicesForRaysXRight-layerPixelIndicesForRaysXLeft)==(numLayerPixelsPerRayPixelX-1)) ) ~= 0) || ...
                            (sum( 1 - ((layerPixelIndicesForRaysYDown-layerPixelIndicesForRaysYUp)==(numLayerPixelsPerRayPixelY-1)) ) ~= 0)
                            error('Dimensions are wrong!');
                        end                                        
                        
                        % construct matrix columns
                        currentIndicesColumnsX  = zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                        currentIndicesColumnsY  = zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                        currentIndicesColumnsZ  = layer + zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);

                        for ky=1:numLayerPixelsPerRayPixelY
                            currentIndicesColumnsY(:,:,ky,:) = repmat( (layerPixelIndicesForRaysYUp')+ky-1, [1 lightFieldResolution(4) 1 numLayerPixelsPerRayPixelX]);
                        end

                        for kx=1:numLayerPixelsPerRayPixelX
                            currentIndicesColumnsX(:,:,:,kx) = repmat( layerPixelIndicesForRaysXLeft+kx-1, [lightFieldResolution(3) 1 numLayerPixelsPerRayPixelY 1]);
                        end

                        % remember all valid indices
                        validIndices = ones([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                        validIndices(currentIndicesColumnsY<1)                  = 0;
                        validIndices(currentIndicesColumnsY>layerResolution(1)) = 0;
                        validIndices(currentIndicesColumnsX<1)                  = 0;
                        validIndices(currentIndicesColumnsX>layerResolution(2)) = 0;
                        
                        % update number of non-zero matrix elements
                        numNonzeroElements = numNonzeroElements + sum(validIndices(:));
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % shift x positions of rays to next slice
                    rayPositionsX = rayPositionsX - layerDistance*vx;
                    rayPositionsY = rayPositionsY - layerDistance*vy;
                    
                    if basisFunctionType > 0
                        rayPositionsXUpper = rayPositionsXUpper - layerDistance*vxUpper;
                        rayPositionsYUpper = rayPositionsYUpper - layerDistance*vyUpper;
                    end

                end
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
    
    % initialize vectors that will store the non-negative matrix indices and values
    if bLargeScale
        % memory matters - only store datatypes that we actually need
        indexIVector = zeros([numNonzeroElements 1]);
        indexJVector = zeros([numNonzeroElements 1]);
        valueSVector = zeros([numNonzeroElements 1]);
        indexCount   = 1;
    else
        % memory doesn't matter
        indexIVector = [];
        indexJVector = [];
        valueSVector = [];
    end
    
    % for all angles in the light field
	for vyIdx=1:lightFieldResolution(1)
        for vxIdx=1:lightFieldResolution(2)     
        
            % update waitbar
            if drawMode == 1
                currentCalls = (vxIdx+(vyIdx-1)*lightFieldResolution(2));
                if bLargeScale
                    currentCalls = currentCalls + lightFieldResolution(1)*lightFieldResolution(2);
                end
                waitbar( currentCalls / numWaitbarCalls );
            end

            % actual angle in v units
            vy = lightFieldAnglesY(vyIdx);
            vx = lightFieldAnglesX(vxIdx);
            
            % if we use box-area integration filtering, use a range of angles
            if (basisFunctionType==1) || (basisFunctionType==2)
                % lower boundary of angular box
                vy = lightFieldAnglesY(vyIdx) - lightFieldAngleStep(1)/2;
                vx = lightFieldAnglesX(vxIdx) - lightFieldAngleStep(2)/2;
                
                % upper boundary of angular box
                vyUpper = lightFieldAnglesY(vyIdx) + lightFieldAngleStep(1)/2;     
                vxUpper = lightFieldAnglesX(vxIdx) + lightFieldAngleStep(2)/2;
                
            % linear angular sampling - full angular step!
            elseif basisFunctionType == 3

                % lower boundary of angular box
                vy = lightFieldAnglesY(vyIdx) - lightFieldAngleStep(1);
                vx = lightFieldAnglesX(vxIdx) - lightFieldAngleStep(2);

                % upper boundary of angular box
                vyUpper = lightFieldAnglesY(vyIdx) + lightFieldAngleStep(1);     
                vxUpper = lightFieldAnglesX(vxIdx) + lightFieldAngleStep(2);                
                
            end
            
            % shift light field pixels to pixel centers of 1st volume slice
            rayPositionsX = lightFieldPixelCentersX - lightFieldOrigin(3)*vx + layerOrigin(3)*vx;
            rayPositionsY = lightFieldPixelCentersY - lightFieldOrigin(3)*vy + layerOrigin(3)*vy;  
            
            % if we use box-area integration filtering, use a range of angles
            if basisFunctionType > 0
                % shift upper light field pixels to pixel centers of 1st layer
                rayPositionsXUpper = lightFieldPixelCentersX - lightFieldOrigin(3)*vxUpper + layerOrigin(3)*vxUpper;
                rayPositionsYUpper = lightFieldPixelCentersY - lightFieldOrigin(3)*vyUpper + layerOrigin(3)*vyUpper;
                                
                % matrix row indices for all light field rays 
                [rayIndicesX rayIndicesY] = meshgrid(1:lightFieldResolution(4), 1:lightFieldResolution(3));
                rayIndicesVY = vyIdx + zeros(size(rayIndicesX));
                rayIndicesVX = vxIdx + zeros(size(rayIndicesY));
                % convert 4D subscipts to matrix indices
                matrixRowIndicesCurrentLFAngle = sub2ind(lightFieldResolution, rayIndicesVY, rayIndicesVX, rayIndicesY, rayIndicesX);                
                clear rayIndicesX rayIndicesY rayIndicesVY rayIndicesVX;                                
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % iterate over all volume slices
            for layer = 1:layerResolution(3)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % set values in propagation matrix            
                if basisFunctionType == 0       % nearest interpolation
                    
                    % basisFunctionIndicesForRaysX = interp1(basisFunctionsX,basisFunctionsXX,rayPositionsX,'linear',0);
                    
                    % get currently effected layer pixels - matrix column indices for each ligh ray                      
                    layerPixelIndicesForRaysX = ceil(lookupScaleX .* (rayPositionsX-layerOrigin(2)));                        
                    layerPixelIndicesForRaysY = ceil(lookupScaleY .* (rayPositionsY-layerOrigin(1)));   
                                        
                    % kick out stuff that's outside
                    layerPixelIndicesForRaysX(layerPixelIndicesForRaysX>layerResolution(2)) = 0;
                    layerPixelIndicesForRaysX(layerPixelIndicesForRaysX<1) = 0;
                    layerPixelIndicesForRaysY(layerPixelIndicesForRaysY>layerResolution(1)) = 0;
                    layerPixelIndicesForRaysY(layerPixelIndicesForRaysY<1) = 0;                                        
                    
                    % convert to matrix row indices - which rays hit some layer pixels
                    validXIndices = find(layerPixelIndicesForRaysX);
                    validYIndices = find(layerPixelIndicesForRaysY);
                    
                    % turn it into a matrix
                    validXIndices = repmat(validXIndices, [numel(validYIndices) 1]);
                    validYIndices = repmat(validYIndices', [1 size(validXIndices,2)]);
                    % angle indices                    
                    validVXIndices = vxIdx + zeros(size(validXIndices));
                    validVYIndices = vyIdx + zeros(size(validXIndices));

                    % convert 4D subscipts to matrix indices
                    matrixRows = sub2ind(lightFieldResolution, validVYIndices(:), validVXIndices(:), validYIndices(:), validXIndices(:));
                    
                    % convert to matrix column indices               
                    layerPixelIndicesForRaysX = layerPixelIndicesForRaysX(layerPixelIndicesForRaysX~=0);
                    layerPixelIndicesForRaysY = layerPixelIndicesForRaysY(layerPixelIndicesForRaysY~=0);
                    validXXIndices = repmat(layerPixelIndicesForRaysX, [numel(layerPixelIndicesForRaysY) 1]);
                    validYYIndices = repmat(layerPixelIndicesForRaysY', [1 size(layerPixelIndicesForRaysX,2)]);                    
                    validZZIndices = layer + zeros(size(validXXIndices));
   
                    % convert 3D subscripts to matrix indices
                    matrixColumns   = sub2ind(layerResolution, validYYIndices(:), validXXIndices(:), validZZIndices(:));                                        

                    if bLargeScale
                        numCurrentNonzeroElements = numel(matrixRows);
                        indexIVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = matrixRows;
                        indexJVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = matrixColumns;
                        valueSVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = ones(size(matrixRows));
                        indexCount = indexCount + numCurrentNonzeroElements ;
                    else
                        indexIVector = cat(1, indexIVector, matrixRows);
                        indexJVector = cat(1, indexJVector, matrixColumns);
                        valueSVector = cat(1, valueSVector, ones(size(matrixRows)));
                    end
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % rectangular region around ray in angle
                
                elseif basisFunctionType == 1   
                
                    % get currently effected layer pixels - matrix column indices for each ligh ray                      
                    layerPixelIndicesForRaysXLeft   = ceil(lookupScaleX .* (rayPositionsXUpper-layerOrigin(2)));                        
                    layerPixelIndicesForRaysYUp     = ceil(lookupScaleY .* (rayPositionsYUpper-layerOrigin(1)));                    
                    layerPixelIndicesForRaysXRight  = ceil(lookupScaleX .* (rayPositionsX-layerOrigin(2)));                        
                    layerPixelIndicesForRaysYDown   = ceil(lookupScaleY .* (rayPositionsY-layerOrigin(1)));                                                            
                                                            
                    % this only works if the number of effected pixels is
                    % the same everywhere
                    numLayerPixelsPerRayPixelX = 1 + layerPixelIndicesForRaysXRight(1)-layerPixelIndicesForRaysXLeft(1);
                    numLayerPixelsPerRayPixelY = 1 + layerPixelIndicesForRaysYDown(1)-layerPixelIndicesForRaysYUp(1);

                    % check if it really is the same everywhere
                    if  (sum( 1 - ((layerPixelIndicesForRaysXRight-layerPixelIndicesForRaysXLeft)==(numLayerPixelsPerRayPixelX-1)) ) ~= 0) || ...
                        (sum( 1 - ((layerPixelIndicesForRaysYDown-layerPixelIndicesForRaysYUp)==(numLayerPixelsPerRayPixelY-1)) ) ~= 0)
                        error('Dimensions are wrong!');
                    end                                        
                        
                    % matrix row indices are easy: just the current light field indices for each layer pixel
                    currentIndicesRows      = repmat(matrixRowIndicesCurrentLFAngle, [1 1 numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    
                    % construct matrix columns
                    currentIndicesColumnsX  = zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    currentIndicesColumnsY  = zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    currentIndicesColumnsZ  = layer + zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    
                    for ky=1:numLayerPixelsPerRayPixelY
                        currentIndicesColumnsY(:,:,ky,:) = repmat( (layerPixelIndicesForRaysYUp')+ky-1, [1 lightFieldResolution(4) 1 numLayerPixelsPerRayPixelX]);
                    end
                    
                    for kx=1:numLayerPixelsPerRayPixelX
                        currentIndicesColumnsX(:,:,:,kx) = repmat( layerPixelIndicesForRaysXLeft+kx-1, [lightFieldResolution(3) 1 numLayerPixelsPerRayPixelY 1]);
                    end
                                        
                    % remember all valid indices
                    validIndices = ones([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    validIndices(currentIndicesColumnsY<1)                  = 0;
                    validIndices(currentIndicesColumnsY>layerResolution(1)) = 0;
                    validIndices(currentIndicesColumnsX<1)                  = 0;
                    validIndices(currentIndicesColumnsX>layerResolution(2)) = 0;
                    
                    % set all invalid indices to 1, otherwise it crashes
                    currentIndicesColumnsX(validIndices==0) = 1;
                    currentIndicesColumnsY(validIndices==0) = 1;
                    
                    % get the matrix indices
                    currentIndicesColumns   = sub2ind(layerResolution, currentIndicesColumnsY, currentIndicesColumnsX, currentIndicesColumnsZ);
                    clear currentIndicesColumnsX currentIndicesColumnsY currentIndicesColumnsZ;
                                        
                    % matrix entries
                    currentIndicesValues    = sum(sum(validIndices,4),3);
                    currentIndicesValues(currentIndicesValues~=0) = 1 ./ currentIndicesValues(currentIndicesValues~=0);
                    currentIndicesValues    = repmat(currentIndicesValues, [1 1 numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);       
                    
                    % vectorize everything
                    validIndices            = validIndices(:);                    
                    currentIndicesRows      = currentIndicesRows(:);
                    currentIndicesColumns   = currentIndicesColumns(:);
                    
                    if bLargeScale
                        numCurrentNonzeroElements = sum(validIndices(:));
                        indexIVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = currentIndicesRows(validIndices==1);
                        indexJVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = currentIndicesColumns(validIndices==1);
                        valueSVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = currentIndicesValues(validIndices==1);
                        indexCount = indexCount + numCurrentNonzeroElements ;
                    else
                        indexIVector = cat(1, indexIVector, currentIndicesRows(validIndices==1));
                        indexJVector = cat(1, indexJVector, currentIndicesColumns(validIndices==1));
                        valueSVector = cat(1, valueSVector, currentIndicesValues(validIndices==1));
                    end                                      
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % rectangular region around ray in space and angle
                
                elseif basisFunctionType == 2  
                    
                    % get currently effected layer pixels - matrix column indices for each ligh ray                      
                    layerPixelIndicesForRaysXLeft   = ceil(lookupScaleX .* (rayPositionsXUpper-areaIntegrationSize(2)/2-layerOrigin(2)));                        
                    layerPixelIndicesForRaysYUp     = ceil(lookupScaleY .* (rayPositionsYUpper-areaIntegrationSize(1)/2-layerOrigin(1)));                    
                    layerPixelIndicesForRaysXRight  = ceil(lookupScaleX .* (rayPositionsX+areaIntegrationSize(2)/2-layerOrigin(2)));                        
                    layerPixelIndicesForRaysYDown   = ceil(lookupScaleY .* (rayPositionsY+areaIntegrationSize(1)/2-layerOrigin(1)));                                                            
                                                            
                    % this only works if the number of effected pixels is
                    % the same everywhere
                    numLayerPixelsPerRayPixelX = 1 + layerPixelIndicesForRaysXRight(1)-layerPixelIndicesForRaysXLeft(1);
                    numLayerPixelsPerRayPixelY = 1 + layerPixelIndicesForRaysYDown(1)-layerPixelIndicesForRaysYUp(1);

                    % check if it really is the same everywhere
                    if  (sum( 1 - ((layerPixelIndicesForRaysXRight-layerPixelIndicesForRaysXLeft)==(numLayerPixelsPerRayPixelX-1)) ) ~= 0) || ...
                        (sum( 1 - ((layerPixelIndicesForRaysYDown-layerPixelIndicesForRaysYUp)==(numLayerPixelsPerRayPixelY-1)) ) ~= 0)
                        error('Dimensions are wrong!');
                    end                                        
                        
                    % matrix row indices are easy: just the current light field indices for each layer pixel
                    currentIndicesRows      = repmat(matrixRowIndicesCurrentLFAngle, [1 1 numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    
                    % construct matrix columns
                    currentIndicesColumnsX  = zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    currentIndicesColumnsY  = zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    currentIndicesColumnsZ  = layer + zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    
                    for ky=1:numLayerPixelsPerRayPixelY
                        currentIndicesColumnsY(:,:,ky,:) = repmat( (layerPixelIndicesForRaysYUp')+ky-1, [1 lightFieldResolution(4) 1 numLayerPixelsPerRayPixelX]);
                    end
                    
                    for kx=1:numLayerPixelsPerRayPixelX
                        currentIndicesColumnsX(:,:,:,kx) = repmat( layerPixelIndicesForRaysXLeft+kx-1, [lightFieldResolution(3) 1 numLayerPixelsPerRayPixelY 1]);
                    end
                                        
                    % remember all valid indices
                    validIndices = ones([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    validIndices(currentIndicesColumnsY<1)                  = 0;
                    validIndices(currentIndicesColumnsY>layerResolution(1)) = 0;
                    validIndices(currentIndicesColumnsX<1)                  = 0;
                    validIndices(currentIndicesColumnsX>layerResolution(2)) = 0;
                    
                    % set all invalid indices to 1, otherwise it crashes
                    currentIndicesColumnsX(validIndices==0) = 1;
                    currentIndicesColumnsY(validIndices==0) = 1;
                    
                    % get the matrix indices
                    currentIndicesColumns   = sub2ind(layerResolution, currentIndicesColumnsY, currentIndicesColumnsX, currentIndicesColumnsZ);
                    clear currentIndicesColumnsX currentIndicesColumnsY currentIndicesColumnsZ;
                                        
                    % matrix entries
                    currentIndicesValues    = sum(sum(validIndices,4),3);
                    currentIndicesValues(currentIndicesValues~=0) = 1 ./ currentIndicesValues(currentIndicesValues~=0);
                    currentIndicesValues    = repmat(currentIndicesValues, [1 1 numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);       
                    
                    % vectorize everything
                    validIndices            = validIndices(:);                    
                    currentIndicesRows      = currentIndicesRows(:);
                    currentIndicesColumns   = currentIndicesColumns(:);
                    
                    if bLargeScale
                        numCurrentNonzeroElements = sum(validIndices(:));
                        indexIVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = currentIndicesRows(validIndices==1);
                        indexJVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = currentIndicesColumns(validIndices==1);
                        valueSVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = currentIndicesValues(validIndices==1);
                        indexCount = indexCount + numCurrentNonzeroElements ;
                    else
                        indexIVector = cat(1, indexIVector, currentIndicesRows(validIndices==1));
                        indexJVector = cat(1, indexJVector, currentIndicesColumns(validIndices==1));
                        valueSVector = cat(1, valueSVector, currentIndicesValues(validIndices==1));
                    end      
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % linearly weighted rectangular region around ray in angle
                
                elseif basisFunctionType == 3   
                
                    % get currently effected layer pixels - matrix column indices for each ligh ray                      
                    layerPixelIndicesForRaysXLeft   = ceil(lookupScaleX .* (rayPositionsXUpper-layerOrigin(2)));                        
                    layerPixelIndicesForRaysYUp     = ceil(lookupScaleY .* (rayPositionsYUpper-layerOrigin(1)));                    
                    layerPixelIndicesForRaysXRight  = ceil(lookupScaleX .* (rayPositionsX-layerOrigin(2)));                        
                    layerPixelIndicesForRaysYDown   = ceil(lookupScaleY .* (rayPositionsY-layerOrigin(1)));      
                                                            
                    % this only works if the number of effected pixels is
                    % the same everywhere
                    numLayerPixelsPerRayPixelX = 1 + layerPixelIndicesForRaysXRight(1)-layerPixelIndicesForRaysXLeft(1);
                    numLayerPixelsPerRayPixelY = 1 + layerPixelIndicesForRaysYDown(1)-layerPixelIndicesForRaysYUp(1);

                    % check if it really is the same everywhere
                    if  (sum( 1 - ((layerPixelIndicesForRaysXRight-layerPixelIndicesForRaysXLeft)==(numLayerPixelsPerRayPixelX-1)) ) ~= 0) || ...
                        (sum( 1 - ((layerPixelIndicesForRaysYDown-layerPixelIndicesForRaysYUp)==(numLayerPixelsPerRayPixelY-1)) ) ~= 0)
                        error('Dimensions are wrong!');
                    end                                        
                        
                    % matrix row indices are easy: just the current light field indices for each layer pixel
                    currentIndicesRows      = repmat(matrixRowIndicesCurrentLFAngle, [1 1 numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    
                    % construct matrix columns
                    currentIndicesColumnsX  = zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    currentIndicesColumnsY  = zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    currentIndicesColumnsZ  = layer + zeros([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    
                    for ky=1:numLayerPixelsPerRayPixelY
                        currentIndicesColumnsY(:,:,ky,:) = repmat( (layerPixelIndicesForRaysYUp')+ky-1, [1 lightFieldResolution(4) 1 numLayerPixelsPerRayPixelX]);
                    end
                    
                    for kx=1:numLayerPixelsPerRayPixelX
                        currentIndicesColumnsX(:,:,:,kx) = repmat( layerPixelIndicesForRaysXLeft+kx-1, [lightFieldResolution(3) 1 numLayerPixelsPerRayPixelY 1]);
                    end
                                        
                    % remember all valid indices
                    validIndices = ones([lightFieldResolution(3) lightFieldResolution(4) numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]);
                    validIndices(currentIndicesColumnsY<1)                  = 0;
                    validIndices(currentIndicesColumnsY>layerResolution(1)) = 0;
                    validIndices(currentIndicesColumnsX<1)                  = 0;
                    validIndices(currentIndicesColumnsX>layerResolution(2)) = 0;
                    
                    % set all invalid indices to 1, otherwise it crashes
                    currentIndicesColumnsX(validIndices==0) = 1;
                    currentIndicesColumnsY(validIndices==0) = 1;
                    
                    % get the matrix indices
                    currentIndicesColumns   = sub2ind(layerResolution, currentIndicesColumnsY, currentIndicesColumnsX, currentIndicesColumnsZ);
                    clear currentIndicesColumnsX currentIndicesColumnsY currentIndicesColumnsZ;
                                        
                    
                    % compute the normalized weights for each layer pixel in the region - only need to do that once
                    layerRegionPixelSizeX   = rayPositionsX(1)-rayPositionsXUpper(1);
                    centerX1                = lookupScaleX .* (rayPositionsXUpper(1) + layerRegionPixelSizeX/2) + 0.5;                                        
                    if layerRegionPixelSizeX > 0
                        layerPixelWeightsX  = max( layerRegionPixelSizeX/2 - abs(centerX1-(layerPixelIndicesForRaysXLeft(1):layerPixelIndicesForRaysXRight(1))), 0);                    
                        layerPixelWeightsX      = layerPixelWeightsX ./ max(layerPixelWeightsX(:));
                    else
                        layerPixelWeightsX  = 1;
                    end
                    
                    
                    layerRegionPixelSizeY   = rayPositionsY(1)-rayPositionsYUpper(1);
                    centerY1                = lookupScaleY .* (rayPositionsYUpper(1) + layerRegionPixelSizeY/2) + 0.5;      
                    if layerRegionPixelSizeY > 0
                        layerPixelWeightsY  = max( layerRegionPixelSizeY/2 - abs(centerY1-(layerPixelIndicesForRaysYUp(1):layerPixelIndicesForRaysYDown(1))), 0);
                        layerPixelWeightsY  = layerPixelWeightsY ./ max(layerPixelWeightsY(:));
                    else
                        layerPixelWeightsY  = 1;
                    end
                                        
                    % add weight in x and y
                    layerPixelWeights       = sqrt(repmat(layerPixelWeightsX, [numLayerPixelsPerRayPixelY 1]).^2 + repmat(layerPixelWeightsY', [1 numLayerPixelsPerRayPixelX]).^2);
                    layerPixelWeights       = layerPixelWeights ./ norm(layerPixelWeights(:));                                        
                    
                    currentIndicesValues = repmat( reshape(layerPixelWeights, [1 1 numLayerPixelsPerRayPixelY numLayerPixelsPerRayPixelX]), [lightFieldResolution(3) lightFieldResolution(4) 1 1]);
                                                          
                    % vectorize everything
                    validIndices            = validIndices(:);                    
                    currentIndicesRows      = currentIndicesRows(:);
                    currentIndicesColumns   = currentIndicesColumns(:);
                    currentIndicesValues    = currentIndicesValues(:);
                    
                    if bLargeScale
                        numCurrentNonzeroElements = sum(validIndices(:));
                        indexIVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = currentIndicesRows(validIndices==1);
                        indexJVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = currentIndicesColumns(validIndices==1);
                        valueSVector(indexCount:indexCount+numCurrentNonzeroElements-1 ) = currentIndicesValues(validIndices==1);
                        indexCount = indexCount + numCurrentNonzeroElements ;
                    else
                        indexIVector = cat(1, indexIVector, currentIndicesRows(validIndices==1));
                        indexJVector = cat(1, indexJVector, currentIndicesColumns(validIndices==1));
                        valueSVector = cat(1, valueSVector, currentIndicesValues(validIndices==1));
                    end
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % shift x positions of rays to next slice
                rayPositionsX = rayPositionsX - layerDistance*vx;
                rayPositionsY = rayPositionsY - layerDistance*vy;
                
                if basisFunctionType > 0
                    rayPositionsXUpper = rayPositionsXUpper - layerDistance*vxUpper;
                    rayPositionsYUpper = rayPositionsYUpper - layerDistance*vyUpper;
                end
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % draw stuff if desired
                if drawMode == 2
                    % plot for debug
                    %imagesc(T); % not enough memory for this
                    spy(T);
                    colorbar;
                    title('Propagation Matrix');
                    drawnow;
                end
                
            end
        end
    end  
         
    clear rayPositionsX rayPositionsY matrixRows matrixColumns validXXIndices validYYIndices validZZIndices validVXIndices validVYIndices layerPixelIndicesForRaysX layerPixelIndicesForRaysY lightFieldPixelCentersX lightFieldPixelCentersY;
    
    % construct the sparse matrix and add all non-zeros entries
	T = sparse(indexIVector,indexJVector,valueSVector,prod(lightFieldResolution),prod(layerResolution));
	clear indexIVector indexJVector valueSVector;
    
    if drawMode == 1
        close(h);        
    end
    
end


