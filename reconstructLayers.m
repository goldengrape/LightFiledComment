% 这是readme.txt中所述的快速起始运行程序
clear all;
global T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user parameters

% number of layers
% 用户参数，定义layers的层数和最大迭代次数。
% 我不确定layers是指的观察角度还是说指的用来做光场显示的显示器层数。
% 好像layers是指的用来做光场显示的显示器层数
numLayers                           = 5;

% max iterations for optimization
maxIters                            = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% datapath to light field
% 数据所用的路径，用来产生光场数据的并不是从一个虚拟3D物体用光路追迹算出来的。
% 而是使用了用相机从不同方向上对同一个物体拍照获得的图像数据。
datapath = 'data/';

% load light field
% 光场的数据需要实现使用data目录下的generateLightField.m子程序来计算，
% 计算的结果会保存在LightField4D.mat文件中，文件有45M，所以我不会保存在github上。
lightFieldFilename = [datapath 'LightField4D.mat'];

% 如果光场数据文件还没有产生，就报错。
% 此时需要运行data目录下的generateLightField.m子程序来计算，
% 如果有光场数据文件，就调用之。
if ~exist(lightFieldFilename, 'file')
	error('No light field given and the generic file does not exist!');
end
    
% load everything
load(lightFieldFilename);                
    
% 在LightField4D.mat数据文件中保存着：
% lightField--这是一个5D的数组，大小是(7,7,384,512,3)，也就是7*7个视角的512*384彩色图像
% lightFieldAnglesX
% lightFieldAnglesY
% lightFieldOrigin
% lightFieldResolution=[7,7,7,384,512,3]
% lightFieldSize=[75,100]

 
% set layer resolution here
layerResolution = [lightFieldResolution(3) lightFieldResolution(4) numLayers];

% layer size is size of light field
layerSize       = lightFieldSize; 

% distance between layers
layerDistance = 16.7 / (numLayers-1);

% origin of the volume in world space [y x z] [mm] - center in layers
layerOrigin         = [0 0 0];
depthRange          = layerDistance * (numLayers-1);                
lightFieldOrigin(3) = -depthRange/2; 

% minimum transmission of each layer
minTransmission = 0.001;

% shift the light field origin to the center of the layers' depth range
bLightFieldOriginInLayerCenter      = true;

                                   
% pad layer borders a little bigger
bPadLayerBorders = true;
if bPadLayerBorders    
    % pixel size of the layers and the light field
    lfPixelSize         = lightFieldSize ./ [lightFieldResolution(3) lightFieldResolution(4)];    

    % get maximum angles 
    lfMaxAngle          = [ max(abs(lightFieldAnglesY)) max(abs(lightFieldAnglesX)) ];

    % depth range of layers
    depthRange          = layerDistance * (numLayers-1);
    if bLightFieldOriginInLayerCenter
        depthRange      = depthRange / 2;
    end

    % number of pixels to add on each side of the layers
    numAddPixels        = ceil( lfMaxAngle.*depthRange ./ lfPixelSize );   

    % size of the clamped regions
    addedRegionSize     = numAddPixels .* lfPixelSize;    
    % width of the layers [y x] in mm
    layerSize   = [lightFieldSize(1)+2*addedRegionSize(1) lightFieldSize(2)+2*addedRegionSize(2)];    
    % origin of the volume in world space [y x z] [mm]
    layerOrigin = [lightFieldOrigin(1)-addedRegionSize(1) lightFieldOrigin(2)-addedRegionSize(2) 0];
    % adjust resolution
    layerResolution = [lightFieldResolution(3)+2*numAddPixels(1) lightFieldResolution(4)+2*numAddPixels(2) numLayers];
end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                   
            
numColorChannels = 1;
if numel(lightFieldResolution) > 4
    numColorChannels = lightFieldResolution(5);
end
             
                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-compute sparse propagation matrix
                
basisFunctionType = 0;

% single channel light field resolution
lfResolution = [lightFieldResolution(1) lightFieldResolution(2) lightFieldResolution(3) lightFieldResolution(4)];

% compute the sparse propagation matrix - this will internally generate and populate the global variable T
precomputeSparsePropagationMatrixLayers3D(  lightFieldAnglesY, lightFieldAnglesX, lightFieldSize, lfResolution, lightFieldOrigin,...
                                            layerResolution, layerSize, layerOrigin, layerDistance, ...
                                            basisFunctionType, 1 );
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do it for each color channel                               
                
% monochromatic                                                      
if numColorChannels == 1                                        

    % scale light field into possible dynamic range        
	lightField(lightField<minTransmission) = minTransmission;

    disp('Reconstructing Layers from Grayscale Light Field');

    % reconstruct layers  
    tic;
    [layersRec lightFieldRec] = ...
    computeAttenuationLayersFromLightField4D(lightField, lightFieldResolution, layerResolution, maxIters, minTransmission);

    tt = toc;
    disp(['  Reconstruction took ' num2str(tt) ' secs']);

% RGB
else                    

    % reconstruct each color channel
    for c=1:numColorChannels

        disp(['Processing light field color channel ' num2str(c) ' of ' num2str(numColorChannels)]);                                               

        % load RGB light field
        load(lightFieldFilename, 'lightField');
        lightField(lightField<minTransmission) = minTransmission;
        
        % delete currently unused color channels 
        lightField = lightField(:,:,:,:,c);                                                     

        % reconstruct layers for current channel 
        tic;
        [layersRecTmp lightFieldRecTmp] = ...
         computeAttenuationLayersFromLightField4D(lightField, lfResolution, layerResolution, maxIters, minTransmission);

        tt = toc;
        disp(['  Reconstruction took ' num2str(tt) ' secs']);


        % save to temp file
        save([datapath 'RecTemp_' num2str(c) '.mat'], 'layersRecTmp', 'lightFieldRecTmp');
        clear layersRecTmp lightFieldRecTmp;

    end


    % initialize reconstructed layers and light field
    layersRec       = zeros([layerResolution(1) layerResolution(2) layerResolution(3) numColorChannels]);
    lightFieldRec   = zeros(lightFieldResolution);                    

    % load each color channel
    for c=1:numColorChannels       
        clear lightField;
        % load temp file
        tmpFilename = [datapath 'RecTemp_' num2str(c) '.mat'];
        load(tmpFilename);
        % set color channel in reconstruction
        layersRec(:,:,:,c)          = layersRecTmp;
        lightFieldRec(:,:,:,:,c)    = lightFieldRecTmp;
        % delete temp file
        if isunix
            system(['rm -f ' tmpFilename]);
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove added border pixels                                  

if bPadLayerBorders    
    layersRec           = layersRec(numAddPixels(1)+1:end-numAddPixels(1), numAddPixels(2)+1:end-numAddPixels(2), :, :); 
    layerResolution(1)  = layerResolution(1) - 2*numAddPixels(1);
    layerResolution(2)  = layerResolution(2) - 2*numAddPixels(2);
end

% filename for reconstruction
filename = [datapath 'Reconstruction3D_' num2str(numLayers) 'layers_dist' num2str(layerDistance) '.mat'];

% save data
disp(['Done. Saving data as ' filename]);
if ~exist('lightField', 'var')
	load(lightFieldFilename);
end
save(filename, 'layersRec', 'lightFieldRec', 'lightFieldAnglesY', 'lightFieldAnglesX', 'lightFieldSize', 'lightFieldOrigin', 'minTransmission', 'lightField');                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data

numSubplots = [2 3];
                
% original light field
subplot(numSubplots(1), numSubplots(2), 1);
img = drawLightField4D(lightField);
imshow(img);
title('Original Light Field');
                
% reconstruction                
subplot(numSubplots(1), numSubplots(2), 2);
img = drawLightField4D(lightFieldRec); 
imshow(img);
title('Reconstructed Light Field');                

% central view
subplot(numSubplots(1), numSubplots(2), 3);
IRecCentral = reshape(lightFieldRec( ceil(size(lightFieldRec,1)/2),ceil(size(lightFieldRec,2)/2),:,:,:), [layerResolution(1) layerResolution(2) size(layersRec,4)]);
imshow( IRecCentral );
title('Central View in Full Resolution');

% attenuation layers
subplot(numSubplots(1), numSubplots(2), 4:6);
img = drawAttenuationLayers3D(layersRec);
imshow(img);
title('Attenuation Layers');
                
colormap gray;
drawnow;

