% 这是readme.txt中所述的快速起始运行程序
clear all;
global T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user parameters

% number of layers
% 用户参数，定义layers的层数和最大迭代次数。
% 好像layers层数是指的用来做光场显示的显示器层数
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
% lightFieldResolution=[7,7,384,512,3]
% lightFieldSize=[75,100]，貌似是实体光场显示器的实际尺寸

 
% set layer resolution here
% 设定每一层的分辨率，这里规定与光场数据每个视角的图片分辨率相等，
% 也就是384*512的图片大小，而且要有5层。
% 注意matlab表现矩阵是先行后列，而从图像显示来说，行描述的是y轴方向，列描述的是x轴方向
layerResolution = [lightFieldResolution(3) lightFieldResolution(4) numLayers];

% layer size is size of light field
% 【尚不明确】
layerSize       = lightFieldSize; 

% distance between layers
% 每两层显示器的之间距离，第一块与最后一块之间有16.7mm，有5-1=4个间距，
% 所以每两层之间的间隔是4.1750mm
layerDistance = 16.7 / (numLayers-1);

% origin of the volume in world space [y x z] [mm] - center in layers
% 世界空间中的原点位置，单位是mm
% 【尚不明确】
% 每层原点的位置是[y=0mm,x=0mm,z=0mm],
% 深度范围是16.7mm，也就是第一块与最后一块屏之间的距离
% lightFieldOrigin(3)貌似是只光场中心中z轴原点的位置，是在几块屏幕组成的空间中央。
layerOrigin         = [0 0 0];
depthRange          = layerDistance * (numLayers-1);                
lightFieldOrigin(3) = -depthRange/2; 

% minimum transmission of each layer
% 每一层屏幕的最小透过率，1/1000，感觉是为了避免0值的出现设置的最小量。
minTransmission = 0.001;

% shift the light field origin to the center of the layers' depth range
% 【尚不明确】
% 一个逻辑变量，说明是否把光场的原点放置到深度范围之内
bLightFieldOriginInLayerCenter      = true;

                                   
% pad layer borders a little bigger
% 【尚不明确】
% 貌似是打算扩大一点点每层显示器的边界
bPadLayerBorders = true;
if bPadLayerBorders    
    % pixel size of the layers and the light field
    % 【尚不明确】
    % lightFieldSize=[75,100]，貌似是实体光场显示器的物理尺寸
    % [lightFieldResolution(3) lightFieldResolution(4)]=[384,512]
    % 好像是计算出每个lf像素的大小。
    % 是两个元素的向量，表示y轴方向的像素大小和x轴方向的像素大小
    lfPixelSize         = lightFieldSize ./ [lightFieldResolution(3) lightFieldResolution(4)];    

    % get maximum angles
    % 最大的光场倾斜视角，应该是说最边上的视点到观察点与中央轴之间的夹角。
    lfMaxAngle          = [ max(abs(lightFieldAnglesY)) max(abs(lightFieldAnglesX)) ];

    % depth range of layers 
    % 【尚不明确】
    % 深度范围=每两层之间的距离*（层数-1），不明确为何要一遍一遍定义
    % 如果光场的原点位于各层显示器之内，就把深度范围变成一半，大约是说+-的意思。
    depthRange          = layerDistance * (numLayers-1);
    if bLightFieldOriginInLayerCenter
        depthRange      = depthRange / 2;
    end

    % number of pixels to add on each side of the layers       
    % 【尚不明确】
    % 可能是为了计算出周边的边框之类
    % 边缘增加的像素数量，=最大视角*深度范围/像素宽度，最大视角大约是近似=tan，
    % 所以相当于最远或者最近处的边缘通过深度范围以后增加了多少个像素
    numAddPixels        = ceil( lfMaxAngle.*depthRange ./ lfPixelSize );   

    % size of the clamped regions
    % 边缘增加的物理宽度
    addedRegionSize     = numAddPixels .* lfPixelSize;    
    % width of the layers [y x] in mm
    % 每层显示器的实际物理宽度，[y,x]
    % [y=75+2*边框y，x=100+2*边框x]
    layerSize   = [lightFieldSize(1)+2*addedRegionSize(1) lightFieldSize(2)+2*addedRegionSize(2)];    
    % origin of the volume in world space [y x z] [mm]
    % 【尚不明确】
    layerOrigin = [lightFieldOrigin(1)-addedRegionSize(1) lightFieldOrigin(2)-addedRegionSize(2) 0];
    % adjust resolution
    % 根据边框增加后的像素数量，增加到新的像素数
    % 
    layerResolution = [lightFieldResolution(3)+2*numAddPixels(1) lightFieldResolution(4)+2*numAddPixels(2) numLayers];
end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                   
% 判断是否是彩色图像，需要有几个颜色通道。            
numColorChannels = 1;
if numel(lightFieldResolution) > 4
    numColorChannels = lightFieldResolution(5);
end
             
                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-compute sparse propagation matrix
% 预先计算“稀疏传输矩阵”                
basisFunctionType = 0;

% single channel light field resolution
% 这是单独一个颜色通道的计算
% lfResolution=[7,7,384,512]
lfResolution = [lightFieldResolution(1) lightFieldResolution(2) lightFieldResolution(3) lightFieldResolution(4)];

% compute the sparse propagation matrix - this will internally generate and populate the global variable T
% 调用了子程序
% 传入的参数有
% lightFieldAnglesY, lightFieldAnglesX, 各个视角的参数
% lightFieldSize光场显示区的物理尺寸，lfResolution[7,7,384,512], lightFieldOrigin,...
% layerResolution带边框显示器的像素分辨率, layerSize物理尺寸, layerOrigin, layerDistance, ...
% basisFunctionType=0, 1 =drawMode 1 show progress bar
                                        
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

