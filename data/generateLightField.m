%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Load a series of pre-rendered images that make up different views of a
%   light field and combine them all into the data structures expected by
%   the tomographic solver.
%
%   Gordon Wetzstein [wetzste1@cs.ubc.ca]
%   PSM Lab | University of British Columbia
%   February 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
addpath('../');

% plot the data
bPlotData = true;
if bPlotData
    clf;
end

% path to output light field and folder containing the scene
%scene = 'dice';
scene = 'dragon';
%scene = 'butterfly';
%scene = 'messerschmitt';

% scale the output light field between this value and what was the maximum
% value
minLightFieldValue = 0.01;

% size of the light field [y x] [mm]
lightFieldSize      = [75 100];

% resolution in [#anglesY @anglesX #pixelsY #pixelsX #colorchannels]
lightFieldResolution = [7 7 384 512 3];

% use alpha channel - if this is true BG must be a background image 
bUseAlpha           = true;
gradientStartEnd    = [155 191];
gradientStartEnd    = [155-20 191+20];
BG                  = repmat((gradientStartEnd(1):(gradientStartEnd(2)-gradientStartEnd(1))/(lightFieldResolution(3)-1):gradientStartEnd(2))' ./ 255, [1 lightFieldResolution(4) lightFieldResolution(5)]);
BG                  = flipdim(BG,1);

% path to images
datapath = [scene '/' num2str(lightFieldResolution(1)) 'x' num2str(lightFieldResolution(2)) 'x' num2str(lightFieldResolution(3)) 'x' num2str(lightFieldResolution(4)) '/'];

% field of view in degrees
fov = 10;

% max angle in both ways in v = tan(theta)
maxAngle            = tan(pi*fov/360);
lightFieldAnglesX   = -maxAngle:2*maxAngle/(lightFieldResolution(1)-1):maxAngle;
lightFieldAnglesY   = -maxAngle:2*maxAngle/(lightFieldResolution(2)-1):maxAngle;
                    
% origin of light field in world space [y x z]
lightFieldOrigin    = [0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% basename for files with count < 10
filename = [datapath scene '-'];
if lightFieldResolution(1)*lightFieldResolution(2)>10
   filename = [filename '0'];
end
if lightFieldResolution(1)*lightFieldResolution(2)>100
   filename = [filename '0'];
end
filetype = '.png';

% initialize light field
lightField = zeros(lightFieldResolution);

count = 1;
for ky=1:lightFieldResolution(1)
    for kx=1:lightFieldResolution(2)
        
        currentFilename = [filename num2str(count) filetype];
        
        % load image
        if count > 9
            currentFilename = [filename(1:end-1) num2str(count) filetype];            
        end
        
        % load image
        if ~bUseAlpha
            I = im2double(imread(currentFilename));
        % load image and alpha
        else
            [I map alpha] = imread(currentFilename);
            I       = im2double(I);
            alpha   = im2double(alpha);
        end
        
        % convert to grayscale if desired
        if lightFieldResolution(5) == 1
            I = rgb2gray(I);
        end
        
        % change background using alpha
        if bUseAlpha
            I = (repmat(alpha, [1 1 size(I,3)])) .* I + (1-repmat(alpha, [1 1 size(I,3)])) .* BG;
        end
        
        lightField(ky,kx,:,:,:) = reshape(I, [1 1 lightFieldResolution(3) lightFieldResolution(4) lightFieldResolution(5)]);
        
        if bPlotData
            imshow(I); drawnow;
        end
                
        % increment counter
        count = count + 1;
    end
end

% scale light field
minMaxLF    = [ min(lightField(:)) max(lightField(:)) ];
if min(lightField(:)) < minLightFieldValue
    lightField  = ((lightField-minMaxLF(1)) ./ (minMaxLF(2)-minMaxLF(1))) .* (minMaxLF(2)-minLightFieldValue) + minLightFieldValue;
end

if bPlotData
    img = drawLightField4D(lightField);
    imshow(img);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(['LightField4D.mat'], 'lightField', 'lightFieldAnglesY', 'lightFieldAnglesX', 'lightFieldSize', 'lightFieldResolution', 'lightFieldOrigin');