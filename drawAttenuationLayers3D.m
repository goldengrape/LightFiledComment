%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Draw reconstructed attenuation layers
%
%   Gordon Wetzstein [wetzste1@cs.ubc.ca]
%   PSM Lab | University of British Columbia
%   February 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function image = drawAttenuationLayers3D(layers)
    
    % create image in full resolution
    bFullResolution = true;

    layerResolution = size(layers);

    if numel(layerResolution) < 3
        layerResolution = [layerResolution 1 1];        
    end
    if numel(layerResolution) < 4
        layerResolution = [layerResolution 1];
    end
    
    % number of pixel between each subplot if an image is assembled
    pixelSpacing = 0;
    
    if nargout > 0                
        % number of super-pixels in the image
        if ~bFullResolution
            subplotResolution = floor([ layerResolution(1)/layerResolution(3) layerResolution(2)/layerResolution(3) ]);                        
        else
            subplotResolution = [ layerResolution(1) layerResolution(2) ];                        
        end
        % initialize image
        image = zeros( [    subplotResolution(1) + 2*pixelSpacing... 
                            subplotResolution(2)*layerResolution(3) + (layerResolution(3)+1)*pixelSpacing ...
                            layerResolution(4)] );                        
    end
    
    for k=1:layerResolution(3)
                
        currenLayer = reshape(layers(:,:,k,:), [layerResolution(1) layerResolution(2) layerResolution(4)]);
        
        if nargout == 0
            subplot(1,layerResolution(3),k);
            imagesc(currenLayer, [0 1]);
            set(gca, 'XTick', [], 'YTick', []);
            axis equal;
        else      
            idxSpacing = (k-1)*(subplotResolution(2)+pixelSpacing) + 1 + pixelSpacing;
            image(pixelSpacing+1:end-pixelSpacing,idxSpacing:idxSpacing+subplotResolution(2)-1,:) = imresize(currenLayer, subplotResolution, 'bilinear'); 
        end
    end
end
