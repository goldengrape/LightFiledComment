%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Given a pre-computed projection matrix and a light field, solve the
%   constrained linear equation system with lsqlin.
%
%   Gordon Wetzstein [wetzste1@cs.ubc.ca]
%   PSM Lab | University of British Columbia
%   February 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [layersRec lightFieldRec] = ...
    computeAttenuationLayersFromLightField4D( lightField, lightFieldResolution, layerResolution, maxIters, minTransmission )         
                                            
	% declare matrix a global variable to save memory and time!
    global T;
                                                        
    % check if system is under-determined
    if prod(layerResolution) > prod(lightFieldResolution)
        error('Optimization is under-determined, this is currently not supported!');
    end  
              
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make sure there are no zeros in the light field (will be inf in log
    % light field)
    lightField(lightField<eps) = eps;
            
    % light field to log light field
    logLightField = log(lightField);       

    % here, the optimization needs to be constrained to transmission values
    % between minTransmission and 1, so the log values need to be < 0
    lb      = zeros([prod(layerResolution) 1]) + log(minTransmission);
    ub      = zeros([prod(layerResolution) 1]);

    % initial guess - use either zeros or the backprojected solution (zeros
    % seem to work much better)
    x0 = zeros(layerResolution)-log(eps);    
    
    % check if everything is finite
    if sum(~isfinite(logLightField(:))) > 0
        error('Log Light Field is not Finite!');
    end        

    % Sparse matrix for preconditioning.
    % This matrix is a required input for the solver;
    % preconditioning is not really being used in this example
    Jinfo = speye(prod(lightFieldResolution), prod(layerResolution));
           
    % use sparse matrix instead of function handles
    fun = @(Jinfo,Y,flag) projectLayersMatrix( Y, flag );

    % function handle
    Afun = @(x,flag) projectLayersMatrix( x, flag );
        
	% set optmimization options
    options = optimset( 'Display', 'final', 'MaxIter', maxIters, 'JacobMult', fun);                      
                 
    % run lsqlin
    [layersRec,resnorm,residual,exitflag,output] = ...
        lsqlin(Jinfo,logLightField(:),[],[],[],[],lb,ub,x0(:),options);  

    if output.iterations < maxIters
        disp(['  stopped after ' num2str(output.iterations) ' iters, of ' num2str(maxIters) ' max']);
    end
                        
    % reshape result
    lightFieldRec   = fun(Jinfo, layersRec, 1);
    lightFieldRec   = reshape(lightFieldRec, lightFieldResolution);
    layersRec       = reshape(layersRec, layerResolution);  
    
    lightFieldRec   = exp(lightFieldRec); 
    layersRec       = exp(layersRec);                       	                
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use explicit matrix (sparse)

function result = projectLayersMatrix( Y, flag ) 

    % declare T a global variable
    global T;

    % forward projection
    if flag > 0                  
        result = T*Y;                                       
    % back-projection
    elseif flag < 0 
          result = T'*Y;                                                                                                         
    % forward and then back-projection
    else
      result = T'*(T*Y);        
    end    
end

