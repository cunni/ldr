%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% test affine solver
%
% This script gives some options for playing around with the error
% it generates reasonably known dynamics (known plus noise), and then
% it fits these dynamics with our model and calculates errors.  It also
% uses checkgrad to calculate the gradient and objective calculation.
% Different noise levels, different data sizes, and so on are all used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
numPoints = 1000;
rs = 122;8;6;
rand('seed',rs);
randn('state',rs);


% choose sizes and noises
n = [ 2:20 ];
epsNoise = [0.001 ; 0.01 ; 0.1; 0.4; 1];

% in this outer loop, we pick an M and generate noisy data.
% we then run the regression for testing purposes.
for i = 1 : length(n)
    for j = 1 : length(epsNoise)
        
        M = randn(n(i));
        % skew-symmetrize it so we know roughly what the right M is
        Mk = 0.5*(M - M');
        % make a random offset
        y = randn(1,n(i));
    
        % draw data from this matrix
        X = randn( numPoints , n(i) );
        dX = X*Mk + epsNoise(j)*randn(size(X)) + repmat(y,size(X,1),1);
    
        % now run the regression.
        [MskewAff, yAff] = skewSymAffineRegress( dX , X );
    
        % now check the least squares.  Note that here we are subtracting
        % off the known y, which we would not have access to in practice.
        % But this makes for a good test of the regression...
        Mun = X\(dX- repmat(y,size(dX,1),1));
 
        % the following error should be small
        symError = norm(MskewAff - Mun, 'fro');
        fprintf('n = %d; The skew - unconstrained error (normalized by n) should very loosely on the order of %g. It is: %0.4f.\n', n(i) , epsNoise(j)^2 , symError/n(i));
        
        % the following error should be small
        offsetError = norm(y - yAff);
        fprintf('n = %d; The offset error should on the order of 0. It is: %0.4f.\n', n(i), offsetError);
        
        % NOTE: you should also run checkgrad always to make sure
        % the gradient is calculating correctly.  
        fprintf('n = %d; The checkgrad error (analytical vs perturbed gradient) should be on the order of 0 for any reasonable matrix. It is: %g.\n', n(i) , checkgrad('skewSymAffineLSeval',[reshapeSkew(Mk);y'] + 0*randn, 1e-8, dX , X ) );
 
        
    end    
end
