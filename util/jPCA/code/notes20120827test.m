%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2012
%
% test symmetric solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
numPoints = 1000;
rs = 8;6;
rand('seed',rs);
randn('state',rs);


% choose sizes and noises
n = [ 2:20 ];
epsNoise = [0.001 ; 0.01 ; 0.1; 0.4; 1];

% in this outer loop, we pick an M and generate noisy data.
% we then run the symmetric regression for testing purposes.
for i = 1 : length(n)
    for j = 1 : length(epsNoise)
        
        M = randn(n(i));
        % symmetrize it so we know roughly what the right M is
        Ms = 0.5*(M + M');
    
        % draw data from this matrix
        X = randn( numPoints , n(i) );
        dX = X*Ms + epsNoise(j)*randn(size(X));
    
        % now run the regression.
        Msym = symRegress( dX , X );
    
        % now check the least squares.
        Mun = X\dX;
 
        % the following error should be small
        symError = norm(Msym - Mun, 'fro');
        fprintf('n = %d; The sym - unconstrained error (normalized by n) should very loosely on the order of %g. It is: %0.4f.\n', n(i) , epsNoise(j)^2 , symError/n(i));
    
        % the diagonal error should be 0
        symErrorDiag = norm(diag(Msym) - diag(Mun));
        fprintf('n = %d; The sym - unconstrained DIAGONAL error (normalized by n) should on the order of 0. It is: %g.\n', n(i) , symErrorDiag/n(i));
    
        % NOTE: you should also run checkgrad always to make sure
        % the gradient is calculating correctly.  
        fprintf('n = %d; The checkgrad error (analytical vs perturbed gradient) should be on the order of 0 for any reasonable matrix. It is: %g.\n', n(i) , checkgrad('symLSeval',reshapeSym(Ms), 1e-8, dX , X ) );
 
    end    
end
