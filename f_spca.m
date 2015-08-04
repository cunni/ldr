%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% sparse_pca.m
%
% approximates the sparse PCA objective
% with the standard LASSO type L1 relaxation:
% minimize 1/2*trace(Q'*A*Q) + alpha*sum(sum(abs(Q))) 
% subject to Q'*Q = I (Q a point on the Stiefel(n,p) manifold
%
% This is written as a proof of concept for the generic 
% linear dimensionality reduction solver in Cunningham
% and Ghahramani 2013.
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Q , fQ ] = sparse_pca( A , p , lambda)

    % first check for symmetry of A 
    if norm(A-A','fro') > 1e-8
        fprintf('warning: expecting a symmetric A... correcting...\n');
        A = A'*A;
    end
    % check alpha... this regulates sparsity
    if nargin < 3 || isempty(lambda)
        lambda = 0.1;
    end
    % specify dimensionality
    if nargin < 2 || isempty(p)
        p = 2;
    end
    % note that there are some silly values of lambda.  If lambda > pn,
    % then there is no constraint, as sum(sum(abs(Q)) < pn by def.  
    % In fact, the largest sum is n/sqrt(n) for any single column (to make
    % a column with squared norm, so p*n/sqrt(n) is the minimum meaningful
    % constraint.  Wait, this definition is not meaningful for the lagrange
    % form of this problem, just for the linear inequality form... IGNORE

    % use svd to initialize.
    [U,S,V] = svd(A + 1*randn(size(A)),0);
    % now Q is a Stiefel point...
    Q = U(:,1:p);
    
    % should perform some error checking here to verify that Q is a stiefel point
    if norm( Q'*Q - eye(p) , 'fro' ) > 1e-4
        % oops
        fprintf('warning: initial point Q is not a Stiefel point. Correcting...\n');
        Q = project_stiefel( Q );
    end

    % now we solve a sparsity problem.
    for i = 0:8
        % alpha is the smoothness term of our approximation to the abs() fn
        alpha = 2^i;
        fprintf('------------Run %d of barrier method-----------\n',i);
        % now run the stiefel...
        [ Q , fQ ] = minimize_stiefel_sd( 'sparse_pca_L1CT' , Q , 1000 , A , lambda , alpha );
        % now the Q is the next central path iterate
        fprintf('-----------------------------------------------\n');
        fprintf('Central path objective: %3.6g.\n',fQ);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




