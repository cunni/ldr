%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2011
%
% symRegress
%
% This function does least squares regression between a matrix dX and X.
% It finds a matrix M such that dX = XM (as close as possible in the 
% least squares sense).  Unlike regular least squares regression (M = X\dX)
% this function find M such that M is a symmetric matrix, that is, 
% M = M^T.
%
% Put another way, this is least squares regression over the constrained set
% of symmetric matrices.  
%
% NOTE: Compare to skewSymRegress()
%
% This can be solved by treating M as a vector of the unique elements that
% exist in a symmetric matrix. A symmetric matrix M of size n by n really only
% has n*(n+1)/2 unique entries.  That is, the diagonal and the
% upper/lower triangle is the transpose of the lower/upper.  So,
% we can just think of such a matrix as a vector x of size n(n+1)/2.
%
% Corresponding to this change in M, we would have to change X to be quite
% big and quite redundant.  That can be done, but an easier and 
% faster and more stable thing to do is to use an iterative solver
% that takes a function and gradient evaluation. 
% This iterative procedure is numerically accurate, etc, etc.
% So, this allows us to never make a big X matrix, and we just have to
% reshape M between vector and symmetric matrix form as appropriate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = symRegress(dX,X)
    
    % check that dX and X are appropriately sized.
    if ~isequal(size(dX),size(X))
        fprintf('ERROR: dX and X are not matched in size... a symmetric matrix can not result.\n');
        keyboard
    end
    if size(dX,2) > 20 || size(X,2) > 20 || size(dX,1) < 20 || size(X,1) < 20
        fprintf('ERROR (likely): dX and X should be ct by k, where ct is big and k is small, like 6.  Currently your ct is small or k is big.\n');
        keyboard
    end
    if size(dX,1) < size(dX,2) || size(X,1) < size(X,2)
        fprintf('ERROR (almost definitely): dX and/or X are fat, not skinny.  This will result in a subrank solution.\n');
        keyboard
    end
   
    % initialize m0 somehow. It does not matter, as this function is
    % convex... ie there is provable one unique minimizer.  But a good
    % initialization will help the optimization converge faster in theory
    % (though in practice this is so fast anyway that it is unnoticed).
    M0 = X\dX; % the non-symmetric matrix.
    M0k = 0.5*(M0 + M0'); % the symmetric matrix component... 
    m0 = reshapeSym(M0k);
    m0 = zeros(size(m0));
    % m0 = 100*randn(size(m0)); 
    % any of the above are fine.
    
    %%%%%%%
    % the following call does all the work
    %%%%%%%
    
    % just call minimize.m with the appropriate function...
    [m, fM, i] = minimize( m0 , 'symLSeval' , 1000, dX , X );
    
    
    % check to make sure that nothing was funky.
    if i > 500
        fprintf('Warning: more than 500 line searches were required for minimize to complete.  It should complete much more quickly.  Check the code or the conditioning of the matrices.\n');
        keyboard
    end
    
    % return the matrix
    M = reshapeSym(m);
    
end
