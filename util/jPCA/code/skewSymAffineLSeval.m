%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% skewSymAfifneLSeval.m
%
% This is the same as skewSymAffineLSeval but it adds an affine term for
% dX = XM + repmat(y) (that is, dx(t) = x(t)*M + y ). The incremental
% changes are entirely to pad the solution with that extra affine term y.
% see 2013 notes p62 for derivation of this incremental bit.
% 
% This function evaluates the least squares function and derivative.
% The derivative is taken with respect to the vector m, which is the vector
% of the unique values in the skew-symmetric matrix M.
%
% A typical least squares regression of Az = b is z = A\b.  
% That is equivalent to an optimization: minimize norm(Az-b)^2.
%
% In matrix form, if we are doing a least squares regression of AM = B,
% that becomes: minimize norm(AM-B,'fro')^2, where 'fro' means frobenius
% norm.  Put another way, define error E = AM-B, then norm(E,'fro') is the
% same as norm(E(:)).
%
% Here, we want to evaluate the objective AM-B, where we call A 'X' and B
% 'dX'.  That is, we want to solve: minimize norm(dX - XM,'fro')^2.
% However, we want to constrain our solutions to just those M that are
% skew-symmetric, namely M = -M^T.  
%
% So, instead of just using M = X\dX, we here evaluate the objective and
% the derivative with respect to the unique elements of M (the vector m),
% for use in an iterative minimizer.
%
% See notes p80-82 for the derivation of the derivative.  
%%%%%%%%%%%%%%%%%%%%%%%%

function [f, df] = skewSymAffineLSeval(z , dX , X)

% since this function is internal, we do very little error checking.  Also
% the helper functions and internal functions should throw errors if any of
% these shapes are wrong.

%%%%%%%%%%%%%
% partition the variable
%%%%%%%%%%%%%
m = z(1 : end-size(X,2));
y = z(end-size(X,2)+1:end)';
   

%%%%%%%%%%%%%
% Evaluate objective function
%%%%%%%%%%%%%
f = norm( dX - X*reshapeSkew(m) - repmat(y,size(dX,1),1) , 'fro')^2;

%%%%%%%%%%%%%
% Evaluate derivative
%%%%%%%%%%%%%
% the projected gradient for M
D = (dX - X*reshapeSkew(m) - repmat(y,size(dX,1),1) )'*X;
dfm = 2*reshapeSkew( D - D' );
% the gradient (unconstrained) for y
dfy = 2*( size(dX,1)*y + sum(X*reshapeSkew(m)) - sum(dX) );
% the total gradient
df = [dfm; dfy'];
