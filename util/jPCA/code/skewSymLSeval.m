%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2010
%
% skewSymLSeval.m
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

function [f, df] = skewSymLSeval(m , dX , X)

% since this function is internal, we do very little error checking.  Also
% the helper functions and internal functions should throw errors if any of
% these shapes are wrong.

%%%%%%%%%%%%%
% Evaluate objective function
%%%%%%%%%%%%%

f = norm( dX - X*reshapeSkew(m) , 'fro')^2;

%%%%%%%%%%%%%
% Evaluate derivative
%%%%%%%%%%%%%
D = (dX - X*reshapeSkew(m))'*X;

df = 2*reshapeSkew( D - D' );