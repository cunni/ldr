%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% f_pca.m
%
% function for a simple eigenvalue problem
% This evaluates the objective and gradient
% for testing basic pca.
%
%%%%%%%%%%%%%%%%%%%%

function [ f , gradf ] = f_pca( Q , X )

    % if needing to check gradients...
    %Q = reshape(Q,5,2);

    % the function
    %f = 0.5*norm( X - Q*Q'*X , 'fro' )^2;
    f = 0.5*trace( ( X - Q*Q'*X )'*( X - Q*Q'*X ) );
    
    % the grad in Q
    gradf = - 2*X*X'*Q;
    % IMPORTANT NOTE: Because we are here treating this objective in an
    % unconstrained fashion, we must take the unconstrained gradient of the
    % objective function.  This means that there is a quartic term in Q,
    % since Q'*Q is not necessarily I in the derivative. 
    gradf = gradf + Q*Q'*X*X'*Q + X*X'*Q*Q'*Q;
    
    % if needing to check gradients...
    %gradf = gradf(:);