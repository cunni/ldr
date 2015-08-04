%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% orth_procrustes.m
%
% function for evaluating orthogonal procrustes
% problem of the form 1/2||A*Q - B||^2.
% This evaluates the objective and gradient
% for testing various stiefel manifold solvers.
%
%%%%%%%%%%%%%%%%%%%%

function [ f , Df ] = orth_procrustes( Q , A , B );

    % the function
    f = 0.5*norm((A*Q - B),'fro')^2;
    
    % the grad in Q
    Df = A'*A*Q - A'*B;
    
    