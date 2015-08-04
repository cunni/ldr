%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% get_results.m
%
% This code is a simple helper function to ensure all results get structed
% in the same way for comparability
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ results ] = get_results( method , fQ_traditional , fQ_stiefel , Q_traditional , Q_stiefel , d , r )

    % objective
    results.name = method;
    results.f.traditional = fQ_traditional;
    results.f.stiefel = fQ_stiefel;
    results.Q.traditional = Q_traditional;
    results.Q.stiefel = Q_stiefel;
    results.d = d;
    results.r = r;
    % the singular values show how different the subspaces are
    results.sv = svd([ Q_traditional , Q_stiefel ]);
    results.subspace = subspace( Q_traditional , Q_stiefel );
    % CCA is a normed objective (a correlation), so in that case we just
    % report the values
    if isequal(method,'CCA') || isequal(method,'MAF')
        results.normed_diff = (fQ_stiefel - fQ_traditional);
    else
        results.normed_diff = (fQ_stiefel - fQ_traditional)/abs(fQ_traditional);
    end
    % print message
    fprintf('Stiefel found a minimum %2.4f smaller than traditional %s.\n', (fQ_traditional - fQ_stiefel) , method );
