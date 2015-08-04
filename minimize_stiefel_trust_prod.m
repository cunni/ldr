%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% minimize_stiefel_trust_prod.m
%
% This code performs minimization of a function f
% over the product of two Stiefel manifolds.  For CCA only.
%
% It uses the excellent manopt package www.manopt.org .
%
% Inputs:
% - f: a function mapping {point on the stiefel
% manifold} --> reals
% - This function should also produce gradf: the Euclidean gradient of f at Q.
% 
% Outputs:
% - Q: the optimal Stiefel point for f
% - fQ: f evaluated at Q
%
% This function is constrained to real valued Q, but it
% can be extended to complex.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Qx , Qy , fQ , info ] = minimize_stiefel_trust_prod( f , Qx , Qy , length , Sxx , Syy , Sxy )

    % set convergence tolerance
    conv_tol = 1e-8;
    converged = 0;
    
    % check inputs
    if isempty(length)
        length = 200;
    end
    % understand the problem size
    [ dx , r ] = size(Qx);
    [ dy , r ] = size(Qy);

    % evaluate beginning fQ
    [fQ_0 , ~] = f_cca( Qx , Qy , [] , [] , Sxx , Syy , Sxy );
    
    % create manopt problem structure.
    Mstruct.A = stiefelfactory(dx,r);
    Mstruct.B = stiefelfactory(dy,r);
    problem.M = productmanifold(Mstruct);

    %[problem.cost , problem.egrad]  = @(M) feval(f , M, varargin{:});
    problem.cost = @(M) get_output( 'f_cca' , 1 , M.A , M.B , [] , [] , Sxx , Syy , Sxy ); 
    problem.egrad = @(M) struct('A', get_output( 'f_cca' , 2 , M.A , M.B , [] , [] , Sxx , Syy , Sxy ) , 'B' , get_output( 'f_cca' , 2 , M.B , M.A , [] , [] , Syy , Sxx , Sxy' ));
    %problem.egrad = @(M) feval(f , M, varargin{:});
    
    % check the gradient consistency.
    % checkgradient(problem);
    
    % Issue a call to a solver. A random initial guess will be chosen and
    % default options are selected except for the ones we specify here.
    %options.Delta_bar = 8*sqrt(p);
    options.verbosity = 0;
    options.tolCost = conv_tol;
    warning('off', 'manopt:getHessian:approx');
    [Q, fQ , info, options] = trustregions(problem , struct('A',Qx,'B',Qy) , options );
    Qx = Q.A;
    Qy = Q.B;

    fprintf('(Stiefel manifold trust region opt completed on %s in %d iterations, objective optimized from %g to %g)\n',f,info(end).iter, fQ_0, fQ);
    %
end    