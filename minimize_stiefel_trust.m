%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% minimize_stiefel_trust.m
%
% This code performs minimization of a function f
% over the Stiefel manifold.
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

function [ Q , fQ , info ] = minimize_stiefel_trust( f , Q , length , varargin )

    % set convergence tolerance
    conv_tol = 1e-8;
    converged = 0;
    
    % check inputs
    if isempty(length)
        length = 200;
    end
    % understand the problem size
    [ d , r ] = size(Q);
    % should perform some error checking here to verify that Q is a stiefel or Grassman point    
    if norm( Q'*Q - eye(r) , 'fro' ) > 1e-4
        % oops
        fprintf('WARNING: initial point Q is not a Stiefel or Grassman point. Correcting...\n');
        Q = project_stiefel( Q );
    end

    % evaluate beginning fQ
    [fQ_0 , ~] = feval(f, Q, varargin{:});
    
    % create manopt problem structure.
    Manif = stiefelfactory(d,r);
    problem.M = Manif;

    %[problem.cost , problem.egrad]  = @(M) feval(f , M, varargin{:});
    problem.cost = @(M) get_output( f , 1 , M , varargin{:} ); 
    problem.egrad = @(M) get_output( f , 2 , M , varargin{:} );
    %problem.egrad = @(M) feval(f , M, varargin{:});
    
    % check the gradient consistency.
    % checkgradient(problem);
    
    % Issue a call to a solver. A random initial guess will be chosen and
    % default options are selected except for the ones we specify here.
    %options.Delta_bar = 8*sqrt(p);
    options.verbosity = 0;
    options.tolCost = conv_tol;
    warning('off', 'manopt:getHessian:approx');
    [Q, fQ , info, options] = trustregions(problem , Q , options );
        
    fprintf('(Stiefel manifold trust region opt completed on %s in %d iterations, objective optimized from %g to %g)\n',f,info(end).iter, fQ_0, fQ);
    %
end    