%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% minimize_stiefel.m
%
% This code performs minimization of a function f
% over the Stiefel manifold.
%
% It uses Algorithm 15 of Manton's 2002 "Optimization
% Algorithms exploiting unitary constraints".
%
% Inputs:
% - f: a function mapping {point on the stiefel
% manifold} --> reals
% - gradf: the Euclidean gradient of f at Q.
% 
% Outputs:
% - Q: the optimal Stiefel point for f
% - fQ: f evaluated at Q
%
% This function is constrained to real valued Q, but it
% can be extended to complex (see Manton 2002).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Q , fQ , info ] = minimize_stiefel_sd( f , Q , length , varargin )

    % set convergence tolerance
    conv_tol = 1e-8;
    converged = 0;
    
    % check inputs
    if isempty(length)
        length = 1000;
    end
    % understand the problem size
    [ d , r ] = size(Q);
    % should perform some error checking here to verify that Q is a stiefel point    
    if norm( Q'*Q - eye(r) , 'fro' ) > 1e-4
        % oops
        fprintf('WARNING: initial point Q is not a Stiefel point. Correcting...\n');
        Q = project_stiefel( Q );
    end

    % evaluate beginning fQ
    [fQ_0 , ~] = feval(f, Q, varargin{:});

    i = 0;
    % main iterative loop
    while ~converged
        % start timer
        t0 = tic;
        % evaluate the function and its derivative
        [fQ , gradfQ] = feval(f, Q, varargin{:});
        % compute the descent direction
        Z = -1*( gradfQ - Q*gradfQ'*Q );
        % compute canonical inner product
        stiefip = trace(Z'*(eye(d) - 0.5*Q*Q')*Z);
        % check for convergence
        if stiefip < conv_tol || i >= length
            % then stop
            converged = 1;
        else
            % now do the step size search over the manifold (Armijo rule)
            % note that this linesearch is in Euclidean flows, not geodesic
            % flows.  Namely, it takes additive steps and projects back to
            % the manifold, rather than taking steps in the manifold as in
            % Abrudan 2008 IEEE Sig Proc (which is only for unitary). 
            gam = 1;
            % first do the increasing gamma search
            while ( fQ - feval(f, project_stiefel(Q + 2*gam*Z) , varargin{:}) ) >= gam*stiefip
                gam = 2*gam;
            end
            % now do the decreasing gamma search
            while ( fQ - feval(f, project_stiefel(Q + gam*Z) , varargin{:}) ) < 0.5*gam*stiefip
                gam = 0.5*gam;
                if gam < 1e-8
                    % then we say the method has converged
                    gam = 0;
                    converged = 1;
                    break
                end
            end
            % now gam is correct, update
            Q = project_stiefel( Q + gam*Z );
            % iterate
            if mod(i,50)==0
                %fprintf('Iteration %d: step size= %0.2f , objective = %3.4f , lognorm descent = %3.2f.\r',i,log(gam), fQ, log(stiefip));
            end
            % fill info struct... to align with manopt, time is
            % cumulative...
            if i > 0
                info(i+1).time = toc(t0) + info(i).time;
            else
                info(i+1).time = toc(t0);
            end
            info(i+1).iter = i;
            info(i+1).cost = fQ;
            info(i+1).gradnorm = norm(gradfQ);
            % iterate
            i = i + 1;
            
        end        
    end
    % now the stepped Q is correct and final.  return Q and fQ
    t0 = tic;
    [fQ,gradfQ] = feval( f , Q , varargin{:} );
    % fill info struct
    if i > 0
        info(i+1).time = toc(t0) + info(i).time;
    else
        info(i+1).time = toc(t0);
    end
    info(i+1).iter = i;
    info(i+1).cost = fQ;
    info(i+1).gradnorm = norm(gradfQ);

    fprintf('(Stiefel manifold steepest descent completed on %s in %d iterations, objective optimized from %g to %g)\n',f,i, fQ_0, fQ);
    %
        
end    