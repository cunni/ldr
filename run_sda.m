%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% run_sda.m
%
% This code runs SDA by the chosen method.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ Q , fQ , D ] = run_sda( Xdot , X , r , method )

    %%%%%%%%%
    % input checking
    %%%%%%%%%
    if nargin < 4 || isempty(method)
        % default stiefel
        method = 'stiefel';
    end
    if nargin < 3 || isempty(r)
        % default the standard choice
        r = 2;
    end
    if mod(r,2)~=0
        % must be even
        fprintf('skew dynamics return only planes... r must be even');
        keyboard;
    end
    if nargin < 2 || isempty(Xdot) || isempty(X)
        % unacceptable
        fprintf('ERROR: Bad X and/or Xdot inputs to run_lda.')
        keyboard
    end
    
    
    %%%%%%%%%
    % SDA in the specified manner
    %%%%%%%%%
    switch(lower(method))
        
        case 'heuristic'
            % SDA in the 'jPCA' manner
            
            %%%%%%%%%%%%%%%
            % find the top plane of the best fitting dynamical system using jPCA
            %%%%%%%%%%%%%%%
            % fprintf('Running jPCA...\n');
            path(path,'util/jPCA/code');
            [ Mskewsym , error_Mskew ] = skewSymRegress( Xdot' , X' );
            % recall that skewsymRegress solves dX = XM, so we give it
            % Xdot' = X' M , so the returned M is actually M' in the
            % formulation at this level, since here we care about 
            % Xdot - MX
            Msk = Mskewsym';
            [U,V] = eig(Msk);
            % take the top r dimensions
            [ ~ , indmax ] = sort(abs(diag(V)),'descend');
            indmax = indmax(1:r);
            % make real bases out of complex conjugates
            % a nice place to see this spelled out rigorously is gamal
            % elsayed normal system paper
            for rr = 1 : r/2 
                Q( : , 2*rr - 1 ) = (1/sqrt(2))*sqrt(-1)*( U(:,indmax(2*rr-1)) - U(:,indmax(2*rr)) );
                Q( : , 2*rr ) = (1/sqrt(2))*( U(:,indmax(2*rr-1)) + U(:,indmax(2*rr)) );
                D(2*rr - 1:2*rr, 2*rr-1:2*rr) = [ 0 imag(V(indmax(2*rr-1),indmax(2*rr-1))); imag(V(indmax(2*rr),indmax(2*rr))) 0 ]; 
            end
            % now everything should be real
            % test verify
            errnorm = norm( Q*D*Q' - U(:,1:r)*V(1:r,1:r)*U(:,1:r)' , 'fro' );
            %keyboard
            
            % and the objective
            fQ = f_sda( Q , Xdot , X , D );
            
        case 'stiefel'
            % SDA in the stiefel manner... we need a loop here to alternate
            % optimization of D and Q. 
            
            % also, for thoroughness, local optima do seem to exist.
            % Starting with the heuristic method always produces a
            % reasonable result, but a handful of random restarts can
            % improve life...
            
            % ADD RANDOM RESTARTS
            num_restarts = 10;

            [ Q, fQ_check , D] = run_sda( Xdot , X , r , 'heuristic' );
            % note that for the first iteration, we want Q to start at
            % the heuristic solution.
            fQ_best = fQ_check;
            Q_best = Q;
            D_best = D;
            
 
            % now loop
            for nr = 1 : num_restarts
                % decrease num_restarts remaining
                if nr ==1
                    % do no initialization, so that we will start with the
                    % heuristic solution
                else
                    % initial random guess for Q
                    Q = project_stiefel(randn(size(X,1),r));
                    % random guess for D
                    DD = 4*rand(r/2);
                    for rr = 1 : r/2
                        D(2*rr - 1:2*rr, 2*rr-1:2*rr) = [ 0 -DD(rr); DD(rr) 0 ];
                    end
                end
                
                
                % empirically, the starting point does not much matter.  What
                % it does do is seem to help the Q1 Q2 axes be the same,
                % instead of permuted as Q2 Q1 (of course, as a Stiefel method,
                % these are equivalent).  This is of no consequence except for
                % consistency when comparing heuristic and stiefel results.
                %[ Q_jpca , fQ_jpca , D_jpca ] = run_sda( Xdot , X , r , 'heuristic' );
                % [ Q , fQ_jpca , D ] = run_sda( Xdot , X , r , 'heuristic' );
                %rng(11110)
                %Q = project_stiefel(Q + randn(size(X,1),r) );
                %D = D*1.2;
                
                % now alternate with D and Q until convergence
                conv_tol = 1e-4;
                converged = 0;
                loop_length = 100;
                i = 0;
                % evaluate starting difference
                [ Q_check , fQ_check , D_check ] = run_sda(  Xdot , X , r , 'heuristic' );
                % initialize with heuristic? (else random from above)
                %Q = Q_check;
                %D = D_check;
                fQ = f_sda( Q , Xdot , X , D );
                fprintf('==Outer SDA run (random restart %d of %d); diff vs heuristic (%g), iter %d=====\n',nr,num_restarts,fQ - fQ_check,i)
                % main iterative loop
                fQ_prev = fQ;
                while ~converged
                    % check for convergence
                    if i >= loop_length
                        % then stop
                        converged = 1;
                    else
                        %fprintf('\n---------jpca obj: (%g)-----------\n',fQ_jpca)
                        %fprintf('\n--pre step-------(%g)-----------\n',f_sda( Q , Xdot , X , D ) )
                        % take a Q step
                        [ Q , fQ ] = minimize_stiefel_sd( 'f_sda' , Q , [] , Xdot , X , D );
                        % fprintf('\n--after Q-------(%g)-----------\n',f_sda( Q , Xdot , X , D ) )
                        % take a D step
                        % can be done in closed form...
                        for rr = 1 : r/2
                            A = Q'*Xdot;
                            B = Q'*X;
                            % see notes p14 for tedious and uninteresting derivation
                            Dz = (B(1,:)*A(2,:)' - B(2,:)*A(1,:)')/(B(1,:)*B(1,:)' + B(2,:)*B(2,:)');
                            D(2*rr - 1:2*rr, 2*rr-1:2*rr) = [ 0 -Dz; Dz 0 ];
                        end
                    end
                    % reevaluate
                    fQ = f_sda( Q , Xdot , X , D );
                    %fprintf('\n--after D-------(%g)-----------\n',f_sda( Q , Xdot , X , D ) )
                    %[ Q_check , fQ_check , D_check ] = run_sda(  Xdot , X , r , 'heuristic' );
                    fprintf('--Diff vs heuristic----(%g)-----------\n',fQ - fQ_check)
                    % check convergence and iterate
                    i = i + 1;
                    if i >= loop_length || norm((fQ - fQ_prev)/fQ_prev) < conv_tol
                        % converged
                        converged = 1;
                    else
                        fQ_prev = fQ;
                    end
                end
                
                % now check if this restart is better than previous...
                fQ_diff_list(nr) = fQ - fQ_check;
                if fQ < fQ_best
                    % then update
                    fQ_best = fQ;
                    Q_best = Q;
                    D_best = D;
                end
                
            end
            fQ = fQ_best;
            Q = Q_best;
            D = D_best;
            
            %hist(fQ_diff_list);
            %keyboard;
            
            
            
        otherwise
            fprintf('Your method is not implemented.  Try stiefel or heuristic.');
            error()
    end
    
    
            