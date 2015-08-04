%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% run_cca.m
%
% This code runs CCA by the chosen method.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ Qx , Qy , fQ , info ] = run_cca(  X , Y , r , method , Qx_0 , Qy_0 )

    %%%%%%%%%
    % input checking
    %%%%%%%%%
    if nargin < 6 || isempty(Qy_0)
        % initial random guess
        Qy_0 = project_stiefel(randn(size(Y,1),r));
    end
    if nargin < 5 || isempty(Qx_0)
        % initial random guess
        Qx_0 = project_stiefel(randn(size(X,1),r));
    end
    if nargin < 4 || isempty(method)
        % default stiefel
        method = 'stiefel';
    end
    if nargin < 3 || isempty(r)
        % default the standard choice
        r = 2;
    end
    if nargin < 2 || isempty(Y) || isempty(X)
        % unacceptable
        error('ERROR: Bad datasets X and/or Y.')
    end
    if size(X,2)~=size(Y,2)
        % unacceptable
        error('ERROR: X and Y must have the same number of data points (columns)');
    end
    if size(X,1) < r || size(Y,1) < r
        % unacceptable
        error('ERROR: X and Y must not have lower dimension than r');
    end
    
    
    %%%%%%%%%
    % Make CCA matrices, regardless of the method
    %%%%%%%%%
    Xm = X - repmat(mean(X,2),1,size(X,2));
    Ym = Y - repmat(mean(Y,2),1,size(Y,2));
    Sxy = Xm*Ym';
    Sxx = Xm*Xm';
    Syy = Ym*Ym';
    
            
    
    %%%%%%%%%
    % CCA in the specified manner
    %%%%%%%%%
    switch(lower(method))
        
        case 'heuristic'
            
            % start timer
            t0 = tic;

            % CCA in the 'trace of quotient' classical manner
            [Ux,Vx] = eig(Sxx\(Sxy*(Syy\(Sxy'))));
            [Uy,Vy] = eig(Syy\((Sxy')*(Sxx\Sxy)));
            
            [Ux,Vx] = eig(inv(Sxx)*Sxy*inv(Syy)*(Sxy'));
            [Uy,Vy] = eig(Syy\((Sxy')*(Sxx\Sxy)));
            
            % the above methods are seemingly quite messy
            % we use Muirhead formulation as it is a bit more stable than
            % others.
            
            %%%% Muirhead option (exactly)
            isqSxx = Sxx^(-1/2);
            isqSyy = Syy^(-1/2);
            [H,P,Q] = svd(isqSxx*Sxy*isqSyy);
            H = H';
            Q = Q';
            % such that H'PQ = isq...
            Lx = H*isqSxx;
            Ly = Q*isqSyy;
            
            %%%% MATLAB built-in option
            [Ax,Ay,R,Yx,Yy] = canoncorr(X',Y');
            
            %%%% generic option
            [Cx,Dx] = eig(isqSxx*Sxy*inv(Syy)*Sxy'*isqSxx);
            Ux = isqSxx*Cx';

            %%%%
            % now compare Lx, Ly, Ux, ...
            % Lx and Ax' should be the same, or at least are when things
            % are well conditioned... and up to a sign. 
            % Lx is such that Lx*X is the low d proj... so it is the ROWS
            % of Lx that form the projection...
            % when poorly conditioned, these are not the same...
            Ux_muirhead = Lx';
            Uy_muirhead = Ly';
            Ux_matlab = Ax;
            Ux_generic = Ux;
            % debug
            % Ux_muirhead./Ux_matlab
            
            % assign the good methodology
            Ux = Ux_muirhead;
            Uy = Uy_muirhead;
            Vx = P;
            Vy = P;
            
            % check to make sure these are real and correct if not
            Ux = get_real_basis( Ux , 'Ux' );
            Uy = get_real_basis( Uy , 'Uy' );
            
            % take the top r dimensions
            [ ~ , indmax ] = sort(abs(diag(Vx)),'descend');
            Qx = Ux(:,indmax(1:r));
            [ ~ , indmax ] = sort(abs(diag(Vy)),'descend');
            Qy = Uy(:,indmax(1:r));
            % these will not be stiefel points, so correct that:
            Qx = project_stiefel(Qx);
            Qy = project_stiefel(Qy);
            
            
            % now evaluate the objective
            fQ = f_cca( Qx , Qy , [] , [] , Sxx , Syy , Sxy );
            
            % or do it directly so we can go easier on traditional cca
            % note the use of the absolute value here.  This is to give
            % traditional CCA as much of an advantage as possible.
            % Sometimes the orthogonalization procedure (project_stiefel)
            % can point in the wrong direction, resulting in anticorrelated
            % points.  That is "correct" in terms of the strictest
            % definition of the problem, but it is needlessly punitive.
            % Instead we take the absolute value so CCA is given a helping
            % hand.  So that is (1) do nothing (TOO PUNITIVE); (2) make
            % values as positive as possible (A BIT TOO GENEROUS, but
            % conservative); or another choice is (3) look for the best
            % basis amongst the r dimensional space spanned by Ux (and
            % similar for Uy).  This is a fine idea, but in and of itself
            % it requires a stiefel-type method (general orthogonal group
            % O^{r \times r}), so it would require just as much machinery.
            
            %trXY = trace(abs(Qx'*Sxy*Qy));
            %trXX = trace(Qx'*Sxx*Qx);
            %trYY = trace(Qy'*Syy*Qy);
            %fQ = - trXY / sqrt(trXX*trYY) ;
            
            %keyboard
            % info 
            info(1).time = toc(t0);
            info(1).iter = 1;
            
        case 'stiefel'
            % CCA in the simplest stiefel manner... we need a loop here to alternate
            % optimization of Qx and Qy, since they are both stiefel points

            % we will also do some number of random restarts.
            num_restarts = 1;
            [ Qx , Qy , fQ_check , info ] = run_cca(  X , Y , r , 'heuristic' );
            % note that for the first iteration, we want Qx,Qy to start at
            % the heuristic solution.
            fQ_best = fQ_check;
            Qx_best = Qx;
            Qy_best = Qy;
            
            % now loop
            for nr = 1 : num_restarts
                % decrease num_restarts remaining
                if nr ==1
                    % do no initialization, so that we will start with the
                    % heuristic solution
                else
                    % initial random guess for Q
                    Qx = project_stiefel(randn(size(X,1),r));
                    Qy = project_stiefel(randn(size(Y,1),r));
                end
                
                % now alternate with Qx and Qy until convergence
                conv_tol = 1e-8;
                converged = 0;
                loop_length = 100;
                i = 0;
                lastind=0;
                % evaluate baseline
                fQ = f_cca( Qx , Qy , [] , [] , Sxx , Syy , Sxy );
                fprintf('==Outer CCA run (random restart %d of %d); diff vs heuristic (%g), iter %d=====\n',nr,num_restarts,fQ - fQ_check,i)
                
                %fprintf('--Diff vs heuristic, iter %d----(%g)-----------\n',i,fQ - fQ_check)
                % main iterative loop
                fQ_prev = fQ;
                while ~converged
                    % check for convergence
                    if i >= loop_length
                        % then stop
                        converged = 1;
                    else
                        %fprintf('\n---------heuristic obj: (%g)-----------\n',fQ_check)
                        %fprintf('\n--pre step-------(%g)-----------\n',f_cca( Qx , Qy , [] , [] , Sxx , Syy , Sxy ) )
                        % take a Qx step
                        [ Qx , fQ , info_x] = minimize_stiefel_sd( 'f_cca' , Qx , [] , Qy , [] , [] , Sxx , Syy , Sxy );
                        % fprintf('\n--after Qx-------(%g)-----------\n',f_cca(Qx , Qy , [] , [] , Sxx , Syy , Sxy ) )
                        % take a Qx step
                        % NOTE: there is a bit of monkey business here, as
                        % minimize() only iterates over the first input
                        % argument of f_cca (which is "Qx").  So, here we
                        % reverse all the arguments so that f_cca is operating
                        % with Qy first.  Conveniently, f_cca is invariant to
                        % this ordering, so the answer is correct.
                        [ Qy , fQ , info_y] = minimize_stiefel_sd( 'f_cca' , Qy , [] , Qx , [] , [] , Syy , Sxx , Sxy' );
                    end
                    % reevaluate
                    fQ = f_cca( Qx , Qy , [] , [] , Sxx , Syy , Sxy );
                    % check convergence and iterate
                    i = i + 1;
                    % fprintf('\n--after Qx-------(%g)-----------\n',f_cca(Qx , Qy , [] , [] , Sxx , Syy , Sxy ) )
                    fprintf('--Diff vs heuristic, iter %d----(%g)-----------\n',i,fQ - fQ_check)
                    if i >= loop_length || norm((fQ - fQ_prev)/fQ_prev) < conv_tol
                        % converged
                        converged = 1;
                    else
                        fQ_prev = fQ;
                    end
                    
                    % extend info
                    
                    % extend info by concatenating info_x and info_y... just in
                    % time and iteration.
                    if lastind > 0
                        for j = 1 : length(info_x)
                            info(lastind+j).iter = info(lastind).iter + j;
                            info(lastind+j).time = info(lastind).time + info_x(j).time;
                            info(lastind+j).cost = info_x(j).cost;
                        end
                    else
                        for j = 1 : length(info_x)
                            info(lastind+j).iter = j;
                            info(lastind+j).time = info_x(j).time;
                            info(lastind+j).cost = info_x(j).cost;
                        end
                    end
                    lastind = lastind + length(info_x);
                    for j = 1 : length(info_y)
                        info(lastind+j).iter = info(lastind).iter + j;
                        info(lastind+j).time = info(lastind).time + info_y(j).time;
                        info(lastind+j).cost = info_y(j).cost;
                    end
                    lastind = lastind + length(info_y);
                end
                
                % now check if this restart is better than previous...
                fQ_diff_list(nr) = fQ - fQ_check;
                if fQ < fQ_best
                    % then update
                    fQ_best = fQ;
                    Qx_best = Qx;
                    Qy_best = Qy;
                end
                
            end
            
            fQ = fQ_best;
            Qx = Qx_best;
            Qy = Qy_best;
            fQ_diff_list;
            
            
        % note: the trust region methods do not use random restarts, as we
        % found them to be unnecessary in CCA generally (it should be
        % similarly switched off above with num_restarts = 1).
        case 'stiefel_trust'
            % CCA in the stiefel manner with trust regions manopt
            % implementation. 

            [ Qx , Qy , fQ_check , info ] = run_cca(  X , Y , r , 'heuristic' );
            % now alternate with Qx and Qy until convergence
            conv_tol = 1e-8;
            converged = 0;
            loop_length = 100;
            i = 0;
            lastind = 0; % for info struct
            % evaluate baseline
            fQ = f_cca( Qx , Qy , [] , [] , Sxx , Syy , Sxy );
            fprintf('==Outer CCA run (random restart %d of %d); diff vs heuristic (%g), iter %d=====\n',1,1,fQ - fQ_check,i)
            
            fQ_prev = fQ;
            while ~converged
                % check for convergence
                if i >= loop_length
                    % then stop
                    converged = 1;
                else
                    %fprintf('\n---------heuristic obj: (%g)-----------\n',fQ_check)
                    %fprintf('\n--pre step-------(%g)-----------\n',f_cca( Qx , Qy , [] , [] , Sxx , Syy , Sxy ) )
                    % take a Qx step
                    [ Qx , fQ , info_x] = minimize_stiefel_trust( 'f_cca' , Qx , [] , Qy , [] , [] , Sxx , Syy , Sxy );
                    % fprintf('\n--after Qx-------(%g)-----------\n',f_cca(Qx , Qy , [] , [] , Sxx , Syy , Sxy ) )
                    % take a Qx step
                    % NOTE: there is a bit of monkey business here, as
                    % minimize() only iterates over the first input
                    % argument of f_cca (which is "Qx").  So, here we
                    % reverse all the arguments so that f_cca is operating
                    % with Qy first.  Conveniently, f_cca is invariant to
                    % this ordering, so the answer is correct.
                    [ Qy , fQ , info_y] = minimize_stiefel_trust( 'f_cca' , Qy , [] , Qx , [] , [] , Syy , Sxx , Sxy' );
                end
                % reevaluate
                fQ = f_cca( Qx , Qy , [] , [] , Sxx , Syy , Sxy );
                % check convergence and iterate
                i = i + 1;
                % fprintf('\n--after Qx-------(%g)-----------\n',f_cca(Qx , Qy , [] , [] , Sxx , Syy , Sxy ) )
                fprintf('--Diff vs heuristic, iter %d----(%g)-----------\n',i,fQ - fQ_check)
                if i >= loop_length || norm((fQ - fQ_prev)/fQ_prev) < conv_tol
                    % converged
                    converged = 1;
                else
                    fQ_prev = fQ;
                end
                % extend info by concatenating info_x and info_y... just in
                % time and iteration.
                if lastind > 0
                    for j = 1 : length(info_x)
                        info(lastind+j).iter = info(lastind).iter + j;
                        info(lastind+j).time = info(lastind).time + info_x(j).time;
                        info(lastind+j).cost = info_x(j).cost;
                    end
                else
                    for j = 1 : length(info_x)
                        info(lastind+j).iter = j;
                        info(lastind+j).time = info_x(j).time;
                        info(lastind+j).cost = info_x(j).cost;
                    end
                end
                lastind = lastind + length(info_x);
                for j = 1 : length(info_y)
                    info(lastind+j).iter = info(lastind).iter + j;
                    info(lastind+j).time = info(lastind).time + info_y(j).time;
                    info(lastind+j).cost = info_y(j).cost;
                end
                lastind = lastind + length(info_y);
          
            end

        case 'stiefel_trust_prod'
            % CCA in the stiefel manner with trust regions manopt
            % implementation.  

            [ Qx , Qy , fQ_check , info ] = run_cca(  X , Y , r , 'heuristic' );
            % now alternate with Qx and Qy until convergence
            conv_tol = 1e-8;
            converged = 0;
            loop_length = 100;
            i = 0;
            % evaluate baseline
            fQ_prev = f_cca( Qx , Qy , [] , [] , Sxx , Syy , Sxy );
            % optimize over the product manifold.
            [ Qx , Qy , fQ , info ] = minimize_stiefel_trust_prod('f_cca' , Qx , Qy , [] , Sxx , Syy , Sxy );
        
        otherwise
            fprintf('Your method is not implemented.  Try stiefel or heuristic.');
            error()
    end
    
    if any(any(~isreal(Qx))) || any(any(~isreal(Qy)))
        fprintf('Complex Q from %s CCA method!',method);
        keyboard
    end
    
    
            