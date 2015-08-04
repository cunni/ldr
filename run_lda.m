%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% run_lda.m
%
% This code runs LDA by the chosen method.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ Q , fQ , info ] = run_lda( X , X_labels , r , method , Q_0 )

    %%%%%%%%%
    % input checking
    %%%%%%%%%
    if nargin < 5 || isempty(Q_0)
        % initial random guess
        Q_0 = project_stiefel(randn(size(X,1),r));
    end
    if nargin < 4 || isempty(method)
        % default stiefel
        method = 'stiefel';
    end
    if nargin < 3 || isempty(r)
        % default the standard choice
        r = length(unique(X_labels)) - 1;
    end
    if nargin < 2 || isempty(X_labels) || isempty(X)
        % unacceptable
        fprintf('ERROR: Bad X and/or X_labels inputs to run_lda.')
        keyboard
    end
    
    %%%%%%%%%
    % form Sb and Sw for LDA, regardless of method
    %%%%%%%%%
    % find the number of classes
    [ classes , ~ , class_ind ] = unique(X_labels);
    % find means
    mu = zeros(size(X,1), length(classes));
    mu_all = mean(X,2);
    for c = 1 : length(classes)
        % find means
        mu(:,c) = mean( X(:, X_labels==classes(c) ) , 2);
    end
    % loop over all data
    Sb = zeros(size(X,1));
    Sw = zeros(size(X,1));
    for n = 1 : size(X,2)
        % calc Sb
        Sb = Sb + ( mu(:,class_ind(n)) - mu_all )*( mu(:,class_ind(n)) - mu_all )';
        % calc Sw
        Sw = Sw + ( X(:,n) - mu(:,class_ind(n)) )*( X(:,n) - mu(:,class_ind(n)) )';
    end
    Sw = Sw/size(X,2);
    Sb = Sb/size(X,2);
    
    %save('SwSb.mat','Sw','Sb');
    
    %%%%%%%%%
    % LDA in the specified manner
    %%%%%%%%%
    switch(lower(method))
        
        case 'heuristic'
            % LDA in the heuristic manner
            % start timer
            t0 = tic;

            % do svd to get the bases in the d space
            [ U , ~ , ~ ] = svd(Sw\Sb);
            % the result is the first r columns of U
            Q = U(:,1:r);
            % and the objective
            fQ = f_lda( Q , [] , [] , Sb , Sw );
            %
            % info 
            info(1).time = toc(t0);
            info(1).iter = 1;

        case 'stiefel'
            % LDA in the stiefel manner with random restarts
            
            % run the method
            [ Q , fQ , info ] = minimize_stiefel_sd( 'f_lda' , Q_0 , [] , [] , [] , Sb , Sw );

            % note that there seems to be no local issues with LDA, so we
            % do not actually need random restarts. 
            num_restarts = 1; % effectively deprecates code below
            if num_restarts > 1
                [ Q, fQ_check ] = run_lda( X , X_labels , r , 'heuristic' );
                % note that for the first iteration, we want Q to start at
                % the heuristic solution.
                fQ_best = fQ_check;
                Q_best = Q;
                
                % now loop
                for nr = 1 : num_restarts
                    % decrease num_restarts remaining
                    if nr ==1
                        % do no initialization, so that we will start with the
                        % heuristic solution
                    else
                        % initial random guess for Q
                        Q = project_stiefel(randn(size(X,1),r));
                    end
                    
                    % run the method
                    [ Q , fQ , info ] = minimize_stiefel_sd( 'f_lda' , Q , [] , [] , [] , Sb , Sw );
                    
                    % now check if this restart is better than previous...
                    fQ_diff_list(nr) = fQ - fQ_check;
                    if fQ < fQ_best
                        % then update
                        fQ_best = fQ;
                        Q_best = Q;
                    end
                    
                end
                fQ = fQ_best;
                Q = Q_best;
            end
            
        % note: the trust region methods do not use random restarts, as we
        % found them to be unnecessary in LDA generally (it should be
        % similarly switched off above with num_restarts = 1).
        case 'stiefel_trust'
            % LDA in the stiefel manner with trust regions manopt
            % implementation
            
            [ Q , fQ , info ] = minimize_stiefel_trust( 'f_lda' ,  Q_0 , [] , [] , [] , Sb , Sw );
            
        case 'grassmann_trust'
            % LDA in the grassmann manner with trust regions manopt
            % implementation
            
            [ Q , fQ , info ] = minimize_grassmann_trust( 'f_lda' ,  Q_0 , [] , [] , [] , Sb , Sw );

        case 'stiefel_mosd'
            % LDA in the stiefel manner with steepest descent manopt
            % implementation
            
            [ Q , fQ , info ] = minimize_stiefel_mosd( 'f_lda' ,  Q_0 , [] , [] , [] , Sb , Sw );
            
        case 'grassmann_mosd'
            % LDA in the grassmann manner with steepest descent manopt
            % implementation
            
            [ Q , fQ , info ] = minimize_grassmann_mosd( 'f_lda' ,  Q_0 , [] , [] , [] , Sb , Sw );

            
        otherwise
            fprintf('Your method is not implemented.  Try stiefel or heuristic.');
            error()
    end
    
            
            