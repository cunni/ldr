%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% run_pca.m
%
% This code runs PCA by the chosen method.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ Q , fQ , info ] = run_pca( X , r , method , Q_0 )

    %%%%%%%%%
    % input checking
    %%%%%%%%%
    if nargin < 4 || isempty(Q_0)
        % initial random guess
        Q_0 = project_stiefel(randn(size(X,1),r));
    end
    if nargin < 3 || isempty(method)
        % default stiefel
        method = 'stiefel';
    end
    

    %%%%%%%%%
    % PCA in the specified manner
    %%%%%%%%%
    switch(lower(method))
        
        case 'heuristic'
            % PCA in the heuristic manner
            % start timer
            t0 = tic;
            
            % do svd to get the bases in the d space
            [U,~,~] = svd(X);
            % the result is the first r columns of U
            Q = U(:,1:r);
            % and the objective
            fQ = f_pca( Q , X );
            % info 
            info(1).time = toc(t0);
            info(1).iter = 1;
            info(1).cost = fQ;
            
        case 'stiefel'
            % PCA in the stiefel sd manner
            
            % debug
            % [U,~,~] = svd(X);
            % the result is the first r columns of U
            %Q_0 = U(:,1:r);
            %Q_0 = project_stiefel(Q_0 + 0.00*randn(size(Q_0)));
            
            [ Q , fQ , info ] = minimize_stiefel_sd( 'f_pca' , Q_0 , [] , X );
            
        case 'stiefel_trust'
            % PCA in the stiefel manner with trust regions manopt
            % implementation
            
            [ Q , fQ , info ] = minimize_stiefel_trust( 'f_pca' , Q_0 , [] , X );
            
        case 'grassmann_trust'
            % PCA in the grassmann manner with trust regions manopt
            % implementation
            
            [ Q , fQ , info ] = minimize_grassmann_trust( 'f_pca' , Q_0 , [] , X );
            
        case 'stiefel_mosd'
            % PCA in the stiefel manner with manopt steepest descent
            % implementation
            
            [ Q , fQ , info ] = minimize_stiefel_mosd( 'f_pca' , Q_0 , [] , X );
            
        case 'grassmann_mosd'
            % PCA in the grassmann manner with manopt steepest descent
            % implementation
            
            [ Q , fQ , info ] = minimize_grassmann_mosd( 'f_pca' , Q_0 , [] , X );
            
        otherwise
            error('Your method is not implemented.  Try stiefel or heuristic.');
            
    end
    
            
            