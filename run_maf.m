%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% run_maf.m
%
% This code runs MAF by the chosen method.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ Q , fQ , info ] = run_maf( X , r , method , Q_0 )

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
    if nargin < 2 || isempty(r)
        % default the standard choice
        r = 2;
    end
    if nargin < 1 || isempty(X)
        % unacceptable
        fprintf('ERROR: Bad X inputs to run_maf.')
        keyboard
    end
    
    %%%%%%%%%
    % form S and St for MAF, regardless of method
    %%%%%%%%%
    % note that there is a presumed temporal structure in this data, as is
    % a requirement for MAF
    % first mean center the data
    Xm = X - repmat(mean(X,2),1,size(X,2));
    % now form the covariance
    S = (1/size(Xm,2))*Xm*Xm' + 0.0*eye(size(X,1));
    % now form the time lag covariance, and symmetrize it.
    St = (1/(size(Xm,2)-1))*Xm(:,2:end)*Xm(:,1:end-1)' + 0*eye(size(X,1));
    St = 0.5*(St + St');
    
    %%%%%%%%%
    % MAF in the specified manner
    %%%%%%%%%
    switch(lower(method))
        
        case 'heuristic'
            % MAF in the heuristic manner
            % start timer
            t0 = tic;

            % do svd to get the bases in the d space
            [ U , ~ , ~ ] = svd(S\St);
            % the result is the first r columns of U
            Q = U(:,1:r);
            % and the objective
            fQ = f_maf( Q , []  , S , St );
            % info 
            info(1).time = toc(t0);
            info(1).iter = 1;
            
        case 'stiefel'
            % MAF in the stiefel manner
            % random restarts do not seem particularly useful here.  Some
            % local optima do exist, but MAF gives a good initialization
            % thus we use that here.
            [ Q_0 , fQ_0 ] = run_maf( X , r , 'heuristic');
            
            [ Q , fQ , info] = minimize_stiefel_sd( 'f_maf' , Q_0 , [] , [] , S , St );
            
        case 'stiefel_trust'

            % better to initialize with heuristic
            [ Q_0 , fQ_0 ] = run_maf( X , r , 'heuristic');
            
            [ Q , fQ , info ] = minimize_stiefel_trust( 'f_maf' , Q_0 , [] , [] , S , St );
            
        case 'grassmann_trust'

            % better to initialize with heuristic
            [ Q_0 , fQ_0 ] = run_maf( X , r , 'heuristic');
            
            [ Q , fQ , info ] = minimize_grassmann_trust( 'f_maf' , Q_0 , [] , [] , S , St );
            
        case 'stiefel_mosd'
            
            % better to initialize with heuristic
            [ Q_0 , fQ_0 ] = run_maf( X , r , 'heuristic');
            
            [ Q , fQ , info ] = minimize_stiefel_mosd( 'f_maf' , Q_0 , [] , [] , S , St );
            
        case 'grassmann_mosd'

            % better to initialize with heuristic
            [ Q_0 , fQ_0 ] = run_maf( X , r , 'heuristic');
            
            [ Q , fQ , info ] = minimize_grassmann_mosd( 'f_maf' , Q_0 , [] , [] , S , St );
            
        otherwise
            fprintf('Your method is not implemented.  Try stiefel or heuristic.');
            error()
    end
    
            
            