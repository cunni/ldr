%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% test_pca.m
%
% This code performs PCA by the stiefel manifold methods and compares 
% those results to the known true answer.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ results ] = test_pca( d , r , parms )
    
    %%%%%%%%%%%%
    % setup
    %%%%%%%%%%%%
    parms = verify_test_parms( d , r , parms );
    % set random seed
    rng(parms.randseed);

    %%%%%%%%%%%%
    % Generate data for PCA
    %%%%%%%%%%%%
    num_points = 2000;
    X = randn( d , num_points );
    % scale them to make things interesting
    X = diag(exprnd(2 , d , 1 ))*X;    
    % mean center the data
    X = X - repmat(mean(X,2),1,size(X,2));
    
    %%%%%%%%%%%%
    % RUN PCA 
    %%%%%%%%%%%%
    % initial random guess to make all methods comparable
    Q_0 = project_stiefel(randn(size(X,1),r));

    optims = {'heuristic','grassmann_trust'};
    for i = 1 : length(optims)
        % run in specified manner
        [ res_arg(i).Q , res_arg(i).fQ , res_arg(i).info ] = run_pca( X , r , optims{i} , Q_0 );
        % add to greater struct
        res_arg(i).optim = optims{i};
    end

    [ results ] = get_results_all( 'PCA' , d , r , res_arg );
    

    
