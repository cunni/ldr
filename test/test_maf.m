%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% test_maf.m
%
% This code performs MAF by the stiefel manifold methods and compares 
% those results to the heuristic solver.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ results ] = test_maf( d , r , parms)
    
    %%%%%%%%%%%%
    % setup
    %%%%%%%%%%%%
    parms = verify_test_parms( d , r , parms );
    % set random seed
    rng(parms.randseed);

    
    %%%%%%%%%%%%
    % Generate data for MAF
    %%%%%%%%%%%%
    num_points = 10*d;
    ddyn = d;
    % draw from spline basis
    for dd = 1 : ddyn
        rp = randperm(num_points);
        X(dd,:) = interp1( [0 sort(rp(1:2)) num_points+1] , randn(1,4), [1:num_points] ,'pchip');
    end
    X(ddyn+1:d,:) = 0.01*randn( max(0,d - ddyn) , num_points);
    % scale and rotate them to make things interesting in terms of covariance
    % also add noise for realism (makes things more conservative...)
    X = 1/sqrt(d)*rand(d)*X + 1/sqrt(d)*rand(d)*randn(d , num_points);  
    % now mean center the data
    X = X - repmat(mean(X,2),1,size(X,2));

    %%%%%%%%%%%%
    % RUN MAF
    %%%%%%%%%%%%
    % initial random guess to make all methods comparable
    Q_0 = project_stiefel(randn(size(X,1),r));
    
    optims = {'heuristic','grassmann_trust'};
    for i = 1 : length(optims)
        % run in specified manner
        [ res_arg(i).Q , res_arg(i).fQ , res_arg(i).info ] = run_maf( X , r , optims{i} , Q_0 );
        % add to greater struct
        res_arg(i).optim = optims{i};
    end

    [ results ] = get_results_all( 'MAF' , d , r , res_arg );

    
    %%%%%%%%%%%%
    % Figures
    %%%%%%%%%%%%
    if parms.show_fig
        % plot reduced r dimensions
        figure;
        htmp = redgreencmap(3, 'interpolation', 'linear');
        t = [1:num_points];
        
        for rr = 1 : r
            subplot(1,r,rr);
            hold on;
            set(gca,'linewidth',2,'fontsize',22);
            set(gca,'xtick',[-inf inf]);
            set(gca,'ytick',[-inf inf]);
            axis off
            
            plot( t, X( rr , : ) , 'color', htmp(1,:) , 'linewidth',2 );
            plot( t, res_arg(1).Q(:,rr)'*X ,'color', htmp(2,:) , 'linewidth',2 );
            plot( t, res_arg(2).Q(:,rr)'*X , 'color', htmp(3,:) , 'linewidth',2 );
            legend(gca, 'data',optims{1},optims{2});
        end
        
        % save as appropriate
        if parms.save_fig
            print(gcf,'-depsc','results/fig_maf.eps');
        end

    end
    
