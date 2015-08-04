%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% test_cca.m
%
% This code performs CCA by the stiefel manifold methods and compares 
% those results to the heuristic solver.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ results ] = test_cca( dx , dy , r , parms )
    
    %%%%%%%%%%%%
    % setup
     %%%%%%%%%%%%
     % run twice since we have two data sets in CCA
    parms = verify_test_parms( dx , r , parms );
    parms = verify_test_parms( dy , r , parms );
    % set random seed
    if isfield(parms,'randseed')
        rng(parms.randseed);
    end
    %
    if isfield(parms,'num_points')
        num_points = parms.num_points;
    else
        num_points = 20*dx;
    end
    %
    if isfield(parms,'dz')
        dz = parms.dz;
    else
        dz = ceil(0.5*max(dx,dy));
    end
    %
    if isfield(parms,'dw')
        dw = parms.dw;
    else
        dw = ceil(0.4*max(dx,dy));
    end
    
    
    %%%%%%%%%%%%
    % Generate data for CCA
    %%%%%%%%%%%%
    % two latent datasets just for flexibility if one wishes to change
    % them, but they are effectively one data set.
    Z = randn(dz , num_points);
    W = randn(dw , num_points);
    % noise power
    Q = 0.0002;
    Wx = randn(dx,dw)*W;
    Wy = randn(dy,dw)*W;
    Zx = randn(dx,dz)*Z;
    Zy = randn(dy,dz)*Z;
    
    X = Zx + Wx + Q*rand(dx)*randn(dx,num_points);
    Y = Zy + Wy + Q*rand(dy)*randn(dy,num_points);
    
    % now mean center the data
    X = X - repmat(mean(X,2),1,size(X,2));
    Y = Y - repmat(mean(Y,2),1,size(Y,2));

    %%%%%%%%%%%%%%%
    % Run CCA methods
    %%%%%%%%%%%%%%%
    % initial random guess to make all methods comparable
    % note that CCA begins with running heuristic in any case, so the
    % random starting points do not matter
    Qx_0 = project_stiefel(randn(size(X,1),r));
    Qy_0 = project_stiefel(randn(size(Y,1),r));
    
    optims = {'heuristic','stiefel_trust_prod'};%,'stiefel_trust'};
    for i = 1 : length(optims)
        % run in specified manner
        [ res_arg(i).Q.x , res_arg(i).Q.y , res_arg(i).fQ , res_arg(i).info ] = run_cca( X , Y , r , optims{i} , Qx_0 , Qy_0 );
        % add to greater struct
        res_arg(i).optim = optims{i};
    end

    [ results ] = get_results_all( 'CCA' , struct('x',dx,'y',dy) , r , res_arg );
    
    %%%%%%%%%%%%
    % Figures
    %%%%%%%%%%%%
    if parms.show_fig
        % plot reduced r dimensions
        figure;
            for kk = 1 : 2
                % plot data dimensions
                subplot(1,2,kk);
                hold on;
                set(gca,'linewidth',2,'fontsize',22);
                set(gca,'xtick',[-inf inf]);
                set(gca,'ytick',[-inf inf]);
                axis off
                % choose data for plot
                switch kk
                    case 1
                        X1 = res_arg(1).Q.x(:,1)'*X;
                        X2 = res_arg(1).Q.x(:,2)'*X;
                        Y1 = res_arg(1).Q.y(:,1)'*Y;
                        Y2 = res_arg(1).Q.y(:,2)'*Y;
                        title('Traditional CCA');
                    case 2
                        X1 = res_arg(2).Q.x(:,1)'*X;
                        X2 = res_arg(2).Q.x(:,2)'*X;
                        Y1 = res_arg(2).Q.y(:,1)'*Y;
                        Y2 = res_arg(2).Q.y(:,2)'*Y;
                        title('Stiefel CCA');
                end
                % debug
                realness = [isreal(X) isreal(Y) isreal(res_arg(1).Q.x) isreal(res_arg(1).Q.y) isreal(res_arg(2).Q.x) isreal(res_arg(2).Q.y) ]
                if any(~realness)
                    keyboard
                end
                % plot data points
                plot(X1,X2,'r.','markersize',14);
                plot(Y1,Y2,'b.','markersize',14);
                
            end
            
        end
        
        % save as appropriate
        if parms.save_fig
            print(gcf,'-depsc','results/fig_cca.eps');
        end

    end
    
