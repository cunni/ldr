%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% test_lda.m
%
% This code performs LDA by the stiefel manifold methods and compares 
% those results to the heuristic solver.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ results ] = test_lda( d , r , parms)
    
    %%%%%%%%%%%%
    % setup
    %%%%%%%%%%%%
    parms = verify_test_parms( d , r , parms );
    % set random seed
    rng(parms.randseed);
    
    
    %%%%%%%%%%%%
    % make LDA data
    %%%%%%%%%%%%
    % note that the number of classes in LDA essentially specify the
    % dimensionality desired, since you are looking for a hyperplane that
    % can separate those classes
    obvious_counterexample = 0;
    
    if obvious_counterexample
        % use Yan and Tang 2006
        d = 3;
        r = 2;
        num_classes = 4;
        Sb0 = diag([10, 100, 2]);
        Sw0 = diag([ 1 , 20, 1]);        
        % draw data from this model
        points_per_class = 10000;
        % set the total number of data points
        n = num_classes*points_per_class;
        % initialize the data matrix
        X = zeros( d , n );
        % now loop through the classes and make data
        % first make cluster means that will result in Sb0
        % note that a regular polyhedron is an easy way to do this.
        mu0 = 2*sqrt(Sb0)*(1/sqrt(2)*[ 1 -1 0 0 ; 0 0 1 -1 ; -1/sqrt(2) -1/sqrt(2) 1/sqrt(2) 1/sqrt(2) ]);
        
        for c = 1 : num_classes
            % draw points according to Sb and Sw above
            X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) = sqrt(Sw0)*randn( d , points_per_class ) + repmat( mu0(:,c) , 1 , points_per_class) ;
        end
        % now make labels
        X_labels = repmat( [ 1 : num_classes ] , points_per_class , 1 );
        X_labels = X_labels(:);
        
    else
        % just draw data sensibly
        points_per_class = 2000;
        cloud_spread = 5/d;
        % number of dimensions - 1 is the dim of the separating hyperplane, and
        % thus we will have one more class than reduced dimensions r.
        num_classes = r + 1;
        % set the total number of data points
        n = num_classes*points_per_class;
        % initialize the data matrix
        X = zeros( d , n );
        % now loop through the classes and make data
        % for each class, randomly draw ppc points in d dim space.  
        % Then, those points are spread a random direction away from zero,
        % where the distance of that direction is
        % norm(cloud_spread*randn(d,1)).  Note that this is the
        % cloud_spread times d times the st dev of these rvs (with
        % sufficiently large d)... so cloud_spread*d*1.  Thus cloud_spread
        % should collapse like 1/sqrt(d) or 1/d, to control the spread...
        for c = 1 : num_classes
            % draw points with random eccentricity and random mean... use
            % an exp distribution to keep the eccentricity from getting
            % worse and worse conditioned, as would happen (more) with
            % randn(d).
            K = diag(exprnd(5,d,1))*project_stiefel(randn(d));
            X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) = K*randn( d , points_per_class ) + cloud_spread*repmat(randn( d , 1) , 1 , points_per_class) ;
        end
        % now make labels
        X_labels = repmat( [ 1 : num_classes ] , points_per_class , 1 );
        X_labels = X_labels(:);
        
    end
    

    %%%%%%%%%%%%
    % RUN LDA 
    %%%%%%%%%%%%
    % initial random guess to make all methods comparable
    Q_0 = project_stiefel(randn(size(X,1),r));
    
    optims = {'heuristic','grassmann_trust'};
    for i = 1 : length(optims)
        % run in specified manner
        [ res_arg(i).Q , res_arg(i).fQ , res_arg(i).info ] = run_lda( X , X_labels , r , optims{i} , Q_0 );
        % add to greater struct
        res_arg(i).optim = optims{i};
    end

    [ results ] = get_results_all( 'LDA' , d , r , res_arg );


    
    %%%%%%%%%%%%
    % Figures
    %%%%%%%%%%%%
    if parms.show_fig
        % plot first two data dimensions
%         figure;
%         subplot(1,3,1);
%         hold on;
%         set(gca,'linewidth',2,'fontsize',22);
%         set(gca,'xtick',[-inf inf]);
%         set(gca,'ytick',[-inf inf]);
%         axis off
        if mod(num_classes,2)==0
            htmp = redgreencmap(num_classes+1, 'interpolation', 'linear');
            htmp = [htmp(1:num_classes/2,:) ; htmp(num_classes/2+2:end,:)];
        elseif num_classes ==3
            htmp = [ 1 0 0 ; 0 0 0 ; 0 0 1 ];
        else
            htmp = redgreencmap(num_classes, 'interpolation', 'linear');
        end
%         for c = 1 : num_classes
%             %plot3( X( 1 , (c - 1)*points_per_class + 1 : c*points_per_class ) , X( 2 , (c - 1)*points_per_class + 1 : c*points_per_class ) , X( 3 , (c - 1)*points_per_class + 1 : c*points_per_class ) , 'o' , 'color', htmp(c,:) , 'markersize',2 );
%             plot( X( 1 , (c - 1)*points_per_class + 1 : c*points_per_class ) , X( 2 , (c - 1)*points_per_class + 1 : c*points_per_class ) , 'o' , 'color', htmp(c,:) , 'markersize',2 );
%         end
%         title('Random Proj.');
        
        figure
        ha = tight_subplot(1, 2, [.01 .03],[.04 .0],[.01 .01]);
        ms = 3;
        
        % plot projections
        % verify sensible action... plot first 2 dimensions of LDA
        axes(ha(1));
        hold on;
        set(gca,'linewidth',2,'fontsize',18);
        set(gca,'xtick',[-inf inf]);
        set(gca,'ytick',[-inf inf]);
        %axis off
        for c = 1 : num_classes
            %plot3( res_arg(2).Q(:,1:2)'*X( 1 , (c - 1)*points_per_class + 1 : c*points_per_class ) , res_arg(2).Q(:,1:2)'*X( 2 , (c - 1)*points_per_class + 1 : c*points_per_class ) , X( 3 , (c - 1)*points_per_class + 1 : c*points_per_class ) , 'o' , 'color', htmp(c,:) , 'markersize',2 );
            if r > 2
                plot3( res_arg(1).Q(:,1)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , res_arg(1).Q(:,2)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , res_arg(1).Q(:,3)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , 'o' , 'color', htmp(c,:) , 'markersize',2 );
            elseif r==2
                plot( res_arg(1).Q(:,1)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , res_arg(1).Q(:,2)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , 'o' , 'color', htmp(c,:) , 'markersize',ms,'linewidth',ms);
            else
                plot( res_arg(1).Q(:,1)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , res_arg(1).Q(:,1)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) ,'o' , 'color', htmp(c,:) , 'markersize',ms,'linewidth',ms);
            end
        end
        %title('Traditional LDA');
        title('Heuristic LDA');
        
        % verify sensible action... plot first 2 dimensions of LDA
        axes(ha(2));
        hold on;
        set(gca,'xtick',[-inf inf]);
        set(gca,'ytick',[-inf inf]);
        set(gca,'linewidth',2,'fontsize',18);
        for c = 1 : num_classes
            if r > 2
                plot3( res_arg(2).Q(:,1)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , res_arg(2).Q(:,2)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , res_arg(2).Q(:,3)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , 'o' , 'color', htmp(c,:) , 'markersize',2 );
            elseif r == 2
                plot( res_arg(2).Q(:,1)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , res_arg(2).Q(:,2)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , 'o' , 'color', htmp(c,:) , 'markersize',ms,'linewidth',ms);
            else
                plot( res_arg(2).Q(:,1)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , res_arg(2).Q(:,1)'*X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) , 'o' , 'color', htmp(c,:) , 'markersize',ms,'linewidth',ms);
            end
        end
        title('Stiefel LDA');
        
        suptitle(sprintf('(improvement of %0.2f)',results.normed_diff(2)));
        
        % save as appropriate
        if parms.save_fig
            print(gcf,'-depsc',sprintf('results/fig_lda_%d.eps',parms.randseed));
        end

    end

