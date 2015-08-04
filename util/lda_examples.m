% make LDA expository example...

% note:  
% 
%     %%%%%%%%%%%%
%     % make LDA data
%     %%%%%%%%%%%%
%     % note that the number of classes in LDA essentially specify the
%     % dimensionality desired, since you are looking for a hyperplane that
%     % can separate those classes
%     obvious_counterexample = 0;
%     
%     if obvious_counterexample
%         % use Yan and Tang 2006
%         d = 3;
%         r = 2;
%         num_classes = 4;
%         Sb0 = diag([10, 100, 2]);
%         Sw0 = diag([ 1 , 20, 1]);        
%         % draw data from this model
%         points_per_class = 10000;
%         % set the total number of data points
%         n = num_classes*points_per_class;
%         % initialize the data matrix
%         X = zeros( d , n );
%         % now loop through the classes and make data
%         % first make cluster means that will result in Sb0
%         % note that a regular polyhedron is an easy way to do this.
%         mu0 = 2*sqrt(Sb0)*(1/sqrt(2)*[ 1 -1 0 0 ; 0 0 1 -1 ; -1/sqrt(2) -1/sqrt(2) 1/sqrt(2) 1/sqrt(2) ]);
%         
%         for c = 1 : num_classes
%             % draw points according to Sb and Sw above
%             %X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) = sqrt(Sw0)*randn( d , points_per_class ) + repmat( num_classes/(num_classes-1)*sqrt(Sb0(:,c)) + (1/(num_classes-1))*sqrt(diag(Sb0)) , 1 , points_per_class) ;
%             X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) = sqrt(Sw0)*randn( d , points_per_class ) + repmat( mu0(:,c) , 1 , points_per_class) ;
%         end
%         % now make labels
%         X_labels = repmat( [ 1 : num_classes ] , points_per_class , 1 );
%         X_labels = X_labels(:);
%         
%     else
%         % just draw data sensibly
%         points_per_class = 1000;
%         cloud_spread = 5;
%         % number of dimensions - 1 is the dim of the separating hyperplane, and
%         % thus we will have one more class than reduced dimensions r.
%         num_classes = r + 1;
%         % set the total number of data points
%         n = num_classes*points_per_class;
%         % initialize the data matrix
%         X = zeros( d , n );
%         % now loop through the classes and make data
%         % for each class, randomly draw ppc points in d dim space.  Then
%         % give those eccentricity that is Wishart distributed.  So each
%         % cluster is conditionally normally distributed, given a Wishart
%         % covariance.  
%         % Then, those points are spread a random direction away from zero,
%         % where the distance of that direction is
%         % norm(cloud_spread*randn(d,1)).  Note that this is the
%         % cloud_spread times d times the st dev of these rvs (with
%         % sufficiently large d)... so cloud_spread*d*1.  Thus cloud_spread
%         % should collapse like 1/d...
%         for c = 1 : num_classes
%             % draw points with random eccentricity and random mean... use
%             % an exp distribution to keep the eccentricity from getting
%             % worse and worse conditioned, as would happen (more) with
%             % randn(d).
%             K = diag(exprnd(5,d,1))*project_stiefel(randn(d));
%             X( : , (c - 1)*points_per_class + 1 : c*points_per_class ) = K*randn( d , points_per_class ) + cloud_spread*repmat(randn( d , 1) , 1 , points_per_class) ;
%         end
%         % now make labels
%         X_labels = repmat( [ 1 : num_classes ] , points_per_class , 1 );
%         X_labels = X_labels(:);
%         
%     end
%    


for i = [ 30 27 31 36 ] ;
%for i = [201:251]
    test_lda(3,2, struct('show_fig',1,'save_fig',0,'randseed',i));
    
end






