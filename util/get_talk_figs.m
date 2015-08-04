%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% Byron M Yu
% 2014
%
% a bit of code adapted from the code of Figure 4 in nn_review project
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;
clear all;

save_figs = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate data points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_points = 100;
rng(0);
% possibly load ../results/nice_rng_states.mat
rho = 0.4;
R = [ 1 rho rho ; rho 1 rho ; rho rho 1 ];
C = R * diag([1 1 1]) * R';
X = mvnrnd([0 0 0], C, num_points)'; % 3 x num_points
X = R*rand(3,num_points);
X = X - repmat(mean(X,2),1,size(X,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do PCA
[U,D] = svd(X*X');
W = U(:,1:2);

% plot
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_pca_3d_no_w.eps');
f3 = plot_data_3d( X , [] , [] , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_pca_3d.eps');
f3 = plot_data_3d( X , W , [] , [] , parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_pca_2d.eps','vec_plot',0);
f2 = plot_data_2d( W'*X , 'w' , [] , [] , parms);


% make a PCA video
% random starting point
W0 = project_stiefel(randn(3,2));
% we will mix this with the result

% video calls
fname_vid = '../results/video_pca_2w.avi';
pca_vo = VideoWriter(fname_vid);
pca_vo.FrameRate = 10;
open(pca_vo);
total_time = 5; % seconds
total_frames = total_time*pca_vo.FrameRate;

parms = struct('save_fig',0);
for i = 1:total_frames
    f3 = plot_data_3d( X , project_stiefel( i/total_frames*W + (total_frames-i)/total_frames*W0) , [] , [], parms );
    pause(0.2)
    writeVideo(pca_vo , getframe(f3) );
    %close gcf
end
close(pca_vo);

%%%%%%%%%%%%%%%%
% PCA objective plot
%%%%%%%%%%%%%%%%
% for plotting purposes we have to do this on two d data
X2= X([1,3],:);
w = [-.6;0.8];
% plot the data
parms = struct('save_fig',save_figs,'fname_fig',sprintf('../results/obj_pca_data_a.eps'),'vec_plot',0);
f3 = plot_data_2d( X2 , 'x' , [] , [], parms);

% plot surf
parms = struct('save_fig',save_figs,'plot_arc',0,'fname_fig','../results/obj_pca_a.tiff');
f = plot_objective( X2 ,[], @f_pca , parms );
%plot surf + constraint
parms = struct('save_fig',save_figs,'plot_arc',1,'fname_fig','../results/obj_pca_b.tiff');
f = plot_objective( X2 ,[], @f_pca , parms );
% plot surf + constraint + point
step = struct('plot_step',0,'plot_pt',1,'label','w','num',0,'pt',w,'grad',[0.2; 0.4]);
parms = struct('save_fig',save_figs,'plot_arc',1,'fname_fig','../results/obj_pca_c.tiff');
f = plot_objective( X2 ,step, @f_pca , parms );


% that also produces a new figure of the data...
parms = struct('save_fig',save_figs,'fname_fig',sprintf('../results/obj_pca_data_c.eps'),'vec_plot',1,'vec',w,'vec_label','w');
f3 = plot_data_2d( X2 , 'x' , [] , [], parms);
    
% now iterate a few times...
stepsize = 0.08;
for i = 0:2:10
    
    
    % calc, step, plot
    % plot surf + constraint + step
    [ff , gradf ] = f_pca( w , X2 );
    step = struct('plot_step',1,'plot_pt',1,'label','w','num',i,'pt',w,'grad',-stepsize*gradf);
    parms = struct('save_fig',save_figs,'plot_arc',1,'fname_fig',sprintf('../results/obj_pca_%d.tiff',i));
    f = plot_objective( X2 ,step, @f_pca , parms );
    
    % step
    w = w - stepsize*gradf;
    % that also produces a new figure of the data...
    parms = struct('save_fig',save_figs,'fname_fig',sprintf('../results/obj_pca_data_%d.eps',i),'vec_plot',1,'vec',w,'vec_label','w');
    f3 = plot_data_2d( X2 , 'x' , [] , [], parms);
    
    % snap back correction term
    snap = (w)/norm(w) - w;
    step = struct('plot_step',1,'plot_pt',1,'label','w','num',i+1,'pt',w,'grad',snap);
    parms = struct('save_fig',save_figs,'plot_arc',1,'fname_fig',sprintf('../results/obj_pca_%d.tiff',i+1));
    f = plot_objective( X2 ,step, @f_pca , parms );
    
    % step
    w = w + snap;    
    % that also produces a new figure of the data...
    parms = struct('save_fig',save_figs,'fname_fig',sprintf('../results/obj_pca_data_%d.eps',i+1),'vec_plot',1,'vec',w,'vec_label','w');
    f3 = plot_data_2d( X2 , 'x' , [] , [], parms);    
    
end



% save a colorbar just for visualization
figure;
h = colorbar
axis off
set(h , 'fontsize', 24)
set(h , 'Ytick' , [-inf inf])
set(h , 'Yticklabel', {'',''})
if save_figs
    print(gcf, '-depsc','../results/colorbar.eps');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LDA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make three classes
% pick a random basis for splitting these groups
C = project_stiefel(rand(3,2));
data_labels = zeros(size(X,2),1);
data_colors = zeros(3,num_points);
for i = 1 : size(X,2)
   if (C(:,1)'*X(:,i) > 0)
       data_labels(i) = 1;
       data_colors(:,i) = [0;0;0.5];
   elseif (C(:,2)'*X(:,i) > -0.2)*(C(:,1)'*X(:,i) < 0)
       data_labels(i) = 2;
       data_colors(:,i) = [0;0;0];
   else
       data_labels(i) = 3;
       data_colors(:,i) = [0;0;1];
   end
end
% run LDA
[Q , fQ ] = run_lda( X , data_labels , 2 , 'heuristic');

% plot
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_lda_3d_no_w.eps');
f3 = plot_data_3d( X , [] , data_colors , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_lda_3d.eps');
f3 = plot_data_3d( X , Q , data_colors , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_lda_2d.eps','vec_plot',0);
f2 = plot_data_2d( Q'*X , 'w' , data_colors , [], parms);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add labels for demixed example
% dPCA should find 'th' as one axis, and orthogonal axis as the
% other axis
th = -0.1 * pi;
R = [cos(th) -sin(th) 0 ; sin(th) cos(th) 0 ; 0 0 1 ];
y = R*X ;  % y indicates 2 attributes of each point x
% Make attributes discrete
[ct1, att1] = histc(y(1,:), [-100 -.5 0 .5 100]);
[ct2, att2] = histc(y(3,:), [-100 -.5 0 .5 100]);
color_strength = [0 0.3 0.7 1];
point_size  = [4 10 20 30];

% make data_colors
for i = 1 : num_points
    data_colors(:,i) = color_strength(att1(i))*[0;0;1];
    data_sizes(i) = point_size(att2(i));
end

% Run DPCA a la Machens 2010
covsum1 = zeros(3);
covsum2 = zeros(3);
for i = 1:4
  idx1 = (att1 == i);
  covsum1 = covsum1 + X(:,idx1)*X(:,idx1)' * sum(idx1); 

  idx2 = (att2 == i);
  covsum2 = covsum2 + X(:,idx2)*X(:,idx2)' * sum(idx2); 
end
covsum1 = covsum1 / num_points;
covsum2 = covsum2 / num_points;
[V, D] = eig(covsum1 - covsum2);
W = V(:,[1,3]);

% plot data
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_dca_3d_no_w.eps');
f3 = plot_data_3d( X , [] , data_colors , data_sizes, parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_dca_3d.eps');
f3 = plot_data_3d( X , W , data_colors , data_sizes, parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_dca_2d.eps','vec_plot',0);
f2 = plot_data_2d( W'*X , 'w' , data_colors , data_sizes, parms);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_figs = 0;
p = struct('show_fig',0,'save_fig',0,'randseed',36,'num_points',100,'dz',2,'dw',2);
res = test_cca( 3 , 3 , 2 , p );


% sort X and Y and color them so people can see what is what
XYsort = sortrows([res.X;res.Y]')';
X = XYsort(1:3,:);
Y = XYsort(4:6,:);

data_colorsX = [ linspace(0.1,1,size(X,2)); zeros(1,size(X,2)); zeros(1,size(X,2)) ];
data_colorsY = [ zeros(1,size(Y,2)); zeros(1,size(Y,2)) ; linspace(0.1,1,size(Y,2)) ];

% plot data
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_X_3d.eps','vec_plot',0);
f3x = plot_data_3d( X , [] , data_colorsX , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_Y_3d.eps','vec_plot',0);
f3y = plot_data_3d( Y , [] , data_colorsY , [], parms);


% plot data heuristic
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_X_3d_heur.eps','vec_plot',0);
f3x = plot_data_3d( X , res.Qx.traditional , data_colorsX , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_Y_3d_heur.eps','vec_plot',0);
f3y = plot_data_3d( Y , res.Qy.traditional , data_colorsY , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_X_2d_heur.eps','vec_plot',0);
f2x = plot_data_2d( res.Qx.traditional'*X , 'w' , data_colorsX , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_Y_2d_heur.eps','vec_plot',0);
f2y = plot_data_2d( res.Qy.traditional'*Y , 'w' , data_colorsY , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_XY_2d_heur.eps','vec_plot',0);
f2both = plot_data_2d( [res.Qx.traditional'*X res.Qy.traditional'*Y] , 'w' , [data_colorsX data_colorsY] , [], parms);

% plot data stiefel
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_X_3d_stief.eps','vec_plot',0);
f3x = plot_data_3d( X , res.Qx.stiefel , data_colorsX , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_Y_3d_stief.eps','vec_plot',0);
f3y = plot_data_3d( Y , res.Qy.stiefel , data_colorsY , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_X_2d_stief.eps','vec_plot',0);
f2x = plot_data_2d( res.Qx.stiefel'*X , 'w' , data_colorsX , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_Y_2d_stief.eps','vec_plot',0);
f2y = plot_data_2d( res.Qy.stiefel'*Y , 'w' , data_colorsY , [], parms);
parms = struct('save_fig',save_figs,'fname_fig','../results/simple_cca_XY_2d_stief.eps','vec_plot',0);
f2both = plot_data_2d( [res.Qx.stiefel'*X res.Qy.stiefel'*Y] , 'w' , [data_colorsX data_colorsY] , [], parms);



% make three CCA video
% optimize, moving from the heuristic to the stiefel result

% video calls
fname_vid = '../results/video_cca_X.avi';
pca_vo = VideoWriter(fname_vid);
pca_vo.FrameRate = 10;
open(pca_vo);
total_time = 5; % seconds
total_frames = total_time*pca_vo.FrameRate;

parms = struct('save_fig',0);
for i = 1:total_frames
    f3 = plot_data_3d( X , project_stiefel( i/total_frames*res.Qx.stiefel + (total_frames-i)/total_frames*res.Qx.traditional) , data_colorsX , [], parms );
    pause(0.2)
    writeVideo(pca_vo , getframe(f3) );
    close gcf
end
close(pca_vo);


% video calls
fname_vid = '../results/video_cca_Y.avi';
pca_vo = VideoWriter(fname_vid);
pca_vo.FrameRate = 10;
open(pca_vo);
total_time = 5; % seconds
total_frames = total_time*pca_vo.FrameRate;

parms = struct('save_fig',0);
for i = 1:total_frames
    f3 = plot_data_3d( X , project_stiefel( i/total_frames*res.Qy.stiefel + (total_frames-i)/total_frames*res.Qy.traditional) , data_colorsY , [], parms );
    pause(0.2)
    writeVideo(pca_vo , getframe(f3) );
    close gcf
end
close(pca_vo);



% video calls
fname_vid = '../results/video_cca_XY.avi';
pca_vo = VideoWriter(fname_vid);
pca_vo.FrameRate = 10;
open(pca_vo);
total_time = 5; % seconds
total_frames = total_time*pca_vo.FrameRate;

parms = struct('save_fig',0,'vec_plot',0);
for i = 1:total_frames
    f2both = plot_data_2d( [project_stiefel( i/total_frames*res.Qx.stiefel + (total_frames-i)/total_frames*res.Qx.traditional)'*X project_stiefel( i/total_frames*res.Qy.stiefel + (total_frames-i)/total_frames*res.Qy.traditional)'*Y] , 'w' , [data_colorsX data_colorsY] , [], parms);
    set(gcf, 'position', [200 200 600 400]);
    pause(0.2)
    writeVideo(pca_vo , getframe(f2both) );
    close gcf
end
close(pca_vo);




%{
p = struct('show_fig',1,'save_fig',0,'randseed',12,'num_points',100, 'dz',2,'dw',2)


p = struct('show_fig',1,'save_fig',0,'randseed',15,'num_points',100, 'dz',2,'dw',2)

res = test_cca( 3 , 3 , 2 , p );

for i = 31 : 70
    p = struct('show_fig',1,'save_fig',0,'randseed',i,'num_points',100, 'dz',2,'dw',2);
    res = test_cca( 3 , 3 , 2 , p );
    title(sprintf('%d|%0.2f|%0.2f',i,res.f.traditional, res.f.stiefel));
end

% 36, 37, 66, 59, 50, 44
%}

