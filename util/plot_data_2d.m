%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 
% 2014
%
% plot for COSYNE 2014 talk
%%%%%%%%%%%%%%%%%%%%%%%%


function [f] = plot_data_2d( X , axis_label , data_colors , data_sizes , parms)

%%%%%%%%%%%%%%%%%
% input checking
%%%%%%%%%%%%%%%%%
if nargin < 5 || isempty(parms)
    parms.save_fig = 0;
    parms.plot_vec = 0;
end
if nargin < 4 || isempty(data_sizes)
    % one size for each data point
    data_sizes = 8*ones(1,size(X,2));
end
if nargin < 3 || isempty(data_colors)
    % make a simple full blue data point for each data
    data_colors = repmat([0 0 1]' , 1 , size(X,2) );
end

% some data normalization to make it fit better in the plot (scale and
% position only, no linear transformations)
% ensure data centering, since the axes will be plotted wrt 0
X = X - repmat(mean(X,2),1,size(X,2));
% now normalize the data by the largest dimension so everything appears in
% the unit box
maxptlen = max( sqrt(sum(X.*X,1)) );
X = 1.4*X / maxptlen;

%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%
% some figure parameters
f = figure;
hold on;
axis off
set(gcf,'color','w');
set(f, 'position', [200 200 600 600]);
lw = 2;
xbump = 0.1;
axbump = 0.2;
xlim([-1 1]);
ylim([-1 1]);

% plot cardinal axes
for i = 1 : 2
    plot( -[1 i~=1] + axbump , -[1 i~=2] + axbump , 'k' , 'linewidth', lw );
    text('string', sprintf('$$\\mathbf{%s}_%d$$',axis_label,i) , 'position', -[i~=1 i~=2]' + axbump + xbump*[i==1 i==2]' , 'fontsize', 24, 'interpreter', 'latex');
end

% scatter data 
for i = 1:size(X,2)
  plot(X(1,i), X(2,i), 'o', 'markerfacecolor', data_colors(:,i) , 'markersize', data_sizes(1,i) , 'markeredgecolor', 'k');
end

% if the parms carry an extra vector, plot it also.
if parms.vec_plot
    plot( [ 0 parms.vec(1) ] , [ 0 parms.vec(2)] , 'k' , 'linewidth', lw );
    text('string', sprintf('$$\\mathbf{%s}_%d$$',parms.vec_label,1) , 'position', parms.vec + xbump*parms.vec , 'fontsize', 24, 'interpreter', 'latex');
end

% save as determined
if parms.save_fig
    % presumes a full pathed parms.fname_fig
    print(gcf,'-depsc',parms.fname_fig);
end


