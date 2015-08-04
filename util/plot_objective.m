%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 
% 2014
%
% plot for COSYNE 2014 talk
%%%%%%%%%%%%%%%%%%%%%%%%


function [f] = plot_objective( X , step , f_handle , parms)

%%%%%%%%%%%%%%%%%
% input checking
%%%%%%%%%%%%%%%%%
if nargin < 4 || isempty(parms)
    parms.save_fig = 0;
    parms.plot_arc = 1;
end
if nargin < 3 || isempty(f_handle)
    % one size for each data point
    f_handle = @f_pca;
end
if nargin < 2 || isempty(step)
    % then plot without a step
    step.plot_step = 0;
    step.plot_pt = 0;
    step.label = 'w';
end
% plotting looks nice at 22, saving as tiff looks better at 14. 
fontsize = 11;
lw = 2;
mksize = 6;

% here we want to plot the objective for PCA with one dimension
outlim = 1.15;
[W1,W2] = meshgrid( [-outlim:.04:outlim] , [0:.02:outlim] );
F = zeros(size(W1));
for i= 1:size(W1,1)
    for j=1:size(W1,2)
        F(i,j) = f_pca( [W1(i,j);W2(i,j)] , X );
    end
end

f = figure;
set(f, 'position', [200 200 1200 500]);
hold on;
surf(W1,W2,F);
shading interp
set( gca, 'xtick' , [-1 1]);
set( gca, 'xlim' , 0.95*[-outlim outlim]);
set( gca, 'ytick' , [0 1]);
set( gca, 'ylim' , 0.95*[0 outlim]);
set( gca, 'fontsize', fontsize );
xlabel(sprintf('$$\\mathbf{%s}_{1,1}$$',step.label),'fontsize',fontsize,'interpreter','latex');
ylabh = get(gca,'ylabel');
set(ylabh,'Rotation',0)
ylabel(sprintf('$$\\mathbf{%s}_{1,2}$$',step.label),'fontsize',fontsize,'interpreter','latex');
set(ylabh,'Position',get(ylabh,'Position') - [.1 0 0])
view([0 90])

% bigz ensures that these things appear atop the surf
bigz = 100+max(max(F));

if parms.plot_arc
    arcx = [-1:.005:1];
    plot3( arcx , sqrt(1 - arcx.^2) , bigz*ones(size(arcx)), 'k' , 'linewidth',lw );
end


if step.plot_step
    % presumes there is a pt, a grad, and a label
    % first plot white arrow
    plot3( [ step.pt(1) step.pt(1)+step.grad(1) ] , [ step.pt(2) step.pt(2)+step.grad(2) ] , [bigz bigz]+10, 'k' , 'linewidth', lw+1 );
    plot3( [ step.pt(1) step.pt(1)+step.grad(1) ] , [ step.pt(2) step.pt(2)+step.grad(2) ] , [bigz bigz]+10, 'w' , 'linewidth', lw );
    % then plot prev and next
    plot3( [ step.pt(1) step.pt(1)+step.grad(1) ] , [ step.pt(2) step.pt(2)+step.grad(2) ] , [bigz bigz]+20, 'o' , 'markerfacecolor', [1 1 1] , 'markersize', mksize , 'markeredgecolor', 'k');
    % labels
    text('string', sprintf('$$\\mathbf{%s}^{(%d)}_1$$',step.label,step.num) , 'position', [step.pt + -0.1*step.grad; bigz+30] , 'fontsize', fontsize, 'interpreter', 'latex');
    text('string', sprintf('$$\\mathbf{%s}^{(%d)}_1$$',step.label,step.num+1) , 'position', [step.pt + 1.1*step.grad; bigz+30] , 'fontsize', fontsize, 'interpreter', 'latex');    
elseif step.plot_pt
    % presumes there is a pt and a label
    % then plot point
    plot3( [ step.pt(1)  ] , [ step.pt(2) ] , [bigz bigz]+20, 'o' , 'markerfacecolor', [1 1 1] , 'markersize', mksize , 'markeredgecolor', 'k');
    % labels
    text('string', sprintf('$$\\mathbf{%s}^{(%d)}_1$$',step.label,step.num) , 'position', [step.pt + -0.1*step.grad; bigz+30] , 'fontsize', fontsize, 'interpreter', 'latex');
end

    
if parms.save_fig
    % presumes a full pathed parms.fname_fig
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 2])
    print(gcf,'-dtiff',parms.fname_fig);
end
    

