% test_plot_data_3d

rng(0);

A = randn(3);
X = A*randn(3,100);

W = project_stiefel(randn(3,2));



data_colors = [];

data_sizes = [];

% sort X and assign colors
X = sortrows(X')';
data_colors = [ linspace(0.1,1,size(X,2)); zeros(1,size(X,2)); zeros(1,size(X,2)) ];

f = plot_data_3d( X , W , data_colors , data_sizes )



% sort X and assign colors
X = sortrows(X')';
data_colors = [ zeros(1,size(X,2)); zeros(1,size(X,2)) ; linspace(0.1,1,size(X,2)) ];

f = plot_data_3d( X , W , data_colors , data_sizes )
