function plot_ax(lw)

% Helper function to plot axes

ogn = [-5 -5];
plot([ogn(1) ogn(1)+5], [ogn(2) ogn(2)], 'k', 'linewidth', lw);
plot([ogn(1) ogn(1)], [ogn(2) ogn(2)+5], 'k', 'linewidth', lw);

text('string', '$$ r_1 $$', 'position', ogn + [4 -1], 'fontsize', 24, ...
     'interpreter', 'latex');
text('string', '$$ r_2 $$', 'position', ogn + [-1.2 4], 'fontsize', 24, ...
     'interpreter', 'latex');

%text('string', '$$ s_1 $$', 'position', [4 1.8], 'fontsize', 24, ...
%     'interpreter', 'latex');
%text('string', '$$ s_2 $$', 'position', [-4 3.8], 'fontsize', 24, ...
%     'interpreter', 'latex');
