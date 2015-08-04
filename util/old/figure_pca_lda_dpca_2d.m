close all;
clear all;

f = figure;
set(f, 'position', [200 200 1500 1000]);
lw = 1.5;
ax = [-6 5 -6 5];
edg = -5:5;

%---------------------
% Generate data points
nPts = 100;
randn('state', 2);
th = pi/6 + pi/2;
R = [cos(th) -sin(th); sin(th) cos(th)];
C = R * diag([1.5 4]) * R';
X = mvnrnd([0 0], C, nPts); % nPts x 2
X = bsxfun(@minus, X, mean(X)); % center data

%-------------------
% PCA
mu = mean(X);
[U, D] = eig(cov(X, 1));
u1 = -U(:,2);
u2 = -U(:,1);

subplot(2,3,1);
hold on;
PC1 = [[0 0]' 5*u1]; 
PC2 = [[0 0]' 5*u2]; 
plot(PC1(1,:), PC1(2,:), 'k', 'linewidth', lw);
plot(PC2(1,:), PC2(2,:), 'k', 'linewidth', lw);

for j = 1:100
  plot(X(j,1), X(j,2), 'o', 'markerfacecolor', 'k', 'markersize', ...
       8, 'markeredgecolor', 'k');
end

title('PCA', 'fontsize', 24);
plot_ax(lw);
axis equal;
axis(ax);
axis off;

subplot(2,3,4);
Xproj = X * u1;
N = histc(Xproj, edg);
b = bar(edg+0.5, N, 'stacked');
set(b, 'facecolor', [0 0 0]);
axis([ax(1) ax(2) -10 70]);
axis off;
text('string', '$$ s_1 $$', 'position', [0 -5], 'fontsize', 24, ...
     'interpreter', 'latex', 'horizontalalignment', 'center');

%-------------------
% LDA
clsMask = (X(:,1)<0);

cls(1).pts = X(clsMask, :);
cls(2).pts = X(~clsMask, :);

for i = 1:2
  cls(i).mu = mean(cls(i).pts);
  cls(i).C = cov(cls(i).pts, 1); 
end

Chat = (cls(1).C + cls(2).C) / 2;
gmn = (cls(1).mu + cls(2).mu) / 2;

% Center data based on grand mean
for i = 1:2
  cls(i).pts = bsxfun(@minus, cls(i).pts, gmn);
  cls(i).mu = cls(i).mu - gmn;
end

w = inv(Chat) * (cls(2).mu - cls(1).mu)';
w = w / norm(w);

subplot(2,3,2);
hold on;

cls(1).col = 0;
cls(2).col = 1;

for i = 1:2 
  for j = 1:size(cls(i).pts, 1) % this loop is needed to make dots opaque
    plot(cls(i).pts(j,1), cls(i).pts(j,2), 'o', 'markerfacecolor', ...
         cls(i).col * [1 1 1], 'markersize', 8, 'markeredgecolor', ...
         'k'); 
  end
end

plot([0 5*w(1)], [0 5*w(2)], 'k', 'linewidth', lw);
plot([0 -5*w(2)], [0 5*w(1)], 'k', 'linewidth', lw);

title('LDA', 'fontsize', 24);
plot_ax(lw);
axis equal;
axis(ax);
axis off;

subplot(2,3,5);
hold on;
for i = 1:2
  % Project points onto w
  cls(i).prj = cls(i).pts * w;
  cls(i).N = histc(cls(i).prj, edg);
end

b = bar(edg+0.5, [cls(1).N cls(2).N], 'stacked');
for i=1:2
  set(b(i), 'facecolor', cls(i).col * [1 1 1]);
end

axis([ax(1) ax(2) -10 70]);
axis off;
text('string', '$$ s_1 $$', 'position', [0 -5], 'fontsize', 24, ...
     'interpreter', 'latex', 'horizontalalignment', 'center');


%----------------
% dPCA

% dPCA should find 'th' as one axis, and orthogonal axis as the
% other axis
th = -0.4 * pi;
R = [cos(th) -sin(th); sin(th) cos(th)];
y = X * R';  % y indicates 2 attributes of each point x
% Make attributes discrete
[ct1, att1] = histc(y(:,1), [-100 -2 0 2 100]);
[ct2, att2] = histc(y(:,2), [-100 -1.5 0 1.5 100]);

subplot(2,3,3);
hold on;

col = [0 0.3 0.7 1];
sz  = [4 10 20 30];

for i = 1:nPts
  plot(X(i,1), X(i,2), 'o', 'markerfacecolor', col(att1(i)) * ...
       [1 1 1], 'markersize', sz(att2(i)), 'markeredgecolor', 'k');
end

% Apply dPCA
covsum1 = zeros(2);
covsum2 = zeros(2);
for i = 1:4
  idx1 = (att1 == i);
  covsum1 = covsum1 + cov(X(idx1,:), 1) * sum(idx1); 

  idx2 = (att2 == i);
  covsum2 = covsum2 + cov(X(idx2,:), 1) * sum(idx2); 
end
covsum1 = covsum1 / nPts;
covsum2 = covsum2 / nPts;

[V, D] = eig(covsum1 - covsum2);
%V = -V;

plot(V(1,1) * [0 5], V(2,1) * [0 5], 'k', 'linewidth', lw);
plot(V(1,2) * [0 5], V(2,2) * [0 5], 'k', 'linewidth', lw);

title('Demixed', 'fontsize', 24);
plot_ax(lw);
axis equal;
axis(ax);
axis off;

subplot(2,3,6);
hold on;

% Project points onto w
prj = X * V(:,1);

Nall = [];
for i = 1:4
  N = histc(prj(att1==i), edg);
  Nall = [Nall N];
end

b = bar(edg+0.5, Nall, 'stacked');
for i=1:4
  set(b(i), 'facecolor', col(i) * [1 1 1]);
end

axis([ax(1) ax(2) -10 70]);
axis off;
text('string', '$$ s_1 $$', 'position', [0 -5], 'fontsize', 24, ...
     'interpreter', 'latex', 'horizontalalignment', 'center');


% Check dPCA answer
if 0
  pts2 = [0 0 10; 0 10 0];
  th = -th;
  R = [cos(th) -sin(th); sin(th) cos(th)];
  y2 = R * pts2;
  plot(y2(1,1:2), y2(2,1:2));
  plot(y2(1,[1 3]), y2(2, [1 3]));
end

