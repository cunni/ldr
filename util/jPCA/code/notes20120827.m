%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2012
%
% test constrained solvers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%
% Example 1:
% helical toy data to show difference between PCA and jPCA
%%%%%%%%%
rs = 8;6;
rand('seed',rs);
randn('state',rs);
saveFigs = 1;
%M = [0 -0.5; 0.5 0];
theta = pi/8;
M = [cos(theta) -sin(theta); sin(theta) cos(theta)];
M = M-eye(2);
M = [[M; [ 0 0]] [0;0;00.01]];

%%%%%%%%%%%%%%%
%%%%% draw data from this matrix %%%%%%%%%%
%%%%%%%%%%%%%%%
numPoints = 3;
numTimes = 80;
% lets have them start out on the axis a la the data we see in jPCA (hence
% rand not randn)
X(:,:,1) = rand(3,numPoints);
X(3,:,1) = X(3,:,1) + 1;
for i = 2 : numTimes
    X(:,:,i) = (M + eye(3))*X(:,:,i-1) + 0.05*randn(3,numPoints);
end

% dump the 3rd dimension and readd it linearly.
for j = 1: numPoints
    X(3,j,:) = [1:numTimes];
end

%%%%%%
% PCA
%%%%%%%%
Xpca = reshape(X,size(X,1),size(X,2)*size(X,3));
%
% now fit and plot the covariance
mu = mean(Xpca,2);
Sigma = 1/size(Xpca,2)*((Xpca*Xpca')) - mu*mu';
[Us,Vs] = eig(Sigma);
XpcaRed = Us(:,2:3)'*Xpca;

%%%%%%%
% unconstrained dynamical PCA
%%%%%%%
dX = X(:,:,2:end)-X(:,:,1:end-1);
Xm1 = X(:,:,1:end-1);
dXr = reshape(dX,size(dX,1),size(dX,2)*size(dX,3));
Xm1r = reshape(Xm1,size(Xm1,1),size(Xm1,2)*size(Xm1,3));
% now do the fit
Mun = (Xm1r'\dXr')';
[Uj,Vj] = eig(Mun);
UUj = Uj*Uj';
XdpcaRed = UUj(:,1:2)'*Xpca;
XdpcaRed = reshape(XdpcaRed,2,size(X,2),size(X,3));

%{
%%%%%%%
% skew symmetric dynamical PCA
%%%%%%%
dX = X(:,:,2:end)-X(:,:,1:end-1);
Xm1 = X(:,:,1:end-1);
dXr = reshape(dX,size(dX,1),size(dX,2)*size(dX,3));
Xm1r = reshape(Xm1,size(Xm1,1),size(Xm1,2)*size(Xm1,3));
% now do the fit
%Mfit = (Xm1r'\dXr')';
Mskew = skewSymRegress(dXr',Xm1r');
[Uj,Vj] = eig(Mskew);
UUj = Uj*Uj';
XjpcaRed = UUj(:,1:2)'*Xpca;
XjpcaRed = reshape(XjpcaRed,2,size(X,2),size(X,3));
%}

%%%%%%%
% symmetric dynamical PCA
%%%%%%%
dX = X(:,:,2:end)-X(:,:,1:end-1);
Xm1 = X(:,:,1:end-1);
dXr = reshape(dX,size(dX,1),size(dX,2)*size(dX,3));
Xm1r = reshape(Xm1,size(Xm1,1),size(Xm1,2)*size(Xm1,3));
% now do the fit
%Mfit = (Xm1r'\dXr')';
Msym = symRegress(dXr',Xm1r');
[Uj,Vj] = eig(Msym);
UUj = Uj*Uj';
XspcaRed = UUj(:,1:2)'*Xpca;
XspcaRed = reshape(XspcaRed,2,size(X,2),size(X,3));



