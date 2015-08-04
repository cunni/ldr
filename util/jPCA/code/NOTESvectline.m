% run vectlineJPC
clear all
close all

saveFigs = 1;

% I do not really know what this syms stuff means but it is what I want for this
% figure at late notice.
syms x y;
F = [x y];

% create a matrix M;
M = 10*randn(2);
% a nice clean one...
%M = [5.2655 6.0014; -2.6025 5.9393];

% better for simulating data...
%M = [2.2655 6.0014; -2.6025 1.9393];
M = [0.2 0.3; 0.3 0.02] + [0 -0.5; 0.5 0];



Msym = 0.5*(M + M');
Mskew = 0.5*(M - M');

% plot the total field
figure;
vectlineJPC(F*M',[x y], [-1,1,-1,1]);
axis off

if saveFigs
    figname = sprintf('M');
    print(gcf,'-dpsc2',sprintf('../written/cosyne2011/cosyne2011poster/%s.ps',figname));
end

% plot the symmetric field
figure;
vectlineJPC(F*Msym',[x y], [-1,1,-1,1]);
axis off

if saveFigs
    figname = sprintf('Msym');
    print(gcf,'-dpsc2',sprintf('../written/cosyne2011/cosyne2011poster/%s.ps',figname));
end

% plot the skew-symmetric field
figure;
vectlineJPC(F*Mskew',[x y], [-1,1,-1,1]);
axis off

if saveFigs
    figname = sprintf('Mskew');
    print(gcf,'-dpsc2',sprintf('../written/cosyne2011/cosyne2011poster/%s.ps',figname));
end

%%%%%%%%%%%%%%%
%%%%% draw data from this matrix %%%%%%%%%%
%%%%%%%%%%%%%%%
numPoints = 40;
numTimes = 11;
% lets have them start out on the axis a la the data we see in jPCA (hence
% rand not randn)
X(:,:,1) = rand(2,numPoints);
for i = 2 : numTimes
    X(:,:,i) = (M + eye(2))*X(:,:,i-1) + 0.1*randn(2,numPoints);
end


%%%%%%
% PCA
%%%%%%%%
Xpca = reshape(X,size(X,1),size(X,2)*size(X,3));
%
figure
scatter(Xpca(1,:),Xpca(2,:),'b.','linewidth',3);
axis off
%
if saveFigs
    figname = sprintf('dataNoTime');
    print(gcf,'-dpsc2',sprintf('../written/cosyne2011/cosyne2011poster/%s.ps',figname));
end
% now fit and plot the covariance
mu = mean(Xpca,2);
Sigma = 1/size(Xpca,2)*((Xpca*Xpca')) - mu*mu';
% plot the fit covariance
figure
plotcov2JPC(mu,Sigma);
axis off
if saveFigs
    figname = sprintf('dataNoTimePCA');
    print(gcf,'-dpsc2',sprintf('../written/cosyne2011/cosyne2011poster/%s.ps',figname));
end

%%%%%%%%%%%%%%%%%    
% dynamical PCA
%%%%%%%%%%%%%%%%%%
figure
hold on
%rgb = [[1:numPoints]',[numPoints:-1:1]',0*[1:numPoints]']/numPoints;
hf = [[numPoints:-2:1]' ; zeros(length([1:2:numPoints]),1)]/numPoints;
rgb = [ hf flipud(hf) 0*hf];
for i = 1: numPoints
    Xthispoint = squeeze(X(:,i,:));
    plot(Xthispoint(1,:),Xthispoint(2,:),'linewidth',2,'color',rgb(i,:));
end
hold off;
axis off

if saveFigs
    figname = sprintf('dataTime');
    print(gcf,'-dpsc2',sprintf('../written/cosyne2011/cosyne2011poster/%s.ps',figname));
end

dX = X(:,:,2:end)-X(:,:,1:end-1);
Xm1 = X(:,:,1:end-1);
dXr = reshape(dX,size(dX,1),size(dX,2)*size(dX,3));
Xm1r = reshape(Xm1,size(Xm1,1),size(Xm1,2)*size(Xm1,3));
% now do the fit
Mfit = (Xm1r'\dXr')';

% plot the fit field
figure;
vectlineJPC(F*Mfit',[x y], [-1,1,-1,1]);
axis off

if saveFigs
    figname = sprintf('dataTimedynPCA');
    print(gcf,'-dpsc2',sprintf('../written/cosyne2011/cosyne2011poster/%s.ps',figname));
end


    