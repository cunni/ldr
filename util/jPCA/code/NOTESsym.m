%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2010
% 
% notes for testing symPCA.
%
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first, a simple test to show that the algorithm works quickly and is
% different than jPCA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

fprintf('\n\n======Small, simple example=========\n\n');

X0 = randn(100,4); % 100 is ct, 6 is k.
dX = X0(2:end,:) - X0(1:end-1,:);
X = X0(1:end-1,:);

% do jPCA1
fprintf('_____________symPCA1___________\n');
Mpre = X\dX;
tic; M1 = 0.5*( Mpre + Mpre' ); toc;
M1

% do jPCA2
fprintf('_____________symPCA2___________\n');
tic; M2 = symRegress(dX,X); toc;
M2

fprintf('====================================\n');
fprintf('\n\n(Press a key to continue...)\n\n')
pause



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, more data... and they seem to increasingly produce the same answer...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

fprintf('\n\n===More data, simple example=====\n\n');

X0 = randn(10000,4); % 100 is ct, 6 is k.
dX = X0(2:end,:) - X0(1:end-1,:);
X = X0(1:end-1,:);

% do jPCA1
fprintf('_____________symPCA1___________\n');
Mpre = X\dX;
tic; M1 = 0.5*( Mpre + Mpre' ); toc;
M1

% do jPCA2
fprintf('_____________symPCA2___________\n');
tic; M2 = symRegress(dX,X); toc;
M2

fprintf('====================================\n');
fprintf('\n\n(Press a key to continue...)\n\n')
pause



