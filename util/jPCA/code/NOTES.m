%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2010
% 
% notes for testing jPCA2 and the differences with jPCA original.
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
fprintf('_____________jPCA1___________\n');
Mpre = X\dX;
tic; M1 = 0.5*( Mpre - Mpre' ); toc;
M1

% do jPCA2
fprintf('_____________jPCA2___________\n');
tic; M2 = skewSymRegress(dX,X); toc;
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
fprintf('_____________jPCA1___________\n');
Mpre = X\dX;
tic; M1 = 0.5*( Mpre - Mpre' ); toc;
M1

% do jPCA2
fprintf('_____________jPCA2___________\n');
tic; M2 = skewSymRegress(dX,X); toc;
M2

fprintf('====================================\n');
fprintf('\n\n(Press a key to continue...)\n\n')
pause







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, huge data... and they seem to increasingly produce the same answer...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

fprintf('\n\n===Huge data, simple example=====\n\n');

X0 = randn(100000,4); % 100 is ct, 6 is k.
C = [ones(4,1) [zeros(1,3); eye(3)]];
X0 = X0*C;

dX = X0(2:end,:) - X0(1:end-1,:);
X = X0(1:end-1,:);

% do jPCA1
fprintf('_____________jPCA1___________\n');
Mpre = X\dX;
tic; M1 = 0.5*( Mpre - Mpre' ); toc;
M1

% do jPCA2
fprintf('_____________jPCA2___________\n');
tic; M2 = skewSymRegress(dX,X); toc;
M2

fprintf('====================================\n');
fprintf('\n\n(Press a key to continue...)\n\n')
pause






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is a bit more interesting of a toy data set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

fprintf('\n\n===More interesting data example=====\n\n');
fprintf('(There is just oscillation in the first two dimensions)\n');

X0(:,1) = sin([1:10000]'); % oscillatory
X0(:,2) = cos([1:10000]'); % oscillatory
X0(:,3) = log([1:10000]') + tan([1:10000]'); % just something not oscillatory
X0(:,4) = randn(10000,1);

dX = X0(2:end,:) - X0(1:end-1,:);
X = X0(1:end-1,:);

% do jPCA1
fprintf('_____________jPCA1___________\n');
Mpre = X\dX;
tic; M1 = 0.5*( Mpre - Mpre' ); toc;
M1

% do jPCA2
fprintf('_____________jPCA2___________\n');
tic; M2 = skewSymRegress(dX,X); toc;
M2

fprintf('NOTE how this example differentiates M1 and M2.  M2 seems to crush the non-oscillatory components much better than does jPCA1.\n');

fprintf('====================================\n');
fprintf('\n\n(Press a key to continue...)\n\n')
pause






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is a bit more interesting of a toy data set. huge data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

fprintf('\n\n===huge data, more interesting data example=====\n\n');
fprintf('(There is just oscillation in the first two dimensions)\n');

X0(:,1) = sin([1:100000]'); % oscillatory
X0(:,2) = cos([1:100000]'); % oscillatory
X0(:,3) = log([1:100000]') + tan([1:100000]'); % just something not oscillatory
X0(:,4) = randn(100000,1);

dX = X0(2:end,:) - X0(1:end-1,:);
X = X0(1:end-1,:);

% do jPCA1
fprintf('_____________jPCA1___________\n');
Mpre = X\dX;
tic; M1 = 0.5*( Mpre - Mpre' ); toc;
M1

% do jPCA2
fprintf('_____________jPCA2___________\n');
tic; M2 = skewSymRegress(dX,X); toc;
M2

fprintf('NOTE how this example differentiates M1 and M2 LESS with MORE data.\n');

fprintf('====================================\n');
fprintf('\n\n(Press a key to continue...)\n\n')
pause






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% same data, but projected to mixed dimensions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

fprintf('\n\n===More interesting data example=====\n\n');
fprintf('(There is just oscillation in the first two dimensions)\n');

X0(:,1) = sin([1:10000]'); % oscillatory
X0(:,2) = cos([1:10000]'); % oscillatory
X0(:,3) = log([1:10000]') + tan([1:10000]'); % just something not oscillatory
X0(:,4) = randn(10000,1);

% project/rotate
P = rand(4);
P = P./repmat(sum(P),4,1);
X0 = X0*rand(4);

dX = X0(2:end,:) - X0(1:end-1,:);
X = X0(1:end-1,:);

% do jPCA1
fprintf('_____________jPCA1___________\n');
Mpre = X\dX;
tic; M1 = 0.5*( Mpre - Mpre' ); toc;
M1

% do jPCA2
fprintf('_____________jPCA2___________\n');
tic; M2 = skewSymRegress(dX,X); toc;
M2

fprintf('NOTE here that the expansive big expanding part of original column 3 washes out a lot of the simpler rotations.\n');

fprintf('This is sort of interesting... eigenvalue decomps of the jPCA1 and jPCA2 approaches...\n\n')

[V1, D1] = eig(M1)

[V2, D2] = eig(M2)

fprintf('\n\n It seems that jPCA2 does isolate the two rotations better (smaller eigs, but they are more distinct from the trailing noise eigs).\n');

fprintf('====================================\n');
