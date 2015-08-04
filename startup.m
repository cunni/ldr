%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generic linear dimensionality reduction
% Copyright (C) 2014 John P. Cunningham (see full notice in README)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P. Cunningham, 2014
% This is a startup file which initializes Matlab's settings to
% allow for use of this package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path(sprintf('%s/util', pwd), path) % Adds the subdirectories of the root folder to the path, allowing us to call functions from them.
path(sprintf('%s/test', pwd), path)
path(sprintf('%s/results', pwd), path)
path(sprintf('%s/notes', pwd), path)
path(sprintf('%s/data', pwd), path)
% add manopt subdir
path(sprintf('%s/util/manopt', pwd), path)
cd util/manopt
importmanopt
cd ../..
% now to verify, we can follow with 
% run('util/manopt/checkinstall/basicexample.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code checks for the relevant MEX files (such as .mexa64
% or .mexglx, depending on the machine architecture), and it creates the
% mex file if it can not find the right one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the mex file if necessary.
% this template project has no MEX, but if relevant, there is good
% functionality in startup_checkMEX (see toeplitz project)
% startup_checkMEX('decomp', 'logdetToeplitzFastZohar') % Toeplitz log determinant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This simply clears all variables and closes all windows opened by Matlab.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all % Clears all variables.
close all % Closes all windows opened by Matlab.