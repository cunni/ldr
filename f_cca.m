%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% f_cca.m
%
% function for a simple eigenvalue problem
% This evaluates the objective and gradient
% for testing CCA.
%
% Important: please note this only takes the gradient wrt Qx... if you want
% the gradient wrt Qy, note that the function is symmetric, so just call
% the function with all X and Y reversed.
%%%%%%%%%%%%%%%%%%%%

function [ f , gradf ] = f_cca( Qx , Qy , X , Y , Sxx , Syy , Sxy )

    % ideally, the summary matrices Sxx etc will have been computed
    % externally.  if not, we do it here, but it is inefficient
    if nargin < 7 || isempty(Sxx) || isempty(Syy) || isempty(Sxy)
        %%%%%%%%%
        % Make CCA matrices, regardless of the method
        %%%%%%%%%
        Xm = X - repmat(mean(X,2),1,size(X,2));
        Ym = Y - repmat(mean(Y,2),1,size(Y,2));
        Sxy = Xm*Ym';
        Sxx = Xm*Xm';
        Syy = Ym*Ym';
        
    end
    
    
    %%%%%%%%%%
    % now eval
    %%%%%%%%%%
    % both singletons
    trXY = trace(Qx'*Sxy*Qy);
    trXX = trace(Qx'*Sxx*Qx);
    trYY = trace(Qy'*Syy*Qy);
        
    % the function to *minimize*
    f = - trXY / sqrt(trXX*trYY) ;
    
    % the grad in Qx ... usual quotient rule
    gradf = - ( trXX*Sxy*Qy - trXY*Sxx*Qx ) / ( trYY * trXX^(3/2) );
    
    