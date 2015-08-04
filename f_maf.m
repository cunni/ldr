%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% f_maf.m
%
% function for a simple eigenvalue problem
% This evaluates the objective and gradient
% for testing basic MAF.
%
%%%%%%%%%%%%%%%%%%%%

function [ f , gradf ] = f_maf( Q , X , S , St )

    % ideally, the summary matrices S and St will have been computed
    % externally.  if not, we do it here, but it is inefficient
    if nargin < 4 || isempty(S) || isempty(St)
        %%%%%%%%%
        % form S and St
        %%%%%%%%%
        % note that there is a presumed temporal structure in this data, as is
        % a requirement for MAF
        % first mean center the data
        Xm = X - repmat(mean(X,2),1,size(X,2));
        % now form the covariance
        S = (1/size(Xm,2))*Xm*Xm';
        % now form the time lag covariance, and symmetrize
        St = (1/(size(Xm,2)-1))*Xm(:,2:end)*Xm(:,1:end-1)';
        St = 0.5*(St + St');
    end
    
    %%%%%%%%%%
    % now eval
    %%%%%%%%%%
    % both singletons
    trS = trace(Q'*S*Q);
    trSt = trace(Q'*St*Q);
    
    % the function to *minimize*
    f = - trSt / trS ;
    
    % the grad in Q ... usual quotient rule
    gradf =  - ( (2*trS)*St*Q - (2*trSt)*S*Q )./(trS^2);
    
    