%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% f_lda.m
%
% function for a simple eigenvalue problem
% This evaluates the objective and gradient
% for testing basic LDA.
%
%%%%%%%%%%%%%%%%%%%%

function [ f , gradf ] = f_lda( Q , X , X_labels , Sb , Sw)

    % ideally, the summary matrices Sb and Sw will have been computed
    % externally.  if not, we do it here, but it is inefficient
    if nargin < 5 || isempty(Sb) || isempty(Sw)
        %%%%%%%%%
        % form Sb and Sw
        %%%%%%%%%
        % find the number of classes
        [ classes , ~ , class_ind ] = unique(X_labels);
        % find means
        mu = zeros(size(X,1), length(classes));
        mu_all = mean(X,2);
        for c = 1 : length(classes)
            % find means
            mu(:,c) = mean( X(:, X_labels==classes(c) ) , 2);
        end
        % loop over all data
        Sb = zeros(size(X,1));
        Sw = zeros(size(X,1));
        for n = 1 : size(X,2)
            % calc Sb
            Sb = Sb + ( mu(:,class_ind(n)) - mu_all )*( mu(:,class_ind(n)) - mu_all )';
            % calc Sw
            Sw = Sw + ( X(:,n) - mu(:,class_ind(n)) )*( X(:,n) - mu(:,class_ind(n)) )';
        end
    end
        
    %%%%%%%%%%
    % now eval
    %%%%%%%%%%
    % both singletons
    trB = trace(Q'*Sb*Q);
    trW = trace(Q'*Sw*Q);
    
    % the function to *minimize*
    f = - trB / trW ;
    
    % the grad in Q ... usual quotient rule
    gradf =  - ( (2*trW)*Sb*Q - (2*trB)*Sw*Q )./(trW^2);
    
    