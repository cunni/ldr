%%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% sparse_pca_L1CT.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the sparse_PCA objective with smooth_L1CT approximation to the
% L1 penalty (as in Schmidt et al 2009, section 3.3)
function [ f , Df ] = sparse_pca_L1CT( Q , A , lambda, alpha)

    reshape_out = 0;
    % reshape as necessary
    if size(Q,2) ==1
        Q = reshape(Q,size(A,1),length(Q)/size(A,1));
        reshape_out = 1;
    end
    
    % calculate the barrier terms.
    [ fa , Dfa ] = smooth_L1CT( Q , alpha );

    % the function
    f = 0.5*trace(Q'*A*Q) + lambda*fa ;
    
    % the grad in Q
    Df = A*Q + lambda*Dfa;
        
    if reshape_out ==1
        Df = Df(:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the smooth_L1CT approximation to the L1 penalty as in Schmidt
function [ fa , Dfa ] = smooth_L1CT( Q , alpha )
    
    % calculate log_one_plus_exp
    [ lope_minus , d_lope_minus ] = log_one_plus_exp(-alpha*Q);
    [ lope_plus , d_lope_plus ] = log_one_plus_exp(alpha*Q);
    
    % the function 
    fa = (1/alpha)*sum(sum( lope_minus + lope_plus ));
    
    % the grad in Q
    Dfa = ( d_lope_plus - d_lope_minus );
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

