%%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% skew dynamical PCA
%
% for Cunningham and Ghahramani 2014. 
% This calculates the objective and gradient of ||dX-XM|| over Q,
% where M = QDQ', for a specified D.
%
% This is for validation of the Stiefel method.  A problem specific method
% is also available and, though computationally better, may have worse optima.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ f , gradf ] = f_sda( Q , A , B , D )
    % reshape
    % Q = reshape(q,size(B,1),size(D,1));
    %{
    % D specifies the spectrum and must be 2x2 block diagonal, where the
    % blocks are each skew symmetric.
    if norm(D + D','fro') > 1e-5
        fprintf('ERROR: D is not block diag with skew-symm blocks.\n');
    end
    % should perform some error checking here to verify that Q is a stiefel point    
    if norm( Q'*Q - eye(size(D)) , 'fro' ) > 1e-4
        % oops
        fprintf('WARNING: initial point Q is not a Stiefel point. Correcting...\n');
        Q = project_stiefel( Q );
    end
    %}
    
    % evaluate objective 
    f = trace( (A - Q*D*Q'*B)'*(A - Q*D*Q'*B) ); %( A - (Q*D*Q')*B , 'fro' )^2;
    
    % IMPORTANT NOTE: Because we are here treating this objective in an
    % unconstrained fashion, we must take the unconstrained gradient of the
    % objective function.  This means that there is a quartic term in Q,
    % since Q'*Q is not necessarily I in the derivative.  In practice this
    % should not matter, but if you are checking gradients, it will.
    % WRONG: evaluate derivative when you know that Q'Q=I
    %Df = 2*( B*B'*Q*D'*D - A*B'*Q*D' - B*A'*Q*D );
    %Df = Df(:);
    % CORRECT: the quartic version
    gradf = -2*(  A*B'*Q*D' + B*A'*Q*D );
    gradf = gradf - (D*Q'*Q*D*Q'*B*B' + D'*Q'*B*B'*Q*D'*Q' + D*Q'*B*B'*Q*D*Q' + D'*Q'*Q*D'*Q'*B*B')';
    % gradf = gradf(:);
end
