%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2014 John P. Cunningham
% a modification of bklinesearchG.m, from Gamaleldin F. Elsayed (with
% permission and thanks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P. Cunningham and Gamaleldin F. Elsayed
% 
% linesearch_geodesic.m
%
% function that determines the next step size of the algorithm using armijo
% track line search along geodesic paths in the Stiefel manifold. 
% To calculate the step size the function needs the following inputs:
% -The function handle to evaluate the objective (incl. its arguments)
% -The current algorithm point M (a matrix of orthonormal columns, namely 
% a point on the Stiefel manifold)
% -the projection of the gradient onto the tangent space of M
% -current step size
% -back track line search parameter beta.
%
% Inputs:
%    - f: function handle
%    - M: current Stiefel point matrix 
%    - g: change direction.
%    - stp: current step size
%    - alph: comparing parameter (usually related to movement direction)
% Outputs: 
%    - M_k: next geodesic step.
%    - stp: stp size used to generate the next geodesic step. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ M_k, stp ] = linesearch_geodesic( f , M , g , stp , alph )

% stp=1; % works best in terms of accuracy but slower (CG is more
% vulnerable)
fQ=fn_eval(X_dot,X, Q,D);
P_k=expm(stp*g);
W_k=P_k*P_k; 
Q_k=W_k*Q;

% % if (fQ-fn_eval(X_dot,X, P_k*Q,D)<stp./2*alph)
% %    stp=stp/2^2; 
% % end

while (fQ-fn_eval(X_dot,X, Q_k,D)>=stp*alph && stp<100)
    P_k=W_k;
    W_k=P_k*P_k;
    Q_k=W_k*Q;
    stp=2*stp;
end
Q_k=P_k*Q;

while (fQ-fn_eval(X_dot,X, Q_k,D)<stp./2*alph && stp>=eps)
   stp=stp/2;
   P_k=expm(stp*g);
   Q_k=P_k*Q;
end
% Q_k=expm(stp*g)*Q;
end