%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% helper function to project a point onto the Stiefel manifold
%%%%%%%%%%%%%%%%%%%%%%%%

function Q = project_stiefel( Q )
    % this is a well known projection, see eg Manton or Fan and Hoffman 1955
    [U,S,V] = svd(Q,0);
    Q = U*V';
end

    
