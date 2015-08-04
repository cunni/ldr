%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% get_real_basis.m
%
%%%%%%%%%%%%%%%%%%%%%%%%

function Q = get_real_basis( U , mat_name)
    % CCA and other methods can sometimes produce complex conjugate eigenvectors
    if nargin < 2 || isempty(mat_name)
        mat_name = 'U';
    end
    % initialize
    Q = U;
    % now check
    if any(any(~isreal(U)))
        fprintf('Complex %s found in heuristic method... correcting\n',mat_name);
        skip_next = 0;
        for i = 1 : size(U,2)
            if any(~isreal(U(:,i))) && ~skip_next
                % then this pair needs correcting... they will be
                % in complex conjugate pairs since Sx... is real
                q1 = 0.5*(U(:,i) + U(:,i+1));
                q2 = 0.5*sqrt(-1)*(U(:,i) - U(:,i+1));
                Q(:,i) = q1/norm(q1);
                Q(:,i+1) = q2/norm(q2);
                % now skip ahead
                skip_next = 1;
            else
                % undo the skip ahead
                skip_next = 0;
            end
        end
    end
