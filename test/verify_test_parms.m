%%%%%%%%%%%%%%%%%%
% John P. Cunningham
% 2014
%
% verify_test_parms.m
%
% verifies testing parameters
%
%%%%%%%%%%%%%%%%%%

function parms = verify_test_parms( d , r , parms )

    % key input args
    if nargin < 2 || isempty(r)
        % unacceptable
        error('ERROR: input r is required (r < d).');
    end
    % sensible input args
    if nargin < 1 || isempty(d)
        % unacceptable
        error('ERROR: input d is required (with r < d).');
    end
    % check sizes
    if d <= r
        % unacceptable
        error('ERROR: r must be less than d.');
    end
    if nargin < 3 || isempty(parms)
        % figure options
        parms.show_fig = 1;
        parms.save_fig = 0;
        parms.randseed = 0; %'shuffle';
    else
        if ~isfield(parms,'randseed');
            parms.randseed = 'shuffle';
        end
    end
