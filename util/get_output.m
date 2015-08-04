%%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% this helper function pulls the output_num 'th output argument from 
% the function call.  Very useful for dealing with function handles and
% multiple argument functions.
% note that this version expects f to be a string of the function name or a
% function handle to that name: 'foo' or @foo.  Other choices will have
% unexpected behavior.
% Adapted from http://stackoverflow.com/questions/3096281/skipping-outputs-with-anonymous-function-in-matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ varargout ] = get_output( f , output_number , varargin )
    % create output structure; default 1.
    varargout = cell( max(output_number) , 1 );
    % call f
    if isstr(f)
        [varargout{:}] = feval( f , varargin{:} );
    else % presumes f is a function handle
        [varargout{:}] = f( varargin{:} );
    end
    % assign correct output
    varargout = varargout(output_number);
end
