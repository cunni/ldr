%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2010
%
% reshapeSkew.m
%
% jPCA2 support function
%
% A skew-symmetric matrix xMat of size n by n really only
% has n*(n-1)/2 unique entries.  That is, the diagonal is 0, and the
% upper/lower triangle is the negative transpose of the lower/upper.  So,
% we can just think of such a matrix as a vector x of size n(n-1)/2.  This
% function reshapes such a vector into the appropriate skew-symmetric
% matrix.  
% 
% The required ordering in x is row-minor, namely that xMat(1,1) = 0,
% xMat(2,1) = x(1), xMat(3,1) = x(2), and so on.
%
% This function goes either from vector to matrix or vice versa, depending
% on what the argument x is.
%
% In short, this function just reindexes a vector to a matrix or vice
% versa.
%%%%%%%%%%%%%%%%%%%%%%%%


function [ Z ] = reshapeSkew( x )

    % this reshapes a n(n-1)/2 vector to a n by n skew symmetric matrix, or vice versa.
    % First we must check if x is a matrix or a vector.
    if isvector(x)
        % then we are making a matrix
        
        % first get the size of the appropriate matrix
        % this should be n(n-1)/2 entries.
        % this is the positive root
        n = (1 + sqrt(1 + 8*length(x)))/2;
        % error check
        if n~=round(n) % if not an integer
            % this is a bad argument
            fprintf('ERROR... the size of the x vector prevents it from being shaped into a skew symmetric matrix.\n');
            keyboard;
        end
        
        % now make the matrix
        % initialize the return matrix
        Z = zeros(n);
        % and the marker index
        indMark = 1;
        
        for i = 1 : n-1
            % add the elements as appropriate.
            Z(i+1:end,i) = x(indMark:indMark+(n-i)-1);
            % now update the index Marker
            indMark = indMark + (n-i);
        end
        
        % now add the skew symmetric part
        Z = Z - Z';
        
    else
        % then we are making a vector from a matrix (note that the 
        % standard convention of lower case being a vector and upper case
        % being a matrix is now reversed).
        
        % first check that everything is appropriately sized and skew
        % symmetric
        if size(x) ~= size(x')
            % this is not symmetric
            fprintf('ERROR... the matrix x is not square, let alone skew-symmetric.\n');
            keyboard;
        end
        % now check for skew symmetry
        if abs(norm(x + x'))>1e-8
            % this is not skew symmetric.
            fprintf('ERROR... the matrix x is not skew-symmetric.\n');
            keyboard;
        end
        % everything is ok, so take the size
        n = size(x,1);
       
        
        % now make the vector Z
        indMark = 1;
        for i = 1 : n-1
            % add the elements into a column vector as appropriate.
            Z(indMark:indMark+(n-i)-1,1) = x(i+1:end,i);
            % now update the index Marker
            indMark = indMark + (n-i);
        end
        
    end
                
end