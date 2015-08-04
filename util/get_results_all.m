%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2014
%
% get_results_all.m
%
% This code is a simple helper function to ensure all results get structed
% in the same way for comparability.  Modified to include multiple
% optimization algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ results ] = get_results_all( method , d , r , res_arg )

    % objective
    results.name = method;
    results.d = d;
    results.r = r;
    % set f and Q and info and so on
    results.f = [];
    results.Q = struct();
    results.info = struct();
    results.time = struct();
    results.iter = struct();
    for i = 1 : length(res_arg)
        results.f = setfield( results.f , res_arg(i).optim , res_arg(i).fQ );
        results.Q = setfield( results.Q , res_arg(i).optim , res_arg(i).Q );
        results.info = setfield( results.info , res_arg(i).optim , res_arg(i).info );
        results.time = setfield( results.time , res_arg(i).optim , max([res_arg(i).info.time]) );
        results.iter = setfield( results.iter , res_arg(i).optim , length([res_arg(i).info.iter]) );
    end
    % 
    
    for i = 1 : length(res_arg)
        results.optim_order{i} = res_arg(i).optim;
        for j = 1 : i-1 %: length(res_arg)
            % the singular values show how different the subspaces are...
            % will not work with CCA
            %results.sv{i,j} = svd([ res_arg(i).Q , res_arg(j).Q ]);
            %results.subspace(i,j) = subspace( res_arg(i).Q , res_arg(j).Q );
            % norm diff.
            if isequal(method,'CCA') || isequal(method,'MAF')
                results.normed_diff(i,j) = (res_arg(j).fQ - res_arg(i).fQ);
            else
                results.normed_diff(i,j) = (res_arg(j).fQ - res_arg(i).fQ)/abs(res_arg(i).fQ);
            end
            % print message
            fprintf('%s: found a minimum %2.4f smaller using %s (vs %s).\n', method, results.normed_diff(i,j) , res_arg(i).optim , res_arg(j).optim );
        end
    end