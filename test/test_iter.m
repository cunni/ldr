%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2014 John P. Cunningham
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P. Cunningham
% 
% test_iter.m
%
% tests the iteration count of several optimization methods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = test_iter( run_meth , show_fig , save_fig , num_runs )

    %%%%%%%%%%
    % check inputs
    %%%%%%%%%%
    if nargin < 4 || isempty(num_runs)
        num_runs = 10;
    end
    if nargin < 3 || isempty(save_fig)
        save_fig = 0;
    end
    if nargin < 2 || isempty(show_fig)
        show_fig = 1;
    end
    if nargin < 1 || isempty(run_meth)
        run_meth = 0;
    end
    
    
    %%%%%%%%%%
    % run methods
    %%%%%%%%%%
    if run_meth
        
        d = 100;
        r = 10;
        
        for i = 1 : num_runs
            
            parms = struct('show_fig',0,'save_fig',0,'randseed',i);
            % test PCA
            [ result(i) ] = test_pca( d , r , parms );
            
        end
        
        
        %%%%%%%%%%%
        % compile averages
        %%%%%%%%%%%
        % times
        rt = [result.time];
        fnrt = fieldnames(rt);
        meantimes = struct();
        for i = 1 : length(fnrt)
            meantimes = setfield(meantimes, fnrt{i} , mean([rt.(fnrt{i})]) );
        end
        % iterations
        rt = [result.iter];
        fnrt = fieldnames(rt);
        meaniters = struct();
        for i = 1 : length(fnrt)
            meaniters = setfield(meaniters, fnrt{i} , mean([rt.(fnrt{i})]) );
        end
        
        save('results/test_convergence_pca.mat', 'd' , 'r' , 'result' , 'meantimes' , 'meaniters');
    
    %%%%%%%%%%%%
    % load existing runs
    %%%%%%%%%%%%
    else
        load('results/test_convergence_pca.mat');
    
    end
    
    %%%%%%%%%%%%
    % pull average for each method
    %%%%%%%%%%%%
    % setup
    % skip heuristic... but iterate over other opt algos.
    optim_avg = struct();
    for i = 2 : length(result(1).optim_order)
        % stringname
        optim = result(1).optim_order{i};
        % now find max number of iterations associated with each method
        max_iter_optim = 0;
        for j = 1 : length(result)
            max_iter_optim = max( max_iter_optim , length(getfield( result(j).info , optim )) );
        end
        % now create the optim_avg struct
        optim_avg(i).name = optim;
        optim_avg(i).fmfstar = zeros( max_iter_optim , 1 );
        optim_avg(i).time = zeros( max_iter_optim , 1 );
        optim_avg(i).normalizer = zeros( max_iter_optim , 1 );
        optim_avg(i).timescatter = [];
        optim_avg(i).gapscatter = [];
        total_times = [];
        % now fill it
        for j = 1 : length(result)
            this_length = length( getfield( result(j).info , optim ) ) ;
            this_info = getfield( result(j).info , optim);
            ti_cost = [this_info.cost]';
            optim_avg(i).fmfstar(1:this_length) = optim_avg(i).fmfstar(1:this_length) + max( eps , ti_cost(1:end) - result(j).f.heuristic );
            optim_avg(i).normalizer(1:this_length) = optim_avg(i).normalizer(1:this_length) + 1;
            % now normalize for the empirical avg
            optim_avg(i).fmfstar = optim_avg(i).fmfstar ./ optim_avg(i).normalizer;
            total_times(j) = getfield(result(j).time,optim);
        end
        % find median index
        [stt, tt_ind] = sort(total_times);
        exemplar = tt_ind( floor(length(total_times)/2));
        this_info = getfield(result(exemplar).info , optim );
        optim_avg(i).time = [this_info.time]';
        optim_avg(i).fmfstar_time = max(eps, [this_info.cost]' - result(exemplar).f.heuristic);
    end
    
    
    %%%%%%%%%%%%
    % plot
    %%%%%%%%%%%%
    if show_fig
        
        lw=2;
        wt = 0;
        z = 255;
        colR = [204 204 153 153 204 204]/z - wt/255;
        colG = [153 178 204 178 153 178]/z - wt/255;
        colB = [204 153 153 178 153 046]/z - wt/255;
        col = [colR; colG; colB];
        col = col(:,[1 3 4 5 2 6]); % reorder for pref.
        
        figure;
        for i = 3 : length(optim_avg)
            semilogy( gca ,  optim_avg(i).time , optim_avg(i).fmfstar_time , 'linewidth', lw , 'color' , col(:,i-2) );
        hold on;
        end
        set(gca,'linewidth',2,'fontsize',22);
        ylabel(sprintf('optimality gap'));
        xlabel(sprintf('time (s)'));
        set(gca,'ytick',10.^([-10:2:4]) );
        set(gca,'ylim',[1e-11 1e4]);
        set(gca,'xtick',[0 : .25 : 1.3 ]);
        set(gca,'xlim',[0 1.32]);
        legend( gca , 'Stiefel SD' , 'Stiefel TR' , 'Grassmann SD' , 'Grassmann TR' );

        figure;
        for i = 3 : length(optim_avg)
            semilogy( gca ,  [1:length(optim_avg(i).fmfstar)]' , optim_avg(i).fmfstar , 'linewidth', lw , 'color' , col(:,i-2) );
            hold on;
        end
        set(gca,'linewidth',2,'fontsize',22);
        ylabel(sprintf('average optimality gap'));
        xlabel(sprintf('iteration'));
        set(gca,'ytick',10.^([-8:2:2]) );
        set(gca,'ylim',[1e-9 1e3]);
        set(gca,'xtick', [1 [20:20:100] ]);
        set(gca,'xlim',[1 115]);
        legend( gca , 'Stiefel SD' , 'Stiefel TR' , 'Grassmann SD' , 'Grassmann TR' );

    
    end

end

    
    
     