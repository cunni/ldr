%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2014 John P. Cunningham
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% John P. Cunningham
% 
% test_project.m
%
% This function aggregates and calls all testing functionality.  Thus it 
% should be able to be called to call all relevant unit tests, which will
% ensure that the codebase is working properly (to the best of our 
% reasonable testing ability).  
%
% if the (optional) input show_fig is enabled, this will produce the
% summary figure probably shown in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = test_project( run_meth , test_name , show_fig , save_fig )

    %%%%%%%%%%
    % check inputs
    %%%%%%%%%%
    if nargin < 4 || isempty(save_fig)
        save_fig = 0;
    end
    if nargin < 3 || isempty(show_fig)
        show_fig = 1;
    end
    if nargin < 2 || isempty(test_name)
        test_name = 'd_sweep';
    end
    if nargin < 1 || isempty(run_meth)
        run_meth = 1;
    end
    
    %%%%%%%%%%
    % run methods
    %%%%%%%%%%
    if run_meth
        % now choose the data based on the test of interest
        switch lower(test_name)
            case 'perf_basic_2'
                % DEPRECATED... works but won't be plotted as previous
                % the 2d case
                d = [ [3:22] ];
                r = 2*ones(size(d));
                num_runs = length(d);
            case 'perf_basic_4'
                % DEPRECATED... works but won't be plotted as previous
                % the 4d case
                d = [ [5:24] ];
                r = 4*ones(size(d));
                num_runs = length(d);
            case 'd_sweep'
                % for sweeping d and fixed r
                num_repeats = 1;
                d = repmat([4 8 16 32 64 128 256 512 1024], num_repeats , 1 );
                d = d(:);
                r = 3*ones(size(d));
                num_runs = length(d);
            case 'r_sweep'
                % for sweeping r and fixed d
                num_repeats = 20;
                r = repmat([1 2 5 10 20 40 80], num_repeats , 1 );
                r = r(:);
                d = 100*ones(size(r));
                num_runs = length(d);
            case 'dr_sweep'
                % for sweeping d and fixed r
                num_repeats = 20;
                d1 = repmat([4 8 16 32 64 128 256 512 1024], num_repeats , 1 );
                d1 = d1(:);
                r1 = 3*ones(size(d1));
                % for sweeping r and fixed d
                num_repeats = 10;
                r2 = repmat([1 2 5 10 20 40 80], num_repeats , 1 );
                r2 = r2(:);
                d2 = 100*ones(size(r2));
                % now together
                r = [r1;r2];
                d = [d1;d2];
                num_runs = length(d);
            otherwise
                fprintf('this test not implemented here');
                keyboard
        end
        % now run the methods
        for i = 1 : num_runs
            fprintf('\n\n\n---------------Iter %d---------------\n\n\n',i);
            parms = struct('show_fig',0,'save_fig',0,'randseed',i);
            % test PCA
            [ r_pca(i) ] = test_pca( d(i) , r(i) , parms );
            % test LDA
            [ r_lda(i) ] = test_lda( d(i) , r(i) , parms );
            % test CCA
            [ r_cca(i) ] = test_cca( d(i) , d(i) , r(i) , parms );
            % test MAF
            [ r_maf(i) ] = test_maf( d(i) , r(i) , parms );
            % save...
            save(sprintf('results/test_project_%s.mat',test_name));
        end
       
        % if dr_sweep, split here into d_sweep and r_sweep for posterity.
        if isequal(test_name,'dr_sweep')
            % then use d1 and d2 to split it properly and save accordingly.
            d_full = d;
            r_full = r;
            num_runs_full = num_runs;
            r_lda_full = r_lda;
            r_maf_full = r_maf;
            r_pca_full = r_pca;
            r_cca_full = r_cca;
            % first save d_sweep
            d = d1;
            r = r1;
            r_lda = r_lda_full(1:length(d1));
            r_pca = r_pca_full(1:length(d1));
            r_maf = r_maf_full(1:length(d1));
            r_cca = r_cca_full(1:length(d1));
            num_runs = length(d1);
            save('results/test_project_d_sweep.mat','d','d1','r','r1','r_lda','r_pca','r_maf','r_cca','num_runs','num_repeats');
            % save r_sweep
            d = d2;
            r = r2;
            r_lda = r_lda_full(length(d1)+1:end);
            r_pca = r_pca_full(length(d1)+1:end);
            r_maf = r_maf_full(length(d1)+1:end);
            r_cca = r_cca_full(length(d1)+1:end);
            num_runs = length(d2);
            save('results/test_project_r_sweep.mat','d','d2','r','r2','r_lda','r_pca','r_maf','r_cca','num_runs','num_repeats');
        end
        
    %%%%%%%%%
    % load existing
    %%%%%%%%%
    else
        load(sprintf('results/test_project_%s.mat',test_name));
    end
    
    
    %%%%%%%%%
    % plot results
    %%%%%%%%%
    if show_fig 
        % now choose the plot based on the test of interest
        switch lower(test_name)
            case {'perf_basic_2','perf_basic_4'}

                % DEPRECATED
                % careful, the data structures have changed shapes...
                % so this code is deprecated
                
                method_names = { 'PCA' ; 'LDA' ; 'CCA' ; 'MAF' ; 'SDA'};
                R = [ [r_pca.normed_diff]' [r_lda.normed_diff]' [r_cca.normed_diff]' [r_maf.normed_diff]' [r_sda.normed_diff]' ];

                % histogram of performance
                figure;
                hold on;
                Rimp = -R;
                set(gca,'linewidth',2,'fontsize',22);
                set(gca,'ytick',[0.0 0.1 0.2 0.3 0.4 0.5 1.0]);
                set(gca,'ylim',[-0.02 0.42]);
                ylabel(sprintf('Normalized improvement'));
                xpoints = [5 : 10 : 5 + 10*(size(Rimp,2)-1) ];
                set(gca,'xtick',xpoints,'xticklabel', method_names );
                %axis off
                pctpts = prctile( Rimp , [2.5 50 97.5] , 1 );
                %
                dot_pt = median(Rimp,1);
                dpcheck = dot_pt - pctpts(2,:)
                topBar = pctpts(3,:);
                botBar = pctpts(1,:);
                %
                plot( xpoints , dot_pt , 'k.', 'markersize',44);
                h= errorbar( xpoints , dot_pt , botBar - dot_pt , topBar - dot_pt , 'k.', 'linewidth',2);
                errorbar_tick(h,10);
                % plot the 0 line
                plot([0 50],[0 0 ] , 'k--','linewidth',1)
                %
                if save_fig
                    print(gcf , '-depsc', sprintf('results/test_project_%s.eps',test_name))
                end
                
            case {'d_sweep','r_sweep'}
                % for sweeping d and fixed r or vice versa...
                % plot improvement at convergence, runtime of method, and number of iterations
                % plot those in distribution.
                
                % first make the relevant distributions
                if isequal(test_name,'d_sweep')
                    [pts,pts_first,~] = unique(d);
                else
                    [pts,pts_first,~] = unique(r);
                end
                r_summary = struct('improvement',[],'time',[],'iter',[]);
                r_summary.inds = pts;
                r_summary.test_name = test_name;
                pct_lims = [25 50 75]; % percentile 
                
                r_all = {r_pca , r_lda , r_cca , r_maf };
                % choice of optimization methods
                opt_name = {'grassmann_trust', 'grassmann_trust' , 'stiefel_trust_prod', 'grassmann_trust' };
                   
                for i = 1 :length(pts)
                    % the relevant indices for this choice of d (or r)
                    inds = [pts_first(i):pts_first(i)+num_repeats-1];
                    
                    for j = 1 : length(r_all)
                        r_this = r_all{j};
                   
                        % runtime
                        rt = [r_this(inds).time];
                        tmp_pts = [];
                        for k = 1 : length(rt)
                            tmp_pts(k) = getfield(rt(k),opt_name{j});
                        end
                        pctpts = prctile( tmp_pts , pct_lims );
                        r_summary(j).time(i).median = pctpts(2);
                        r_summary(j).time(i).low = pctpts(1);
                        r_summary(j).time(i).high = pctpts(3);
                        
                        % iter
                        rt = [r_this(inds).iter];
                        tmp_pts = [];
                        for k = 1 : length(rt)
                            tmp_pts(k) = getfield(rt(k),opt_name{j});
                        end
                        pctpts = prctile( tmp_pts , pct_lims );
                        r_summary(j).iter(i).median = pctpts(2);
                        r_summary(j).iter(i).low = pctpts(1);
                        r_summary(j).iter(i).high = pctpts(3);
                        
                        % normed_diff
                        rt = [r_this(inds).normed_diff];
                        pctpts = prctile( rt(2,:) , pct_lims );
                        r_summary(j).improvement(i).median = pctpts(2);
                        r_summary(j).improvement(i).low = pctpts(1);
                        r_summary(j).improvement(i).high = pctpts(3);
                        
                    end

                end
                                    
                % figure preliminaries
                lw=2;
                wt = 0;
                z = 255;
                colR = [204 204 153 153 204 204]/z - wt/255;
                colG = [153 178 204 178 153 178]/z - wt/255;
                colB = [204 153 153 178 153 046]/z - wt/255;
                col = [colR; colG; colB];
                col = col(:,[1 3 4 5 2 6]); % reorder for pref.
                r_xtick = [1 10 100];
                r_xlim = [0.9 110];
                d_xtick = [10 100 1000];
                d_xlim = [3.5 1200];
                
                %%%%%%%%%
                % now time figure
                %%%%%%%%%
                figure;
                for i = 1 : length(r_summary)
                    loglog( gca ,  r_summary(1).inds' , [r_summary(i).time.median] , 'linewidth', lw , 'color' , col(:,i) );
                    hold on;
                end
                set(gca,'linewidth',2,'fontsize',22);
                % x properties
                if isequal(test_name,'d_sweep')
                    xlabel(sprintf('data dimensionality (d)'));
                    set(gca,'xtick',d_xtick);
                    set(gca,'xlim',d_xlim);
                else
                    xlabel(sprintf('projected dimensionality (r)'));
                    set(gca,'xtick',r_xtick);
                    set(gca,'xlim',r_xlim);
                end
                % put numbers in non exp format
                new_XTickLabel = get(gca,'xtick');
                set(gca,'XTickLabel',new_XTickLabel);
                % y properties
                ylabel(sprintf('runtime (seconds)'));
                set(gca,'ytick',[10.^[-3:1:2]]);
                set(gca,'ylim',[1e-2 480]);
                if isequal(test_name,'d_sweep')
                    legend( gca , 'PCA' , 'LDA' , 'CCA', 'MAF' ,'Location','NorthWest');
                end
                for i = 1 : length(r_summary)
                    h= errorbar( r_summary(1).inds' , [r_summary(i).time.median] , [r_summary(i).time.low] - [r_summary(i).time.median] , [r_summary(i).time.high] - [r_summary(i).time.median] , 'linewidth', lw , 'color' , col(:,i) );
                end
                %
                if save_fig
                    print(gcf , '-depsc', sprintf('results/test_project_time_%s.eps',test_name))
                end

                %%%%%%%%%
                % now iter figure
                %%%%%%%%%
                figure;
                for i = 1 : length(r_summary)
                    semilogx( gca ,  r_summary(1).inds' , [r_summary(i).iter.median] , 'linewidth', lw , 'color' , col(:,i) );
                    hold on;
                end
                set(gca,'linewidth',2,'fontsize',22);
                % x properties
                if isequal(test_name,'d_sweep')
                    xlabel(sprintf('data dimensionality (d)'));
                    set(gca,'xtick',d_xtick);
                    set(gca,'xlim',d_xlim);
                else
                    xlabel(sprintf('projected dimensionality (r)'));
                    set(gca,'xtick',r_xtick);
                    set(gca,'xlim',r_xlim);
                end
                % put numbers in non exp format
                new_XTickLabel = get(gca,'xtick');
                set(gca,'XTickLabel',new_XTickLabel);
                % y properties
                ylabel(sprintf('number of iterations'));
                set(gca,'ytick',[5 10 15 20 25]);
                set(gca,'ylim',[4 28]);
                if isequal(test_name,'d_sweep')
                    legend( gca , 'PCA' , 'LDA' , 'CCA', 'MAF' ,'Location','NorthWest');
                end
                
                for i = 1 : length(r_summary)
                    h= errorbar( r_summary(1).inds' , [r_summary(i).iter.median] , [r_summary(i).iter.low] - [r_summary(i).iter.median] , [r_summary(i).iter.high] - [r_summary(i).iter.median] , 'linewidth', lw , 'color' , col(:,i) );
                end
                %
                if save_fig
                    print(gcf , '-depsc', sprintf('results/test_project_iter_%s.eps',test_name))
                end
                
                %%%%%%%%%
                % now improvement figure
                %%%%%%%%%
                figure;
                for i = 1 : length(r_summary)
                    semilogx( gca ,  r_summary(1).inds' , [r_summary(i).improvement.median] , 'linewidth', lw , 'color' , col(:,i) );
                    hold on;
                end
                set(gca,'linewidth',2,'fontsize',22);
                % x properties
                if isequal(test_name,'d_sweep')
                    xlabel(sprintf('data dimensionality (d)'));
                    set(gca,'xtick',d_xtick);
                    set(gca,'xlim',d_xlim);
                else
                    xlabel(sprintf('projected dimensionality (r)'));
                    set(gca,'xtick',r_xtick);
                    set(gca,'xlim',r_xlim);
                end
                % put numbers in non exp format
                new_XTickLabel = get(gca,'xtick');
                set(gca,'XTickLabel',new_XTickLabel);
                % y properties
                ylabel(sprintf('normalized improvement'));
                set(gca,'ytick',[0 .05 .10 .15 .2 .25 .3 .4 .5]);
                set(gca,'ylim',[-0.02 0.28]);
                if isequal(test_name,'d_sweep')
                    legend( gca , 'PCA' , 'LDA' , 'CCA', 'MAF' ,'Location','NorthWest');
                end
                for i = 1 : length(r_summary)
                    h= errorbar( r_summary(1).inds' , [r_summary(i).improvement.median] , [r_summary(i).improvement.low] - [r_summary(i).improvement.median] , [r_summary(i).improvement.high] - [r_summary(i).improvement.median] , 'linewidth', lw , 'color' , col(:,i) );
                end
                %
                % plot the 0 line
                semilogx([0.0001 5000],[0 0 ] , 'k--','linewidth',1)
                
                if save_fig
                    print(gcf , '-depsc', sprintf('results/test_project_improvement_%s.eps',test_name))
                end
 
                
            case 'dr_sweep'
                fprintf('if you ran a dr sweep, those first need to be separated...\n');
                fprintf('This should have been done above... Please call the plot on r_sweep and d_sweep individually.\n');
                keyboard;

            otherwise
                fprintf('this test not implemented here');
                keyboard
        end
        
        

    end
        
end

    
    
     