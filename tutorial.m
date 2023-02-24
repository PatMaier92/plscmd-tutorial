clear; close all; clc; format long;
addpath(genpath(pwd)) % add subfolder functions to path

%% Tutorial for Partial Least Square Correlation (PLSC) Analysis 
% requires the plscmd toolbox (http://www.rotman-baycrest.on.ca/pls).

%--------------------------------------------------------------------------
% LOAD DATA 
%--------------------------------------------------------------------------
load('data/demo_data.mat');
demo_data = cellfun(@str2num, m); clear m;

data_cell = { demo_data demo_data demo_data demo_data};
data_names = { 'demo_m3' 'demo_m3' 'demo_m5' 'demo_m5'}; 
analysis_method = { 3 3 5 5 }; 
n_conditions = { 1 1 1 2 }; 
by_group = { 0 1 1 1 }; 
design_contrasts = { [] [] ... % n rows = n conditions x n groups, sorted as condition in group; n colums =  n desired contrasts
    [1 -1]' [1 -1 1 -1; 1 -1 0 0; 0 0 1 -1; 1 1 -1 -1; 1 0 -1 0; 0 1 0 -1]' }; 

for i=1:numel(data_cell)
    %--------------------------------------------------------------------------
    % CONFIGURATION
    %-------------------------------------------------------------------------- 
    data = data_cell{i};
    file_name = data_names{i};
    
    % behavioral output data
    plsinput.y = data(:,4);
    
    % behavioral explanatory data
    plsinput.X = data(:,5:size(data,2));
    
    % z-standardization
    plsinput.y = zscore(plsinput.y,0,1);
    plsinput.X = zscore(plsinput.X,0,1);

    % cfg settings
    cfg.pls = [];
    cfg.pls.num_perm = 500; % number of permutations
    cfg.pls.num_boot = 500; % number of bootstrap tests
    cfg.pls.clim     = 95; % confidence interval level

    % analysis method: 3=regular behavior PLS; 5=non-rotated behavior PLS
    cfg.pls.method   = analysis_method{i}; 
    
    % add behavioral output data
    cfg.pls.stacked_behavdata = plsinput.y;
    
    % number of conditions 
    n_con = n_conditions{i}; 
    if cfg.pls.method==5
        cfg.pls.stacked_designdata=design_contrasts{i};
    end 

    % number of subjects and data preparation
    if by_group{i}==0
        n_subj = size(plsinput.y,1) / n_con;
        datamat1_allgroups = plsinput.X;
    elseif by_group{i}==1
        n_subj = histc(data(:,2),unique(data(:,2))) / n_con;
        datamat1_group1 = plsinput.X(1:n_subj(1),:);
        datamat1_group2 = plsinput.X(n_subj(1)+1:end,:);
    end   
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % RUN PLS
    %--------------------------------------------------------------------------
    % input arguments: data, number of subjects, number of conditions, specific settings
    if by_group{i}==0
        plsres = pls_analysis({ datamat1_allgroups }, n_subj, n_con, cfg.pls);
    elseif by_group{i}==1
        plsres = pls_analysis({ datamat1_group1, datamat1_group2 }, n_subj, n_con, cfg.pls);
    end
    
    % save([path, 'PLSC_full_results_', file_name, '_m', int2str(cfg.pls.method),'.mat'],'plsres');
    
    clear n_subj n_con cfg datamat*; 
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % SAVE OUTPUT
    % -------------------------------------------------------------------------
    for LV_n=1:numel(plsres.perm_result.sprob)
        
        % Latent variable significance (should be < than 0.05)
        p=plsres.perm_result.sprob(LV_n); 
        
        % Latent variable weights
        % Bootstrap ratios (BSR) (should be < -1.96 or > +1.96)
        % correlations with behavioral data and standard error 
        BSR = plsres.boot_result.compare_u(:,LV_n);
        if size(plsres.datamatcorrs_lst{1,1},1)==1 % TBD CHECK for different methods/conditions
            cor = plsres.datamatcorrs_lst{1,1}'; 
        else
            cor = plsres.datamatcorrs_lst{1,1}(LV_n,:)'; 
        end
        se = plsres.boot_result.u_se(:,LV_n); % u = plsres.u(:,LV_n); 
               
        % BSR plots 
        if plsres.method==3
            
            figure; subplot(1,2,1);
            n_dim=numel(cor);
            bar(cor(:)./se,'k'); hold on;
            set(gca,'xticklabels',{'behavior1','behavior2','behavior3','behavior4','behavior5', 'behavior6'});
            box off; grid on;
            lh = line([0,size(BSR,1)+1],[2,2]);
            set(lh, 'color','r','linestyle','--');
            lh = line([0,size(BSR,1)+1],[-2,-2]);
            set(lh, 'color','r','linestyle','--');
            ylim([-12 12]);
            title(['BSR for LV with p-value=', num2str(round(p,3))]);
            hold off;
            
            subplot(1,2,2);
            bar(cor(:),'k'); hold on;
            set(gca,'xticklabels',{'behavior1','behavior2','behavior3','behavior4','behavior5', 'behavior6'});            box off; grid on;
            for nd = 1:n_dim
                lh1 = line([nd,nd],[cor(nd)+1.96*se(nd),cor(nd)-1.96*se(nd)]);
                set(lh1, 'color','r');
            end
            ylim([-1 1]);
            title('Mean Correlation +- 1.96*SE');
            hold off;
            clear n_dim lh*; 

        elseif plsres.method==5 
       
            figure; subplot(1,2,1);
            bar(BSR(:),'k'); hold on;
            set(gca,'xticklabels',{'behavior1','behavior2','behavior3','behavior4','behavior5', 'behavior6'});
            box off; grid on;
            lh = line([0,size(BSR,1)+1],[2,2]);
            set(lh, 'color','r','linestyle','--');
            lh = line([0,size(BSR,1)+1],[-2,-2]);
            set(lh, 'color','r','linestyle','--');
            ylim([-12 12]);
            title(['BSR for LV with p-value=', num2str(round(p,3))]);
            hold off;
            
            subplot(1,2,2);
            n_dim = size(plsres.boot_result.orig_corr,1);
            bar(1:n_dim, plsres.boot_result.orig_corr(:,LV_n),'k'); hold on;
            box off; grid on; grid minor;
            for nd = 1:n_dim
                lh1 = line([nd,nd],[plsres.boot_result.llcorr_adj(nd,LV_n), plsres.boot_result.ulcorr_adj(nd,LV_n)]);
                set(lh1, 'color','r');
            end
            xlim([0,n_dim+1]); ylim([-1,1]);
            title('Correlation LV with memory');
            hold off;
            clear n_dim nd lh*; 
            
        end
        % saveas(p,'plot.jpeg')
        
    end 
    
    clear p BSR cor u se; 
    
    % Latent profile scores
    % indicates an individual's expression of the profile
    if plsres.method==3
        
        group = data(:,2);
        figure;
        gscatter(plsres.usc(:,LV_n), plsinput.y, group);
        xlabel(upper('LV profile score'),'fontweight','bold');
        ylabel(upper('outcome'),'fontweight','bold');
        [R,P]=corrcoef(plsres.usc(:,LV_n), plsinput.y, 'rows', 'complete');
        title(strcat('r=',num2str(R(2,1)),', p=',num2str(P(2,1))));
        clear R P group;
        % saveas(p,'plot.jpeg')
        
    end

    clear data file_name plsinput plsres; 
    %--------------------------------------------------------------------------
end

clear all; 