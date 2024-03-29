# plscmd-tutorial
Tutorial on the Partial Least Square Correlation (PLSC) with the plscmd toolbox (http://www.rotman-baycrest.on.ca/pls) in Matlab. 

## General information on PLSC
PLSC is a multivariate correlation method (Krishnan, Williams, McIntosh & Abdi, 2011) to analyze associations between two sets of indicators. The correlation matrix between the two sets is decomposed into three components using singular value decomposition (SVD), based on which one latent variable (LV) is extracted in a least-square sense. The LV represents a distinct profile of indicators that has the strongest relation to the outcome. One can also look at the LV weights for each indicator, representing the degree and direction (positive or negative association) to which they contribute to the LV. Last, but not least one can calculate a latent profile score (LPS) for each individual, by multiplying the LV weights with the original data matrix (higher values correspond to stronger expression of the profile by the individual) which can then be correlated with other variables of interest (example: correlation of age-related latent profile scores with memory to see if a more mature behavior results in better memory). 

We focus on behavior PLSC here i.e., correlations between behavior and brain data, or two types of behavioral data (for other variants like task, seed or multi-block PLSC, check Krishnan et al., 2011) and we explore data-driven and confirmatory versions of behavior PLSC. If you use brain data, you can use the plsgui toolbox, which is pretty well-documented on http://www.rotman-baycrest.on.ca/pls. However, if you want to use any other type of data, you will need to use the plscmd toolbox, which is not so well documented. This tutorial was inspired by code from Muelroth et al. (2020; available on https://osf.io/w76f3/) and Keresztes et al. (2017).

First, we need to download the plscmd package from http://www.rotman-baycrest.on.ca/pls and add it to our path in Matlab. 

```
addpath(genpath(pwd))

load(['data.mat']);
```

The data needs to be in format "subject in group in condition".

| id | group | condition | outcome| predictor1 | predictor2 | predictor3 | predictor4 | predictor5 | predictor6 |
| -- | ----- | --------- | ------ | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| 1  | 1     | 1         | 0.64   | 0.126      | 10.5       | 5          | 0.97       | 105.7      | 77.6       |
| 1  | 1     | 2         | 0.56   | 0.567      | 10.9       | 3          | 0.53       | 202.6      | 89.0       |
| 2  | 1     | 1         | 0.23   | 0.767      | 11.9       | 7          | 0.98       | 254.8      | 93.4       |
| 2  | 1     | 2         | 0.54   | 0.300      | 15.5       | 2          | 0.33       | 155.9      | 55.6       |
| 3  | 1     | 1         | 0.76   | 0.150      | 19.3       | 9          | 0.88       | 148.7      | 72.8       |
| 3  | 1     | 2         | 0.18   | 0.337      | 13.7       | 5          | 0.45       | 122.0      | 43.9       |
| 15 | 2     | 1         | 0.99   | 0.836      | 18.3       | 0          | 0.76       | 200.7      | 59.1       |
| 15 | 2     | 2         | 0.17   | 0.237      | 11.8       | 8          | 0.29       | 227.3      | 51.8       |
| 16 | 2     | 1         | 0.85   | 0.378      | 13.2       | 1          | 0.46       | 176.9      | 69.4       |
| 16 | 2     | 2         | 0.38   | 0.217      | 18.0       | 2          | 0.83       | 134.0      | 41.8       |
| 17 | 2     | 1         | 0.66   | 0.786      | 12.5       | 5          | 0.26       | 127.8      | 49.0       |
| 17 | 2     | 2         | 0.27   | 0.378      | 12.0       | 4          | 0.96       | 245.3      | 78.1       |

If you have missing values, you need to either exclude cases or use an imputation method to replace the missing values. 

## Regular behavior PLSC
This approach is data-driven. 

### Data preparation
```
% outcome behavioral data
plsinput.y = data(:,4);

% explanatory behavioral data or brain data 
plsinput.X = data(:,5:size(data,2));
    
% z-standardization
plsinput.y = zscore(plsinput.y,0,1);
plsinput.X = zscore(plsinput.X,0,1);
```

### Running the analysis
```
% configuration 
cfg.pls = [];
cfg.pls.num_perm 	= 500;  	% number of permutations
cfg.pls.num_boot 	= 500; 		% number of bootstrap tests
cfg.pls.clim     	= 95; 		% confidence interval level
cfg.pls.stacked_behavdata = plsinput.y;	% outcome data 
cfg.pls.method   	= 3; 		% regular behavior PLS
    
% number of conditions
n_con = 1; % condition
    
% number of subjects and data preparation
n_subj = size(plsinput.y,1) / n_con;
datamat1_allgroups = plsinput.X;
    
% run plsc 
% input arguments: data, number of subjects, number of conditions, specific settings
plsres = pls_analysis({ datamat1_allgroups }, n_subj, n_con, cfg.pls);

save(['PLSC_full_results.mat'], 'plsres');
```

### Interpreting the output 

#### Latent variable (LV) 
The LV represents a distinct profile of indicators that has the strongest relation to the outcome.
Significance is evaluated with permutations and here set to < 0.05.

```
% LV significance (only interpret LVs < 0.05)
p=plsres.perm_result.sprob; 

% here: LV 1 was significant
LV_n = 1;
```

#### Latent variable weights (LV weights) 
The LV weights represent the degree and direction to which they contribute to the LV.
Significance is evaluated with bootstrap ratios (BSR).

```
% Bootstrap ratios for LV weights (significant if < -1.96 or > +1.96)
BSR = plsres.boot_result.compare_u(:,LV_n);

figure; subplot(1,2,1);
n_dim=numel(cor);
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

% Alternatively one can look at the correlations with behavioral data and standard error
cor = plsres.datamatcorrs_lst{1,1}'; % u = plsres.u(:,LV_n); 
se = plsres.boot_result.u_se(:,LV_n); 

subplot(1,2,2);
bar(cor(:),'k'); hold on;
set(gca,'xticklabels',{'behavior1','behavior2','behavior3','behavior4','behavior5', 'behavior6'});            
box off; grid on;
for nd = 1:n_dim
	lh1 = line([nd,nd],[cor(nd)+1.96*se(nd),cor(nd)-1.96*se(nd)]);
	set(lh1, 'color','r');
end
ylim([-1 1]);
title('Mean Correlation +- 1.96*SE');
hold off;
clear n_dim lh*; 
```

#### Latent profile scores (LPS)
LPS indicate an individual's expression of the profile.

```
group = data(:,2);
figure;
gscatter(plsres.usc(:,LV_n), plsinput.y, group);
xlabel(upper('LV profile score'),'fontweight','bold');
ylabel(upper('outcome'),'fontweight','bold');
[R,P]=corrcoef(plsres.usc(:,LV_n), plsinput.y, 'rows', 'complete');
title(strcat('r=',num2str(R(2,1)),', p=',num2str(P(2,1))));
clear R P group; 
```

The LPS values can be used for further in-depth correlation analyses (not shown). 

## Non-rotated (constraint) behavior PLSC
This approach is confirmatory and uses pre-defined contrasts. 

### Data preparation
```
% outcome behavioral data
plsinput.y = data(:,4);

% explanatory behavioral data or brain data 
plsinput.X = data(:,5:size(data,2));
    
% z-standardization
plsinput.y = zscore(plsinput.y,0,1);
plsinput.X = zscore(plsinput.X,0,1);
```

### Running the analysis 
```
% configuration
% configuration 
cfg.pls = [];
cfg.pls.num_perm 	= 500;  	% number of permutations
cfg.pls.num_boot 	= 500; 		% number of bootstrap tests
cfg.pls.clim     	= 95; 		% confidence interval level
cfg.pls.stacked_behavdata = plsinput.y;	% outcome data 
cfg.pls.method   	= 5; 		% non-rotated behavior PLS
    
% number of conditions  
n_con = 2; % condition
cfg.pls.stacked_designdata=[1 -1 1 -1; 1 1 -1 -1]'; % contrasts 
% the data are subjects in group in condition and thus the first contrast compares condition 1 and condition 2 (across groups), whereas the second contrast compares group 1 and group 2 (across conditions). One could also more contrasts to test interactions effects. 
    
% number of subjects and data preparation
n_subj = histc(data(:,2),unique(data(:,2))) / n_con;
datamat1_group1 = plsinput.X(1:n_subj(1),:);
datamat1_group2 = plsinput.X(n_subj(1)+1:end,:);
    
% run plsc 
% input arguments: data, number of subjects, number of conditions, specific settings
plsres = pls_analysis({ datamat1_group1, datamat1_group2 }, n_subj, n_con, cfg.pls);

save(['PLSC_full_results.mat'], 'plsres');
```

### Interpreting the output 

#### Latent variable (LV) 
The LV represents a distinct profile of indicators that has the strongest relation to the outcome.
Significance is evaluated with permutations and here set to < 0.05.
Here, a significant LV mean that the corresponding pre-defined contrast was significant (e.g. group 1 and group 2 differ across condition 1 and condition 2).

```
% LV significance (only interpret LVs < 0.05)
p=plsres.perm_result.sprob; 

% here: contrast 2 ([1 1 -1 -1]) was significant
% this means that group 1 and group 2 differ across condition 1 and condition 2
LV_n = 2;
```

#### Latent variable weights (LV weights) 
The LV weights represent the degree and direction to which they contribute to the LV.
Significance is evaluated with bootstrap ratios (BSR).
Here, we can see which indicators contribute to the significant contrast. 

```
% Bootstrap ratios for LV weights (significant if < -1.96 or > +1.96)
BSR = plsres.boot_result.compare_u(:,LV_n);

figure; subplot(1,2,1);
n_dim=numel(cor);
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
```

#### Latent profile scores (LPS)
Don't make much sense for non-rotated (constraint) behavior PLSC. 
