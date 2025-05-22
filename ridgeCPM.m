function [q_s, r_pearson, r_rank, y, coef_total, coef0_total,...
    all_test_idx] = ridgeCPM(all_edges, all_behav, varargin)
%ridgeCPM Connectome-based predictive modeling using univariate
%feature selection and ridge regression
%
%   [q_s, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total] = ridgeCPM(all_mats, all_vectors, all_behav, thresh, v_alpha, lambda, k, seed)
%
%   Input:      all_edges,          vector of all subjects and their edges
%                                   [edges x subjects]
%
%              
%               all_behav,          behavior of all the subjects
%                                   [subjects x 1]
%
%               thresh (optional),  feature selection threshold (p-value)
%
%               family (optional),  vector of family IDs
%
%
%               v_alpha(optional),  value of the alpha parameter in elastic
%                                   net, default is 1e-6 which makes the
%                                   regression method to be ridge
%                                   regression, v_alpha=1 makes it lasso.
%
%               lambda(optional),   value of the lambda, if not provided,
%                                   cross-validation will be used
%
%               k(optional),        number of folds in k-fold cross
%                                   validation, default is 10
%
%               seed(optional),     random seed, default is 665
%
%               run_cpm(optional),  whether or not to run cpm, default 1
%
%   Output:     q_s,                cross-validated R^2 between predicted
%                                   value with ground truth (all_behav)
%
%               r_pearson,          WRONG! direct pearson correlation
%                                   between predicted value with ground
%                                   truth, only kept here for comparison
%                                   will be removed afterwards
%
%               r_rank,             cross-validated spearman correlation
%                                   between predicted value with ground
%                                   truth (all_behav)
%
%               y,                  predicted value for all the subjects
%
%               coef_total,         regression coefficients of all the edges
%                                   in all the k folds
%
%               coef0_total,        regression intercept in all the k folds
%
%               lambda_total,       penalty parameter chosen at each
%                                   iteration
%
%   Siyuan Gao, Yale University, 2018-2019 
%   https://www.sciencedirect.com/science/article/pii/S1053811919306196?via%3Dihub#bib2
%   Adapted by Abby Greene, 2019
%   Adapted further by Matt Rosenblatt, 2022
%% parse inputs
p=inputParser;
addRequired(p,'all_edges',@isnumeric);
addRequired(p,'all_behav',@isnumeric);
addParameter(p,'thresh',0.05,@isnumeric);
addParameter(p,'family',NaN,@isnumeric);
addParameter(p,'v_alpha',1e-6,@isnumeric);
addParameter(p,'lambda',NaN,@isnumeric);
addParameter(p,'kfolds',10,@isnumeric);
addParameter(p,'seed',1,@isnumeric);
addParameter(p,'run_cpm',1,@isnumeric);

parse(p, all_edges, all_behav, varargin{:});

thresh = p.Results.thresh;
family = p.Results.family;
v_alpha = p.Results.v_alpha;
lambda = p.Results.lambda;
kfolds = p.Results.kfolds;
seed = p.Results.seed;
run_cpm = p.Results.run_cpm;

%% initialization

num_sub_total = size(all_edges,2);
num_edge = size(all_edges,1);
nfams = length(unique(family));

coef_total = zeros(num_edge, kfolds); %store all the coefficients
coef0_total = zeros(1, kfolds); % store all the intercept
lambda_total = zeros(1, kfolds); % store all the lambda

%% main
y = zeros(num_sub_total, 1);
all_test_idx = nan(length(all_behav),1);

rng(seed);
[~,fam_rank] = ismember(family,unique(family)); % adjust family numbers to range from 1 to max number of current families
fold_size = floor(nfams/kfolds);
all_ind = [];
for idx=1:kfolds
    all_ind=[all_ind; idx*ones(fold_size,1)];
end
leftover=mod(nfams,kfolds);
fam_indices=[all_ind; randperm(kfolds,leftover)'];
fam_indices=fam_indices(randperm(length(fam_indices)));

for i_fold = 1 : kfolds
    fprintf('%dth fold\n', i_fold);
    if kfolds == 1
        fam_indices(:) = 0;
    end
    fam_testinds = find(fam_indices==i_fold);
    fam_traininds = find(fam_indices~=i_fold);
    test_idx = ismember(fam_rank,fam_testinds);
    train_idx = ismember(fam_rank,fam_traininds);  
    all_test_idx(test_idx) = i_fold;
    
    train_mats = all_edges(:, train_idx);
    test_mats = all_edges(:, test_idx);
    
    
    train_behav = all_behav;
    train_behav(test_idx) = [];
    
    if run_cpm
        % first step univariate edge selection
        [~, edge_p] = corr(train_mats', train_behav);
        edges_1 = find(edge_p<thresh);
        disp(['edge size = ' num2str(size(edges_1))]);
        
        % build model on TRAIN subs
    
        if isnan(lambda)
            [fit_coef, fit_info] = lasso(train_mats(edges_1, :)', train_behav, 'Alpha',v_alpha, 'CV', 10);
            idxLambda1SE = fit_info.Index1SE;
            coef = fit_coef(:,idxLambda1SE);
            coef0 = fit_info.Intercept(idxLambda1SE);
            lambda_total(i_fold) = fit_info.Lambda(idxLambda1SE);
        elseif length(lambda)==1
            [coef, fit_info] = lasso(train_mats(edges_1, :)', train_behav, 'Alpha',v_alpha, 'Lambda', lambda);
            coef0 = fit_info.Intercept;
        elseif length(lambda)>1
            [fit_coef, fit_info] = lasso(train_mats(edges_1, :)', train_behav, 'Alpha',v_alpha, 'CV', 10, 'Lambda', lambda);
            idxLambda1SE = fit_info.Index1SE;
            coef = fit_coef(:,idxLambda1SE);
            coef0 = fit_info.Intercept(idxLambda1SE);
            lambda_total(i_fold) = fit_info.Lambda(idxLambda1SE);
        end
        
        if kfolds~=1
            % run model on TEST sub with the best lambda parameter  
            y(test_idx) = test_mats(edges_1, :)'*coef+coef0; 
        end
        coef_total(edges_1, i_fold) = coef;
        coef0_total(:, i_fold) = coef0;
    else
        y = nan;
        coef_total = nan;
        coef0_total = nan;
    end
    
end

if run_cpm
    % compare predicted and observed behaviors
    [r_pearson, ~] = corr(y, all_behav);
    [r_rank, ~] = corr(y, all_behav, 'type', 'spearman');
    mse = sum((y - all_behav).^2) / num_sub_total;
    q_s = 1 - mse / var(all_behav, 1);
else
    q_s = nan;
    r_pearson = nan;
    r_rank = nan;
end

end

