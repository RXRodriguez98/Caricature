function [ypred,mean_mdl_sz,all_test_idx,all_mdl,all_feat_loc] = ...
    cpm_classifier_family_cv(X, y, k, pthresh,...
    learner, optimize_hyperparams, seed, family, run_cpm)
% inputs:
% X (connectome) and y (phenotype) data
%
% seed: random seed
%
% per_feat: top % of features to use (e.g. 0.1 for 10%)
%
% learner: svm or logistic supported. typically better
% performance is seen with svm for connectomes
%
% balance: balance the data such that numbers are equal between classes
%
% NOTE: this is for binary classification only

% outputs:
% ypred: prediction values
% acc: accuracy of the prediction

% classes, often 0/1
y_class1 = min(y); y_class2 = max(y);

% Set seed and initialize predicted y variable
rng(seed)
nsub = size(X,2);
nfams = length(unique(family));
[~,fam_rank] = ismember(family,unique(family)); % adjust family numbers to range from 1 to max number of current families
fold_size = floor(nfams/k);
all_ind = [];
for idx=1:k
    all_ind=[all_ind; idx*ones(fold_size,1)];
end
leftover=mod(nfams,k);
fam_indices=[all_ind; randperm(k,leftover)'];
fam_indices=fam_indices(randperm(length(fam_indices)));

ypred = zeros(nsub, 1);
all_feat_loc = cell(1,k);
all_mdl = cell(1,k);
mean_mdl_sz = cell(1,k);
mdl_sz = zeros(1,k);
all_test_idx = nan(nsub,1);
for fold_idx = 1:k
    
    disp(['**** Fold ', num2str(fold_idx), ' ****']) 
    if k == 1
        fam_indices(:) = 0;
    end
    fam_testinds = find(fam_indices==fold_idx);
    fam_traininds = find(fam_indices~=fold_idx);
    testinds = ismember(fam_rank,fam_testinds);
    traininds = ismember(fam_rank,fam_traininds);
    all_test_idx(testinds) = fold_idx;
    if run_cpm
        % get train and test
        Xtrain = X(:,traininds)';
        Xtest = X(:,testinds)';
        ytrain = y(traininds);
        
        % feature selection
        [~, p_all] = ttest2(Xtrain(ytrain==y_class1, :), Xtrain(ytrain==y_class2, :));
        feat_loc = find(p_all<pthresh);
        
        % search for regularization parameters (Note: you can typically remove
        % the optimization and see similar results in much shorter time)
        if optimize_hyperparams
            mdl = fitclinear(Xtrain(:, feat_loc), ytrain, 'Learner', learner, ...
                'Regularization' ,'Ridge', 'OptimizeHyperparameters', {'Lambda'},...
                'HyperparameterOptimizationOptions', struct('Kfold', 5, 'ShowPlots', 0, 'Verbose', 0));
        else
            mdl = fitclinear(Xtrain(:, feat_loc), ytrain, 'Learner', learner, 'Lambda', 0);  
        end
        
        if k ~= 1
            ypred(testinds) = predict(mdl, Xtest(:, feat_loc));
        end
        all_feat_loc{fold_idx} = feat_loc;
        all_mdl{fold_idx} = mdl;
        mdl_sz(fold_idx) = length(feat_loc);
    end
end
if run_cpm
    mean_mdl_sz = mean(mdl_sz);
else
    ypred = nan;
    mean_mdl_sz = nan;
end
end





