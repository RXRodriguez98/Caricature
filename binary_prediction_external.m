function results = ...
    binary_prediction_external(timeseries_train,timeseries_test,...
    eigenvectors,phenotype_train,phenotype_test,pthresh,seeds)
%% builds predictive models for standard and caricatured data vs binary phenotypes

% INPUTS %
% Name: timeseries_train, Data Type: matrix, Size: n_frames1 X n_nodes X
% n_subjects1
% Description: This variable is a 3D matrix of timeseries data for a
% particular task across all subjects, used to train the model.
%
% Name: timeseries_test, Data Type: matrix, Size: n_frames2 X n_nodes X
% n_subjects2
%
% Description: This variable is a 3D matrix of timeseries data used to test
% the model.
%
% Name: eigenvectors, Data Type: matrix, Size: n_nodes X n_nodes
% Description: This variable is a matrix of eigenvectors.
%
% Name: phenotype_train, Data Type: vector, Size: n_subjects1 X 1
% Description: This variable is a vector representing the binary phenotype
% (as 0 or 1) across training subjects.
%
% Name: phenotype_test, Data Type: vector, Size: n_subjects2 X 1
% Description: This variable is a vector representing the binary phenotype
% (as 0 or 1) across test subects.
%
% Name: pthresh, Data Type: int
% Description: This variable is the p-value threshold to be used for
% feature selection.
%
% Name: seeds, Data Type: vector, Size: n_iter X 1
% Description: This variable is a vector of randomization seeds to be used
% when initializing a model.

% OUTPUTS %
% Name: results, Data Type: struct
% Description: This variable stores the predictions made by both the
% standard and caricatured models as well as the subjects used for
% prediction in each iteration.

%% Set eigenvectors to project away from
eigs2remove = 1:5;

%% Get required numbers
n_nodes = size(timeseries_train,2);

%% Check for subjects missing the phenotype
all_subj_train = 1:length(phenotype_train);
subj_no_pheno_train = find(isnan(phenotype_train)); % subjects with no phenotype
phenotype_train(subj_no_pheno_train) = [];

all_subj_test = 1:length(phenotype_test);
subj_no_pheno_test = find(isnan(phenotype_test)); 
phenotype_test(subj_no_pheno_test) = [];
pred_subj_test = all_subj_test;
pred_subj_test(subj_no_pheno_test) = [];

%% Process timeseries
timeseries_train = zscore(timeseries_train,0,1);
timeseries_train = permute(timeseries_train,[3 1 2]);
timeseries_train(subj_no_pheno_train,:,:) = [];

timeseries_test = zscore(timeseries_test,0,1);
timeseries_test = permute(timeseries_test,[3 1 2]);
timeseries_test(subj_no_pheno_test,:,:) = [];

%% Create standard connectomes
mats_train = zeros(n_nodes*(n_nodes-1)/2,length(phenotype_train));
for subj = 1:length(phenotype_train)
    % fisher transform
    curr_mat = atanh(corr(squeeze(timeseries_train(subj,:,:))));
    % vectorize
    mats_train(:,subj) = curr_mat(find(triu(ones(n_nodes),1)));
end

mats_test = zeros(n_nodes*(n_nodes-1)/2,length(phenotype_test));
for subj = 1:length(phenotype_test)
    % fisher transform
    curr_mat = atanh(corr(squeeze(timeseries_test(subj,:,:))));
    % vectorize
    mats_test(:,subj) = curr_mat(find(triu(ones(n_nodes),1)));
end

%% Create caricatured connectomes
mats_car_train_tmp = caricature(permute(timeseries_train,[2 3 1]),eigenvectors,eigs2remove,1);
mats_car_train = zeros(size(mats_train));
for subj = 1:length(phenotype_train)
    curr_mat = mats_car_train_tmp(:,:,subj);
    mats_car_train(:,subj) = curr_mat(find(triu(ones(n_nodes),1)));
end

mats_car_test_tmp = caricature(permute(timeseries_test,[2 3 1]),eigenvectors,eigs2remove,1);
mats_car_test = zeros(size(mats_test));
for subj = 1:length(phenotype_test)
    curr_mat = mats_car_test_tmp(:,:,subj);
    mats_car_test(:,subj) = curr_mat(find(triu(ones(n_nodes),1)));
end

%% Initialize results
results = struct;
results.standard.pred = nan(length(all_subj_test),length(seeds));
results.caricatured.pred = results.standard.pred;
results.predicted_subjects = pred_subj_test;

%% Iterate through seeds and perform CPM
for seed_ind = 1:length(seeds)
    seed = seeds(seed_ind);
    rng(seed);
    train_inds = sort(randsample(1:length(phenotype_train),floor(0.9*length(phenotype_train))));
    % standard train
    [~,~,~,curr_standard_mdl,curr_standard_feat] = ...
        cpm_classifier_family_cv(mats_train(:,train_inds),...
        phenotype_train(train_inds),1,pthresh,'svm',0,seed,...
        1:length(train_inds),1);
    % caricatured train
    [~,~,~,curr_car_mdl,curr_car_feat] = ...
        cpm_classifier_family_cv(mats_car_train(:,train_inds),...
        phenotype_train(train_inds),1,pthresh,'svm',0,seed,...
        1:length(train_inds),1);

    % standard test
    curr_pred = nan(length(all_subj_test),1);
    curr_pred(pred_subj_test) = ...
        predict(curr_standard_mdl{1},mats_test(curr_standard_feat{1},:)');
    % caricatured test
    curr_pred_car = nan(length(all_subj_test),1);
    curr_pred_car(pred_subj_test) = ...
        predict(curr_car_mdl{1},mats_car_test(curr_car_feat{1},:)');
    
    results.standard.pred(:,seed_ind) = curr_pred;
    results.caricatured.pred(:,seed_ind) = curr_pred_car;
end

end