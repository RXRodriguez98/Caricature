function results = ...
    continuous_prediction(timeseries,eigenvectors,family,phenotype,pthresh,seeds)
%% builds predictive models for standard and caricatured data vs continuous phenotypes

% INPUTS %
% Name: timeseries, Data Type: matrix, Size: n_frames X n_nodes X
% n_subjects
% Description: This variable is a 3D matrix of timeseries data for a
% particular task across all subjects.
%
% Name: eigenvectors, Data Type: matrix, Size: n_nodes X n_nodes
% Description: This variable is a matrix of eigenvectors.
%
% Name: family, Data Type: vector, Size: n_subjects X 1
% Description: This variable indicates which family each subject belongs
% to.
%
% Name: phenotype, Data Type: vector, Size: n_subjects X 1
% Description: This variable is a vector representing the phenotype across 
% subjects.
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
n_frames = size(timeseries,1);
n_nodes = size(timeseries,2);

%% Check for subjects missing the phenotype
all_subj = 1:length(phenotype);
subj_no_pheno = find(isnan(phenotype)); % subjects with no phenotype
phenotype(subj_no_pheno) = [];
family(subj_no_pheno) = [];
pred_subj = all_subj;
pred_subj(subj_no_pheno) = [];

%% Process timeseries
timeseries = zscore(timeseries,0,1);
timeseries = permute(timeseries,[3 1 2]);

%% Create standard connectomes
mats = zeros(n_nodes*(n_nodes-1)/2,length(all_subj));
for subj = 1:length(all_subj)
    % fisher transform
    curr_mat = atanh(corr(squeeze(timeseries(subj,:,:))));
    % vectorize
    mats(:,subj) = curr_mat(find(triu(ones(n_nodes),1)));
end

%% Create caricatured connectomes
mats_car_tmp = caricature(permute(timeseries,[2 3 1]),eigenvectors,eigs2remove,1);
mats_car = zeros(size(mats));
for subj = 1:length(all_subj)
    curr_mat = mats_car_tmp(:,:,subj);
    mats_car(:,subj) = curr_mat(find(triu(ones(n_nodes),1)));
end

%% Initialize results
results = struct;
results.standard.pred = nan(length(all_subj),length(seeds));
results.caricatured.pred = results.standard.pred;
results.predicted_subjects = pred_subj;

%% Iterate through seeds and perform CPM
for seed_ind = 1:length(seeds)
    seed = seeds(seed_ind);
    % standard
    curr_pred = nan(length(all_subj),1);
    [~,~,~,curr_pred(pred_subj)] = ridgeCPM(mats(:,pred_subj),phenotype,...
        'thresh',pthresh,'k',10,'family',family,'seed',seed);
    % caricatured
    curr_pred_car = nan(length(all_subj),1);
    [~,~,~,curr_pred_car(pred_subj)] = ridgeCPM(mats_car(:,pred_subj),phenotype,...
        'thresh',pthresh,'k',10,'family',family,'seed',seed);
    
    results.standard.pred(:,seed_ind) = curr_pred;
    results.caricatured.pred(:,seed_ind) = curr_pred_car;
end

end