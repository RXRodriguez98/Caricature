function results = ...
    binary_prediction_subsample_nulls(timeseries,eigenvectors,subsample_inds,...
    family,phenotype,pthresh)
%% builds null predictive models for standard and caricatured data vs binary phenotypes
%% for use when timeseries and eigenvectors come from the same subjects
%% calls binary_prediction_nulls.m

% INPUTS %
% Name: timeseries, Data Type: matrix, Size: n_frames X n_nodes X
% n_subjects
% Description: This variable is a 3D matrix of timeseries data for a
% particular task across all subjects.
%
% Name: eigenvectors, Data Type: matrix, Size: n_nodes X n_nodes X n_iter
% Description: This variable is a matrix of eigenvectors.
%
% Name: subsample_inds, Data Type: cell, Size: n_iter X 1
% Description: This variable indicates which subjects were used to create
% the manifold for each iteration.
%
% Name: family, Data Type: vector, Size: n_subjects X 1
% Description: This variable indicates which family each subject belongs
% to.
%
% Name: phenotype, Data Type: vector, Size: n_subjects X 1
% Description: This variable is a vector representing the binary phenotype
% (as 0 or 1) across subjects.
%
% Name: pthresh, Data Type: int
% Description: This variable is the p-value threshold to be used for
% feature selection.

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
n_iter = size(eigenvectors,3);

% Check for subjects missing the phenotype
all_subj = 1:length(phenotype);
subj_no_pheno = find(isnan(phenotype));
pred_subj = all_subj;
pred_subj(subj_no_pheno) = [];

%% Process timeseries
timeseries = zscore(timeseries,0,1);
timeseries = permute(timeseries,[3 1 2]);

%% Initialize results
results = struct;
results.standard.pred = nan(length(all_subj),n_iter);
results.standard.null_phenotype = nan(length(all_subj),n_iter);
results.caricatured.pred = results.standard.pred;
results.caricatured.null_phenotype = results.standard.null_phenotype;
results.predicted_subjects = cell(n_iter,1);

%% Iterate through seeds and call binary_prediction
for seed_ind = 1:n_iter
    curr_subj = setdiff(pred_subj,subsample_inds{seed_ind});
    curr_timeseries = timeseries(curr_subj,:,:);
    curr_timeseries = permute(curr_timeseries,[2 3 1]);
    curr_eigenvectors = eigenvectors(:,:,seed_ind);
    curr_family = family(curr_subj);
    curr_phenotype = phenotype(curr_subj);
    curr_results = binary_prediction_nulls(curr_timeseries,curr_eigenvectors,...
        curr_family,curr_phenotype,pthresh,seed_ind);

    results.standard.pred(curr_subj,seed_ind) = curr_results.standard.pred;
    results.standard.null_phenotype(curr_subj,seed_ind) = ...
        curr_results.standard.null_phenotype;
    results.caricatured.pred(curr_subj,seed_ind) = curr_results.caricatured.pred;
    results.caricatured.null_phenotype(curr_subj,seed_ind) = ...
        curr_results.caricatured.null_phenotype;
    results.predicted_subjects{seed_ind} = curr_subj;
end

end