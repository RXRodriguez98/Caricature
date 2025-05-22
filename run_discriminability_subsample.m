function results = ...
    run_discriminability_subsample(timeseries,eigenvectors,subsample_inds)
% INPUTS %
% Name: timeseries, Data Type: matrix, Size: n_frames X n_nodes X
% n_subjects X n_runs
% Description: This variable is a 4D matrix of timeseries data for a
% particular task across all subjects.
%
% Name: eigenvectors, Data Type: matrix, Size: n_nodes X n_nodes X n_iter
% Description: This variable is a matrix of eigenvectors.
%
% Name: subsample_inds, Data Type: cell, Size: n_iter X 1
% Description: This variable indicates which subjects were used to create
% the manifold for each iteration.

%% Set eigenvectors to project away from
eigs2remove = 1:5;

%% Process timeseries
timeseries = zscore(timeseries,0,1);
timeseries = permute(timeseries,[3 1 2 4]);
timeseries = timeseries(:,1:176,:,:);

%% Get required numbers
n_subjects = size(timeseries,1);
n_nodes = size(timeseries,3);
n_runs = size(timeseries,4);
n_iter = length(subsample_inds);

%% Create standard connectomes
conn = zeros(n_nodes,n_nodes,n_subjects,n_runs);

for subj = 1:n_subjects
    for run = 1:n_runs
        conn(:,:,subj,run) = atanh(corr(squeeze(timeseries(subj,:,:,run))));
    end
end

%% Initialize results
results = struct;
results.standard = nan(n_iter,1);
results.caricatured = nan(n_iter,1);

%% Run discriminability
for seed_ind = 1:n_iter
    curr_conn = conn;
    curr_conn(:,:,subsample_inds{seed_ind},:) = [];
    % create caricatured data
    curr_ts = timeseries;
    curr_ts(subsample_inds{seed_ind},:,:,:) = [];
    conn1_car = caricature(permute(curr_ts(:,:,:,1),[2 3 1]),eigenvectors(:,:,seed_ind),eigs2remove,1);
    conn2_car = caricature(permute(curr_ts(:,:,:,2),[2 3 1]),eigenvectors(:,:,seed_ind),eigs2remove,1);
    curr_conn_car = cat(4,conn1_car,conn2_car);
    
    % discriminability
    % standard
    results.standard(seed_ind) = discriminability(curr_conn,'correlation');
    % caricatured
    results.caricatured(seed_ind) = discriminability(curr_conn_car,'correlation');
end
end