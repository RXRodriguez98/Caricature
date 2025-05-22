function results = ...
    icc_subsample(timeseries1,timeseries2,eigenvectors,subsample_inds)
% INPUTS %
% Name: timeseries1, Data Type: matrix, Size: n_frames X n_nodes X
% n_subjects
% Description: This variable is a 3D matrix of timeseries data for a
% particular task across all subjects. Scan # 1.
%
% Name: timeseries2, Data Type: matrix, Size: n_frames X n_nodes X
% n_subjects
% Description: This variable is a 3D matrix of timeseries data for a
% particular task across all subjects. Scan # 2.
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
timeseries1 = zscore(timeseries1,0,1);
timeseries1 = permute(timeseries1,[3 1 2]);

timeseries2 = zscore(timeseries2,0,1);
timeseries2 = permute(timeseries2,[3 1 2]);

%% Get required numbers
n_subjects = size(timeseries1,1);
n_nodes = size(timeseries1,3);
n_iter = length(subsample_inds);

%% Create standard connectomes
conn1 = zeros(n_nodes,n_nodes,n_subjects);
conn2 = conn1;

for subj = 1:n_subjects
    conn1(:,:,subj) = atanh(corr(squeeze(timeseries1(subj,:,:))));
    conn2(:,:,subj) = atanh(corr(squeeze(timeseries2(subj,:,:))));
end

%% Initialize results
results = struct;
results.standard = nan(n_nodes,n_nodes,n_iter);
results.caricatured = results.standard;

%% Iterate through n_iter and perform icc
for seed_ind = 1:n_iter
    % set ftbl for this seed
    ftbl = [sort([1:(n_subjects-length(subsample_inds{seed_ind})) 1:(n_subjects-length(subsample_inds{seed_ind}))])' ...
        repmat([1 2]',n_subjects-length(subsample_inds{seed_ind}),1)]; 
    curr_conn1 = conn1;
    curr_conn1(:,:,subsample_inds{seed_ind}) = [];
    curr_conn2 = conn2;
    curr_conn2(:,:,subsample_inds{seed_ind}) = [];
    all_conn = zeros(n_nodes,n_nodes,2*size(curr_conn1,3));
    all_conn(:,:,1:2:end) = curr_conn1;
    all_conn(:,:,2:2:end) = curr_conn2;
    % create caricatured data
    curr_ts1 = timeseries1;
    curr_ts1(subsample_inds{seed_ind},:,:) = [];
    curr_conn1_car = caricature(permute(curr_ts1,[2 3 1]),eigenvectors(:,:,seed_ind),eigs2remove,1);
    curr_ts2 = timeseries2;
    curr_ts2(subsample_inds{seed_ind},:,:) = [];
    curr_conn2_car = caricature(permute(curr_ts2,[2 3 1]),eigenvectors(:,:,seed_ind),eigs2remove,1);
    all_conn_car = zeros(n_nodes,n_nodes,2*size(curr_conn1_car,3));
    all_conn_car(:,:,1:2:end) = curr_conn1_car;
    all_conn_car(:,:,2:2:end) = curr_conn2_car;
    
    % run icc
    [v1,v2,ve] = calculate_variance_comps(all_conn,ftbl);
    icc = calculate_icc_from_variance(v1,v2,ve);
    results.standard(:,:,seed_ind) = icc;
    [v1_car,v2_car,ve_car] = calculate_variance_comps(all_conn_car,ftbl);
    icc_car = calculate_icc_from_variance(v1_car,v2_car,ve_car);
    results.caricatured(:,:,seed_ind) = icc_car;

end

end