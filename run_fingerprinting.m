function results = ...
    run_fingerprinting(timeseries1,timeseries2,eigenvectors)
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
% Name: eigenvectors, Data Type: matrix, Size: n_nodes X n_nodes
% Description: This variable is a matrix of eigenvectors.

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

%% Create standard connectomes
conn1 = zeros(n_nodes,n_nodes,n_subjects);
conn2 = conn1;

for subj = 1:n_subjects
    conn1(:,:,subj) = atanh(corr(squeeze(timeseries1(subj,:,:))));
    conn2(:,:,subj) = atanh(corr(squeeze(timeseries2(subj,:,:))));
end

%% Create caricatured connectomes
conn1_car = caricature(permute(timeseries1,[2 3 1]),eigenvectors,eigs2remove,1);
conn2_car = caricature(permute(timeseries2,[2 3 1]),eigenvectors,eigs2remove,1);

%% Run Fingerprinting
results = struct;

% standard
[curr_accs,curr_within,curr_between,curr_counts_LR,curr_counts_RL] = ...
    fingerprinting(conn1,conn2);
results.standard.acc = mean(curr_accs);
results.standard.within = curr_within;
results.standard.between = curr_between;
results.standard.LR_id_count = curr_counts_LR;
results.standard.RL_id_count = curr_counts_RL;

% caricatured 
[curr_accs,curr_within,curr_between,curr_counts_LR,curr_counts_RL] = ...
    fingerprinting(conn1_car,conn2_car);
results.caricatured.acc = mean(curr_accs);
results.caricatured.within = curr_within;
results.caricatured.between = curr_between;
results.caricatured.LR_id_count = curr_counts_LR;
results.caricatured.RL_id_count = curr_counts_RL;

%% run stats
[p,car_null_accs,standard_null_accs] = ...
    fingerprinting_stat(cat(4,conn1,conn2),...
    cat(4,conn1_car,conn2_car),1000);
results.p_rest = 4*p; % correct for multiple comps (one-tail and two tests)
results.standard_null_accs = standard_null_accs;
results.caricatured_null_accs = car_null_accs;

end