function results = ...
    icc(timeseries1,timeseries2,eigenvectors)
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

%% Create data table for ICC
ftbl = [sort([1:n_subjects 1:n_subjects])' ...
    repmat([1 2]',n_subjects,1)];
all_conn = zeros(n_nodes,n_nodes,2*n_subjects);
all_conn(:,:,1:2:end) = conn1;
all_conn(:,:,2:2:end) = conn2;

%% Create caricatured connectomes
conn1_car = caricature(permute(timeseries1,[2 3 1]),eigenvectors,eigs2remove,1);
conn2_car = caricature(permute(timeseries2,[2 3 1]),eigenvectors,eigs2remove,1);

all_conn_car = zeros(size(all_conn));
all_conn_car(:,:,1:2:end) = conn1_car;
all_conn_car(:,:,2:2:end) = conn2_car;

%% Run ICC
results = struct;

% standard
[v1,v2,ve] = calculate_variance_comps(all_conn,ftbl);
icc = calculate_icc_from_variance(v1,v2,ve);
results.standard = icc;

% caricatured
[v1_car,v2_car,ve_car] = calculate_variance_comps(all_conn_car,ftbl);
icc_car = calculate_icc_from_variance(v1_car,v2_car,ve_car);
results.caricatured = icc_car;

end