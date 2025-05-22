function results = ...
    icc_3_factors(timeseries,eigenvectors)
% INPUTS %
% Name: timeseries, Data Type: matrix, Size: n_frames X n_nodes X
% n_subjects X n_sessions X n_runs
% Description: This variable is a 5D matrix of timeseries data for a
% particular task across all subjects.
%
% Name: eigenvectors, Data Type: matrix, Size: n_nodes X n_nodes
% Description: This variable is a matrix of eigenvectors.

%% Set eigenvectors to project away from
eigs2remove = 1:5;

%% Process timeseries
timeseries = zscore(timeseries,0,1);
timeseries = permute(timeseries,[3 1 2 4 5]);

%% Get required numbers
n_subjects = size(timeseries,1);
n_nodes = size(timeseries,3);
n_sessions = size(timeseries,4);
n_runs = size(timeseries,5);

%% Create standard connectomes and data table
conn = zeros(n_nodes,n_nodes,n_subjects*n_sessions*n_runs);
ftbl = zeros(n_subj*n_sessions*n_runs,3);

index = 1;
for subj = 1:n_subjects
    for sess = 1:n_sessions
        for run = 1:n_runs
            conn(:,:,index) = atanh(corr(squeeze(timeseries(subj,:,:,sess,run))));
            ftbl(index,:) = [subj sess run];
            index = index + 1;
        end
    end
end

%% Create caricatured connectomes
conn_car = zeros(size(conn));

index = 1;
for subj = 1:n_subjects
    for sess = 1:n_sessions
        for run = 1:n_runs
            curr_ts = timeseries(subj,:,:,sess,run);
            conn_car(:,:,index) = caricature(permute(curr_ts,[2 3 1]),...
                eigenvectors,eigs2remove,1);
            index = index + 1;
        end
    end
end

%% Run ICC
results = struct;

% standard
[v1,v2,v3,v1x2,v1x3,v2x3,ve] = calculate_variance_comps_3_factors(conn,ftbl);
icc = calculate_icc_from_variance_3_factors(v1,v2,v3,v1x2,v1x3,v2x3,ve);
results.standard = icc;

% caricatured
[v1,v2,v3,v1x2,v1x3,v2x3,ve] = calculate_variance_comps_3_factors(conn_car,ftbl);
icc_car = calculate_icc_from_variance_3_factors(v1,v2,v3,v1x2,v1x3,v2x3,ve);
results.caricatured = icc_car;

end