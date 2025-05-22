function results = ...
    separability(timeseries,eigenvectors)
% INPUTS %
% Name: timeseries, Data Type: matrix, Size: n_frames X n_nodes X
% n_subjects X n_runs
% Description: This variable is a 4D matrix of timeseries data for a
% particular task across all subjects.
%
% Name: eigenvectors, Data Type: matrix, Size: n_nodes X n_nodes
% Description: This variable is a matrix of eigenvectors.

%% Set eigenvectors to project away from
eigs2remove = 1:5;

%% Process timeseries
timeseries = zscore(timeseries,0,1);
timeseries = permute(timeseries,[3 1 2 4]);

%% Get required numbers
n_subjects = size(timeseries,1);
n_nodes = size(timeseries,3);
n_runs = size(timeseries,4);

%% Create standard connectomes
conn = zeros(n_nodes,n_nodes,n_subjects,n_runs);

for subj = 1:n_subjects
    for run = 1:n_runs
        conn(:,:,subj,run) = atanh(corr(squeeze(timeseries(subj,:,:,run))));
    end
end

%% Create caricatured connectomes
conn_car = zeros(size(conn));

for subj = 1:n_subjects
    for run = 1:n_runs
        curr_ts = timeseries(subj,:,:,run);
        conn_car(:,:,subj,run) = caricature(permute(curr_ts,[2 3 1]),...
            eigenvectors,eigs2remove,1);
    end
end

%% Run separability
results = struct;

% standard
[results.standard.sep,results.standard.within,...
    results.standard.between] = separation(conn);

% caricatured
[results.caricatured.sep,results.caricatured.within,...
    results.caricatured.between] = separation(conn_car);

%% run stats
[p,s1_null,s2_null] = ...
    separation_stat(conn_car,conn,1000);

results.p = 4*p; % correct for multiple comps (two tails and two tests)
results.s1_null = s1_null;
results.s2_null = s2_null;

end