function [eigenvectors,variance_explained] = create_task_manifold(timeseries)
%% concatenates time-series for all subjects and tasks and performs
%% eigendecomposition on the corresponding covariance matrix

% INPUTS %
% Name: timeseries, Data Type: cell, Size: n_tasks X 1
% Description: This variable is a cell of length equal to the number of
% tasks (n_tasks). In each cell entry i is a 3D matrix for which the 
% columns are the nodes (n_nodes), the rows are the frames
% (length of task i), and the third dimension represents the subjects 
% (n_subjects). 

% OUTPUTS %
% Name: eigenvectors, Data Type: matrix, Size: n_nodes X n_nodes
% Description: This variable is a square matrix of size equal to the
% number of nodes in the parcellation atlas. Each column represents the
% eigenvectors in order of variance explained, and each row represents how
% that node contributes to the eigenvector.
%
% Name: variance_explained, Data Type: vector, Size: n_nodes X 1
% Description: This variable is a vector of length equal to the number of
% nodes in the parcellation atlas. Each entry represents the amount of
% variance explained in timeseries by the corresponding eigenvector.

%% Set basic values
n_tasks = length(timeseries);
n_subjects = size(timeseries{1},3);
n_nodes = size(timeseries{1},2);

%% Z-score data and concatenate data in time
full_timeseries = [];
for task_index = 1:n_tasks
    current_task = timeseries{task_index};
    current_n_frames = size(current_task,1);
    % z-score
    current_task_z = zscore(current_task,0,1);
    % create full timeseries for the current task
    full_timeseries_current_task = zeros(n_subjects*current_n_frames,n_nodes);
    for subject_index = 1:n_subjects
        full_timeseries_current_task(((subject_index-1)*current_n_frames+1):...
            (subject_index*current_n_frames),:) = current_task_z(:,:,subject_index);
    end
    % concatenate
    full_timeseries = cat(1,full_timeseries,full_timeseries_current_task);
end

%% Perform Eigendecomposition
% calculate covariance matrix
covariance_mat = cov(full_timeseries,1);
% eigendecomposition
[eigenvectors,variance_explained] = eig(covariance_mat);
% sort eigenvectors in descending order of explained variance
[variance_explained,sort_inds] = sort(variance_explained(find(eye(n_nodes))),'descend');
eigenvectors = eigenvectors(:,sort_inds);

end