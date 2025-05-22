function [all_eigenvectors,all_variance_explained,all_subsample_inds] = ...
    create_task_manifold_subsample(timeseries,family_ID,n_iter)
% INPUTS %
% Name: timeseries, Data Type: cell, Size: n_tasks X 1
% Description: This variable is a cell of length equal to the number of
% tasks (n_tasks). In each cell entry i is a 3D matrix for which the 
% columns are the nodes (n_nodes), the rows are the frames
% (length of task i), and the third dimension represents the subjects 
% (n_subjects). 
%
% Name: family_ID, Data Type: vector, Size: n_subjects X 1
% Description: This variable is a vector of length equal to the number of
% subjects (n_subjects). In each entry is a number indicating which family 
% each subject belongs to.
%
% Name: n_iter, Data Type: int
% Description: Number of subsamples to perform when creating task
% manifolds.

% OUTPUTS %
% Name: all_eigenvectors, Data Type: matrix, Size: n_nodes X n_nodes X
% n_iter
% Description: This variable is a square matrix of size equal to the
% number of nodes in the parcellation atlas. Each column represents the
% eigenvectors in order of variance explained, and each row represents how
% that node contributes to the eigenvector.
%
% Name: all_variance_explained, Data Type: matrix, Size: n_nodes X n_iter
% Description: This variable is a matrix of dimensions equal to the number of
% nodes in the parcellation atlas and the nubmer of iterations. Each entry 
% represents the amount of variance explained in timeseries by the 
% corresponding eigenvector.
%
% Name: subsample_inds, Data Type: cell, Size: n_iter X 1
% Description: This variable stores the indices used for each subsample to
% calculate the manifold.

%% Set basic values
n_tasks = length(timeseries);
n_nodes = size(timeseries{1},2);

%% Initialize outputs
all_eigenvectors = zeros(n_nodes,n_nodes,n_iter);
all_variance_explained = zeros(n_nodes,n_iter);
all_subsample_inds = cell(n_iter,1);

%% Iterate through subsamples
family_size = 10;
for seed = 1:n_iter
    rng(seed);
    families = randsample(unique(family_ID),family_size);
    subsample_inds = find(ismember(family_ID,families));
    full_timeseries = [];
    for task_index = 1:n_tasks
        current_task = timeseries{task_index};
        current_task = current_task(:,:,subsample_inds);
        current_n_frames = size(current_task,1);
        current_n_subjects = size(current_task,3);
        % z-score
        current_task_z = zscore(current_task,0,1);
        % create full timeseries for the current task
        full_timeseries_current_task = zeros(current_n_subjects*current_n_frames,n_nodes);
        for subject_index = 1:current_n_subjects
            full_timeseries_current_task(((subject_index-1)*current_n_frames+1):...
                (subject_index*current_n_frames),:) = current_task_z(:,:,subject_index);
        end
        % concatenate
        full_timeseries = cat(1,full_timeseries,full_timeseries_current_task);
    end

    % perform eigendecompostion
    % calculate covariance matrix
    covariance_mat = cov(full_timeseries,1);
    % eigendecomposition
    [eigenvectors,variance_explained] = eig(covariance_mat);
    % sort eigenvectors in descending order of explained variance
    [variance_explained,sort_inds] = sort(variance_explained(find(eye(n_nodes))),'descend');
    eigenvectors = eigenvectors(:,sort_inds);
    all_eigenvectors(:,:,seed) = eigenvectors;
    all_variance_explained(:,seed) = variance_explained;
    all_subsample_inds{seed} = subsample_inds;
end

end