function discr = discriminability(mats,metric)
%% Computes the discriminability of a set of matrices with multiple measures
%% per subject

%% NOTE: Currently only supports the case in which all subjects have the same
%% number of measures

% INPUTS %
% Name: mats, Data Type: 4D matrix, Size: n_nodes X n_nodes X n_subjects X n_measures
%   OR        Data Type: 3D matrix, Size: n_edges X n_subjects X n_measures
% Description: This variable is a set of connectomes across subjects with
% multiple measures. When 4-dimensional, the first two dimensions represent
% a typical connectome which can be indexed by the third and fourth
% dimensions (subjects and measures, respectively). When 3-dimensional, the
% first dimension represents a vectorized connectomes which can be indexed
% by the second and third dimensions (subjects and measures, respectively).
%
% Name: metric, Data Type: string
% Description: This variable is a string indicating the distance metric
% desired for the descriminability calculation. For example, 'correlation'
% or 'euclidean' could be used here.

% OUTPUTS %
% Name: discr, Data Type: double
% Description: This variable is the discriminability of the given dataset.

%% Convert connectomes to 2D
if length(size(mats)) == 4 % check if mats is 4D
    n_subjects = size(mats,3);
    n_nodes = size(mats,1);
    n_edges = n_nodes*(n_nodes-1)/2;
    n_measures = size(mats,4);
    
    % convert connectomes to vectors
    vectors = zeros(n_edges,n_measures*n_subjects);
    for subject_index = 1:n_subjects
        for measure_index = 1:n_measures
            curr_mat = mats(:,:,subject_index,measure_index);
            curr_vec = curr_mat(find(triu(ones(n_nodes),1)));
            vectors(:,(subject_index-1)*n_measures+measure_index) = curr_vec;
        end
    end
else % if mats is 3D
    n_edges = size(mats,1);
    n_subjects = size(mats,2);
    n_measures = size(mats,3);

    % convert connectomes to vectors
    vectors = zeros(n_edges,n_measures*n_subjects);
    for subject_index = 1:n_subjects
        curr_vec = squeeze(mats(:,subject_index,:));
        vectors(:,((subject_index-1)*n_measures+1):(subject_index*n_measures)) = curr_vec;
    end
end

%% Create pairwise distance matrix
dist_mat = pdist2(vectors',vectors',metric);

%% Calculate discriminability
%% iterate through rows and add to count
count = 0; % number of times a within-subject comparison is closer than a between-subject comparison for the same scan
for subject_index = 1:n_subjects
    dist_inds = (subject_index-1)*n_measures + (1:n_measures);
    for row = dist_inds
        full_row = dist_mat(row,:);
        cross_scan_inds = dist_inds;
        cross_scan_inds(dist_inds==row) = [];
        for cross_scan_ind = cross_scan_inds
            curr_comparison = full_row(cross_scan_ind);
            other_comparisons = full_row;
            other_comparisons(dist_inds) = [];
            count = count + length(find(other_comparisons>=curr_comparison));
        end
    end
end
discr = count/(n_subjects*n_measures*(n_subjects*n_measures-n_measures)*(n_measures-1));


end