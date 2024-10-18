function caricatures = ...
    caricature(timeseries,eigenvectors,eigenvectors2remove,fisher)
%% creates caricatured connectomes

% INPUTS %
% Name: timeseries, Data Type: matrix, Size: n_frames X n_nodes X n_subjects
% Description: This variable is a 3D matrix of time-series where the rows
% are frames (n_frames), the columns are nodes (n_nodes), and the third
% dimension represents the subjects (n_subjects). Note that the matrix 
% should already be z-scored along the time dimension. If you are inputting
% a truncated time-series, then the z-scoring should be done prior to the 
% truncation.
%
% Name: eigenvectors, Data Type: matrix, Size: n_nodes X n_nodes
% Description: This variable is a square matrix of size equal to the
% number of nodes in the parcellation atlas. Each column represents the
% eigenvectorss in order, and each row represents how that node contributes
% to the eigenvector.
%
% Name: eigenvectors2remove, Data Type: matrix, Size: 1 X n_eigenvectors2remove
% Description: This variable is a vector indicating which eigenvectors
% should be projected away in the caricaturing.
%
% Name: fisher, Data Type: binary
% Description: This variable indicates whether to Fisher transform the 
% caricatured connectomes.

% OUTPUTS %
% Name: caricatures, Data Type: matrix, Size: n_nodes X n_nodes X n_subjects
% Description: This variable is a 3D matrix of caricatured connectomes.

%% Set basic values
n_subjects = size(timeseries,3);
n_frames = size(timeseries,1);
n_nodes = size(timeseries,2);

%% Caricature
eigenvectors(:,eigenvectors2remove) = [];
projection_matrix = eigenvectors * eigenvectors';
% permute timeseries to by n_frames X n_subjects X n_nodes
timeseries = permute(timeseries,[1 3 2]);
projected_ts = reshape(timeseries,[],size(timeseries,3)) * projection_matrix;
caricatures = zeros(n_nodes,n_nodes,n_subjects);
for subject_index = 1:n_subjects
    caricatures(:,:,subject_index) = corr(projected_ts(((subject_index-1)*n_frames+1):(subject_index*n_frames),:));
end

%% Fisher transform
if fisher
    caricatures = atanh(caricatures);
end

end