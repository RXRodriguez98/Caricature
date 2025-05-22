function [sep,within,between] = separation(mats)
% takes sets of connectomes that are matched and performs
% separation analysis
% metric is either 'euclidean' or 'correlation'

% INPUTS %
% mats -- nodeXnodeXsubjXn_measures matrix
% metric -- 'correlation' or 'euclidean'

%%
n_subj = size(mats,3);
n_nodes = size(mats,1);
n_edges = n_nodes*(n_nodes-1)/2;
n_measures = size(mats,4);

%% convert connectomes to vectors
vecs = zeros(n_edges,n_measures*n_subj);
for subj = 1:n_subj
    for meas = 1:n_measures
        curr_mat = mats(:,:,subj,meas);
        curr_vec = curr_mat(find(triu(ones(n_nodes),1)));
        vecs(:,(subj-1)*n_measures+meas) = curr_vec;
    end
end

%% create subjectwise similarity matrix
% also find within and between-subject similarity
corr_mat = corr(vecs);
within_mask = zeros(size(corr_mat));
for block = 1:n_subj
    within_mask((n_measures*(block-1)+1):(n_measures*block),...
        (n_measures*(block-1)+1):(n_measures*block)) = 1;
end
within_mask(find(tril(ones(size(corr_mat))))) = nan;
within = corr_mat(within_mask==1);
between = corr_mat(within_mask==0);

%% iterate through scans and add to count
count = 0;
for subj = 1:n_subj
    subj_row_inds = (subj-1)*n_measures + (1:n_measures);
    for row = subj_row_inds
        full_row = corr_mat(row,:);
        within_subj_inds = subj_row_inds;
        within_subj_inds(subj_row_inds==row) = [];
        between_subj_inds = setdiff(1:(n_measures*n_subj),subj_row_inds);
        min_within = min(full_row(within_subj_inds));
        max_between = max(full_row(between_subj_inds));
        if min_within > max_between
            count = count + 1;
        end
    end
end

%% create stat
sep = 100*count/(n_measures*n_subj);

end