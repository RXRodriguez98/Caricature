function [acc,within,between,mats1_id_count,mats2_id_count] = ...
    fingerprinting(mats1,mats2)
% takes two sets of connectomes that are matched and performs
% fingerprinting
% mats could be vectors then no conversion needed (vec is edges X subj)
% acc is vector with [database is mats1, database is mats2]
% mats#_id_count is vector showing how many times each subject in mats# was
%   chosen as the ID when mats# was the database

if length(size(mats1)) == 3
    n_subj = size(mats1,3);
    n_nodes = size(mats1,1);
    n_edges = n_nodes*(n_nodes-1)/2;

    % convert connectomes to vectors
    mats1_vecs = zeros(n_edges,n_subj);
    mats2_vecs = mats1_vecs;
    for subj = 1:n_subj
        curr1 = mats1(:,:,subj);
        mats1_vecs(:,subj) = curr1(find(triu(ones(n_nodes),1)));
        curr2 = mats2(:,:,subj);
        mats2_vecs(:,subj) = curr2(find(triu(ones(n_nodes),1)));
    end
    % create subjectwise correlation matrix
    similarity = corr(mats1_vecs,mats2_vecs);
else
    n_subj = size(mats1,2);
    similarity = corr(mats1,mats2);
end

% within vs between subject similarity
within = diag(similarity);
diag_inds = eye(n_subj);
between = similarity(~diag_inds);

% fingerprinting
prediction_of_group2 = zeros(n_subj,1);
prediction_of_group1 = zeros(n_subj,1);
predicted_ids = cell(n_subj,2);
for subj = 1:n_subj
    % group1 subj is target; group2 is database
    predicted_ids{subj,1} = find(similarity(subj,:)==max(similarity(subj,:)));
    group1_pred = predicted_ids{subj,1} == subj;
    % group2 subj is target; group1 is database
    predicted_ids{subj,2} = find(similarity(:,subj)==max(similarity(:,subj)));
    group2_pred = predicted_ids{subj,2} == subj;
    % ensures if multiple ids are selected that it is counted as incorrect
    prediction_of_group1(subj) = min(group1_pred);
    prediction_of_group2(subj) = min(group2_pred);
end
acc = 100*[mean(prediction_of_group2) mean(prediction_of_group1)];

% gets counts of how many times each subject was selected as the id
mats1_id_count = zeros(n_subj,1);
mats2_id_count = zeros(n_subj,1);
for subj = 1:n_subj
    mats1_id_count(predicted_ids{subj,2}) = mats1_id_count(predicted_ids{subj,2}) + 1;
    mats2_id_count(predicted_ids{subj,1}) = mats2_id_count(predicted_ids{subj,1}) + 1;
end
end