function [p,null_discs1,null_discs2] = discriminability_stat(mats1,mats2,metric,nperm)
% takes two datasets of repeated measures connectomes and does a
% permutation test to determine if one is more discriminable

% INPUTS %
% mats1 -- nodeXnodeXsubjXn_measures matrix
% mats2 -- nodeXnodeXsubjXn_measures matrix
% metric -- 'correlation' or 'euclidean'
% nperm -- number of permutations to use

% mats1 and mats2 need subjects in same order

% OUTPUTS %
% p -- pvalue associated with the one-tailed test that mats1 is more
%      discriminable than mats2
% null_discs1 -- distribution of discriminabilities for null_mats1
% null_discs2 -- distribution of discriminabilities for null_mats2

%% numbers
n_subj = size(mats1,3);

%% convert to vectors
vecs1 = zeros(35778,n_subj,size(mats1,4));
vecs2 = vecs1;
for subj = 1:n_subj
    for meas = 1:size(mats1,4)
        curr1 = mats1(:,:,subj,meas);
        vecs1(:,subj,meas) = curr1(find(triu(ones(268),1)));
        curr2 = mats2(:,:,subj,meas);
        vecs2(:,subj,meas) = curr2(find(triu(ones(268),1)));
    end
end

%% get true discriminability
d1 = discriminability(vecs1,metric);
d2 = discriminability(vecs2,metric);

%% run permutations and create distribution of discriminabilities and differences
null_discs1 = zeros(nperm,1);
null_discs2 = null_discs1;
for perm = 1:nperm
    perm
    tic
    % determine which scans to switch
%     idx = randsample([0 1],n_subj,true);
    null_mats1 = vecs1;
    null_mats2 = vecs2;
    for scan = 1:size(mats1,4)
        idx = randsample([0 1],n_subj,true);
        null_mats1(:,idx==1,scan) = vecs2(:,idx==1,scan);
        null_mats2(:,idx==1,scan) = vecs1(:,idx==1,scan);
    end
    % create null mats
%     null_mats1 = cat(3,mats1(:,:,idx==0,:),mats2(:,:,idx==1,:));
%     null_mats2 = cat(3,mats1(:,:,idx==1,:),mats2(:,:,idx==0,:));
    % compute null discs
    null_discs1(perm) = discriminability(null_mats1,metric);
    null_discs2(perm) = discriminability(null_mats2,metric);
    toc
end
% compute differences
null_dist = null_discs1-null_discs2;
figure; histogram(null_dist);

%% calculate p_val
true_diff = d1-d2;
p = 1-length(find(true_diff>null_dist))/length(null_dist);