function [p,null_accs1,null_accs2] = fingerprinting_stat(mats1,mats2,nperm)
% takes two datasets of 2-measure connectomes and does a permutation 
% test to determine if one has higher fingerprinting accuracy

% INPUTS %
% mats1 -- nodeXnodeXsubjX2 matrix
% mats2 -- nodeXnodeXsubjX2 matrix
% nperm -- number of permutations to use

% mats1 and mats2 need subjects in same order

% OUTPUTS %
% p -- pvalue associated with the one-tailed test that mats1 has a higher
%      fingerprinting accuracy than mats2
% null_accs1 -- distribution of accuracies for null_mats1
% null_accs2 -- distribution of accuracies for null_mats2

%% numbers
n_subj = size(mats1,3);

%% convert matrices to vectors
vecs1 = zeros(35778,n_subj,2);
vecs2 = zeros(35778,n_subj,2);
for subj = 1:n_subj
    curr_1_1 = mats1(:,:,subj,1);
    curr_1_2 = mats1(:,:,subj,2);
    vecs1(:,subj,1) = curr_1_1(find(triu(ones(268),1)));
    vecs1(:,subj,2) = curr_1_2(find(triu(ones(268),1)));

    curr_2_1 = mats2(:,:,subj,1);
    curr_2_2 = mats2(:,:,subj,2);
    vecs2(:,subj,1) = curr_2_1(find(triu(ones(268),1)));
    vecs2(:,subj,2) = curr_2_2(find(triu(ones(268),1)));
end

%% get true fingerprinting accuracy
% accs1 = fingerprinting(mats1(:,:,:,1),mats1(:,:,:,2));
accs1 = fingerprinting(vecs1(:,:,1),vecs1(:,:,2));
accs1 = mean(accs1);
% accs2 = fingerprinting(mats2(:,:,:,1),mats2(:,:,:,2));
accs2 = fingerprinting(vecs2(:,:,1),vecs2(:,:,2));
accs2 = mean(accs2);

%% run permutations and create distribution of discriminabilities and differences
% change shape of mats1 and mats2
null_accs1 = zeros(nperm,1);
null_accs2 = null_accs1;
for perm = 1:nperm
    perm
%     % determine which subjects to switch
%     idx = randsample([0 1],n_subj,true);
%     % create null mats
%     null_mats1 = cat(3,mats1(:,:,idx==0,:),mats2(:,:,idx==1,:));
%     null_mats2 = cat(3,mats1(:,:,idx==1,:),mats2(:,:,idx==0,:));
    % determine which scans to switch
    idx = randsample([0 1],n_subj*2,true);
    idx1 = idx(1:661);
    idx2 = idx(662:end);
    null_vecs1 = vecs1;
    null_vecs1(:,idx1==1,1) = vecs2(:,idx1==1,1);
    null_vecs1(:,idx2==1,2) = vecs2(:,idx2==1,2);
    null_vecs2 = vecs2;
    null_vecs2(:,idx1==1,1) = vecs1(:,idx1==1,1);
    null_vecs2(:,idx2==1,2) = vecs1(:,idx2==1,2);
%     null_mats1 = mats1;
%     null_mats1(:,:,idx1==1,1) = mats2(:,:,idx1==1,1);
%     null_mats1(:,:,idx2==1,2) = mats2(:,:,idx2==1,2);
%     null_mats2 = mats2;
%     null_mats2(:,:,idx1==1,1) = mats1(:,:,idx1==1,1);
%     null_mats2(:,:,idx2==1,2) = mats1(:,:,idx2==1,2);
    % compute null fingerprinting accuracies
%     curr_accs1 = fingerprinting(null_mats1(:,:,:,1),null_mats1(:,:,:,2));
    curr_accs1 = fingerprinting(null_vecs1(:,:,1),null_vecs1(:,:,2));
    null_accs1(perm) = mean(curr_accs1);
%     curr_accs2 = fingerprinting(null_mats2(:,:,:,1),null_mats2(:,:,:,2));
    curr_accs2 = fingerprinting(null_vecs2(:,:,1),null_vecs2(:,:,2));
    null_accs2(perm) = mean(curr_accs2);
end
% compute differences
null_dist = null_accs1-null_accs2;

%% calculate p_val
true_diff = accs1-accs2;
p = 1-length(find(true_diff>null_dist))/length(null_dist);

end