function [p,null_seps1,null_seps2] = separation_stat(mats1,mats2,nperm)
% takes two datasets of repeated measures connectomes and does a
% permutation test to determine if one is more separable

% INPUTS %
% mats1 -- nodeXnodeXsubjXn_measures matrix
% mats2 -- nodeXnodeXsubjXn_measures matrix
% nperm -- number of permutations to use

% mats1 and mats2 need subjects in same order

% OUTPUTS %
% p -- pvalue associated with the one-tailed test that mats1 is more
%      discriminable than mats2
% null_seps1 -- distribution of separations for null_mats1
% null_seps2 -- distribution of separation for null_mats2

%% numbers
n_subj = size(mats1,3);

%% get true separability
s1 = separation(mats1);
s2 = separation(mats2);

%% run permutations and create distribution of separabilities and differences
null_seps1 = zeros(nperm,1);
null_seps2 = null_seps1;
for perm = 1:nperm
    perm
    tic
    % determine which scans to switch
%     idx = randsample([0 1],n_subj,true);
    null_mats1 = mats1;
    null_mats2 = mats2;
    for scan = 1:size(mats1,4)
        idx = randsample([0 1],n_subj,true);
        null_mats1(:,:,idx==1,scan) = mats2(:,:,idx==1,scan);
        null_mats2(:,:,idx==1,scan) = mats1(:,:,idx==1,scan);
    end
%     for sess = 1:4
%         idx = randsample([0 1],n_subj,true);
%         null_mats1(:,:,idx==1,(4*(sess-1)+1):(sess*4)) = ...
%             mats2(:,:,idx==1,(4*(sess-1)+1):(sess*4));
%         null_mats2(:,:,idx==1,(4*(sess-1)+1):(sess*4)) = ...
%             mats1(:,:,idx==1,(4*(sess-1)+1):(sess*4));
%     end
    % create null mats
%     null_mats1 = cat(3,mats1(:,:,idx==0,:),mats2(:,:,idx==1,:));
%     null_mats2 = cat(3,mats1(:,:,idx==1,:),mats2(:,:,idx==0,:));
    % compute null discs
    null_seps1(perm) = separation(null_mats1);
    null_seps2(perm) = separation(null_mats2);
    toc
end
% compute differences
null_dist = null_seps1-null_seps2;

%% calculate p_val
true_diff = s1-s2;
p = 1-length(find(true_diff>null_dist))/length(null_dist);

end