function icc = calculate_icc_from_variance(subj_var,measure_var,error_var)
%% takes variance components from one-way repeated measures ANOVA and returns
%% the icc for each edge
% INPUTS %
% subj_var -- edges X n matrix (n is number of iterations)
% measure_var -- edges X n matrix
% error_var -- edges X n matrix

% OUTPUT %
% icc -- icc for each edge and iteration

%% set negatives to zero
subj_var(subj_var<0) = 0;
measure_var(measure_var<0) = 0;
error_var(error_var<0) = 0;

%% compute ICC
icc = subj_var./(subj_var + measure_var + error_var);

end