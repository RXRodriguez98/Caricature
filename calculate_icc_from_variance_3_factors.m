function icc = calculate_icc_from_variance_3_factors(subj_var,session_var,...
    run_var,subjXsession_var,subjXrun_var,sessionXrun_var,error_var)
%% takes variance components from nested repeated measures ANOVA and returns
%% the icc for each edge
% INPUTS %
% subj_var -- edges X n matrix (n is number of iterations)
% session_var -- edges X n matrix
% run_var -- edges X n matrix
% subjXsession_var -- edges X n matrix
% subjXrun_var -- edges X n matrix
% sessionXrun_var -- edges X n matrix
% error_var -- edges X n matrix

% OUTPUT %
% icc -- icc for each edge and iteration

%% set negatives to zeros
subj_var(subj_var<0) = 0;
session_var(session_var<0) = 0;
run_var(run_var<0) = 0;
subjXsession_var(subjXsession_var<0) = 0;
subjXrun_var(subjXrun_var<0) = 0;
sessionXrun_var(sessionXrun_var<0) = 0;
error_var(error_var<0) = 0;

%% compute icc
icc = subj_var./(subj_var+session_var+run_var+subjXsession_var+...
    subjXrun_var+sessionXrun_var+error_var);

end