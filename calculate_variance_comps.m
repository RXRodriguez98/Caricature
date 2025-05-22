function [var_subj,var_measure,var_error] = calculate_variance_comps(full_conn,ftbl)
%% Computes variance components for one-way repeated measures ANOVA
% INPUTS %
% full_conn -- all connectomes used (268X268X(num_measures*num_subj))
% ftbl -- factor table where first column is subject and second is measure number

% OUTPUTS % 
% var_subj -- variance due to subject
% var_measure -- variance due to scan "session" or "run"
% var_error -- variance due to error

%% number of subjects and measures
n_subj = length(unique(ftbl(:,1)));
n_measure = length(unique(ftbl(:,2)));

%% degrees of freedom
df_subj = n_subj-1;
df_measure = n_measure-1;
df_error = df_subj*df_measure;

%% calculate relevant means
grand_mean = mean(full_conn,3);
subj_means = zeros(268,268,n_subj);
for i = 1:n_subj
    subj_means(:,:,i) = mean(full_conn(:,:,ftbl(:,1)==i),3);
end
measure_means = zeros(268,268,n_measure);
for i = 1:n_measure
    measure_means(:,:,i) = mean(full_conn(:,:,ftbl(:,2)==i),3);
end

%% sums of squares
SS_total = sum((full_conn-grand_mean).^2,3);
SS_subj = sum((subj_means-grand_mean).^2,3)*n_measure;
SS_measure = sum((measure_means-grand_mean).^2,3)*n_subj;
SS_error = SS_total - SS_subj - SS_measure;

%% mean squares
MS_subj = SS_subj/df_subj;
MS_measure = SS_measure/df_measure;
MS_error = SS_error/df_error;

%% variance estimates
var_error = MS_error;
var_subj = (MS_subj-var_error)/n_measure;
var_measure = (MS_measure-var_error)/n_subj;

end