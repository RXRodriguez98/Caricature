function [var_subj,var_session,var_run,var_subjXsession,...
    var_subjXrun,var_sessionXrun,var_error] = ...
    calculate_variance_comps_3_factors(full_conn,ftbl)
%% Computes variance components for one-way repeated measures ANOVA
% INPUTS %
% full_conn -- all connectomes used (268X268X(num_session*num_subj*n_runs))
% ftbl -- factor table where first column is subject, second is session,
% and third is run

% OUTPUTS % 
% var_subj -- variance due to subject
% var_session -- variance due to session
% var_run -- variance due to run
% var_{factor_i}X{factor_j} -- variance due to interaction between factors
% var_error -- variance due to error

%% number of subjects, sessions, and runs
n_subj = length(unique(ftbl(:,1)));
n_sessions = length(unique(ftbl(:,2)));
n_runs = length(unique(ftbl(:,3)));

%% degrees of freedom
df_subj = n_subj-1;
df_session = n_sessions-1;
df_run = n_runs-1;
df_subjXsession = df_subj*df_session;
df_subjXrun = df_subj*df_run;
df_sessionXrun = df_session*df_run;
df_error = df_subj*df_session*df_run;

%% calculate relevant means
grand_mean = mean(full_conn,3);
subj_means = zeros(268,268,n_subj);
for i = 1:n_subj
    subj_means(:,:,i) = mean(full_conn(:,:,ftbl(:,1)==i),3);
end
session_means = zeros(268,268,n_sessions);
for i = 1:n_sessions
    session_means(:,:,i) = mean(full_conn(:,:,ftbl(:,2)==i),3);
end
run_means = zeros(268,268,n_runs);
for i = 1:n_runs
    run_means(:,:,i) = mean(full_conn(:,:,ftbl(:,3)==i),3);
end
subjXsession_rows = unique(ftbl(:,1:2),'rows');
subjXsession_means = zeros(268,268,n_subj*n_sessions);
for i = 1:(n_subj*n_sessions)
    subjXsession_means(:,:,i) = ...
        mean(full_conn(:,:,find(ismember(ftbl(:,1:2),subjXsession_rows(i,:),'rows'))),3);
end
subjXrun_rows = unique(ftbl(:,[1 3]),'rows');
subjXrun_means = zeros(268,268,n_subj*n_runs);
for i = 1:(n_subj*n_runs)
    subjXrun_means(:,:,i) = ...
        mean(full_conn(:,:,find(ismember(ftbl(:,[1 3]),subjXrun_rows(i,:),'rows'))),3);
end
sessionXrun_rows = unique(ftbl(:,2:3),'rows');
sessionXrun_means = zeros(268,268,n_sessions*n_runs);
for i = 1:(n_sessions*n_runs)
    sessionXrun_means(:,:,i) = ...
        mean(full_conn(:,:,find(ismember(ftbl(:,2:3),sessionXrun_rows(i,:),'rows'))),3);
end

%% sums of squares
SS_total = sum((full_conn-grand_mean).^2,3);
SS_subj = sum((subj_means-grand_mean).^2,3)*n_sessions*n_runs;
SS_session = sum((session_means-grand_mean).^2,3)*n_subj*n_runs;
SS_run = sum((run_means-grand_mean).^2,3)*n_subj*n_sessions;
% interactions
SS_subjXsession = zeros(268);
for i = 1:(n_subj*n_sessions)
    SS_subjXsession = SS_subjXsession + ...
        (subjXsession_means(:,:,i)-subj_means(:,:,subjXsession_rows(i,1))-...
        session_means(:,:,subjXsession_rows(i,2))+grand_mean).^2;
end
SS_subjXsession = SS_subjXsession*n_runs;
SS_subjXrun = zeros(268);
for i = 1:(n_subj*n_runs)
    SS_subjXrun = SS_subjXrun + ...
        (subjXrun_means(:,:,i)-subj_means(:,:,subjXrun_rows(i,1))-...
        run_means(:,:,subjXrun_rows(i,2))+grand_mean).^2;
end
SS_subjXrun = SS_subjXrun*n_sessions;
SS_sessionXrun = zeros(268);
for i = 1:(n_sessions*n_runs)
    SS_sessionXrun = SS_sessionXrun + ...
        (sessionXrun_means(:,:,i)-session_means(:,:,sessionXrun_rows(i,1))-...
        run_means(:,:,sessionXrun_rows(i,2))+grand_mean).^2;
end
SS_sessionXrun = SS_sessionXrun*n_subj;

% error
SS_error = SS_total - SS_subj - SS_session - SS_run - SS_subjXsession - ...
    SS_subjXrun - SS_sessionXrun;

%% mean squares
MS_subj = SS_subj/df_subj;
MS_session = SS_session/df_session;
MS_run = SS_run/df_run;
MS_subjXsession = SS_subjXsession/df_subjXsession;
MS_subjXrun = SS_subjXrun/df_subjXrun;
MS_sessionXrun = SS_sessionXrun/df_sessionXrun;
MS_error = SS_error/df_error;

%% variance estimates
var_error = MS_error;
var_subjXsession = (MS_subjXsession-var_error)/n_runs;
var_subjXrun = (MS_subjXrun-var_error)/n_sessions;
var_sessionXrun = (MS_sessionXrun-var_error)/n_subj;
var_subj = (MS_subj-var_error-n_runs*var_subjXsession-...
    n_sessions*var_subjXrun)/(n_sessions*n_runs);
var_session = (MS_session-var_error-n_runs*var_subjXsession-...
    n_subj*var_sessionXrun)/(n_subj*n_runs);
var_run = (MS_run-var_error-n_sessions*var_subjXrun-...
    n_subj*var_sessionXrun)/(n_subj*n_sessions);