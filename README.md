# Caricature Repository
The functions in this repository are used to caricature resting-state fMRI time series data and perform subsequent analyses on the resulting connectomes. Note that discriminability.m is based on code that can be originally found at https://github.com/ebridge2/Discriminability/tree/master, fingerprinting.m and separation.m are based on code originally available at https://github.com/SNeuroble/fingerprinting, the functions calculate_icc_from_variance.m, calculate_icc_from_variance_3_factors.m, calculate_variance_comps.m, and calculate_variance_comps_3_factors.m are based on code originally found at https://github.com/SNeuroble/Multifactor_ICC, cpm_classifier_family_cv.m is an updated code based on that found in https://github.com/YaleMRRC/CPM/tree/master/matlab/func, and ridgeCPM.m is an updated code based on that found in https://github.com/mattrosenblatt7/trust_connectomes/tree/main/utils.

# Examples
1) Creating a task manifold using the CNP dataset  
    i) Load the timeseries data for the CNP task data as a cell with length equal to the number of tasks. Each element in the cell should be a timeseries 
       corresponding to a particular task with dimensions n_frames X n_nodes X n_subjects.  
   ii) Run the create_task_manifold.m function with the timeseries cell as the input.  
  iii) The outputs will be the eigenvector matrix and the variance explained by each eigenvector.  
  
2) Predicting age in the HCP dataset using the CNP-derived task manifold for caricaturing  
    i) Load the resting-state HCP data (size: n_frames X n_nodes X n_subjects), the eigenvectors from the CNP manifold, the family structure of the HCP subjects (can 
       set to a vector 1:n_subjects if not interested), and the age vector.  
   ii) Run the continuous_prediction.m function with these as inputs as well as a p-value threshold for feature selection and seeds for randomization. We set these 
       values to 0.05 and 1:1000 respectively.  
  iii) The output is a structure with predictions for both the standard and caricatured models.

3) Running fingerprinting in the HCP dataset using the CNP-derived task manifold for caricaturing  
    i) Load the LR phase-encoded resting-state HCP data as timeseries1 (size: n_frames X n_nodes X n_subjects), the RL phase-encoded resting-state HCP data as 
       timeseries2 (size: n_frames X n_nodes X n_subjects), and the eigenvectors from the CNP manifold.  
   ii) Run the run_fingerprinting.m function with these as inputs.  
  iii) The output is a structure with fingerprinting accuracy, within-participant similarities, between-participant similarities, and statistics information based 
       on permutation testing.
