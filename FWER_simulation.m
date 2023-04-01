% This is a port of R code published at https://egap.org/resource/10-things-to-know-about-multiple-comparisons/
% retrieved on July 9 of 2021
% by Juan Ignacio Mendoza Garay, at the University of Jyväskylä

% Set a seed:
rng(1,'twister')

% Generate correlated outcomes:
% ( Outcomes are unrelated to treatment and all null hypotheses are true)
n = 1000; % observations
k = 100;  % outcomes
r = .7;   % expected correlation
s = 1;
this_sigma = ones( k ) * s * r ;
this_sigma(1:k+1:end) = s * ones(1,k);
outcomes = mvnrnd(zeros(1,k),this_sigma,n);

% Complete Random Assignment:
treatment = logical( mod( reshape(randperm(n), 1, n), 2 ));

% Conduct k hypothesis tests:
[~,p_obs]  = splitapply( @ttest , outcomes(treatment,:) , outcomes(~treatment,:) , 1:k );

% Simulate under the sharp null:
n_iters = 1000;
many_ps = NaN(k,n_iters);
for i_iter = 1:n_iters

   treatment_sim = logical( mod( reshape(randperm(n), 1, n), 2 ));
   
   for i_outcome = 1:k
       
       [~,many_ps(i_outcome,i_iter)] = ttest( outcomes(treatment_sim,i_outcome) , outcomes(~treatment_sim,i_outcome) );
   end
end

% Obtain the Type I error rate for a series of thresholds:
n_thr = 1000;
thresholds = 0 : 0.05 / (n_thr-1) : 0.05 ;
type_I_rate = NaN(n_thr,1);
for i_thr = 1:n_thr
    counts = zeros(n_iters,1);
   for i_iter = 1:n_iters 
       counts(i_iter) = any( many_ps(:,i_iter) <= thresholds(i_thr) );
   end
   type_I_rate(i_thr) = mean(counts);
end

% Find the largest threshold that yields an alpha type I error rate:
target_p_value = max( thresholds( type_I_rate <= 0.05 ) )

% Apply target p_value to observed p_values:
sig_simulated = p_obs <= target_p_value;

% Compare to raw p-values
sig = p_obs <= 0.05;
