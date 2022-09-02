%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                                   BINSEQSI                                   %
%                             MONTECARLO SIMULATION                            %
%                                                                              %
%                                    May 2020                                  %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program has been tested with Matlab R2015a 
% and requires the function 'binseqsi'.

% ==============================================================================
% Instructions:

% If necessary, edit the parameters marked with an arrow like this: <---
% Run the thing, close your eyes and hope for the best.

% ==============================================================================
clc

   n_iter = 10000;  % <--- iterations (use qith caution)
        L = 1000;   % <--- length of binary sequences
disp_mode = 0;      % <--- 1 = update display at each iteration, 0 = display at the end

% ------------------------------------------------------------------------------

if ~exist('this_dir','var')
    this_dir = fileparts(matlab.desktop.editor.getActiveFilename);
    addpath(fileparts(this_dir))
end

   n_bins = 100;
lengths.a = NaN(1,n_iter);
lengths.b = lengths.a;

S = NaN(1,n_iter);
c = S;
p = S;
l = L-1;
n_sp = 3;

for i_iter = 1:n_sel_iter
   
   a = unique( ceil( rand(1, ceil(rand * l) ) * l ) ) ;
   
   b = unique( ceil( rand(1, ceil(rand * l) ) * l ) ) ;
   
   [S(i_iter), ~, c(i_iter), p(i_iter)] = binseqsi( a,b,L,0 );
   
   if disp_mode || (i_iter == n_iter)

       subplot(n_sp,1,1)
       hist(S,n_bins);
       xlabel('S')
       xlim([0,1])
       title(sprintf('Histograms at iteration %g of %g',i_iter,n_iter))
       
       subplot(n_sp,1,2)
       hist(c,n_bins);
       xlabel('c')
       xlim([0,1])
       
       subplot(n_sp,1,3)
       hist(f,n_bins);
       xlabel('f')
       xlim([0,1])
       
       drawnow
   end
end

[N,X] = hist(S,n_bins);
[~,i_max_N] = max(N);
dist_pk = X(i_max_N);
fprintf('\ndistribution peak = %.2g\n',dist_pk)

save([this_dir,'/binseqsi_montecarlo'],'n_iter','L','S')

%% ------------------------------------------------------------------------------

S_obs = 0.66; % <--- Observed S

n_ge = sum( S >= S_obs );
p_val = ( n_ge + 1 ) / ( n_iter + 1 );

fprintf('\n p-value = %.2g\n',p_val)

%% ------------------------------------------------------------------------------

p_val = 0.05; % <--- target upper p-value

step = 0.001;
n_ge = (( n_iter + 1 ) * p_val ) - 1;

S_obs = max(S) + step;
n_ge_1 = 0;

while n_ge_1 < n_ge

    S_obs = S_obs - step;
    n_ge_1 = sum( S >= S_obs );
end

fprintf('\n Observed S for target p-value = %.2g\n',S_obs)
