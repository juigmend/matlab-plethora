%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                                 NOVELTY SCORE                                %
%                           OF A ONE-DIMENSIONAL SIGNAL                        %
%                      COMPARISON OF METHODS FULL AND FAST                     %
%                                                                              %
%                                 October, 2018                                %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tested on Octave 4 and Matlab 2015a

%===============================================================================

% DESCRIPTION:

%   Produce a signal with some randomness, then a self-similarity matrix of it.
%   The matrix is bounded, to reduce computing time.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE THIS:

%   Then convolve a Gaussian-tapered checkerboard kernel along the diagonal of 
%   the self-similarity matrix. The result is a novelty score.
%   (c.f. Foote & Cooper, 2003)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUCTIONS:

%   1) Modify the values indicated by an arrow (<---)
%   2) Run the script and enjoy.

%%%=============================================================================

clc
clear
close all

%%%=============================================================================

   signal_type ='sqr'; % <--- signal carrier type ('sin' or 'sqr')
 signal_length = 99;   % <--- signal length
 signal_period = 33;   % <--- signal period
  kernel_width = 11;   % <--- kernels' width
       full_sw = 1;    % <--- 0 = upper and lower triangles; 0 = upper triangle

%%%-----------------------------------------------------------------------------
% MAKE TEST SIGNAL

if strcmp(signal_type,'sin')
  x = [1:signal_length]; % grid
  signal = sin(x*pi*2/(signal_period) + randn(1)/8) + (sin(x*pi/4*randn(1)*2)/6+randn(1)/100);
elseif strcmp(signal_type,'sqr')
  signal = ones(1,signal_period);
  signal(ceil(signal_period/2):end) = -1;
  signal = repmat(signal,1,ceil(signal_length/signal_period));
  signal = signal(1:signal_length) + randn(1,signal_length)/4;
end
signal = signal .* [1:signal_length]/signal_length/2; 

%%%-----------------------------------------------------------------------------
% SET KERNELS' PARAMETERS

gaussian_height = 1/(sqrt(2*pi)*kernel_width);
gaussian_variance = (kernel_width/8)^2;	

%%%-----------------------------------------------------------------------------
% NOVELTY SCORE - METHOD 1: FULL (Foote & Cooper, 2003)

tic

selfsim_matrix = zeros(signal_length); % initialise self-similarity matrix

% boundary's ending indexes:
boundary = kernel_width - 1;
end_bounds = [1:signal_length]; % increasing indexes
end_bounds(end-boundary+1:end) = ...
  end_bounds(end-boundary+1:end) - [1:boundary]; % at the end the available space shrinks
end_bounds = end_bounds + boundary; % shift for length of boundary

% compute distances;
for row = 1:signal_length
  for col = row+1:end_bounds(row) % start at diagonal + 1; end at corresponding index
    selfsim_matrix(row,col) = abs( signal(row) - signal(col) ); % Euclidean distance
  end
end

if full_sw
  selfsim_matrix = selfsim_matrix + flipud(rot90(selfsim_matrix)); % upper and lower triangles
end

% add inner-margin-average padding to beginning and ending of similarity matrix's diagonal
if rem(kernel_width,2)
  kernel_width = kernel_width + 1; % has to be even to make the checkerboard (method 1)
end
length_selfsim_matrix_avpad = length(signal) + kernel_width;
margin_selfsim_matrix_avpad = fix(kernel_width/2);
selfsim_matrix_avpad = zeros(length_selfsim_matrix_avpad,length_selfsim_matrix_avpad);
beginning_average = selfsim_matrix(1:margin_selfsim_matrix_avpad,...
  1:margin_selfsim_matrix_avpad);
beginning_average = mean(beginning_average(:));
selfsim_matrix_avpad(1:length_selfsim_matrix_avpad,1:length_selfsim_matrix_avpad) = ...
  beginning_average;
ending_average = selfsim_matrix( length(signal) - margin_selfsim_matrix_avpad : length(signal),...
  length(signal) - margin_selfsim_matrix_avpad : length(signal) );
ending_average = mean(ending_average(:));
selfsim_matrix_avpad(1:length_selfsim_matrix_avpad,1:length_selfsim_matrix_avpad) = ...
  ending_average;
selfsim_matrix_avpad( margin_selfsim_matrix_avpad+1:length(signal)+margin_selfsim_matrix_avpad,...
  margin_selfsim_matrix_avpad+1:length(signal)+margin_selfsim_matrix_avpad ) = ...
  selfsim_matrix;
	
% Gaussian Checkerboard Kernel:
[mesh_2D_x, mesh_2D_y] = meshgrid((-(kernel_width-1)/2):((kernel_width-1)/2),...
  (-(kernel_width-1)/2):((kernel_width-1)/2));
kernel_2D.gauss = exp( -(mesh_2D_x.^2) / (2*gaussian_variance) ...
  - (mesh_2D_y.^2) / (2*gaussian_variance));
kernel_2D.gauss = kernel_2D.gauss .* gaussian_height;
kernel_2D.gausscb = kernel_2D.gauss .* kron([-1,1;1,-1],ones(kernel_width/2));

% convolve Gaussian Checkerboard kernel along diagonal of similarity matrix:
nov_full.novelty_padded = zeros(1,length_selfsim_matrix_avpad); % init
for i = [1:length(signal)]
  window_start = i;
  window_end = window_start + kernel_width - 1;
  function_window = selfsim_matrix_avpad(window_start:window_end,window_start:window_end);
  function_output = sum(sum(function_window.*kernel_2D.gausscb));
  nov_full.novelty_padded(1,i + margin_selfsim_matrix_avpad) = function_output;
end

% output of same length (get rid of outer margins):
nov_full.novelty_same = nov_full.novelty_padded(margin_selfsim_matrix_avpad + 1 :...
  end - margin_selfsim_matrix_avpad );
 
% normalise it, don't criticise it: 
nov_full.novelty_same_norm = nov_full.novelty_same - min(min(nov_full.novelty_same));
nov_full.novelty_same_norm = nov_full.novelty_same_norm / max(max(nov_full.novelty_same_norm));

timer_1 = toc;

%%%-----------------------------------------------------------------------------
% NOVELTY SCORE - METHOD 2: FAST

tic

if rem(kernel_width,2) == 0
  kernel_width = kernel_width - 1; % has to be uneven so that it is symmetric (method 2)
end

% 1D Gaussian kernel:
kernel_1D_x =((-(kernel_width-1)/2):((kernel_width-1)/2));
kernel_1D.gauss = exp( -(kernel_1D_x.^2) / (2*gaussian_variance));
kernel_1D.gauss = kernel_1D.gauss .* gaussian_height;
kernel_1D.gausscb = kernel_1D.gauss;
kernel_1D.gausscb(ceil(kernel_width/2):end) = -( kernel_1D.gauss(ceil(kernel_width/2):end) ); % make second half of kernel negative

% compute novelty score:
nov_fast.novelty_same = conv(signal,kernel_1D.gausscb,'same'); % convolve test signal with kernel
nov_fast.novelty_same = abs(nov_fast.novelty_same); % absolute value
nov_fast.novelty_same_norm = nov_fast.novelty_same - min(nov_fast.novelty_same(:)); % normalise it
nov_fast.novelty_same_norm = nov_fast.novelty_same_norm./ max(nov_fast.novelty_same_norm(:)); % don't criticise it

timer_2 = toc;

%%% ----------------------------------------------------------------------------
% COMPUTE CORRELATION COEFFICIENT:

R = corrcoef(nov_full.novelty_same_norm, nov_fast.novelty_same_norm);

%%%-----------------------------------------------------------------------------
% VISUALISATION

figure;
subplot(14,1,1)
plot(signal)
title('signal')
subplot(14,1,3:10)
imagesc(selfsim_matrix)
title('self-similarity matrix')
subplot(14,1,12)
plot(nov_full.novelty_same_norm)
title('full novelty score')
subplot(14,1,14)
plot(nov_fast.novelty_same_norm)
title('fast novelty score')

%%% ----------------------------------------------------------------------------
% DISPLAY CORRELATION COEFFICIENT AND TIME:

disp('----------------------------------------------------------')
disp('Computing times of Novelty Scores:')
disp(' ')
disp(['METHOD 1 = ',num2str(timer_1),' sec.'])
disp(['METHOD 2 = ',num2str(timer_2),' sec.'])
disp(' ')
disp('----------------------------------------------------------')
disp('Correlation coefficient:')
disp(' ')
disp(['R = ',num2str(R(1,2))])
disp(' ')
disp('----------------------------------------------------------')