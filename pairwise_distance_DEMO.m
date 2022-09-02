%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%            PAIRWISE DISTANCE BETWEEN TWO MULTI-DIMENSIONAL SIGNALS           %
%                                                                              %
%                                September, 2018                               %
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

%   Produce two multi-dimensional signals where columns are dimensions and rows
%   are samples. Then compute the Euclidean distance of each point of one signal
%   to each point of the other signal.

%   Three methods are compared:
%      METHOD 1: The 'pdist2' function from the Statistics package.
%      METHOD 2: Embedded code (not a function) using nested loops.
%      METHOD 3: Embedded code (not a function) using vectorisation.
%      METHOD 4: Custom-made 'pdist_2' function (uses vectorisation).

% INSTRUCTIONS:

%   1) Modify the values indicated by an arrow (<---)
%   2) Run the script and enjoy.

%===============================================================================

clc
clear
close all

%===============================================================================

 signal_length = 99;   % <--- signal length
 signal_period = 33;   % <--- signal period
    dimensions =  3;   % <--- dimensions
 
% only for METHOD 2: 
      boundary = 0;    % <--- amount of distances from each value (0 = no bounds)

%-------------------------------------------------------------------------------
% Make multi-dimensional testing signals

x = [1:signal_length]'; % time grid
x = repmat(x,1,dimensions); % time grid for each dimension
ramp = x/signal_length/2; % incrementing ramp
sinusoid = sin(x*pi*2/(signal_period)); % plain sinusoid
signals = zeros(signal_length,dimensions,2); % initialise signals' matrix
for i = 1:2
  signals(:,:,i) = sinusoid + randn(signal_length,dimensions); % plain sinusoid with noise
  signals(:,:,i) = signals(:,:,i) .* ramp; % multiply with ramp
end

%% -------------------------------------------------------------------------------
% METHOD 1: distance matrix using function 'pdist2'

tic
distance_matrix_1 = pdist2(signals(:,:,1),signals(:,:,2));
t_d_1 = toc;

%% -------------------------------------------------------------------------------
% METHOD 2: distance matrix using nested loops

distance_matrix_2 = zeros(signal_length); % initialise distance matrix

% boundary's ending indexes:
if boundary == 0
  boundary = signal_length;
end
end_bounds = [1:signal_length]; % increasing indexes
end_bounds(end-boundary+1:end) = ...
  end_bounds(end-boundary+1:end) - [1:boundary]; % at the end the available space shrinks
end_bounds = end_bounds + boundary; % shift for length of boundary

% compute Euclidean distances;
tic
for row = 1:signal_length
  for col = row:end_bounds(row)
    distance_matrix_2(row,col) = sqrt( sum( ( signals(row,:,1) - signals(col,:,2) ).^2 ));
    distance_matrix_2(col,row) = sqrt( sum( ( signals(col,:,1) - signals(row,:,2) ).^2 ));
  end
end
t_d_2 = toc;

%% -------------------------------------------------------------------------------
% METHOD 3: distance matrix using vectorisation
%   This method considers the property of the square of a substraction:
%         (x-y)^2 = x^2 + y^2 - 2*x*y
%   This has to be repeated for all vectors (broadcasting).

% compute Euclidean distances;
tic 
sumsq_1 = sum(signals(:,:,1).^2, 2);
sumsq_2 = sum(signals(:,:,2).^2, 2);
sumsq_1_2 = repmat(sumsq_1,1,signal_length) + repmat(sumsq_2',signal_length,1); % broadcasting
distance_matrix_3 = sumsq_1_2 - 2 * signals(:,:,1) * (signals(:,:,2)');
distance_matrix_3 = sqrt(distance_matrix_3);
t_d_3 = toc;

%% -------------------------------------------------------------------------------
% METHOD 4: distance matrix using function 'pdist_2'

tic
distance_matrix_4 = pdist_2(signals(:,:,1),signals(:,:,2));
t_d_4 = toc;

%% -------------------------------------------------------------------------------
% Display data and results

% report computing times:
disp('Computing times for distance matrices')
disp(sprintf('  METHOD 1 (pdist2)   = %d s.',t_d_1))
disp(sprintf('  METHOD 2 (loops)    = %d s.',t_d_2))
disp(sprintf('  METHOD 3 (vectors)  = %d s.',t_d_3))
disp(sprintf('  METHOD 4 (pdist_2)  = %d s.',t_d_4))

% report comparison:
d_1_2 = sum(sum(abs(distance_matrix_1 - distance_matrix_2)));
d_1_3 = sum(sum(abs(distance_matrix_1 - distance_matrix_3)));
d_1_4 = sum(sum(abs(distance_matrix_1 - distance_matrix_4)));
disp(sprintf('\nSum of the absolute differences between distance matrices'))
disp(sprintf('  1 and 2 = %d',d_1_2))
disp(sprintf('  1 and 3 = %d',d_1_3))
disp(sprintf('  1 and 4 = %d',d_1_4))

% visualise:
subplot(4,1,1)
plot(signals(:,:,1))
title('signal 1')
subplot(4,1,2)
plot(signals(:,:,2))
title('signal 2')
subplot(2,4,5)
imagesc(distance_matrix_1)
title('pdist2')
subplot(2,4,6)
imagesc(distance_matrix_2)
title('loops')
subplot(2,4,7)
imagesc(distance_matrix_3)
title('vectors')
subplot(2,4,8)
imagesc(distance_matrix_4)
title('pdist_2','interpreter','none')