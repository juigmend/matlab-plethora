%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                      GAUSSIAN KERNEL DENSITY ESTIMATION                      %
%                                                                              %
%                                 DEMONSTRATION                                %
%                                                                              %
%                                  28 May 2018                                 %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Description:
%   Computes kernel density estimation with a Gaussian kernel
%   (Silverman, 1986).

% Instructions:
%   1) Edit the variables indicated with an arrow (<---)
%   2) Run the script and enjoy.

%%==============================================================================

close all
clear
clc

%%==============================================================================
% Define variables:

resfact = 1; % <--- resolution factor

 data_points = resfact * [10,50,54,82,85,87]; % <--- data points
  datalength = resfact * 100;  % <------------------ length of data
      sigma  = resfact * 1; % <--------------------- bandwidth

%  data_points = resfact * [2, 5,5, 8,9]; % <------- data points
%   datalength = resfact * 10;  % <----------------- length of data
%        sigma = resfact * 1/4; % <----------------- bandwidth

%  data_points = resfact * [434, 440,442, 670, 674, 828, 830, 830, 847, 852, 961, 1134, 1132]; % <--- data points
%   datalength = resfact * 1227; % <----------------------------------------------------------------- length of data
%        sigma = resfact * 6; % <-------------------------------------------------------------------- bandwidth

%-------------------------------------------------------------------------------
% METHOD 1: Embedded code (no toolboxes needed)
tic

% If the width is even, add one so that the kernel is odd.
% This makes the kernel to be symmetrical. Also it facilitates to test for peaks
% by taking the difference between each sample, given enough resoultion.

window_width = datalength; % this gives more than enough numerical precision, although theoretically it should be datalength * 2
if rem(window_width,2) ~= 1
    window_width = window_width + 1;
end

%...............................................................................
% Make a Gaussian Kernel (comment/uncomment method 1.1 or 1.2)

% 1.1 embedded code:
n = linspace( -(window_width-1)/2, (window_width-1)/2, window_width )';
kernel = exp( -(n.^2) / (2*sigma^2) ) ; 

% 1.2 'gausswin' function from the Signal Processing Toolbox:
% alpha = (window_width -1) / (2*sigma);
% kernel = gausswin(window_width,alpha);  

%...............................................................................

kernel = ( 1 / ( sqrt( 2 * pi ) * sigma ) ) * kernel; % this makes the kernel to have unit area

% Each data point is treated as the index of a one amongst a vector of zeroes
% The result is a pulse train (a pulse is a one in a vector of zeroes).
% Data points that are the same point will be added.

pulse_train = zeros(1,datalength);
for i = data_points
    pulse_train(i) = pulse_train(i)+1; 
end
kde_embedded = conv(pulse_train,kernel,'same'); % KDE : convolve the Gaussian kernel with the pulse train
kde_embedded = kde_embedded * max(pulse_train) / max(kde_embedded); % rescale
times(1) = toc;

%-------------------------------------------------------------------------------
% METHOD 2: 'ksdensity' function from the Statistics and Machine Learning Toolbox by Mathworks

tic
kde_mathworks = ksdensity(data_points,1:datalength,'bandwidth',sigma,'kernel','normal'); % KDE : Mathworks' implementation
kde_mathworks = kde_mathworks * max(pulse_train) / max(kde_mathworks); % rescale
times(2) = toc;

%-------------------------------------------------------------------------------
% Correlation test and visualisation: 

disp(sprintf('Time method 1 (embedded)  = %g s.',times(1)))
disp(sprintf('Time method 2 (Mathworks) = %g s.',times(2)))
rho = corr(kde_embedded',kde_mathworks');
disp(sprintf('Linear correlation between KDE methods = %f',rho))

figure
plot(kernel,'linewidth',2,'linestyle','-')
xlim([1,window_width])
title('kernel','FontSize',16)
figure('Position',[0,0,1600,300])

hold on
plot(pulse_train,'g')
rep_data_points = repmat(data_points,2,1);
y = repmat([0,max(pulse_train)/10],size(rep_data_points,2),1)';
line(rep_data_points,y,'Color','k','linewidth',4,'linestyle','-');
plot(kde_embedded,'r','linewidth',2,'linestyle','-')
plot(kde_mathworks,'b','linewidth',4,'linestyle','--')
xlim([1,datalength])
title('data points (black), pulse train (green), embedded code KDE (red) and Mathwork''s implementation KDE (blue)',...
    'FontSize',16)

