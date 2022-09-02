%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                   2D Gaussian Tapered Checkerboard Kernel                    %
%                                                                              %
%                                 DEMONSTRATION                                %
%                                                                              %
%                                  29 May 2018                                 %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script has been tested with Matlab R2015a

% ==============================================================================
% Description:

% Produce a Gaussian Kernel, then transform it so that it has positive and 
% negative zones like a checkerboard of 4 squares (c.f. Foote & Cooper, 2003).
% Plot each stage.

% ==============================================================================
% Initialisation:

clc
clear
close all

scrsz = get(groot,'ScreenSize');
fig = figure('Position',[1 scrsz(4)-scrsz(4)*1/2 scrsz(3) scrsz(4)*1/2]); % adjust position

%% % ----------------------------------------------------------------------------
% 2D Gaussian Kernel:

window_width = 30;
[kernel_x, kernel_y] = meshgrid((-(window_width-1)/2):((window_width-1)/2), (-(window_width-1)/2):((window_width-1)/2));
variance = (window_width-1)/2;
height =  1 / ( ( 2 * pi * variance) ) ; % this produces unit volume
gaussian_kernel = height *  exp( - ( (kernel_x.^2) / ( 2 * variance) ) - ( (kernel_y.^2) / ( 2 * variance) ) );

subplot(1,3,1)
surf(gaussian_kernel)
title('2D Gaussian Kernel')

%% % ----------------------------------------------------------------------------
% Checkerboard Matrix:

checkerboard = kron([-1, 1; 1,  -1],ones(window_width/2));
 
subplot(1,3,2)
imagesc(checkerboard)
title('Checkerboard Matrix')

%% % ----------------------------------------------------------------------------
% 2D Gaussian-tapered Checkerboard Kernel:

gaussian_checkerboard_kernel = gaussian_kernel .* checkerboard;

subplot(1,3,3)
surf(gaussian_checkerboard_kernel)
title('2D Gaussian-tapered Checkerboard Kernel')


