%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                          Second-Order Isomporphism                           %
%                                                                              %
%                              Juan Ignacio Mendoza Garay and Yorgos Diapoulis %
%                                   Music Department - University of Jyv?skyl? %
%                                                                  March, 2016 %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script has been tested with Matlab R2015a

% ==============================================================================
% Description:

% This script shows the workings of second order isomorphism.

% ==============================================================================
% Initialisation:

clc
clear
close all

% ==============================================================================

W = 50; % <--- input size of matrix

% ------------------------------------------------------------------------------
% comment/uncomment the options:

%...............................................................................
% original_A = magic(W);
% %...............................................................................
% original_B = kron(...
%               [3 2 1 2 3;
%               3 4 1 1 2;
%               1 2 4 1 5;
%               1 2 4 5 5;
%               12 4 2 1 2]...
%               ,ones(W/5));
%...............................................................................     
[gauss_x gauss_y] = meshgrid((-(W-1)/2):((W-1)/2), (-(W-1)/2):((W-1)/2));
gauss_bell = exp( -((2*pi*gauss_x/W).^2) / 2 - ((2*pi*gauss_y/W).^2) / 2); % gaussian bell

original_A = gauss_bell;
original_B = gauss_bell .* (kron([-1, 1; 1,  -1],ones(W/2))); % gaussian tapered checkerboard kernel (c.f. Foote & Cooper, 2003)
% ------------------------------------------------------------------------------

cov_A = cov(original_A);
corrcoef_cov_A = corrcoef(cov_A);

cov_B = cov(original_B);
corrcoef_cov_B = corrcoef(cov_B);

clc
R = corrcoef(corrcoef_cov_A,corrcoef_cov_B)

% ------------------------------------------------------------------------------

scrsz = get(groot,'ScreenSize');
fig = figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2]); 

subplot(2,3,1)
imagesc(original_A)
axis square
title('original A')
subplot(2,3,2)
imagesc(cov_A)
axis square
title('covariance of A')
subplot(2,3,3)
imagesc(corrcoef_cov_A)
axis square
title('corr. of cov. of A')

subplot(2,3,4)
imagesc(original_B)
axis square
title('original B')
subplot(2,3,5)
imagesc(cov_B)
axis square
title('covariance of B')
subplot(2,3,6)
imagesc(corrcoef_cov_B)
axis square
title('corr. of cov. of B')

