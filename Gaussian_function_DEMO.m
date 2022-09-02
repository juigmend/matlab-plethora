%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                               GAUSSIAN FUNCTION                              %
%                                                                              %
%                                 DEMONSTRATION                                %
%                                                                              %
%                                  27 May 2018                                 %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instructions:
%   1) Edit the variables indicated with an arrow (<---)
%   2) Run the script and enjoy.

%%==============================================================================

close all
clear
clc

%%==============================================================================
  
  window_width = 10; % <--------------- length of window
   start_sigma = 1/sqrt(2*pi);   % <--- kernel's bandwidth starting value
interval_sigma = 0.5/sqrt(2*pi); % <--- interval
  ending_sigma = 1;  % <--------------- kernel's bandwidth ending value
        resfac = 10; % <--------------- resolution

% bandwidth is equal to the standard deviation (sigma)
% sigma^2 = variance
%-------------------------------------------------------------------------------  

if rem(window_width,2) ~= 1
        window_width = window_width + 1;
end
sigma_values = [ start_sigma : interval_sigma : ending_sigma ];
n = linspace( -(window_width-1)/2, (window_width-1)/2, window_width*resfac )';
kernel_h1 = zeros(length(n),1);
kernel_a1 = kernel_h1;
figpos = get(0,'ScreenSize');
figpos(1) = figpos(3)/2;
figpos(3) = figpos(1);
figure('Position',figpos,'Units','Pixel')
subplot(2,1,1)
hold on
xlim([n(1),n(end)])
set(gca,'Fontsize',16);
title('Gaussian Kernel, unit height','Fontsize',16)
subplot(2,1,2)
hold on
xlim([n(1),n(end)])
yticks = [0:0.2:1];
lyt = length(yticks);
yticks(3) = 1/sqrt(2*pi);
yticklabels = cell(1,lyt);
for i = 1:lyt
    yticklabels{i} = num2str(yticks(i));
end
yticklabels{3} = '1/sqrt(2*pi)';
set(gca,'ytick',yticks,'yticklabel',yticklabels,'Fontsize',16)
title('Gaussian Kernel, unit area','Fontsize',16)
std_h_linestyle = ':';
std_h_linewidth = 1;

for i = 1:length(sigma_values)

    alpha = (window_width -1) / (2*sigma_values(i));

    %...............................................................................
    % Make a Gaussian (normal) Kernel:
    
    kernel_h1 = exp( -(n.^2) / (2 * sigma_values(i)^2) ) ; % height = 1
%     kernel_h1 = gausswin(window_width,alpha);     % height = 1 , 'gausswin' is part of the Signal Processing Toolbox
    
    height_a1 = 1 / ( sqrt( 2 * pi ) * sigma_values(i) ) ;
    kernel_a1 = height_a1 * kernel_h1 ; % area = 1
     
    %...............................................................................
    
    if (sigma_values(i) > 0.99) && (sigma_values(i) < 1.01) % enhance std and height lines when sigma is one or very close
        std_h_linestyle = '--';
        std_h_linewidth = 2;
    end
    
    subplot(2,1,1)
    plot(n,kernel_h1,'linewidth',3)
    line(repmat([-sigma_values(i),sigma_values(i)],2,1),repmat([kernel_h1(ceil(resfac*(sigma_values(i)+(window_width)/2))+1);0],1,2),...
        'Color','k','linestyle',std_h_linestyle,'linewidth',std_h_linewidth); % std

    subplot(2,1,2)
    plot(n,kernel_a1,'linewidth',3)
    line(repmat([-sigma_values(i),sigma_values(i)],2,1),repmat([kernel_a1(ceil(resfac*(sigma_values(i)+(window_width)/2))+1);0],1,2),...
        'Color','k','linestyle',std_h_linestyle,'linewidth',std_h_linewidth); % std
    line([n(1),0],[height_a1;height_a1],...
        'Color','k','linestyle',std_h_linestyle,'linewidth',std_h_linewidth); % height
    pause(1)
end


