%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                        Segmentation of an AMFM Signal                        %
%                                                                              %
%                                                                              %
%                                                   Juan Ignacio Mendoza Garay %
%                                                             doctoral student %
%                                   Music Department - University of Jyv?skyl? %
%                                                                January, 2017 %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script has been tested in Matlab R2015a.

% ==============================================================================
% Description

% Segmentation of a signal using the following program:

% 1) Generate a signal by modulating amplitude and frequency of a sine wave.
% 2) Make a self-similarity matrix of the signal.
% 3) Compute a novelty curve form the self-similarity matrix and extract peaks
%    of said curve above a threshold.
% 4) Measure the distance between the segments bounded by the peaks.
% 5) Cluster the segments by their distance.

% ==============================================================================
% Instructions

% Parameters that may be changed for exploration are marked with a long arrow
% like this:

% variable_name = value; <------------------------------------------------------ parameter description  

% ==============================================================================
% Initialisation:

clc
clear
close all
% addpath(genpath('PATH')) % <-------------------------------------------------- add path of scripts or enclosing folder

scrsz = get(groot,'ScreenSize');
figure('Position',[1 0 scrsz(3) scrsz(4)*3/4])
colormap(parula)
signal_parameters = ... % default parameters to generate signal (5 options)
    [80  80   100  50  90;   % FM modulating frequency
    1600 1600 1000 500 800;  % FM carrier frequency
    20   20   10   10  20;   % FM strength 
    40   160  50   25  40;   % AM modulating frequency
    0.5  0.5  1    1   0.5]; % AM strength (between 0 and 1)

% ==============================================================================
%%% ----------------------------------------------------------------------------
% 1) AMFM Signal
%    Generate a signal containing distinct patterns by modulating a
%    sine wave in frequency and amplitude, then plot.

% Use preset signal parameters:
use_preset = 1; % <------------------------------------------------------------- Signal parameters preset to use

% ... or input parameters:
save_to = 5; % <---------------------------------------------------------------- Location where to save these parameters
signal_parameters(1,save_to) = 90;  % <----------------------------------------- FM modulating frequency
signal_parameters(2,save_to) = 800; % <----------------------------------------- FM carrier frequency
signal_parameters(3,save_to) = 20;  % <----------------------------------------- FM strength 
signal_parameters(4,save_to) = 40;  % <----------------------------------------- AM modulating frequency
signal_parameters(5,save_to) = 0.5; % <----------------------------------------- AM strength (between 0 and 1)

sample_time = 0.0001; 
total_time = 0.1-sample_time; 
t = [0:sample_time:total_time]; % time grid

FM_signal = sin(...
    2*pi*signal_parameters(2,use_preset)*t+...
    (signal_parameters(3,use_preset).*sin(2*pi*signal_parameters(1,use_preset)*t))...
    );

AM_FM_signal = (1+signal_parameters(5,use_preset)*sin(2*pi*signal_parameters(4,use_preset)*t)).*FM_signal;

subplot(5,1,1)
plot(t,AM_FM_signal)
title('1) AMFM Signal')

%%% ----------------------------------------------------------------------------
% 2) Segmentation Distance Matrix
%    Generate a distance matrix of even and small-windowed distances 
%    to find boundaries of patterns and display a heatmap.

segmentation_bandwidth = size(AM_FM_signal,2)*0.02; % <--- size of the moving window (samples)
segmentation_hop = segmentation_bandwidth/2; % <------------ hop size (samples)
segmentation_bound = 10; % <------------------------------ boundary for closest DTW distance (percent) 
% Notes:
% - segmentation_bandwidth: it seems that 1% of total length works best.
% - segmentation_hop: it seems that no overlap works best.
% - segmentation_bound is a percentage of the DTW moving window, it restricts the maximum
% closest distance so that computation is faster.

% segmentation_bandwidth = 100; % <------------------------- size of the moving window (samples)
% segmentation_hop = 10; % <------------------------------- hop size (samples)
% segmentation_bound = 10; % <----------------------------- boundary for closest DTW distance (percent)

% segmentation_bandwidth = 100; % <------------------------- size of the moving window (samples)
% segmentation_hop = 5; % <------------------------------- hop size (samples)
% segmentation_bound = 10; % <----------------------------- boundary for closest DTW distance (percent) 

amount_segmentation_windows = fix((size(AM_FM_signal,2)-segmentation_bandwidth)/segmentation_hop + 1);
segmentation_distance_matrix = zeros(amount_segmentation_windows);
fractot = 100/(amount_segmentation_windows-1); 
laps = 1;
clc
disp('----------------------------------------------------------')
disp(fprintf('Progress of segmentation distance matrix = '));
tic
for row = 2:amount_segmentation_windows % This could be done with mcwindow from the Mocap Toolbox? <--- !!!
   a =  AM_FM_signal( segmentation_hop*(row-1)+1 : segmentation_bandwidth+(segmentation_hop*(row-1)) );
   for column = 1:(row-1)
       b = AM_FM_signal( segmentation_hop*(column-1)+1 : segmentation_bandwidth+(segmentation_hop*(column-1)) ); 
       
       %segmentation_distance_matrix(row,column) = pdist([a;b],'spearman');
       % options: euclidean cityblock minkowski chebychev mahalanobis cosine correlation spearman jaccard
       
       %segmentation_distance_matrix(row,column) = (a*b'); % product distance (fast and inaccurate)
       
       segmentation_distance_matrix(row,column) = norm(a-b); % euclidean distance (fast and a bit inaccurate)
       
       %segmentation_distance_matrix(row,column) = dtw(a',b',segmentation_bandwidth*segmentation_bound/100); % dynamic time warping (slow and accurate)

       % Note:
       % DTW and Euclidean give very similar results, in practice thay can be regarded as with no difference.
       % It seems pertinent to use Euclidean as it is faster.
       
       segmentation_distance_matrix(column,row) = segmentation_distance_matrix(row,column);
   end
   display_progress = round(fractot*laps);
   fprintf(repmat('\b',1,numel(mat2str(display_progress))+3))
   fprintf(' %d%%',display_progress);
   fprintf('\n')
   laps = laps+1;
end
timer_1 = toc;
disp(['Segmentation distance matrix computing time = ',num2str(timer_1),' sec.'])
disp('----------------------------------------------------------')

% add a "tail" to the matrix so that scaled it matches the length of the original signal:
segmentation_tail = (segmentation_bandwidth/segmentation_hop)-1;
segmentation_distance_matrix_display = zeros(amount_segmentation_windows,amount_segmentation_windows+segmentation_tail);
segmentation_distance_matrix_display(:,1:amount_segmentation_windows) = segmentation_distance_matrix;
segmentation_distance_matrix_display(:,amount_segmentation_windows+1:end) = repmat(segmentation_distance_matrix(:,end),1,segmentation_tail);

hold off
subplot(5,1,2)
% imagesc(segmentation_distance_matrix)
imagesc(segmentation_distance_matrix_display)
cbarpos = get(gca,'Position');
colorbar('Position', [cbarpos(3)+cbarpos(1)+0.01 cbarpos(2)  cbarpos(3)*0.01  cbarpos(4)*0.97])
title('2) Segmentation Distance Matrix')

%%% ----------------------------------------------------------------------------
% 3) Novelty Curve and Peak Picking
%    Extract novelty curve and peaks (segmentation boundaries), then plot.

tic
%...............................................................................
% Different methods to obtain a Novelty curve from the
% segmentation distance matrix (comment/uncomment):

%%%. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% Crude, dirty, fast and shameless methods:

%novelty_curve = mean(segmentation_distance_matrix); % average (fast)
%novelty_curve = sum(segmentation_distance_matrix); % sum (faster)

%%%. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% Convolve the main diagonal with a Gaussian-tapered checkerboard kernel:
% (c.f. Foote & Cooper, 2003)

% make a 2D Gaussian-tapered checkerboard kernel :
W = 2; % <------------------------------------- kernel width in hops (minimum 2) 
[kernel_x kernel_y] = meshgrid((-(W-1)/2):((W-1)/2), (-(W-1)/2):((W-1)/2));
checkerboard_gauss_kernel = (exp( -((2*pi*kernel_x/W).^2) / 2 - ((2*pi*kernel_y/W).^2) / 2)) ...
    .* (kron([-1, 1; 1,  -1],ones(W/2)));

% % METHOD 1: Sliding kernel along matrix' diagonal
% % add a margin to the original matrix: 
% margin = W/2; % <---- margin (W/2 will cause the kernel to start with is center at the beginning of the original matrix)
% margined_matrix = ...
%     zeros(amount_segmentation_windows + margin*2); % initialise with zeros a bigger matrix (A.K.A. margined matrix) to contain the original matrix
% margined_begin = margin + 1;
% margined_end = margin + amount_segmentation_windows;
% margined_matrix( margined_begin : margined_end , margined_begin : margined_end ) = ...
%     segmentation_distance_matrix; % put the original matrix into the bigger matrix
% 
% % Convolve it, baby!
% end_offset = W - 1;
% novelty_curve = zeros(1, (amount_segmentation_windows + 2* (margin-W) ));
% for i = 1:( amount_segmentation_windows + 2*margin - W )
%     sliding_window_beginning = i;
%     sliding_window_end = end_offset + i ;
%     sliding_window = margined_matrix( sliding_window_beginning : sliding_window_end , sliding_window_beginning : sliding_window_end );
%     novelty_curve(i) = sum(sum( sliding_window .* checkerboard_gauss_kernel ));
% end
% 
% novelty_curve = novelty_curve.^2; % boost

% METHOD 2: Diagonal of kernel-convoluted matrix
convolved_segmentation_matrix = conv2(segmentation_distance_matrix,checkerboard_gauss_kernel); % convolve kernel with segmentation matrix
novelty_curve = diag(convolved_segmentation_matrix); % extract the diagonal
novelty_curve = novelty_curve( fix(W/2) : end - fix(W/2) , 1 )'; % trim borders
novelty_curve = novelty_curve.^2; % boost!

%%%. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

% NOTE:
% The sum and average methods evaluate distinctiveness within the whole signal 
% whereas the kernel method acts locally.

%%%.............................................................................

timer_2 = toc;
disp(['Novelty curve computing time = ',num2str(timer_2),' sec.'])
disp('----------------------------------------------------------')

% add a "tail" to the curve so that scaled it matches the length of the original signal:
novelty_curve_display = zeros(1,amount_segmentation_windows+segmentation_tail);
novelty_curve_display(:,1:amount_segmentation_windows) = novelty_curve;
novelty_curve_display(:,amount_segmentation_windows+1:end) = repmat(novelty_curve(:,end),1,segmentation_tail);

% Retrieve peaks over a threshold:
threshold_factor = 0.1; % <---- peaks threshold factor
peaks_ind = zeros(amount_segmentation_windows,1)';
peaks_ind(2:end-1) = (diff(sign(diff(novelty_curve))) == -2);
peaks_values = novelty_curve.*peaks_ind;
max_peak = max(peaks_values);
selpeaks_bool = (peaks_values >= max_peak * threshold_factor);
selpeaks_curve = novelty_curve.*selpeaks_bool;
selpeaks_values = selpeaks_curve(selpeaks_curve > 0);
selpeaks_ind = find(selpeaks_bool);
selpeaks_ind_scaled = selpeaks_ind*total_time/amount_segmentation_windows;
selpeaks_ind_scaled_display = ...
    selpeaks_ind / (segmentation_tail/amount_segmentation_windows + 1) * total_time/amount_segmentation_windows;

% Plot:
%AM_FM_signal_rescaled = max_peak * (AM_FM_signal + max(AM_FM_signal)) / (2*max(AM_FM_signal));
AM_FM_signal_rescaled = AM_FM_signal + abs(min(AM_FM_signal));
AM_FM_signal_rescaled = (max_peak - min(novelty_curve)) * (AM_FM_signal_rescaled/max(AM_FM_signal_rescaled));
AM_FM_signal_rescaled = AM_FM_signal_rescaled + min(novelty_curve);
hold off
subplot(5,1,3)
plot(t,AM_FM_signal_rescaled,'Color',[0.7 1 1]*0.8) % plot signal rescaled
hold on

segmentation_timehop = segmentation_hop * sample_time;
t_novelty_curve =[segmentation_timehop:segmentation_timehop:total_time+segmentation_timehop];
plot(t_novelty_curve,novelty_curve_display,'b') % plot novelty curve

t_boundaries = repmat(selpeaks_ind_scaled_display,2,1);
y_boundaries = repmat((get(gca,'ylim')),size(find(selpeaks_bool),2),1)';
line(t_boundaries,y_boundaries,'Color',[0.9 0.6 0.6],'linewidth',2,'linestyle','-') % plot segmentation boundaries
title('3) Novelty Curve, Peak Picking and Segmentation of AMFM Signal')

plot(selpeaks_ind_scaled_display,selpeaks_values,'.r','markersize',25) % plot peaks

%%% ----------------------------------------------------------------------------
% 4) Clustering Similarity Matrix
%    Generate a matrix of distance between segments and display a heatmap of the
%    corresponding similarity matrix:

clustering_bound = 10; % <------------------------------ boundary for closest distance (percent)
boundaries_indexes = selpeaks_ind * round(amount_segmentation_windows/segmentation_hop);
amount_segments = size(boundaries_indexes,2)+1;
clustering_distance_matrix = zeros(amount_segments);
fractot = 200/(amount_segments*(amount_segments-1)); 
laps = 1;
vector_counter = 1;

disp(fprintf('Progress of clustering distance matrix = '));
tic

for row = 2:amount_segments
    
    if row == amount_segments
        a_1 = AM_FM_signal( boundaries_indexes(1,row-1)+1 : end );
    else
        a_1 = AM_FM_signal( boundaries_indexes(1,row-1)+1 : boundaries_indexes(1,row) );
    end
    
    for column = 1:(row-1)
        
        if column == 1
            b_1 = AM_FM_signal( 1 : boundaries_indexes(1,column) );
        else
            b_1 = AM_FM_signal( boundaries_indexes(1,column-1)+1 : boundaries_indexes(1,column) );
        end
        
        size_a_1 = size(a_1,2);
        size_b_1 = size(b_1,2);
        
        %if abs(size_a_1 - size_b_1) <= segmentation_hop % windows whose sizes are a segmentation_hop of difference or less (equal size)
            
            clustering_distance_matrix(row,column) = dtw(a_1',b_1',round(size_a_1*clustering_bound/100)); % dynamic time warping (slow and accurate)
            
            % clustering_distance_matrix(row,column) = LCSS; <--- TO-DO !!!
            
        %else % windows whose sizes are more than a segmentation_hop of difference
            %clustering_distance_matrix(row,column) = NaN;
        %end

        clustering_distance_matrix(column,row) = clustering_distance_matrix(row,column);

        display_progress = round(fractot*laps);
        fprintf(repmat('\b',1,numel(mat2str(display_progress))+3))
        fprintf(' %d%%',display_progress);
        fprintf('\n')
        laps = laps+1;
    end
end

timer_3 = toc;
disp(['Clustering distance matrix computing time = ',num2str(timer_3),' sec.'])
disp('----------------------------------------------------------')

clustering_similarity_matrix = 1 - (clustering_distance_matrix / max(max(clustering_distance_matrix)));

subplot(5,1,4)
imagesc(clustering_similarity_matrix)
set(gca,'xlim',[0.5,amount_segments+0.5],'xtick',[1:amount_segments],'Ticklength', [0 0])
cbarpos = get(gca,'Position');
colorbar('Position', [cbarpos(3)+cbarpos(1)+0.01 cbarpos(2)  cbarpos(3)*0.01  cbarpos(4)*0.97])
title('4) Clustering Similarity Matrix')

%%% ----------------------------------------------------------------------------
% 5) Pattern Clustering
%    Find clusters of patterns and display labels of clustered segments 
%    on top of signal.

clustering_similarity_matrix_nonan = clustering_similarity_matrix;
clustering_similarity_matrix_nonan(find(isnan(clustering_similarity_matrix_nonan))) = 0;
Z = linkage(clustering_similarity_matrix_nonan,'ward'); % 

distance_threshold = 2; % <--- threshold for clustering (distance between nodes)
labelled_segments = cluster(Z,'cutoff',distance_threshold,'criterion','distance'); % clustering according to distance

%inconsistency_threshold = 0.4; % <--- threshold for clustering (inconsistency, between 0 and 1)
%labelled_segments = cluster(Z,'cutoff',inconsistency_threshold,'depth',3); % clustering according to inconsistency

hold off
subplot(5,1,5)
plot(t,AM_FM_signal,'Color',[0.7 1 1]*0.8)
hold on

y_limits = (get(gca,'ylim'));
lower_text_position = y_limits(1)+y_limits(1)*1/4;
% t_boundaries = repmat((selpeaks_ind*total_time/amount_segmentation_windows),2,1);
y_boundaries = repmat(-y_limits,size(find(selpeaks_bool),2),1)';
line(t_boundaries,y_boundaries,'Color',[1 0.6 0.6],'linewidth',2,'linestyle','-')
alphabet=('A':'Z');
for i_1 = 1:amount_segments  % display the clusters labelling them with letters
    if i_1 == 1
        text(t_boundaries(1,i_1)/2,0,alphabet(labelled_segments(i_1,1)),...
            'fontsize',20,'Color','k','HorizontalAlignment','Center');
        text(t_boundaries(1,i_1)/2,lower_text_position,num2str(i_1),...
            'fontsize',12,'Color','k','HorizontalAlignment','Center');        
    elseif i_1 == amount_segments
        text(t_boundaries(1,i_1-1)+(t(end)-t_boundaries(1,i_1-1))/2,0,alphabet(labelled_segments(i_1,1)),...
            'fontsize',20,'Color','k','HorizontalAlignment','Center');
        text(t_boundaries(1,i_1-1)+(t(end)-t_boundaries(1,i_1-1))/2,lower_text_position,num2str(i_1),...
            'fontsize',12,'Color','k','HorizontalAlignment','Center');
    else
        text(t_boundaries(1,i_1-1)+(t_boundaries(1,i_1)-t_boundaries(1,i_1-1))/2,0,alphabet(labelled_segments(i_1,1)),...
            'fontsize',20,'Color','k','HorizontalAlignment','Center');
        text(t_boundaries(1,i_1-1)+(t_boundaries(1,i_1)-t_boundaries(1,i_1-1))/2,lower_text_position,num2str(i_1),...
            'fontsize',12,'Color','k','HorizontalAlignment','Center');        
    end
end
set(gca,'xtick',[])
title('5) Pattern Clustering')

%%%.............................................................................
% Display dendrogram:

figure
dendrogram(Z,'colorThreshold','default')
title('Hierarchical Pattern Clustering')
