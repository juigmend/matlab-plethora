%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                             AUTOMATIC SEGMENTATION,                          %
%                  ROTATION AND PRODUCTION OF A MOTION PICTURE                 %  
%                  OF MULTI-MARKER OPTICAL MOTION CAPTURE DATA                 %
%                                                                              %
%                                 October 2018                                 %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program has been tested with:
%   Matlab R2015a
%   Mocap Toolbox v1.5 https://www.jyu.fi/hytk/fi/laitokset/mutku/en/research/materials/mocaptoolbox
 
% ==============================================================================
% Description:

% This program first segments optical motion capture (mocap) data recorded with 
% a Qualisys infrared reflective marker system. Then, a segment is extracted and 
% rotated. Finally, a video file of the rotated segment is produced.

% The input data has to be producced by the Qualisys QTM software. 
% The required files are:
%   - a .TSV file containing 3D marker trajectories data
%   - a .TXT file containing the labels list

% ==============================================================================
% Instructions:

% If necessary, edit the parameters marked with an arrow like this: <---
% Run the cells one by one, following the comments.

% ==============================================================================
% Initialisation:

clc
clear
close all
restoredefaultpath

% ==============================================================================
% Declare toolboxes paths and data files:

% toolboxes_path = ' '; % <------------ toolboxes path
%      data_path = ' '; % <--- folder where mocap files are (.TSV and .TXT)

    mocap_file = 'heaven_six_short.tsv'; % <--- mocap data file
   labels_file = 'heaven_six_short.txt'; % <--- mocap labels file
video_filename = 'heaven_six_short';     % <--- name that the video file will have
 
% ..............................................................................

addpath(genpath(toolboxes_path))
cd(data_path)

%% -----------------------------------------------------------------------------
% IMPORT AND PREPROCESS DATA

rsfreq = 10; % <--- frequency to resample (Hz)

% Import and resample mocap data:
mocap_struct_raw = mcread(mocap_file);
mocap_struct_resampled = mcresample(mocap_struct_raw,rsfreq); % resample

%% -----------------------------------------------------------------------------
% SEGMENTATION

% Similarity matrix: 
mocap_simmat = mcsimmat(mocap_struct_resampled,'euclidean');

%% .............................................................................
% Novelty score:

bw = 10;   % <--- novelty kernel bandwidth  

window_width = round( bw * 8);
if rem(window_width,2) % this ensures that the length of the window is even so that the checkerboard can be produced
    window_width = window_width + 1;
end
[kernel_x, kernel_y] = meshgrid((-(window_width-1)/2):((window_width-1)/2), (-(window_width-1)/2):((window_width-1)/2));
variance = bw^2;
height =  1 / ( ( 2 * pi * variance) ) ; % this produces unit volume
kernel = height *  exp( - ( (kernel_x.^2) / ( 2 * variance) ) - ( (kernel_y.^2) / ( 2 * variance) ) );
kernel = kernel .* (kron([-1, 1; 1,  -1],ones(window_width/2)));

amt_w = length(mocap_simmat);
margin = window_width/2;
beginning = margin + 1;
ending = margin + amt_w;
wis = window_width - 1;
ext = zeros(amt_w + window_width);

% pad beginning and ending corners with averages:
ext(1:window_width,1:window_width) = mean(mean(mocap_simmat(1:window_width/2,1:window_width/2)));
ext(end-window_width:end,end-window_width:end) = mean(mean(mocap_simmat(end-window_width/2:end,end-window_width/2:end)));
ext( beginning : ending , beginning : ending ) = mocap_simmat;
nov = zeros(1, amt_w);

% convolve kernel along diagonal of similarity matrix:
for i_3 = 1:amt_w
    wbeg = i_3;
    wend = wis + i_3 ;
    nov(i_3) = sum(sum(ext( wbeg : wend , wbeg : wend ) .* kernel ));
end

% Normalise it, don't criticise it:
minnov = min(nov);
nov = (nov - minnov) / ( max(nov) - minnov);

%% .............................................................................
% Local maxima:

thr = 0.5; % <--- novelty peaks threshold

peaks_index_all = zeros(1,length(nov));
peaks_index_all(2:end-1) = (diff(sign(diff(nov))) == -2);
peaks_index_all = nov .* peaks_index_all;
peaks_index_selected = (peaks_index_all >= max(peaks_index_all) * thr);
pind = find( peaks_index_selected == 1 );

%% .............................................................................
% Visualise segmentation:

close all
datalength = length(nov);
boundaries_colour = [0.2 0.7 0.2];

% plot similarity matrix:
subplot(2,1,1)
imagesc(mocap_simmat)
set(gca,...
    'yticklabel',[],...
    'xtick',0:10*rsfreq:datalength,...
    'xticklabel',0 : 10 : datalength/rsfreq )
xlabel('time (seconds)','fontsize',12)
title('SELF-SIMILARITY MATRIX (from mocap data)','fontsize',14)

% plot novelty curve:
subplot(2,1,2)
plot(nov)
set(gca,...
    'ylim',[ min(nov) , max(nov) ],'xlim',[1,datalength],...
    'xtick',0:10*rsfreq:datalength,...
    'xticklabel',0 : 10 : datalength/rsfreq )
xlabel('time (seconds)','fontsize',12)
hold on

% plot computed boundaries:
rep_peaks = repmat(pind,2,1);
lines = repmat(get(gca,'ylim'),size(rep_peaks,2),1)';
line(rep_peaks,lines,'Color',boundaries_colour,'linewidth',3,'linestyle','-');
title('COMPUTED BOUNDARIES (novelty peaks from mocap data)','fontsize',14)

% display segment labels:
pbind = [0,pind,datalength];
for i = 1:(length(pind)+1)
    thisx = pbind(i) + (pbind(i+1) - pbind(i))/2;
    text(thisx,0.8,num2str(i),'fontsize',20,'Color',boundaries_colour,'HorizontalAlignment','Center');
end

%% .............................................................................

selected_segment = 3; % <--- select segment to use for the rotating moving picture
         preroll = 2; % <--- pre-roll time (seconds)
        postroll = 2; % <--- post-roll time (seconds)

% Extract segment:
mocap_struct_segment = mctrim(mocap_struct_raw,(pind(selected_segment-1)/rsfreq)-preroll,(pind(selected_segment)/rsfreq)+postroll);

% %% -----------------------------------------------------------------------------
% MOTION PICTURE

start_end_pitchdown = [5, 10];      % <--- start and end of pitch downwards (seconds)
  start_end_pitchup = [13.1, 18.1]; % <--- start and end of pitch upwards (seconds)
start_end_rollright = [9, 14.1];    % <--- start and end of yaw rightwards (seconds)
 
            scrsize = [1280 720];   % <---- screen size
 left_baston_colour = [0.9,0.8,0];  % <---- left baston colour
right_baston_colour = [0,0.8,0.3];  % <---- right baston colour
       trace_length = 0.2; % <------------- trace length (seconds)

% %% .............................................................................

% set parameters for picture:
animpar = mcinitanimpar;
animpar = mccreateconnmatrix(labels_file, animpar);
animpar.scrsize = scrsize;
animpar.msize = ones(27,1) * 10;
animpar.msize(24:27) = 10;
animpar.markercolors = ones(27,3);
animpar.markercolors(24:25,:) = repmat(left_baston_colour,2,1);
animpar.markercolors(26:27,:) = repmat(right_baston_colour,2,1);
animpar.conncolors = ones(length(animpar.conn),3);
animpar.conncolors = animpar.conncolors * 0.7;
animpar.conncolors(34,:) = left_baston_colour;
animpar.conncolors(35,:) = right_baston_colour;
animpar.tracecolors = ones(2,3);
animpar.tracecolors(1,:) = left_baston_colour;
animpar.tracecolors(2,:) = right_baston_colour;
animpar.cwidth = ones(35,1) * 4;
animpar.cwidth(34:35) = 8;
animpar.twidth = 1;
animpar.trm = [25,27];
animpar.trl = trace_length;
animpar.animate = 0;
animpar.fps = 60;
animpar.output = video_filename;
animpar.videoformat = 'mpeg4'; 
animpar.perspective = 1;
animpar.pers.c = [0 -10000 0];
animpar.pers.th = [0 0 0];
animpar.pers.e = [0 -10 0];

% rotate to frontal position:
mocap_struct_segment = mcrotate(mocap_struct_segment,180,[0 0 1]);

% check appearence:
close all
mcplotframe(mocap_struct_segment,1,animpar)

%% .............................................................................

% rotate pitch down-up with a sigmoidal slope:
pitchdown_length = (start_end_pitchdown(2) - start_end_pitchdown(1)) * mocap_struct_segment.freq;
pitchdown_slope = [-2:4/(pitchdown_length-1):2]';
pitchdown_slope = erf(pitchdown_slope); % "error function" (sigmoidal shape)
pitchdown_slope = pitchdown_slope - min(pitchdown_slope);
pitchdown_slope = pitchdown_slope / max(pitchdown_slope);
pitchdown_angles_sloped = pitchdown_slope * 90;
pitchup_length = (start_end_pitchup(2) - start_end_pitchup(1)) * mocap_struct_segment.freq;
pitchup_slope = [-2:4/(pitchdown_length-1):2]';
pitchup_slope = erf(pitchup_slope); % "error function" (sigmoidal shape)
pitchup_slope = pitchup_slope - min(pitchup_slope);
pitchup_slope = pitchup_slope / max(pitchup_slope);
pitchup_slope = flipud(pitchup_slope);
pitchup_angles_sloped = pitchup_slope * 90;
thetas_pitchdownup = zeros(mocap_struct_segment.nFrames,1);
thetas_pitchdownup(start_end_pitchdown(1)*mocap_struct_segment.freq : start_end_pitchdown(1)*mocap_struct_segment.freq+pitchdown_length-1) = pitchdown_angles_sloped;
thetas_pitchdownup(start_end_pitchdown(2)*mocap_struct_segment.freq : start_end_pitchup(1)*mocap_struct_segment.freq + 1) = 90;
thetas_pitchdownup(start_end_pitchup(1)*mocap_struct_segment.freq + 1 : round(start_end_pitchup(1)*mocap_struct_segment.freq+pitchup_length)) = pitchup_angles_sloped;
mocap_struct_pitched = mcrotate(mocap_struct_segment,thetas_pitchdownup,[1 0 0]);

% rotate rolling rightwards:
rollright_interval = 360/((start_end_rollright(2)-start_end_rollright(1))*mocap_struct_segment.freq);
rollright_angles = [rollright_interval:rollright_interval:360]';
thetas_rollright = zeros(mocap_struct_segment.nFrames,1);
thetas_rollright(start_end_rollright(1)*mocap_struct_segment.freq : start_end_rollright(1)*mocap_struct_segment.freq+length(rollright_angles)-1) = rollright_angles;
mocap_struct_rolled = mcrotate(mocap_struct_pitched,thetas_rollright,[0 1 0]);

%% .............................................................................

% make video:
mocap_struct_animate = mocap_struct_rolled;
animpar.animate = 1;
mcanimate(mocap_struct_animate,animpar);
