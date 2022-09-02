%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                        VIDEO PROGRESS BAR AND MARKERS                        %
%                                                                              %
%                                  January 2020                                %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ==============================================================================

% This program has been tested with:
%   Matlab R2015a

% Description:
%   This program produces a video with a progress bar or markers or both.
% The video may be added to the target video using a video editor, 
% by replacement of green or blue pixels.

% Instructions:
%   Edit the values indicated with an arrow like this: <--- (length of the arrow can vary)
% Run the script, close your eyes and hope for the best.

% ==============================================================================
% Initialisation:

clc
clear
close all
restoredefaultpath

% ==============================================================================
% set parameters:

% ..............................................................................

  video_filename = ''; % <--- video file name (empty to use markers filename)
   videos_folder = ''; % <--- path to save videos
markers_filename = ''; % <--- markers file name (empty for no markers or to use all .txt files in markers folder)
  markers_folder = ''; % <--- path for folder with only markers .txt files (empty for no markers)
  
% ..............................................................................

     nframes_out = 1000;            % <--- video number of frames (total time), will be overriden if markers file is loaded
         fps_out = 29.97;           % <--- video frame rate
       frame_dim = [960, 540];      % <--- video width and height (pixels)
         pbar_sw = 1;               % <--- make progress bar (1 = yes ; 0 = no and make only markers)
   pbarmk_height = 20;              % <--- progress bar and markers height (pixels)
        mk_width = 3;               % <--- markers width (pixels)
     pbar_colour = [0.6, 0.3, 0.9]; % <--- progress bar colour (RGB)
       mk_colour = [1, 1, 0];       % <--- markers colour (RGB)
    barbg_colour = [0, 0, 0];       % <--- bar and markers background colour (RGB)
     mk_fontsize = 20;              % <--- marker number font size (0 = no marker numbers)
        r_colour = [1, 1, 1];       % <--- number's enclosing rectangle background colour (RGB)
   canvas_colour = 'green';         % <--- frame background colour ('green' or 'blue', other value will set black)
   
        pv_mk_sw = 0;   % <--- 1 = preview only markers; 0 = preview or save video
   save_video_sw = 1;   % <--- 1 = save video; 0 = preview video
   pv_video_frac = 1;   % <--- fraction of video to preview
        
% The markers text file should be un-formatted .txt, whose first line 
% is the frame rate, the second line is the total number of frames, and the
% following lines are one marker per line, in cronological order, as frames.

% ------------------------------------------------------------------------------
% make MPEG-4 video:

if ~isempty(markers_folder) && ~exist(markers_folder, 'dir')
    error('Specified markers folder doesn''t exist.')
end
 
if isempty(markers_filename)    
    mk_dir = dir([markers_folder,'/*.txt']);
    n_vid = length(mk_dir);
    if n_vid == 0
        error('No markers files in specified markers folder.')
    end
else
    mk_dir = [];
    n_vid = 1;
end

if ~pbar_sw && isempty(markers_folder) 
    error('No markers folder specified and progress bar switch off. Nothing to do.')
end

fprintf('computing...\n');
pbarmk_start = frame_dim(2) - pbarmk_height + 1;
figpos = [50, 50, frame_dim(1), frame_dim(2)];

for i_vid = 1:n_vid
    
    if pv_mk_sw || (~pv_mk_sw && i_vid == 1)
        figure('WindowStyle','normal','Position',figpos);
    end
    
    frame_mov = zeros(frame_dim(2),frame_dim(1),3);
    frame_pbar = zeros(pbarmk_height,frame_dim(1),3);
    if ~isequal(barbg_colour,[0,0,0])
        for i_c = 1:3
            frame_pbar(:,:,i_c) = barbg_colour(i_c);
        end
    end
    if strcmp(canvas_colour,'green')
        frame_mov(1:pbarmk_start-1,:,2) = 1;
    elseif strcmp(canvas_colour,'blue')
        frame_mov(1:pbarmk_start-1,:,3) = 1;
    end
    
    if ~isempty(mk_dir)
        markers_filename = mk_dir(i_vid).name;
    end
    
    if ~isempty(markers_filename) % make markers
        
        [~,mk_name,mk_ext] = fileparts(markers_filename);
        if isempty(mk_ext)
            markers_filename = [mk_name,'.txt'];
        end
        markers_data = importdata([markers_folder,'/',markers_filename]);
        frame_mk = frame_pbar;
        mk_mask = frame_pbar;
        nframes_out = round( (markers_data(2) * fps_out) / markers_data(1) );
        mkpix = round( ( markers_data(3:end) * frame_dim(1) ) / markers_data(2) ); % markers frames to pixels
        mkpix_margin = ceil([NaN; diff(mkpix)] / 2); % half space between markers in pixels
        if max(mkpix) > frame_dim(1)
            error('Marker in file ''%s'' is greater than specified data length. Process stopped.',markers_filename)
        end
        mkpix_start = mkpix - ceil(mk_width/2); % offset for centering
        mkpix_end = mkpix_start + mk_width - 1; % add width
        mkpix_start(mkpix_start < 1) = 1; % trim beginning
        mkpix_start(mkpix_start > (frame_dim(1) - ceil(mk_width/2) - 1) ) = (frame_dim(1) - ceil(mk_width/2) - 1); % adjust ending
        mkpix_end(mkpix_end > frame_dim(1)) = frame_dim(1); % trim ending
        n_mk = length(mkpix);
        for i_m = 1:n_mk
            for i_c = 1:3
                frame_mk(:,mkpix_start(i_m):mkpix_end(i_m),i_c) = mk_colour(i_c);
                mk_mask(:,mkpix_start(i_m):mkpix_end(i_m),i_c) = 1;
            end
        end
        mk_mask = logical(mk_mask);
        frame_mov(pbarmk_start:end,:,:) = frame_mk;
        frame_pbar = frame_mk;
    end
    
    clf;
    axes('position', [0 0 1 1]);

    if ~pv_mk_sw

        start_frames = ones(1,nframes_out+1);
        start_frames(2:end) = ceil((1:nframes_out) * (frame_dim(1) / nframes_out) );
        start_frames(end) = frame_dim(1);
        
        if save_video_sw
            if ~exist(videos_folder, 'dir')
                mkdir(videos_folder)
            end
            if isempty(video_filename)
                if exist('mk_name', 'var')
                    video_filename = mk_name;
                else
                    fprintf(repmat('\b', 1, 13));
                    error('No video file-name or markers specified. Need to specify either.')
                end
            end
            mov_obj = VideoWriter([videos_folder,'/',video_filename],'MPEG-4');
            mov_obj.FrameRate = fps_out;
            open(mov_obj);
            video_filename = [];
        else
            nframes_out = nframes_out * pv_video_frac;
        end
    else
        nframes_out = 0;
    end

    for i_f = 0:nframes_out
        
        if i_f > 0 && pbar_sw
            i_bp = start_frames(i_f) : (start_frames(i_f + 1)-1);
            
            for i_c = 1:3
                frame_pbar(:,i_bp,i_c) = pbar_colour(i_c); % progress bar fill with colour
            end
            if ~isempty(markers_filename) % superimpose markers
                frame_pbar(mk_mask) = frame_mk(mk_mask);
            end
        end
        frame_mov(pbarmk_start:end,:,:) = frame_pbar;

        image(frame_mov)
        axis off
        rpos = zeros(1,4);
        if mk_fontsize > 0 % superimpose numbers
            num_base_vpos = round(pbarmk_start - mk_fontsize / 1.7);
            
            rec_vpos = num_base_vpos - mk_fontsize / 1.6; % rectangle vertical position
            rpos(4) = mk_fontsize * 1.1;                  %     "     height

            for i_m = 1:n_mk 
                
                rpos(1) = mkpix(i_m) - (mk_fontsize / 1.7) * (length(num2str(i_m))^0.25); % rectangle horizontal position
                rpos(3) = rpos(4) * length(num2str(i_m)) / (length(num2str(i_m))^0.7); % rectangle width
                
                if ( mkpix_margin(i_m) <= ( rpos(3) / 2 ) ) && at_base_pos
                    rpos(2) = rec_vpos - rpos(4);
                    num_vpos = num_base_vpos - rpos(4);
                    at_base_pos = 0;
                else
                    rpos(2) = rec_vpos;
                    num_vpos = num_base_vpos;
                    at_base_pos = 1;
                end
                
                rectangle('Position',rpos,'FaceColor',r_colour)
                text(mkpix(i_m)-1,num_vpos,sprintf('%d',i_m),'FontWeight','Bold','Fontsize',mk_fontsize,'HorizontalAlignment','center')
            end
        end
        drawnow
        
        if save_video_sw && ~pv_mk_sw
            writeVideo(mov_obj,getframe(gcf));
        end
    end
    
    if save_video_sw && ~pv_mk_sw
        close(mov_obj);
    end
end

fprintf(repmat('\b', 1, 13));
fprintf('done\n');
beep on
for i_beep = 1:5
    beep
    pause(1)
end