%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                     MAKE FIGURE FROM EXISTING .FIG FILES                     %
%                        INTO VERTICALLY STACKED SUBPLOTS                      %
%                         AND SAVE IT AS A .TIFF PICTURE                       %  
%                                                                              %
%                                 November 2018                                %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program has been tested with:
%   Matlab R2015a

% ==============================================================================
% Description:

% This program takes as input several Matlab-native figures saved in files with 
% extension '.fig', stacking them into a single figure. The resulting figure
% is then saved as a TIFF file with extension '.tif'.

% ==============================================================================
% Instructions:

% Edit the parameters marked with an arrow like this: <---
% Run the program and enjoy.

% ==============================================================================

     fig_path = ' ';   % <--------- folder where the FIG files are and where the resulting figure will be saved
filename_tiff = 'Fig1';      % <--- name of the resulting TIFF file without extension ([] = don't save)
      figdims = [4,3] * 300; % <--- width and height of the picture
  sel_monitor = 1;     % <--------- select a monitor to display the picture
 panel_lbl_sw = 1;     % <--------- display an alphabetic label for each subplot (0 = no, 1 = yes)
    font_mult = 0.75;  % <--------- ratio to rescale fonts
     v_offset = 0.09;  % <--------- ratio to offset vertically
      v_space = 0.022; % <--------- ratio to separate vertically
    v_rescale = 0.7;   % <--------- ratio to rescale vertically
 
fi = 0;
fi = fi + 1;
filenames_fig{fi} = ' '; % <--- name of a FIG file without extension

fi = fi + 1;
filenames_fig{fi} = ' '; % <--- name of a FIG file without extension

fi = fi + 1;
filenames_fig{fi} = ' '; % <--- name of a FIG file without extension

fi = fi + 1;
filenames_fig{fi} = ' '; % <--- name of a FIG file without extension

fi = fi + 1;
filenames_fig{fi} = ' '; % <--- name of a FIG file without extension

fi = fi + 1;
filenames_fig{fi} = ' '; % <--- name of a FIG file without extension

%-------------------------------------------------------------------------------
clc
monpos = get(0,'MonitorPositions');
this_monpos = monpos(sel_monitor,:);
if  sum(figdims > this_monpos(3:4))
    error('The width and height of the picture should not be greater than [%g, %g]',this_monpos(3),this_monpos(4));
end
close all
cd(fig_path)
figs_n = length(filenames_fig);
for i = 1:figs_n
    figs{i} = openfig([filenames_fig{i},'.fig'],'reuse','invisible');
    gcas{i} = gca;
    gcas{i}.FontSize = gcas{i}.FontSize * font_mult;
    gcas{i}.Title.FontSize = gcas{i}.Title.FontSize * font_mult;
    gcas{i}.XLabel.FontSize = gcas{i}.XLabel.FontSize * font_mult;
    gcas{i}.YLabel.FontSize = gcas{i}.YLabel.FontSize * font_mult;
end
figpos = (this_monpos(3:4) - figdims)/2;
figpos = [figpos,figdims];
figsh = figure('Position',figpos);
alphabet = 'abcdefghijklmnopqrstuvwxyz';
this_v_space = 0;
for i = 1:figs_n
    sp_ax{i} = subplot(figs_n,1,i);
    this_v_space = this_v_space + v_space;
    sp_ax{i}.Position(2) = sp_ax{i}.Position(2) + v_offset - this_v_space;
    sp_ax{i}.Position(4) = sp_ax{i}.Position(4) * v_rescale;
    handles{i} = get(gcas{i},'children');
    copyobj(handles{i},sp_ax{i});
    field_names{i} = setdiff(fieldnames(gcas{i}),...
        {'BeingDeleted','Children','CurrentAxes','CurrentCharacter','CurrentObject','CurrentPoint','OuterPosition','Parent',...
        'Position','TightInset','Type'});
    for i_1 = 1:length(field_names{i})
        set(gca, field_names{i}{i_1}, gcas{i}.(field_names{i}{i_1}));
    end
    if panel_lbl_sw
        annotations{i} = annotation('textbox','String',['(',alphabet(i),')'],'EdgeColor','none');
        annotations{i}.FontSize = gcas{i}.Title.FontSize;
        annotations{i}.Position = get(gca,'Position');
        annotations{i}.Position(1) = annotations{i}.Position(1) - annotations{i}.Position(1)/1.3;
        annotations{i}.Position(4) = annotations{i}.Position(4)/1.7;
    end
end
if isempty(filename_tiff) == 0
    formattype = '-dtiff';
    resolution = '-r300';
    figsh.PaperPositionMode = 'auto';
    print(figsh,[filename_tiff,'.tif'],formattype,resolution)
end