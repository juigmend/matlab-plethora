%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                              RASTRIGIN'S FUNCTION                            %
%                                      DEMO                                    %
%                                                                              %
%                                8 February 2018                               %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script has been tested with Matlab R2015a

% ==============================================================================
close all
clear
clc

% Rastrigin's function:
rastrigin = @(XY) (10 * size(XY,2) + sum( (XY .^2 - 10 * cos(2 * pi .* XY)) , 2 ));

range = [-1, 1]; % <--- range to evaluate the function
resol = 100;     % <--- resolution 

figure
colormap copper

%...............................................................................
% make  Rastrigin's function in defined range:
xy_seed = range(1) : ((range(end)-range(1))/resol) : range(end) ;
L = length(xy_seed);
Z = zeros(length(xy_seed));
for rowind = 1:L
    for colind = 1:L
        Z(rowind,colind) = feval(rastrigin,[xy_seed(rowind),xy_seed(colind)]);
    end
end

% plot function's surface:
subplot(2,2,1)
surf(xy_seed,xy_seed,Z)
xlim(range)
ylim(range)
shading interp

% plot function's heat map:
subplot(2,2,2)
imagesc(Z)

%...............................................................................
% Make and plot Rastrigin's function in defined range, 
% given individual random points:

iters = 500; % <--- amount of iterations

itersstr = num2str(iters);
results = zeros(iters,3); % rows are iterations and colums are [X,Y,Z]
colourstep = length(colormap) / max(max(Z));
plotscat = @(x,y,z) (scatter3(x,y,z,20,floor(z*colourstep)+1,'filled'));
maxfunc = max(max(Z));
for i = 1:iters
    results(i,1:2) = (rand(1,2) * range(end) * 2) + range(1);
    results(i,3)  = feval(rastrigin,results(i,1:2));
    
    % plot points from anlge of surface:
    subplot(2,2,3)
    feval(plotscat,results(i,1),results(i,2),results(i,3))
    title(['iteration ',num2str(i),' of ',itersstr])
    xlim(range)
    ylim(range)
    zlim([0, maxfunc])
    hold on
    drawnow
    
    % plot points from anlge of heat map:
    subplot(2,2,4)
    feval(plotscat,results(i,1),results(i,2),results(i,3))
    xlim(range)
    ylim(range)
    zlim([0, maxfunc])
    view(0,0)
    title(['iteration ',num2str(i),' of ',itersstr])
    view(0,90)
    hold on
    drawnow
end



