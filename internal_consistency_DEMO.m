%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                    Consistency of participants'responses                     %
%                                                                              %
%                                   March 2020                                 %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ==============================================================================
% INFORMATION

% This program has been tested with:
%   Matlab R2015a

% Dsecription:
%   Generate artificial data, compute regression models (Partial Least Squares
%   and Stepwise), predict responses and visualise.

% Instructions:
%   Edit the values indicated with an arrow like this: <--- (length of the arrow can vary)
%   Run the script, close your eyes and hope for the best.

% ==============================================================================

% Initialisation:
clc
close all

% ------------------------------------------------------------------------------

         n_q = 16;  % <--- number of questions
         n_p = 30;  % <--- number of responders
length_scale = 127; % <--- length of the measuring scale for each question
         n_l = 130; % <--- disagreement level

      n_iter = 1;   % <--- test iterations to measure fastest algorithm

% ------------------------------------------------------------------------------

fastest_record = zeros(n_iter,4);
alpha = zeros(n_iter,3);
figsize = get(0, 'ScreenSize');
screen_w = figsize(3);
figsize(3:4) = figsize(3:4) * .8;
figsize(1) = screen_w/2 - figsize(3)/2;

for i_iter = 1:n_iter
    
    X = repmat( (rand(n_tr,1) * 100) ,1,n_p); % alpha = 1
    X = X + rand(n_tr,n_p) * n_l; % noise
    
    X = floor( X / max(X(:)) * length_scale );
    
    if n_iter == 1
        figure('Position', figsize)
        
        subplot(2,2,1)
        plot(X)
        xlim([1,n_q])
        axis square
        if n_p < 5
            set(gca,'XTick',[1:n_p])
        end
        xlabel('questions')
        title('raw data')
        
        
        subplot(2,2,2)
        imagesc(X)
        axis square
        if n_p < 5
            set(gca,'XTick',[1:n_p])
        end
        ylabel('questions')
        xlabel('responders')
        title('raw data')
        
        corr_X = corr(X);

        subplot(2,2,3)
        imagesc(corr_X)
        colorbar
        axis square
        title('correlation amongst responders')
        
        corr_X_triu =  triu(corr_X,1);
        mean_corr_X = mean( corr_X_triu( corr_X_triu ~= 0 ) );
        fprintf(' mean corr = %.4g \n\n',mean_corr_X)
        
        dist_X = squareform( pdist(X') );
                
        subplot(2,2,4)
        imagesc(dist_X)
        colorbar
        title('Euclidean distance amongst responders')
        axis square
        
    end
    
    % Cronbach's alpha, algorithm 1:
    tic
    varsum = var( sum( X, 2 )); % variance of all participants
    sumvar = sum( var( X )); % variance of each participant
    alpha(i_iter,1) = ( n_p / (n_p-1) ) * ( (varsum - sumvar) / varsum );
    times(1) = toc;
    
    % Cronbach's alpha, algorithm 2:
    tic
    alpha(i_iter,2) = ( n_p / (n_p - 1) ) * ( 1 - sum(var(X,1)) / var(sum(X,2),1) ) ; % biased variance (w = 1)
    times(2) = toc;
    
    % Cronbach's alpha, algorithm 3:
    tic
    i_triu = triu( true( n_p ),1);
    vc	= cov( X );
    v	= mean( diag( vc ));
    c	= mean( vc( i_triu ));
    alpha(i_iter,3)	= n_p * c / (v + (n_p - 1 ) * c);
    times(3) = toc;
    
    % Standardized Cronbach's alpha:
    tic
    i_triu_2 = triu( true( n_p ),1);
    corr_X_2 = corr(X);
    mean_corr_X_2 = mean( corr_X_2( i_triu_2 ) );
    alpha_std(i_iter) = (n_p * mean_corr_X_2) / (1+(n_p-1)*mean_corr_X_2);
    times(4) = toc;
    
    [~,i_fastest] = min(times);
    fastest_record(i_iter,i_fastest) = 1;
    
    
    
end

disp(' alpha 1  alpha 2  alpha 3 std alpha')
disp(' -----------------------------------')
for i_iter = 1:n_iter
    
    disp(sprintf(' %.4f  ',[alpha(i_iter,:),alpha_std(i_iter)]))
end

if n_iter > 1    
    figure
    imagesc(fastest_record)
    set(gca,'XTick',[1:4])
end

