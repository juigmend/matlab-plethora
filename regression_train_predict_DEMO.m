%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                   Train and Predict using Regression Models                  %
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
%   Generate artificial data, compute regression models, predict responses and 
%   visualise.
%   Available regression methods:
%       Partial Least-Squares using 'plsregress' function (PLS) 
%       Multiple Linear Regression using'regress' function (MLRa)
%       Multiple Linear Regression using'fitlm' function (MLRb)
%       Stepwise Linear Regression using 'stepwisefit' function (SLR)

% Instructions:
%   Edit the values indicated with an arrow like this: <--- (length of the arrow can vary)
%   Run the script, close your eyes and hope for the best.

% ==============================================================================

% Initialisation:
clc
close all

% ------------------------------------------------------------------------------

sel_method = 2;      % <--- regression method(s): 1 = PLS ; 2 = MLRa ; 3 = MLRb ; 3 = SLR

    n_pred = 6;      % <--- number of predictors
     i_nnp = [1,7];  % <--- non-normal predictors (irrelevant for few observations)
       icp = [];     % <--- intercorrelated predictors
       pcr = 4;      % <--- predictors correlated to response
norm_maxit = 10;     % <--- check predictors' normal distribution: 0 = no, 1 = test, N = test and apply power transform up to N times
     n_obs = 16;     % <--- number of observations
    n_meas = 30;     % <--- number of responders
   i_obsrz = [2];    % <--- observations whose responses are mostly zero (some random elements added)
    g_resp = 0;      % <--- 1 = all responses, 0 = average response
  train_pc = 100;    % <--- percentage for training (the remainder will be used as new, unseen data)
   cv_fold = 5;      % <--- cross-validation fold (only PLS)

% ------------------------------------------------------------------------------
   
% X columns are independent variables (predictors):
X = zeros(n_obs,n_pred);
for i_pred = 1:n_pred
    if any(i_pred == i_nnp)
        X(:,i_pred) = rand(n_obs,1);
    else
        X(:,i_pred) = randn(n_obs,1);
        X(:,i_pred) = X(:,i_pred) - min(X(:,i_pred)) + rand(1)*0.1;
    end
end

% correlate predictors:
for i_icp = 2:length(icp)
    X(:,icp(i_icp)) = X(:,icp(1)) * (1 + rand(1,1)/2);
end

if norm_maxit
    
    % test predictors for normality (alpha > 0.05 ):
    % very few observations may render test results irrelevant
    
    fprintf('NORMALITY TEST:\n')
    
    X_notnormal = zeros(1,size(X,2));
    passed_normality_test = 0;
    i_iter = 1;
    boxcox_started = 0;
    ui = [];
    rep_newdata = 0;
    
    while ~passed_normality_test && (i_iter <= norm_maxit)
        
        start_boxcox = 0;
        for i_pred = 1:n_pred
            
            [X_notnormal(i_pred)] = lillietest(X(:,i_pred));
%             [X_notnormal(i_pred)] = kstest(X(:,i_pred));
%             [X_notnormal(i_pred)] = jbtest(X(:,i_pred));
            
            if X_notnormal(i_pred)
                
                if i_iter == 1
                    fprintf('X(%i) might be not normal :O\n',i_pred)
                    start_boxcox = 1;
                end
            end
        end

        passed_normality_test = ~any(X_notnormal);
        i_iter = i_iter + 1;
        
        if norm_maxit > 1 % Box-Cox power transformation
            
            if start_boxcox
                fprintf('computing Box-Cox power transform ...')
                boxcox_started = 1;
            end

            i_all_notnorm = find(X_notnormal);
                        
            for i_notnorm = i_all_notnorm
                
                X(:,i_notnorm) = boxcox(X(:,i_notnorm));
                X(:,i_notnorm) = X(:,i_notnorm) - min(X(:,i_notnorm)) + rand(1)*0.01 ;
            end

            if ~passed_normality_test && isempty(ui)
                fprintf('.')
            end

            if i_iter == norm_maxit && ~passed_normality_test
                
                fprintf('\n')
                fprintf('no normality achieved for X(%i) :(\n',i_all_notnorm)
                ui = input('replace with new data and test again? (y/n)','s');
                
                if strcmp(ui,'y') || strcmp(ui,'Y')
                    
                    rep_newdata = 1;
                    i_iter = 1;
                    for i_notnorm = i_all_notnorm
                        
                        X(:,i_notnorm) = randn(n_obs,1);
                        X(:,i_notnorm) = X(:,i_notnorm) - min(X(:,i_notnorm)) + rand(1)*0.1;
                    end
                end
            end
        end
        
        if passed_normality_test
            if boxcox_started
                if ~rep_newdata
                    disp(' ')
                end
                now_label = sprintf('now ');
            else
                now_label = '';
            end
            fprintf('%sall columns of X may be normal :)\n',now_label)
        end
    end
    
    disp(' ')
end

% dependent variables (responses):
Y = repmat(X(:,pcr),1,n_meas/length(pcr)) + (rand(n_obs,n_meas) .^ 2);

Y(i_obsrz,:) = Y(i_obsrz,:) * 0.1;
Y(i_obsrz, logical(round(rand(n_meas,1))) ) = 0;

% dependent variable (average response):
y(:,1) = mean(Y,2);

% select responses:
if g_resp
    Ys = Y;
else
    Ys = y;
end

% divide data in train / test groups:
if train_pc == 100
    X_train = X;
    X_test  = X;
    Ys_train = Ys;
    Ys_test  = Ys;
else    
    length_train = fix(n_obs * train_pc / 100);
    X_train = X( 1 :  length_train, : );
    X_test  = X( length_train + 1 : end, : );
    Ys_train = Ys( 1 :  length_train, : );
    Ys_test  = Ys( length_train + 1 : end, : );
end

% plotting parameters:
marker_size = 20;
font_size = 14;

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
% plot (X,y)

figsize = get(0, 'ScreenSize');
figsize(4) = figsize(4) / 4;
figure('Position',figsize)
for i_sp = 1:size(X,2)

    subplot(1,size(X,2),i_sp)
    scatter(X(:,i_sp),y,'filled')
    ylabel('y (mean Y)')
    xlabel(sprintf('predictor X(%d)',i_sp))
    box on
    axis square
end

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
% Visualize correlation among predictors

figure
imagesc(corr(X))
title('correlation among predictors (X)','FontSize',font_size)

% ..........................................................................

if any(sel_method == 1)
    
    % ..........................................................................
    % Partial Least-Squares Regression
    
    % train:
    ncomp = 2;
    tic
    
    [~,~,~,~,beta_PLS,pctvar_PLS,~,stats_PLS] = plsregress(X_train,Ys_train,ncomp,'cv',cv_fold);
    
    time_PLS = toc;
    fprintf('time for PLS = %g s. \n',time_PLS)
    
    figure
    bar(beta_PLS)
    set(gca,'XTickLabel',{'int.',[1:size(X,2)]},'FontSize',font_size)
    ylabel(sprintf('\\beta'),'FontSize',font_size*2,'FontWeight','bold')
    xlabel('i_X (predictors)')
    title(sprintf('partial least-squares regression'),'FontSize',font_size)
    
    % predict:
    Y_pred_PLS = [ones(size(X_test ,1),1), X_test] * beta_PLS;
    SST = sum( ( Ys_test - mean( Ys_test(:)) ).^2 ); % total sum of squares
    RSS = sum( ( Ys_test - Y_pred_PLS ).^2 ); % residual sum of squares (same as sum of squared errors, SSE)
    rsquared_PLS = 1 - ( RSS / SST ); % explained variance (R^2)
    % rsquared_PLS = 1 - ( sum(sum(stats_PLS.Yresiduals.^2)) / SST) % alternative explained variance (R^2, only when Ys_train is a vector)
    % rsquared_PLS = corr(Y_pred_PLS,Ys_test)^2 % alternative explained variance (R^2, only when Ys_train is a vector)
    % rsquared_PLS = sum(pctvar_PLS(2,:)) % alternative explained variance (R^2, only for training data)
    adjusted_rsquared_PLS = 1 - (1-rsquared_PLS) * ( (n_obs-1)/(n_obs-n_pred-1) );
    
    PLS_LL = ( -(n_obs/2) * log( 2 * pi * (RSS/n_obs) ) )  - ( ( 1/( 2 * (RSS/n_obs) ) ) *  RSS );
    fprintf(' PLS Log Likelihood  = %g \n', PLS_LL )
    
    figure
    plot(Y_pred_PLS,'ob','MarkerSize',marker_size), hold on , plot(Ys_test,'or','MarkerSize',marker_size)
    set(gca,'XLim',[0,size(X,1)+1],'XTick',[1:size(X,1)])
    ylabel('Y (observations'' values)')
    xlabel('i_Y (observations)')
    set(gca,'FontSize',font_size)
    title(...
            sprintf('Partial Least-Squares Regression\npredicted (blue) and test (red) $R^{2}$ = %.2g , $\\bar{R}^{2}$ = %.2g' , ...
            rsquared_PLS         , ...
            adjusted_rsquared_PLS   ...
            ) , ...
        'FontSize',font_size*1.5,'Interpreter', 'LaTeX')
end

if g_resp % the following methods won't work for a matrix of responses
    return 
end

if any(sel_method == 2)
    
    % ..........................................................................
    % Multiple Linear Regression a

    % train:
    tic
    
    [beta_MLRa,~,~,~,stats_MLRa] = regress(Ys_train,[ones(size(X_train ,1),1), X_train]);
    
    time_MLRa = toc;
    fprintf('time for MLRa = %g s. \n',time_MLRa)

    figure
    bar(beta_MLRa)
    set(gca,'XTickLabel',{'int.',[1:size(X,2)]},'FontSize',font_size)
    ylabel(sprintf('\\beta'),'FontSize',font_size*2,'FontWeight','bold')
    xlabel('i_X (predictors)')
    title(sprintf('multiple linear regression'),'FontSize',font_size)

    % predict:

    Y_pred_MLRa = [ones(size(X_test ,1),1), X_test] * beta_MLRa;
    SST = sum( ( Ys_test - mean( Ys_test(:)) ).^2 ); % total sum of squares
    RSS = sum( ( Ys_test - Y_pred_MLRa ).^2 ); % residual sum of squares (same as sum of squared errors, SSE)
    rsquared_MLRa = 1 - ( RSS / SST ); % explained variance (R^2)
    % rsquared_MLRa = corr(Y_pred_MLRa,Ys_test)^2 % alternative explained variance (R^2, only when Ys_train is a vector)
    % rsquared_MLRa = stats_MLRa(1) % alternative explained variance (R^2, only for training data)
    adjusted_rsquared_MLRa = 1 - (1-rsquared_MLRa) * ( (n_obs-1)/(n_obs-n_pred-1) );
    
    MLRa_LL = ( -(n_obs/2) * log( 2 * pi * (RSS/n_obs) ) )  - ( ( 1/( 2 * (RSS/n_obs) ) ) *  RSS );
    fprintf(' MLRa Log Likelihood  = %g \n', MLRa_LL )
    
    figure
    plot(Y_pred_MLRa,'ob','MarkerSize',marker_size), hold on , plot(Ys_test,'or','MarkerSize',marker_size)
    set(gca,'XLim',[0,size(X,1)+1],'XTick',[1:size(X,1)])
    ylabel('Y (observations'' values)')
    xlabel('i_Y (observations)')
    set(gca,'FontSize',font_size)
    title(...
            sprintf('Multiple Linear Regression a\npredicted (blue) and test (red) $R^{2}$ = %.2g , $\\bar{R}^{2}$ = %.2g, $p_{F}$ = %.2g' , ...
                    rsquared_MLRa          , ...
                    adjusted_rsquared_MLRa , ...
                    stats_MLRa(3)   ...
                    ) , ...
             'FontSize',font_size*1.2,'Interpreter', 'LaTeX')
         
end

if any(sel_method == 3)

    % ..........................................................................
    % Multiple Linear Regression b

    % train:
    tic
        
    MLRb_mdl = fitlm( X_train , Ys_train );
    
    time_MLRb = toc;
    fprintf('time for MLRb = %g s. \n',time_MLRb)

    beta_MLRb = MLRb_mdl.Coefficients.Estimate; 

    disp('MLRb: F-test p-values for each beta coefficient:')
    for i_beta = 1:length(beta_MLRb)
        
        if i_beta == 1
            i_coeff_lbl = ' int.';
        else
            i_coeff_lbl = sprintf('   %i)',i_beta-1);
        end
        
        fprintf('%s %.2g \n',i_coeff_lbl,MLRb_mdl.Coefficients.pValue(i_beta))
    end

    summary_MLRb = anova(MLRb_mdl,'summary');
    stats_MLRb(1,1) = MLRb_mdl.Rsquared.Ordinary; % explained variance ( = R^2)
    stats_MLRb(1,2) = summary_MLRb.F(2); % F-statistic for the model
    stats_MLRb(1,3) = summary_MLRb.pValue(2); % p-value of F-statistic for the model
    
    figure
    bar(beta_MLRb)
    set(gca,'XTickLabel',{'int.',[1:size(X,2)]},'FontSize',font_size)
    ylabel(sprintf('\\beta'),'FontSize',font_size*2,'FontWeight','bold')
    xlabel('i_X (predictors)')
    title(sprintf('multiple linear regression'),'FontSize',font_size)

    % predict:
    Y_pred_MLRb = [ones(size(X_test ,1),1), X_test] * beta_MLRb;
    % Y_pred_MLRb = MLRb_mdl.Fitted; % alternative
    SST = sum( ( Ys_test - mean( Ys_test(:)) ).^2 ); % total sum of squares
    % SST = MLRb_mdl.SST; % alternative
    RSS = sum( ( Ys_test - Y_pred_MLRb ).^2 ); % residual sum of squares (same as sum of squared errors, SSE)
    % RSS = MLRb_mdl.SSE; % alternative
    rsquared_MLRb = 1 - ( RSS / SST ); % explained variance (R^2)
    % rsquared_MLRb = MLRb_mdl.Rsquared.Ordinary; % alternative (R^2, only for training data)
    % rsquared_MLRb = corr(Y_pred_MLRb,Ys_test)^2 % alternative explained variance (R^2, only when Ys_train is a vector)
    % rsquared_MLRb = stats_MLRb(1) % alternative explained variance (R^2, only for training data)
    adjusted_rsquared_MLRb = 1 - (1-rsquared_MLRb) * ( (n_obs-1)/(n_obs-n_pred-1) );
    % adjusted_rsquared_MLRb = MLRb_mdl.Rsquared.Adjusted; % alternative adjusted R squared

    MLRb_LL = ( -(n_obs/2) * log( 2 * pi * (RSS/n_obs) ) )  - ( ( 1/( 2 * (RSS/n_obs) ) ) *  RSS );
    % MLRb_LL = MLRb_mdl.LogLikelihood; % alternative
    fprintf(' MLRb Log Likelihood  = %g \n', MLRb_LL )

    figure
    plot(Y_pred_MLRb,'ob','MarkerSize',marker_size), hold on , plot(Ys_test,'or','MarkerSize',marker_size)
    set(gca,'XLim',[0,size(X,1)+1],'XTick',[1:size(X,1)])
    ylabel('Y (observations'' values)')
    xlabel('i_Y (observations)')
    set(gca,'FontSize',font_size)
    title(...
            sprintf('Multiple Linear Regression\npredicted (blue) and test (red) $R^{2}$ = %.2g , $\\bar{R}^{2}$ = %.2g, $p_{F}$ = %.2g' , ...
                    rsquared_MLRb          , ...
                    adjusted_rsquared_MLRb , ...
                    stats_MLRb(3)            ...
                    ) , ...
             'FontSize',font_size*1.2,'Interpreter', 'LaTeX')
end

if any(sel_method == 4)
    
    % ..........................................................................
    % Stepwise Linear Regression

    display_results = 'off'; % 'off' or 'on'
    if strcmp(display_results,'on')
        fprintf('\nStepwise Linear Regression:\n')
    end

    % train:
    tic

    [beta_SLR,~,pval_SLR,included_SLR,stats_SLR,~,history_SLR] = stepwisefit(X_train,Ys_train,'display',display_results);

    time_SLR = toc;
    fprintf('time for SLR = %g s. \n',time_SLR)

    pval_s_SLR = pval_SLR(included_SLR);
    pval_s_SLR(pval_s_SLR < 0.0001) = 0;
    pval_s_SLR_label = sprintf('%0.3g, ',pval_s_SLR);
    pval_s_SLR_label(end-1) = [];

    if sum(included_SLR) > 1
        pval_s_SLR_label = ['[',pval_s_SLR_label,']'];
    end
    int_beta_SLR = [stats_SLR.intercept; beta_SLR];
    col_included_SLR = [0,included_SLR];

    figure
    hold on
    for i_b = 1:length(int_beta_SLR)
        this_bar = zeros(1,length(int_beta_SLR));
        this_bar(i_b) = int_beta_SLR(i_b);
        if col_included_SLR(i_b)
            this_color = [0.2,0.8,0.2];
        else
            this_color = 'flat';
        end
        bar(this_bar,'FaceColor', this_color) ;
    end

    box on
    set(gca,'XTick',[1:size(X,2)+1])
    set(gca,'XTickLabel',{'int.',[1:size(X,2)]},'FontSize',font_size)
    ylabel(sprintf('\\beta'),'FontSize',font_size*2,'FontWeight','bold')
    xlabel('i_X (predictors)')
    title(sprintf('stepwise linear regression \n included p = %s (green)',pval_s_SLR_label),'fontsize',font_size)

    % predict:
    Y_pred_SLR = X_test * (beta_SLR .* included_SLR') + stats_SLR.intercept;
    SST = sum( ( Ys_test - mean( Ys_test(:)) ).^2 ); % total sum of squares
    RSS = sum( ( Ys_test - Y_pred_SLR ).^2 ); % residual sum of squares (same as sum of squared errors, SSE)
    rsquared_SLR = 1 - ( RSS / SST ); % explained variance (R^2)
    % rsquared_SLR = 1 - ( sum(stats_SLR.yr.^2) / SST) % alternative explained variance (R^2, only when Ys_train is a vector)
    % rsquared_SLR = corr(Y_pred_SLR,Ys_test)^2 % alternative explained variance (R^2)
    % rsquared_SLR = 1 - ( stats_SLR.SSresid / stats_SLR.SStotal ) % alternative explained variance (R^2, only for training data)
    % rsquared_SLR = 1 - (history_SLR.rmse(end).^2/var(Ys_test)) .* ((n_obs-1-history_SLR.df0(end))/(n_obs-1)) % alternative explained variance (R^2)
    adjusted_rsquared_SLR = 1 - (1-rsquared_SLR) * ( (n_obs-1)/(n_obs-n_pred-1) );
    
    SLR_LL = ( -(n_obs/2) * log( 2 * pi * (RSS/n_obs) ) )  - ( ( 1/( 2 * (RSS/n_obs) ) ) *  RSS );
    fprintf(' SLR Log Likelihood  = %g \n', SLR_LL )
    
    figure
    plot(Y_pred_SLR,'ob','MarkerSize',marker_size), hold on , plot(Ys_test,'or','MarkerSize',marker_size)
    set(gca,'XLim',[0,size(X,1)+1],'XTick',[1:size(X,1)])
    ylabel('Y (observations'' values)')
    xlabel('i_Y (observations)')
    set(gca,'FontSize',font_size)
    title(...
            sprintf('Stepwise Linear Regression\npredicted (blue) and test (red) $R^{2}$ = %.2g , $\\bar{R}^{2}$ = %.2g, $p_{F}$ = %.2g' , ...
                    rsquared_SLR          , ...
                    adjusted_rsquared_SLR , ...
                    stats_SLR.pval   ...
                    ) , ...
             'FontSize',font_size*1.2,'Interpreter', 'LaTeX')
end