%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%    LINEAR REGRESSION WITH LEAST ABSOLUTE SHRINKAGE AND SELECTION OPERATOR    %
%                                 DEMONSTRATION                                %
%                                                                              %
%                                   JULY 2020                                  %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ==============================================================================
% INFORMATION

% This program requires:
%   'lasso' function, typically found in Mathworks' Statistics and Machine Learning Toolbox

% This program has been tested with:
%   Matlab R2015a

% Instructions:
%   Edit the values indicated with an arrow like this: <--- (length of the arrow can vary)
%   Run the script, close your eyes and hope for the best.

% ==============================================================================
% Parameters for data:

        info.n_pred       = 6;              % <--- number of predictors
        info.n_obs        = 60;             % <--- number of observations

        info.i_nnp        = [];             % <--- non-normal predictors (irrelevant for few observations)
        info.icp          = [];             % <--- intercorrelated predictors
        info.pcr          = [];             % <--- predictors correlated to response
        info.pcr_strength = [];             % <--- approximated absolute correlation of correlated predictors [0...1]
        info.norm_maxit   = 0;              % <--- test predictors' normal distribution and ask if tests fail: 0 = no, 1 = test, N = test and apply power transform up to N times 
        info.norm_auto    = 0;              % <--- automatically replace non-normality and test again: 0 = no, N = maximum iterations
        info.i_obsrz      = [];             % <--- observations whose responses are mostly zero (some random elements added)

        info.make_data    = 1;              % <--- make new data (0 = no, 1 = yes)

% Parameters for LASSO:

              LASSO.CV    = info.n_obs;     % <--- cross-validation folds
              LASSO.DFmax = 6;              % <--- maximum number of non-zero predictor coefficients in the model
              LASSO.alpha = 1;              % <--- Elastic Net [0...1] ; 0 = ridge , 1 = LASSO

% Parameters for RELAXO:

                RELAXO_sw = 0;              % <--- 1 = do RELAXO , 0 = don't

             RELAXO.CV    = info.n_obs;     % <--- cross-validation folds
             RELAXO.DFmax = 6;              % <--- maximum number of non-zero predictor coefficients in the model
             RELAXO.alpha = 1;              % <--- Elastic Net [0...1] ; 0 = ridge , 1 = LASSO

% Other parameters:

                     n_MC = 1000;           % <--- Monte Carlos
 info.OLS_Rsqr_alpha_pval = [.05,.01,.001]; % <--- upper p-value threshold for OLS Rsqr
           record_results = 1;              % <--- Record results in a table (0 = no, 1 = new table, 2 = add to existing table )

                 disp_tbl = 2;              % <--- display table; 0 = no, 1 = on every iteration, 2 = at the end
               disp_plots = 0;              % <--- display plots; 0 = no, 1 = on every iteration

% Results file name (empty to not save):

   results_filename = sprintf('~/montecarlo_npred%i_DFmax%i_CV%i' , info.n_pred, LASSO.DFmax ,LASSO.CV) ;
%    results_filename = [];

% ------------------------------------------------------------------------------
% Parameters to compare LASSO.DFmax > RELAXO.DFmax:
%   Will override DFmax and other parameters specified above

% lasso_comp = [ 12 , 6 ];     % <--- [ LASSO.DFmax , RELAXO.DFmax ] ; 0 = don't
lasso_comp = 0;     % <--- [ LASSO.DFmax , RELAXO.DFmax ] ; 0 = don't

if any(lasso_comp)
    
    if length(lasso_comp) ~= 2
        
        error('Length of ''lasso_comp'' should be 2.')
        
    elseif diff(lasso_comp) >= 0
        
        error('lasso_comp(1) should be greater than lasso_comp(2).')
    end
    
    info.make_data = 1;
    LASSO.DFmax = lasso_comp(1);
    RELAXO_sw = 1;
    RELAXO.DFmax = lasso_comp(2);
    n_MC = 1;
    record_results = 1;
    disp_tbl = 1;
    disp_plots = 0;
    results_filename = [];
    LASSO.sel_DF = [];
end
   
% ------------------------------------------------------------------------------
% Initialise:
 
clc
if ~exist('scrsz','var')
    scrsz = get(0, 'ScreenSize');
    fig_hpos = [1, floor(scrsz(3)/4), floor((2*scrsz(3))/4),floor((3*scrsz(3))/4)];
    fig_w = floor( 0.98 * (scrsz(3)/4));
    fig_h = floor(0.8 * fig_w);
end

if (record_results == 1) || ~exist('LASSO_results','var')
    
    LASSO_results.info = info;
    LASSO_results.LASSO = LASSO;
    LASSO_results.RELAXO = RELAXO;
    LASSO_results.table = table;
end

do_LASSO = true;
start_time = tic;

while do_LASSO
 
    if any(lasso_comp) 
        
        if LASSO.sel_DF > RELAXO.DFmax
            
            info.make_data = 0;
            LASSO.DFmax = RELAXO.DFmax;    
            record_results = 2;
            disp_tbl = 2;
            do_LASSO = false;
        end
    else
        do_LASSO = false;
    end
        
    for i_MC = 1:n_MC

        close all
        disp(' ')

        if n_MC > 1
            fprintf('Monte Carlo %i / %i \n',i_MC,n_MC);
        end

        % --------------------------------------------------------------------------
        % Make artificial data:

        if info.make_data

            fprintf('Make new data: computing...\n');
            lmsg = 13;

            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
            % Independent variables (predictors):

            X = zeros(info.n_obs,info.n_pred); % rows are observations, cols are variables
            for i_pred = 1:info.n_pred
                if any(i_pred == info.i_nnp)
                    X(:,i_pred) = rand(info.n_obs,1);
                else
                    X(:,i_pred) = randn(info.n_obs,1);
                    X(:,i_pred) = X(:,i_pred) - min(X(:,i_pred)) + rand(1)*0.1;
                end
            end

            % correlate predictors:
            for i_icp = 2:length(info.icp)
                X(:,info.icp(i_icp)) = X(:,info.icp(1)) * (1 + rand(1,1)/2);
            end

            i_nonpcr = 1:info.n_pred;
            i_nonpcr(info.pcr) = [];

            if info.norm_maxit

                % test predictors for normality (alpha > 0.05 ):
                % very few observations may render test results irrelevant

                fprintf('NORMALITY TEST:\n')

                X_notnormal = zeros(1,size(X,2));
                passed_normality_test = 0;
                i_normtest_iter = 1;
                boxcox_started = 0;
                ui = [];
                i_autonorm_iter = 0;
                rep_newdata = 0;

                while ~passed_normality_test && (i_normtest_iter <= info.norm_maxit)

                    start_boxcox = 0;
                    for i_pred = i_nonpcr

                        [X_notnormal(i_pred)] = lillietest(X(:,i_pred));
                        % [X_notnormal(i_pred)] = kstest(X(:,i_pred));
                        % [X_notnormal(i_pred)] = jbtest(X(:,i_pred));

                        if X_notnormal(i_pred)

                            if i_normtest_iter == 1
                                fprintf('X(%i) might be not normal :O\n',i_pred)
                                start_boxcox = 1;
                            end
                        end
                    end

                    passed_normality_test = ~any(X_notnormal);
                    i_normtest_iter = i_normtest_iter + 1;

                    if info.norm_maxit > 1 % Box-Cox power transformation

                        if start_boxcox
                            fprintf('computing Box-Cox power transform ...')
                            boxcox_started = 1;
                            lmsg = 0;
                        end

                        i_all_notnorm = find(X_notnormal);

                        for i_notnorm = i_all_notnorm

                            X(:,i_notnorm) = boxcox(X(:,i_notnorm));
                            X(:,i_notnorm) = X(:,i_notnorm) - min(X(:,i_notnorm)) + rand(1)*0.01 ;
                        end

                        if ~passed_normality_test && isempty(ui)
                            fprintf('.')
                        end

                        if i_normtest_iter == info.norm_maxit && ~passed_normality_test

                            fprintf('\n')
                            fprintf('no normality achieved for X(%i) :(\n',i_all_notnorm)

                            if info.norm_auto
                                fprintf('replacing with new data and testing...\n')
                                ui = 'y';
                                i_autonorm_iter = i_autonorm_iter + 1;
                                if i_autonorm_iter >= info.norm_auto
                                    ui = 'n';
                                    warning('Reached maximum automatic iterations for replacing data and testing for normality.')
                                end
                            else
                                ui = input('replace with new data and test again? (y/n)\n','s');
                            end

                            if strcmp(ui,'y') || strcmp(ui,'Y')

                                rep_newdata = 1;
                                i_normtest_iter = 1;
                                for i_notnorm = i_all_notnorm

                                    X(:,i_notnorm) = randn(info.n_obs,1);
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

            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            % dependent variable (response):

            Y = randn(info.n_obs,1);

            if ~isempty(info.pcr)

                pcr_tolerance = 0.001; % deviation tolerace for correlation with correlated predictors [0...1]
                max_pcr_iter = 100;    % maximum iterations to converge

                if info.pcr_strength == 1
                    info.pcr_strength = info.pcr_strength - pcr_tolerance;
                end

                mmx = minmax(X(:)');

                for i_pcr = 1:length(info.pcr)

                    M = cell(max_pcr_iter,1);
                    pcr_dev = zeros(1,max_pcr_iter);
                    pcr_converged = false;

                    for i_pcr_iter = 1:max_pcr_iter

                        M{i_pcr_iter} = [Y, randn(info.n_obs,1)];
                        M{i_pcr_iter} = M{i_pcr_iter} * chol( [1 , info.pcr_strength; info.pcr_strength , 1] );

                        pcr_dev(i_pcr_iter) = abs( corr( M{i_pcr_iter}(:,1) , M{i_pcr_iter}(:,2) ) -  info.pcr_strength );

                        if pcr_dev(i_pcr_iter) <= pcr_tolerance
                            pcr_converged = true;
                            break
                        end
                    end

                    if pcr_converged
                        X(:,info.pcr(i_pcr)) = M{i_pcr_iter}(:,2);
                    else
                        [~, i_min_pcr_dev] = min(pcr_dev);
                        X(:,info.pcr(i_pcr)) = M{i_min_pcr_dev}(:,2);
                    end

                    X(:,info.pcr(i_pcr)) = X(:,info.pcr(i_pcr)) - min( X(:,info.pcr(i_pcr)) );
                    X(:,info.pcr(i_pcr)) = X(:,info.pcr(i_pcr)) * abs(diff(mmx)) / max( X(:,info.pcr(i_pcr)) );
                    X(:,info.pcr(i_pcr)) = X(:,info.pcr(i_pcr)) + mmx(1);
                end
            end

            Y(info.i_obsrz) = Y(info.i_obsrz) * 0.1;
            Y(info.i_obsrz, logical(round(rand(1,1))) ) = 0;

            SST = sum( (Y - mean( Y )).^2 ); % total sum of squares
        end

        if ~exist('X','var') || ~exist('Y','var')

            error('No data. Please set info.make_data = 1')
        else

            if info.make_data

                fprintf(repmat('\b', 1, lmsg));
                fprintf('done\n');
            else
                fprintf('Using data made earlier. \n');
            end
        end

        % --------------------------------------------------------------------------
        % LASSO:

        fprintf('LASSO: computing...\n');

        if isempty(LASSO.CV)
            [LASSO.B,LASSO.STATS] = lasso(X,Y,'DFmax',LASSO.DFmax,'Alpha',LASSO.alpha);
            [~, LASSO.STATS.IndexMinMSE] = min(LASSO.STATS.MSE);
        else
            [LASSO.B,LASSO.STATS] = lasso(X,Y,'CV',LASSO.CV,'DFmax',LASSO.DFmax,'Alpha',LASSO.alpha);
        end

        LASSO.sel_beta = LASSO.B(:,LASSO.STATS.IndexMinMSE);
        LASSO.i_X_sel = LASSO.sel_beta ~= 0;
        LASSO.X_sel = X(:,LASSO.i_X_sel);
        LASSO.i_model_pred = find(LASSO.i_X_sel);
        LASSO.sel_DF = size(LASSO.X_sel,2);

        LASSO.Y_pred = X * LASSO.sel_beta + LASSO.STATS.Intercept(LASSO.STATS.IndexMinMSE);

        LASSO.SSE = sum( ( Y - LASSO.Y_pred ).^2 );
        LASSO.Rsqr = 1 - LASSO.SSE / SST;

        if disp_plots

            fprintf(repmat('\b', 1, 15));
            fprintf(', %i beta coefficients:\n',LASSO.sel_DF + 1)
            for i_b = [0,LASSO.i_model_pred']

                if i_b == 0
                    fprintf('  %g , constant \n',LASSO.STATS.Intercept(LASSO.STATS.IndexMinMSE) )
                else

                    fprintf('  %g , predictor %i \n',LASSO.sel_beta(i_b),i_b )
                end
            end
            disp(' ')

            figure('Position',[fig_hpos(1),400,fig_w,fig_h])
            scatter(Y,LASSO.Y_pred)
            ylabel('predicted')
            xlabel('observed')
            set(gca,'FontSize',18)
            title(sprintf(' LASSO - {\\it R^2} = %g',round(LASSO.Rsqr,2)))
            drawnow
        else

            fprintf(repmat('\b', 1, 13));
            fprintf('done \n');
        end

        % ..........................................................................
        % Ordinary Least Squares:
        % using LASSO coefficients

        fprintf('OLS: computing...\n');

        lastwarn('');

        [OLS.B,~,~,~,OLS.STATS] = regress( Y , [ ones(info.n_obs,1) , LASSO.X_sel ] ) ; % STATS = [Rsqr, F-statistic, p-val, error variance est.]

        OLS.Y_pred = [ ones(info.n_obs,1) , LASSO.X_sel ] * OLS.B;
        OLS.Rsqr = OLS.STATS(1) ;

        if disp_plots

            [warnMsg, warnId] = lastwarn;
            if isempty(warnMsg)
                fprintf(repmat('\b', 1, 15));
                virgule = ',';
            else
                virgule = 'OLS';
            end

            fprintf('%s beta coefficients for LASSO features (p = %.2g):\n',virgule, OLS.STATS(3))
            disp(OLS.B)

            figure('Position',[fig_hpos(2),400,fig_w,fig_h])
            scatter(Y,OLS.Y_pred)
            ylabel('predicted')
            xlabel('observed')
            set(gca,'FontSize',18)
            title(sprintf(' OLS - {\\it R^2} = %g',round(OLS.Rsqr,2)))
            drawnow
        else

            fprintf(repmat('\b', 1, 13));
            fprintf('done \n');
        end
        
        if RELAXO_sw

            % ..................................................................
            % Relaxed LASSO:
            % 2nd. pass of LASSO using first LASSO coefficients

            fprintf('RELAXO: computing...\n');

            if isempty(LASSO.X_sel)

                RELAXO.sel_DF = 0;
                RELAXO.Rsqr = NaN;
                REOLS.Rsqr = NaN;

                fprintf(repmat('\b', 1, 21));
                warning('No beta coefficients for predictors.')
            else

                if isempty(RELAXO.CV)
                    [RELAXO.B,RELAXO.STATS] = lasso(LASSO.X_sel,Y,'DFmax',RELAXO.DFmax,'Alpha',RELAXO.alpha);
                    [~, RELAXO.STATS.IndexMinMSE] = min(RELAXO.STATS.MSE);
                else
                    [RELAXO.B,RELAXO.STATS] = lasso(LASSO.X_sel,Y,'CV',RELAXO.CV,'DFmax',RELAXO.DFmax,'Alpha',RELAXO.alpha);
                end

                RELAXO.sel_beta = RELAXO.B(:,RELAXO.STATS.IndexMinMSE);
                RELAXO.i_X_sel = RELAXO.sel_beta ~= 0;
                RELAXO.i_model_pred =  LASSO.i_model_pred(RELAXO.i_X_sel);
                RELAXO.X_sel = X(:,RELAXO.i_model_pred);
                RELAXO.sel_DF = size(RELAXO.X_sel,2);

                RELAXO.Y_pred = LASSO.X_sel * RELAXO.sel_beta + RELAXO.STATS.Intercept(RELAXO.STATS.IndexMinMSE);

                RELAXO.SSE = sum( ( Y - RELAXO.Y_pred ).^2 );
                RELAXO.Rsqr = 1 - RELAXO.SSE / SST;

                if disp_plots

                    fprintf(repmat('\b', 1, 15));
                    fprintf(', %i beta coefficients:\n',RELAXO.sel_DF + 1)
                    for i_b = [ 0 , find( RELAXO.i_X_sel )' ]

                        if i_b == 0
                            fprintf('  %g , constant \n',RELAXO.STATS.Intercept(RELAXO.STATS.IndexMinMSE) )
                        else

                            fprintf('  %g , predictor %i \n',RELAXO.sel_beta(i_b),LASSO.i_model_pred(i_b) )
                        end
                    end
                    disp(' ')

                    figure('Position',[fig_hpos(3),400,fig_w,fig_h])
                    scatter(Y,RELAXO.Y_pred)
                    ylabel('predicted')
                    xlabel('observed')
                    set(gca,'FontSize',18)
                    title(sprintf(' RELAXO - {\\it R^2} = %g',round(RELAXO.Rsqr,2)))
                    drawnow
                else

                    fprintf(repmat('\b', 1, 13));
                    fprintf('done \n');
                end

                % ..............................................................
                % Ordinary Least Squares:
                % using RELAXO coefficients

                fprintf('REOLS: computing...\n');

                lastwarn('');

                [REOLS.B,~,~,~,REOLS.STATS] = regress( Y , [ ones(info.n_obs,1) , RELAXO.X_sel ] ) ; % STATS = [Rsqr, F-statistic, p-val, error variance est.]

                REOLS.Y_pred = [ ones(info.n_obs,1) , RELAXO.X_sel ] * REOLS.B;
                REOLS.Rsqr = REOLS.STATS(1) ;

                if disp_plots

                    [warnMsg, warnId] = lastwarn;
                    if isempty(warnMsg)
                        fprintf(repmat('\b', 1, 15));
                        virgule = ',';
                    else
                        virgule = 'REOLS';
                    end

                    fprintf('%s beta coefficients for LASSO features (p = %.2g):\n',virgule, REOLS.STATS(3))
                    disp(REOLS.B)

                    figure('Position',[fig_hpos(4),400,fig_w,fig_h])
                    scatter(Y,REOLS.Y_pred)
                    ylabel('predicted')
                    xlabel('observed')
                    set(gca,'FontSize',18)
                    title(sprintf(' REOLS - {\\it R^2} = %g',round(REOLS.Rsqr,2)))
                    drawnow
                else

                    fprintf(repmat('\b', 1, 13));
                    fprintf('done \n');
                end
            end
            
            RELAXO_CV_DFmax = [RELAXO.CV , RELAXO.DFmax];
        else
            
            RELAXO_CV_DFmax = [NaN , NaN];
            RELAXO.sel_DF = NaN;
            RELAXO.Rsqr = NaN;
            REOLS.Rsqr = NaN; 
        end

        % ......................................................................
        
        if record_results

            LASSO_results.table = [ LASSO_results.table ; { ...
                                                           [LASSO.CV , LASSO.DFmax]   , ...
                                                           RELAXO_CV_DFmax            , ...
                                                           LASSO.sel_DF               , ...
                                                           RELAXO.sel_DF              , ...
                                                           round(LASSO.Rsqr,2)        , ...
                                                           round(OLS.Rsqr,2)          , ...
                                                           round(RELAXO.Rsqr,2)       , ...
                                                           round(REOLS.Rsqr,2)          ...
                                                           } ];    
        end


        if ~isempty(LASSO_results.table)

            if i_MC == 1

                LASSO_results.table.Properties.VariableNames = { ...
                                                               'LASSO_CV_DFmax'  , ...
                                                               'RELAXO_CV_DFmax' , ...
                                                               'LASSO_sel_DF'    , ...
                                                               'RELAXO_sel_DF'   , ...
                                                               'LASSO_Rsqr'      , ...
                                                               'OLS_Rsqr'        , ...
                                                               'RELAXO_Rsqr'     , ...
                                                               'REOLS_Rsqr'        ...
                                                               };
            end

            if ( disp_tbl == 1 ) || ( ( disp_tbl == 2 ) && i_MC == n_MC)

                disp(' ')
                disp(LASSO_results.table)
            end

            if ~isempty(results_filename)
                save(results_filename,'LASSO_results')
            end
        end
    end

    fprintf('OLS Rsqr upper p-values:\n')
    
    max_OLS_Rsqr = max(LASSO_results.table.OLS_Rsqr);
    step = 0.001;
    for i_alpha_pval = 1:length(info.OLS_Rsqr_alpha_pval)
        
        n_ge = (( n_MC + 1 ) * info.OLS_Rsqr_alpha_pval(i_alpha_pval) ) - 1;
        
        Rsqr_obs = max_OLS_Rsqr + step;
        n_ge_1 = 0;
        
        while n_ge_1 < n_ge
            Rsqr_obs = Rsqr_obs - step;
            n_ge_1 = sum( LASSO_results.table.OLS_Rsqr >= Rsqr_obs );
        end
        
        LASSO_results.Rsqr_alpha(i_alpha_pval) = Rsqr_obs;
        fprintf('%g, p < %g \n',LASSO_results.Rsqr_alpha(i_alpha_pval),info.OLS_Rsqr_alpha_pval(i_alpha_pval))
    end
    
    if ~isempty(LASSO_results.table) && ~isempty(results_filename)
        
        save(results_filename,'LASSO_results')
    end
end

% ------------------------------------------------------------------------------
% Report computation time:

seconds_elapsed = toc(start_time);
HH = floor(seconds_elapsed / 3600);
seconds_elapsed = seconds_elapsed - HH * 3600;
MM = floor(seconds_elapsed / 60);
SS = seconds_elapsed - MM * 60;
computation_time = sprintf('completed in %02d:%02d:%02.0f (hours:minutes:seconds)',HH,MM,SS);
fprintf('\n%s\n',computation_time)

