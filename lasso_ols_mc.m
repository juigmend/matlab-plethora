
function varargout = lasso_ols_mc(npred,nobs,CV,DFmax,niter,varargin)

% Compute Montecarlo simulation and upper p-values for distribution of LASSO-OLS regression.
%
% SYNTAX:
%       M = lasso_ols_mc(npred,nobs,CV,DFmax,niter)
%   [A,M] = lasso_ols_mc(npred,nobs,CV,DFmax,niter,'alpha',alpha)
%   [A,M] = lasso_ols_mc(npred,nobs,CV,DFmax,niter,'alpha',alpha,m)
%   [P,M] = lasso_ols_mc(npred,nobs,CV,DFmax,niter,'Rsqr',Rsqr)
%   [P,M] = lasso_ols_mc(npred,nobs,CV,DFmax,niter,'Rsqr',Rsqr,m)
%
% INPUT:
%   npred  : number of predictors
%   nobs   : number of observations
%   CV     : LASSO cross-validation folds
%   DFmax  : LASSO maximum number of predictors
%   niter  : Montecarlo iterations
%   alpha  : optional, scalar or vector of threshods for upper p-values less than alpha
%   Rsqr   : optional, R squared for which to compute an upper p-value 
%   m      : optional, vector of already computed Montecarlo R-squared, will be deduced from niter
%
% OUTPUT:
%   M      : R-squared for all Montecarlo iterations
%   A      : R-squared for upper p-value < alpha
%   P      : upper p-value for given Rsqr (R squared)
%
% VERSION: 9 August 2020
%
% Juan Ignacio Mendoza
% University of Jyv?skyl?

% ------------------------------------------------------------------------------

if ( length(varargin) >= 2 ) ...
        && strcmp(varargin{1},'Rsqr') ...
        && ( ~isnumeric(varargin{2}) || (length(varargin{2}) > 1) )
    error('Rsqr should be a scalar.')
end

fprintf('Montecarlo simulation for LASSO-OLS regression with %i predictors: computing \n',npred)
lmsg = 1;

M = NaN(niter,1);

if length(varargin) == 3
    
    m = varargin{3};
    lm = length(m);
        
    if lm <= niter
    
        M(1:lm) = m;
        i_iter = lm + 1 : niter;
    else
        M(1:end) = m(1:niter);
        i_iter = [];
    end
else
    i_iter = 1:niter;
end


for i_MC = i_iter
    
    fprintf(repmat('\b', 1, lmsg));
    lmsg = fprintf('%i of %i \n', i_MC,niter);
    
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    % random independent variables (predictors):
    
    X = zeros(nobs,npred); % rows are observations, cols are variables
    for i_pred = 1:npred
        
        X(:,i_pred) = rand(nobs,1);
        X(:,i_pred) = X(:,i_pred) - min(X(:,i_pred)) + rand(1)*0.1;
    end
    
    % random dependent variable (response):
    Y = randn(nobs,1);
    
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    % LASSO:
    
    if isempty(CV)
        [LASSO.B,LASSO.STATS] = lasso(X,Y,'DFmax',DFmax,'Alpha',1);
        [~, LASSO.STATS.IndexMinMSE] = min(LASSO.STATS.MSE);
    else
        [LASSO.B,LASSO.STATS] = lasso(X,Y,'CV',CV,'DFmax',DFmax,'Alpha',1);
    end
    
    LASSO.i_X_sel = LASSO.B(:,LASSO.STATS.IndexMinMSE) ~= 0;
    LASSO.X_sel = X(:,LASSO.i_X_sel);
    
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    % Ordinary Least Squares using LASSO coefficients:
    
    [~,~,~,~,OLS.STATS] = regress( Y , [ ones(nobs,1) , LASSO.X_sel ] ) ; % STATS = [Rsqr, F-statistic, p-val, error variance est.]
    
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    % Record results:
    
    M(i_MC) = OLS.STATS(1);
end

if isempty(varargin)

    varargout{1} = M;
else

    if strcmp(varargin{1},'alpha')

        
        alpha = varargin{2};

        A = NaN(1,length(alpha));
        max_OLS_Rsqr = max(M);
        step = min(alpha)/10;

        for i_alpha_pval = 1:length(alpha)

            n_ge = (( niter + 1 ) * alpha(i_alpha_pval) ) - 1;

            Rsqr_obs = max_OLS_Rsqr + step;
            n_ge_1 = 0;

            while n_ge_1 < n_ge
                Rsqr_obs = Rsqr_obs - step;
                n_ge_1 = sum( M >= Rsqr_obs );
            end

            A(i_alpha_pval) = Rsqr_obs;
        end

        varargout{1} = A;
        
    elseif strcmp(varargin{1},'Rsqr')

        Rsqr = varargin{2};
        
        n_ge = sum( M >= Rsqr );
        P = ( n_ge + 1 ) / ( niter + 1 );

        varargout{1} = P;
    end
    
    varargout{2} = M;
end

fprintf(repmat('\b', 1, 10 + lmsg));
disp('done')

end



