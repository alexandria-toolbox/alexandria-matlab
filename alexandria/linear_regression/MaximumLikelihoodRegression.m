classdef MaximumLikelihoodRegression < handle & LinearRegression

    % Maximum likelihood linear regression, developed in section 9.1
    % 
    % Parameters:
    % -----------
    % endogenous : vector of size (n_obs,1)
    %     endogenous variable, defined in (3.9.1)
    % 
    % exogenous : matrix of size (n_obs,n_regressors)
    %     exogenous variables, defined in (3.9.1)
    % 
    % constant : bool, default = true
    %     if true, an intercept is included in the regression
    %
    % trend : bool, default = false
    %     if true, a linear trend is included in the regression
    %
    % quadratic_trend : bool, default = false
    %     if true, a quadratic trend is included in the regression
    % 
    % credibility_level : float, default = 0.95
    %     credibility level (between 0 and 1)
    % 
    % verbose : bool, default = false
    %     if true, displays a progress bar  
    % 
    % 
    % Properties
    % ----------
    % endogenous : vector of size (n_obs,1)
    %     endogenous variable, defined in (3.9.1)
    % 
    % exogenous : matrix of size (n_obs,n_regressors)
    %     exogenous variables, defined in (3.9.1)
    % 
    % constant : bool
    %     if true, an intercept is included in the regression
    %
    % trend : bool
    %     if true, a linear trend is included in the regression
    %
    % quadratic_trend : bool
    %     if true, a quadratic trend is included in the regression
    % 
    % credibility_level : float
    %     credibility level (between 0 and 1)
    % 
    % verbose : bool
    %     if true, displays a progress bar during MCMC algorithms
    %
    % y : vector of size (n,1)
    %     endogenous variable, defined in (3.9.3)
    % 
    % X : matrix of size (n,k)
    %     exogenous variables, defined in (3.9.3)
    %
    % n : int
    %     number of observations, defined in (3.9.1)
    % 
    % k : int
    %     dimension of beta, defined in (3.9.1)
    % 
    % beta : vector of size (k,1)
    %     regression coefficients
    %
    % sigma : float
    %     residual variance, defined in (3.9.1)
    % 
    % estimates_beta : matrix of size (k,4)
    %     estimates for beta
    %     column 1: interval lower bound; column 2: point estimate; 
    %     column 3: interval upper bound; column 4: standard deviation
    %
    % X_hat : matrix of size (m,k)
    %     predictors for the model 
    %
    % m : int
    %     number of predicted observations, defined in (3.10.1)
    %
    % estimates_forecasts : matrix of size (m,3)
    %     estimates for predictions   
    %     column 1: interval lower bound; column 2: point estimate; 
    %     column 3: interval upper bound
    % 
    % estimates_fit : vector of size (n,1)
    %     posterior estimates (median) for in sample-fit
    %
    % estimates_residuals : vector of size (n,1)
    %     posterior estimates (median) for residuals
    %
    % insample_evaluation : structure
    %     in-sample fit evaluation (SSR, R2, ...)
    %
    % forecast_evaluation_criteria : structure
    %     out-of-sample forecast evaluation (RMSE, MAE, ...)
    %
    %
    % Methods
    % ----------
    % estimate
    % forecast
    % fit_and_residuals
    % forecast_evaluation


    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess= protected)
        endogenous
        exogenous
        constant
        trend
        quadratic_trend
        credibility_level
        verbose
        beta
        sigma
        estimates_beta
        X_hat
        m
        estimates_forecasts
        estimates_fit
        estimates_residuals
        insample_evaluation
        forecast_evaluation_criteria
    end    

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    
    
    methods (Access = public)
 
        
        function [self] = MaximumLikelihoodRegression(endogenous, ...
                          exogenous, varargin)
            
            % constructor for the MaximumLikelihoodRegression class
            
            % allow for optional arguments
            parser = inputParser;
            default_constant = true;
            default_trend = false;
            default_quadratic_trend = false; 
            default_credibility_level = 0.95;
            default_verbose = false;
            addRequired(parser,'endogenous');
            addRequired(parser,'exogenous');            
            addParameter(parser, 'constant', default_constant);
            addParameter(parser, 'trend', default_trend);
            addParameter(parser, 'quadratic_trend', default_quadratic_trend);            
            addParameter(parser,'credibility_level', default_credibility_level);
            addParameter(parser,'verbose', default_verbose);
            parse(parser, endogenous, exogenous,varargin{:});
            self.endogenous = endogenous;
            self.exogenous = exogenous;
            self.constant = parser.Results.constant;
            self.trend = parser.Results.trend;
            self.quadratic_trend = parser.Results.quadratic_trend;            
            self.credibility_level = parser.Results.credibility_level;
            self.verbose = parser.Results.verbose;
            % make regressors
            self.make_regressors();
        end
        
        
        function estimate(self)
            
            % estimate()
            % estimates parameters beta and sigma of linear regression model
            %
            % parameters:
            % none
            %
            % returns:
            % none

            % fit to obtain maximum likelihood estimates
            self.fit();
            % obtain posterior estimates for regression parameters
            self.parameter_estimates();
        end        
        
        
        function [estimates_forecasts] = forecast(self, X_hat, credibility_level)
            
            % [estimates_forecasts] = forecast(X_hat, credibility_level)
            % predictions for the linear regression model using (3.10.2)
            %
            % parameters:
            % X_hat : matrix of shape (m, k)
            %     array of predictors
            % credibility_level : float
            %     credibility level for predictions (between 0 and 1)
            %
            % returns:
            % estimates_forecasts : matrix of size (m, 3)
            %     posterior estimates for predictions
            %     column 1: interval lower bound; column 2: median; 
            %     column 3: interval upper bound

            % unpack
            XX = self.XX;
            verbose = self.verbose;
            beta = self.beta;
            sigma = self.sigma;
            n = self.n;
            k = self.k;
            % add constant and trends if included
            X_hat = self.add_intercept_and_trends(X_hat, false);
            % obtain prediction location, scale and degrees of freedom from (3.10.2)
            m = size(X_hat,1);
            location = X_hat * beta;
            scale = sigma * eye(m) + sigma * X_hat * la.invert_spd_matrix(XX) * X_hat';
            s = sqrt(diag(scale));        
            % if verbose, display progress bar
            if verbose
                cu.progress_bar_complete('Predictions:');
            end
            % critical value of Student distribution for credibility level
            df = n - k;
            Z = su.student_icdf((1 + credibility_level) / 2, df);
            % initiate estimate storage; 3 columns: lower bound, median, upper bound
            estimates_forecasts = zeros(m, 3);   
            % fill estimates
            estimates_forecasts(:,1) = location - Z * s;
            estimates_forecasts(:,2) = location;
            estimates_forecasts(:,3) = location + Z * s;
            % save as attributes
            self.X_hat = X_hat;
            self.m = m;
            self.estimates_forecasts = estimates_forecasts;
        end        
        
        
        function fit_and_residuals(self)
            
            % fit_and_residuals()
            % estimates of in-sample fit and regression residuals
            %
            % parameters:
            % none
            %
            % returns:
            % none

            % unpack
            X = self.X;
            y = self.y;
            beta = self.beta;
            sigma = self.sigma;
            k = self.k;
            n = self.n;
            % estimate fits and residuals
            estimates_fit = X * beta;
            estimates_residuals = y - X * beta;
            % estimate in-sample prediction criteria from equation (3.10.8)
            res = estimates_residuals;
            ssr = res' * res;
            tss = (y - mean(y))' * (y - mean(y));
            r2 = 1 - ssr / tss;
            adj_r2 = 1 - (1 - r2) * (n - 1) / (n - k);
            % additionally estimate AIC and BIC from equation (3.10.9)
            aic = 2 * k / n + log(sigma);
            bic = k * log(n) / n + log(sigma);
            insample_evaluation = struct;            
            insample_evaluation.ssr = ssr;
            insample_evaluation.r2 = r2;
            insample_evaluation.adj_r2 = adj_r2;
            insample_evaluation.aic = aic;
            insample_evaluation.bic = bic;
            % save as attributes
            self.estimates_fit = estimates_fit;
            self.estimates_residuals = estimates_residuals;
            self.insample_evaluation = insample_evaluation;
        end
        
        
        function forecast_evaluation(self, y)
            
            % forecast_evaluation(y)
            % forecast evaluation criteria for the linear regression model
            %
            % parameters:
            % y : vector of shape (m,1)
            %     array of realised values for forecast evaluation
            %
            % returns:
            % none
            
            % unpack
            estimates_forecasts = self.estimates_forecasts;
            m = self.m;
            % calculate forecast error
            y_hat = estimates_forecasts(:,2);
            err = y - y_hat;
            % compute forecast evaluation from equation (3.10.11)
            rmse = sqrt(err' * err / m);
            mae = sum(abs(err)) / m;
            mape = 100 * sum(abs(err ./ y)) / m;
            theil_u = sqrt(err' * err) / (sqrt(y' * y) + sqrt(y_hat' * y_hat));
            bias = sum(err) / sum(abs(err));
            forecast_evaluation_criteria = struct;  
            forecast_evaluation_criteria.rmse = rmse;
            forecast_evaluation_criteria.mae = mae;
            forecast_evaluation_criteria.mape = mape;
            forecast_evaluation_criteria.theil_u = theil_u;
            forecast_evaluation_criteria.bias = bias;
            % save as attributes
            self.forecast_evaluation_criteria = forecast_evaluation_criteria;
        end
            
        
    end   
    
        
    methods (Access = protected, Hidden = true)        
        

        function fit(self)
            
            % estimates beta_hat and sigma_hat from (3.9.7)

            % unpack        
            y = self.y;
            X = self.X;
            XX = self.XX;
            Xy = self.Xy;
            n = self.n;
            verbose = self.verbose;
            % estimate beta_hat and sigma_hat
            [beta, sigma] = self.ols_regression(y, X, XX, Xy, n);
            if verbose
                cu.progress_bar_complete('Model parameters:');
            end
            self.beta = beta;
            self.sigma = sigma;
        end
        
        
        function parameter_estimates(self)
            
            % estimates and intervals from Student distribution in (3.9.8)

            % unpack
            beta = self.beta;
            sigma = self.sigma;
            credibility_level = self.credibility_level;
            n = self.n;
            k = self.k;
            XX = self.XX;
            % obtain scale for Student distribution of beta from equation (3.9.8)
            S = la.invert_spd_matrix(XX);
            s = sqrt(sigma * diag(S));
            % critical value of Student distribution for credibility level
            df = n - k;
            Z = su.student_icdf((1+credibility_level)/2, df);
            % initiate storage: 4 columns: lower bound, median, upper bound, standard deviation
            estimates_beta = zeros(k,4);
            % fill estimates
            estimates_beta(:,1) = beta - Z * s;
            estimates_beta(:,2) = beta;
            estimates_beta(:,3) = beta + Z * s;
            estimates_beta(:,4) = sqrt(df / (df - 2)) * s;
            % save as attributes
            self.estimates_beta = estimates_beta;
        end

        
    end
    
    
end

