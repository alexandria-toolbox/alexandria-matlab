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
    % beta_estimates : matrix of size (k,4)
    %     estimates for beta
    %     column 1: point estimate; column 2: interval lower bound; 
    %     column 3: interval upper bound; column 4: standard deviation
    %
    % X_hat : matrix of size (m,k)
    %     predictors for the model 
    %
    % m : int
    %     number of predicted observations, defined in (3.10.1)
    %
    % forecast_estimates : matrix of size (m,3)
    %     estimates for predictions   
    %     column 1: point estimate; column 2: interval lower bound; 
    %     column 3: interval upper bound
    % 
    % fitted_estimates : vector of size (n,1)
    %     posterior estimates (median) for in sample-fit
    %
    % residual_estimates : vector of size (n,1)
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
        beta_estimates
        X_hat
        m
        forecast_estimates
        fitted_estimates
        residual_estimates
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


        function insample_fit(self)
            
            % insample_fit()
            % generates in-sample fit and residuals along with evaluation criteria
            % 
            % parameters:
            % none
            % 
            % returns:
            % none    
    
            % compute fitted and residuals
            self.fitted_and_residual();
            % compute in-sample criteria
            self.insample_criteria();
        end


        function [forecast_estimates] = forecast(self, X_hat, credibility_level)
            
            % [forecast_estimates] = forecast(self, X_hat, credibility_level)
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
            constant = self.constant;
            trend = self.trend;
            quadratic_trend = self.quadratic_trend;
            % add constant and trends if included
            X_hat = ru.add_intercept_and_trends(X_hat, constant, trend, quadratic_trend, n);
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
            forecast_estimates = zeros(m, 3);   
            % fill estimates
            forecast_estimates(:,1) = location;
            forecast_estimates(:,2) = location - Z * s;
            forecast_estimates(:,3) = location + Z * s;
            % save as attributes
            self.X_hat = X_hat;
            self.m = m;
            self.forecast_estimates = forecast_estimates;
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
            forecast_estimates = self.forecast_estimates;
            % obtain regular forecast evaluation criteria
            y_hat = forecast_estimates(:,1);
            forecast_evaluation_criteria = ru.forecast_evaluation_criteria(y_hat, y);
            % save as attributes
            self.forecast_evaluation_criteria = forecast_evaluation_criteria;
        end


    end   
    
        
    methods (Access = protected, Hidden = true)        
        

        function fit(self)
            
            % estimates beta_hat and sigma_hat from (3.9.7)

            % estimate beta_hat and sigma_hat
            [beta, sigma] = self.ols_regression();
            if self.verbose
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
            beta_estimates = zeros(k,4);
            % fill estimates
            beta_estimates(:,1) = beta;
            beta_estimates(:,2) = beta - Z * s;
            beta_estimates(:,3) = beta + Z * s;
            beta_estimates(:,4) = sqrt(df / (df - 2)) * s;
            % save as attributes
            self.beta_estimates = beta_estimates;
        end


        function fitted_and_residual(self)
            
            % in-sample fitted and residuals
        
            [fitted residual] = ru.fitted_and_residual(self.y, self.X, self.beta_estimates(:,1));
            self.fitted_estimates = fitted;
            self.residual_estimates = residual;
        end
        
        
        function insample_criteria(self)
            
            % in-sample fit evaluation criteria
        
            insample_evaluation = ru.ml_insample_evaluation_criteria(self.y, self.residual_estimates, self.n, self.k, self.sigma);
            self.insample_evaluation = insample_evaluation;
        end


    end
    
    
end

