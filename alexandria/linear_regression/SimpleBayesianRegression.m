classdef SimpleBayesianRegression < handle & LinearRegression
    
    % Simplest Bayesian linear regression, developed in section 9.2
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
    % b_exogenous : float or vector of size (n_regressors,1), default = 0
    %     prior mean for regressors
    %
    % V_exogenous : float or vector of size (n_regressors,1), default = 1
    %     prior variance for regressors (positive)
    %
    % b_constant : float, default = 0
    %     prior mean for constant term
    %
    % V_constant : float, default = 1
    %     prior variance for constant (positive)
    % 
    % b_trend : float, default = 0
    %     prior mean for trend
    %
    % V_trend : float, default = 1
    %     prior variance for trend (positive)
    % 
    % b_quadratic_trend : float, default = 0
    %     prior mean for quadratic_trend
    %
    % V_quadratic_trend : float, default = 1
    %     prior variance for quadratic_trend (positive)
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
    % b_exogenous : float or vector of size (n_regressors,1)
    %     prior mean for regressors
    %
    % V_exogenous : float or vector of size (n_regressors,1)
    %     prior variance for regressors (positive)
    %
    % b_constant : float
    %     prior mean for constant term
    %
    % V_constant : float
    %     prior variance for constant (positive)
    % 
    % b_trend : float
    %     prior mean for trend
    %
    % V_trend : float
    %     prior variance for trend (positive)
    % 
    % b_quadratic_trend : float
    %     prior mean for quadratic_trend
    %
    % V_quadratic_trend : float
    %     prior variance for quadratic_trend (positive)
    %
    % b : vector of size (k,1)
    %     prior mean, defined in (3.9.10)
    %
    % V : matrix of size (k,k)
    %     prior variance, defined in (3.9.10)
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
    % sigma : float
    %     residual variance, defined in (3.9.1)
    %
    % b_bar : vector of size (k,1)
    %     posterior mean, defined in (3.9.14)
    %
    % V_bar : matrix of size (k,k)
    %     posterior variance, defined in (3.9.14)    
    % 
    % estimates_beta : matrix of size (k,4)
    %     posterior estimates for beta
    %     column 1: interval lower bound; column 2: median; 
    %     column 3: interval upper bound; column 4: standard deviation
    %
    % X_hat : matrix of size (m,k)
    %     predictors for the model 
    %
    % m : int
    %     number of predicted observations, defined in (3.10.1)
    %
    % estimates_forecasts : matrix of size (m,3)
    %     posterior estimates for predictions   
    %     column 1: interval lower bound; column 2: median; 
    %     column 3: interval upper bound
    % 
    % estimates_fit : vector of size (n,1)
    %     posterior estimates (median) for in sample-fit
    %
    % estimates_residuals : vector of size (n,1)
    %     posterior estimates (median) for residuals
    %
    % insample_evaluation : structure
    %     in-sample fit evaluation (SSR, R2, adj-R2)
    %
    % forecast_evaluation_criteria : structure
    %     out-of-sample forecast evaluation (RMSE, MAE, ...)
    %
    % m_y : float
    %     log10 marginal likelihood
    %
    %
    % Methods
    % ----------
    % estimate
    % forecast
    % fit_and_residuals
    % forecast_evaluation
    % marginal_likelihood
    % optimize_hyperparameters
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
     
    
    properties (GetAccess = public, SetAccess= protected)
        endogenous
        exogenous
        constant
        trend
        quadratic_trend
        b_exogenous
        V_exogenous
        b_constant
        V_constant
        b_trend
        V_trend
        b_quadratic_trend
        V_quadratic_trend        
        b
        V
        credibility_level
        verbose
        sigma
        b_bar
        V_bar
        estimates_beta
        X_hat
        m
        estimates_forecasts
        estimates_fit
        estimates_residuals
        insample_evaluation
        forecast_evaluation_criteria
        m_y
    end
    
    
    properties (GetAccess = private, SetAccess = private)
        inv_sigma_XX
        inv_sigma_Xy
        inv_V
        inv_V_b
        inv_V_bar
        forecast_mean
        forecast_variance
    end
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods (Access = public)
        
        
        function self = SimpleBayesianRegression(endogenous, ...
                exogenous, varargin)
            
            % constructor for the SimpleBayesianRegression class
            
            % allow for optional arguments
            parser = inputParser;
            default_constant = true;
            default_trend = false;
            default_quadratic_trend = false;             
            default_b_exogenous = 0;
            default_V_exogenous = 1;
            default_b_constant = 0;
            default_V_constant = 1;
            default_b_trend = 0;
            default_V_trend = 1;
            default_b_quadratic_trend = 0;
            default_V_quadratic_trend = 1;
            default_credibility_level = 0.95;
            default_verbose = false;
            addRequired(parser, 'endogenous');
            addRequired(parser, 'exogenous');
            addParameter(parser, 'constant', default_constant);
            addParameter(parser, 'trend', default_trend);
            addParameter(parser, 'quadratic_trend', default_quadratic_trend);
            addParameter(parser, 'b_exogenous', default_b_exogenous);
            addParameter(parser, 'V_exogenous', default_V_exogenous);
            addParameter(parser, 'b_constant', default_b_constant);
            addParameter(parser, 'V_constant', default_V_constant);
            addParameter(parser, 'b_trend', default_b_trend);
            addParameter(parser, 'V_trend', default_V_trend);
            addParameter(parser, 'b_quadratic_trend', default_b_quadratic_trend);
            addParameter(parser, 'V_quadratic_trend', default_V_quadratic_trend);            
            addParameter(parser, 'credibility_level', default_credibility_level);
            addParameter(parser, 'verbose', default_verbose);
            parse(parser, endogenous, exogenous,varargin{:});
            self.endogenous = endogenous;
            self.exogenous = exogenous;
            self.constant = parser.Results.constant;
            self.trend = parser.Results.trend;
            self.quadratic_trend = parser.Results.quadratic_trend;
            self.b_exogenous = parser.Results.b_exogenous;
            self.V_exogenous = parser.Results.V_exogenous;
            self.b_constant = parser.Results.b_constant;
            self.V_constant = parser.Results.V_constant;
            self.b_trend = parser.Results.b_trend;
            self.V_trend = parser.Results.V_trend;
            self.b_quadratic_trend = parser.Results.b_quadratic_trend;
            self.V_quadratic_trend = parser.Results.V_quadratic_trend;            
            self.credibility_level = parser.Results.credibility_level;
            self.verbose = parser.Results.verbose;
            % make regressors
            self.make_regressors();
        end
        

        function estimate(self)
            
            % estimate()
            % generates posterior estimates for linear regression model parameters beta and sigma 
            %
            % parameters:
            % none
            %
            % returns:
            % none

            % define prior values
            self.prior();
            % define posterior values
            self.posterior();
            % obtain posterior estimates for regression parameters
            self.parameter_estimates();
        end
        
        
        function [estimates_forecasts] = forecast(self, X_hat, credibility_level)
            
            % [estimates_forecasts] = forecast(X_hat, credibility_level)
            % predictions for the linear regression model using (3.10.4)
            %
            % parameters:
            % X_hat : matrix of shape (m,k)
            %     array of predictors
            % credibility_level : float
            %     credibility level for predictions (between 0 and 1)
            %
            % returns:
            % estimates_forecasts : matrix of size (m,3)
            %     posterior estimates for predictions
            %     column 1: interval lower bound; column 2: median; 
            %     column 3: interval upper bound

            % unpack
            verbose = self.verbose;
            b_bar = self.b_bar;
            sigma = self.sigma;
            V_bar = self.V_bar;    
            % add constant and trends if included
            X_hat = self.add_intercept_and_trends(X_hat, false);
            % obtain forecast mean and variance, defined in (3.10.4)
            m = size(X_hat, 1);
            mean = X_hat * b_bar;
            variance = sigma * eye(m) + X_hat * V_bar * X_hat';
            standard_deviation = sqrt(diag(variance));
            % if verbose, display progress bar
            if verbose
                cu.progress_bar_complete('Predictions:');
            end
            % initiate estimate storage; 3 columns: lower bound, median, upper bound
            estimates_forecasts = zeros(m, 3);
            % critical value of normal distribution for credibility level
            Z = su.normal_icdf((1 + credibility_level) / 2);
            % fill estimates
            estimates_forecasts(:,1) = mean - Z * standard_deviation;
            estimates_forecasts(:,2) = mean;
            estimates_forecasts(:,3) = mean + Z * standard_deviation;
            % save as attributes
            self.X_hat = X_hat;
            self.m = m;
            self.forecast_mean = mean;
            self.forecast_variance = variance;
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
            y = self.y;
            X = self.X;
            beta = self.estimates_beta(:,2);
            k = self.k;
            n = self.n;            
            % get fit and residual estimates
            estimates_fit = X * beta;
            estimates_residuals = y - X * beta;
            % estimate in-sample prediction criteria from equation (3.10.8)
            % estimate in-sample prediction criteria from equation (3.10.8)
            res = estimates_residuals;
            ssr = res' * res;
            tss = (y - mean(y))' * (y - mean(y));
            r2 = 1 - ssr / tss;
            adj_r2 = 1 - (1 - r2) * (n - 1) / (n - k);
            insample_evaluation.ssr = ssr;
            insample_evaluation.r2 = r2;
            insample_evaluation.adj_r2 = adj_r2;
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
            forecast_mean = self.forecast_mean;
            forecast_variance = self.forecast_variance;
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
            % loop over the m predictions
            log_score = zeros(m,1);
            crps = zeros(m,1);
            for i = 1:m
                % get actual, prediction mean, prediction variance
                y_i = y(i);
                mu_i = forecast_mean(i);
                sigma_i = forecast_variance(i,i);
                % get log score from equation (3.10.13)
                [log_pdf, ~] = su.normal_pdf(y_i, mu_i, sigma_i);
                log_score(i) = - log_pdf;
                % get CRPS from equation (3.10.15)
                s_i = sqrt(sigma_i);
                y_tld = (y_i - mu_i) / s_i;
                [~, pdf] = su.normal_pdf(y_tld, 0, 1);
                cdf = su.normal_cdf(y_tld, 0, 1);
                crps(i) = s_i * (y_tld * (2 * cdf - 1) + 2 * pdf - 1 / sqrt(pi));
            end
            log_score = mean(log_score);
            crps = mean(crps);
            forecast_evaluation_criteria = struct;  
            forecast_evaluation_criteria.rmse = rmse;
            forecast_evaluation_criteria.mae = mae;
            forecast_evaluation_criteria.mape = mape;
            forecast_evaluation_criteria.theil_u = theil_u;
            forecast_evaluation_criteria.bias = bias;
            forecast_evaluation_criteria.log_score = log_score;
            forecast_evaluation_criteria.crps = crps;
            % save as attributes
            self.forecast_evaluation_criteria = forecast_evaluation_criteria;
        end
        
        
        function [m_y] = marginal_likelihood(self)
            
            % marginal_likelihood()
            % log10 marginal likelihood, defined in (3.10.20)
            %
            % parameters:
            % none
            %
            % returns:
            % m_y: float
            %     log10 marginal likelihood

            % unpack
            sigma = self.sigma;
            n = self.n;
            y = self.y;
            b = self.b;
            V = self.V;
            inv_sigma_XX = self.inv_sigma_XX;
            inv_V_b = self.inv_V_b;
            inv_V_bar = self.inv_V_bar;
            b_bar = self.b_bar;
            % evaluate the log marginal likelihood from equation (3.10.20)
            term_1 = -(n / 2) * log(2 * pi * sigma);
            term_2 = -0.5 * la.stable_determinant(V * inv_sigma_XX);
            term_3 = -0.5 * (y' * y / sigma + b' * inv_V_b - b_bar' * inv_V_bar * b_bar);
            log_f_y = term_1 + term_2 + term_3;
            % convert to log10
            m_y = log_f_y / log(10);
            % save as attributes
            self.m_y = m_y;
        end
        
        
        function optimize_hyperparameters(self, type)
            
            % optimize_hyperparameters(type)
            % optimize V by maximizing the marginal likelihood
            %
            % parameters
            % type : int 
            %     optimization type (1 or 2)
            %     1 = optimize scalar v; 2 = optimize vector V
            
            % unpack
            constant = self.constant;
            trend = self.trend;
            quadratic_trend = self.quadratic_trend;
            verbose = self.verbose;
            k = self.k;
            % estimate prior elements to get proper estimate of b
            self.prior();
            % optimize: in the simplest case,  V = vI so only scalar v is optimized
            if type == 1
                % initial value of optimizer
                v_0 = 1;
                % bounds for parameter values
                lower_bound = eps;
                upper_bound = 1000;
                % optimize
                [solution, fval, exitflag] = minimize(@self.negative_likelihood_simple_V, ...
                                   [v_0], [], [], [], [], [lower_bound], [upper_bound]);
                % update hyperparameters
                self.V_constant = solution(1);
                self.V_trend = solution(1);
                self.V_quadratic_trend = solution(1);
                self.V_exogenous = solution(1);
            % in the second case, V is diagonal from vector v, so vector v is optimized
            elseif type == 2
                % initial value of optimizer
                v_0 = ones(k, 1);
                % bounds for parameter values
                lower_bound = eps * ones(k, 1);
                upper_bound = 1000 * ones(k, 1);
                % optimize
                [solution, fval, exitflag] = minimize(@self.negative_likelihood_full_V, ...
                                   [v_0], [], [], [], [], [lower_bound], [upper_bound]);
                % update hyperparameters
                self.V_constant = solution(1);
                self.V_trend = solution(1 + constant);
                self.V_quadratic_trend = solution(1 + constant + trend);
                self.V_exogenous = solution(1 + constant + trend + quadratic_trend:end);
            end
            % update V with optimized values
            self.prior();
            % if verbose, display progress bar and success/failure of optimization
            if verbose
                cu.progress_bar_complete('Hyperparameter optimization:');
                cu.optimization_completion(exitflag == 1);
            end
        end
    end
    
    
    methods (Access = protected, Hidden = true)

        
        function make_regressors(self)
            
            % generates sigma defined in (3.9.7)

            % run superclass function
            make_regressors@LinearRegression(self);
            % unpack
            y = self.y;
            X = self.X;
            XX = self.XX;
            Xy = self.Xy;
            n = self.n;
            % obtain maximum likelihood estimates for sigma, and derived hyperparameters
            [~, sigma] = self.ols_regression(y, X, XX, Xy, n);
            inv_sigma_XX = XX / sigma;
            inv_sigma_Xy = Xy / sigma;            
            % save as attributes
            self.sigma = sigma;
            self.inv_sigma_XX = inv_sigma_XX;
            self.inv_sigma_Xy = inv_sigma_Xy;
        end          
        
        
        function prior(self)
            
            % creates prior elements b and V defined in (3.9.10)

            % generate b
            self.generate_b();
            % generate V
            self.generate_V();
        end
        
        
        function generate_b(self)
        
            % creates prior element b
        
            % unpack
            constant = self.constant;
            trend = self.trend;
            quadratic_trend = self.quadratic_trend;
            b_exogenous = self.b_exogenous;
            b_constant = self.b_constant;
            b_trend = self.b_trend;
            b_quadratic_trend = self.b_quadratic_trend;
            n_exogenous = self.n_exogenous;
            % if b_exogenous is a scalar, turn it into a vector replicating the value
            if isscalar(b_exogenous)
                b_exogenous = b_exogenous * ones(n_exogenous, 1);
            end
            b = b_exogenous;
            % if quadratic trend is included, add to prior mean
            if quadratic_trend
                b = [b_quadratic_trend; b];
            end
            % if trend is included, add to prior mean
            if trend
                b = [b_trend; b];
            end
            % if constant is included, add to prior mean
            if constant
                b = [b_constant; b];
            end
            % save as attribute
            self.b = b;
        end 
        
        
        function generate_V(self)
            
            % creates prior element V

            % unpack
            b = self.b;            
            constant = self.constant;
            trend = self.trend;
            quadratic_trend = self.quadratic_trend;
            V_exogenous = self.V_exogenous;
            V_constant = self.V_constant;
            V_trend = self.V_trend;
            V_quadratic_trend = self.V_quadratic_trend;
            n_exogenous = self.n_exogenous;
            % if V_exogenous is a scalar, turn it into a vector replicating the value
            if isscalar(V_exogenous)
                V_exogenous = V_exogenous * ones(n_exogenous, 1);
            end
            V = V_exogenous;
            % if quadratic trend is included, add to prior mean
            if quadratic_trend
                V = [V_quadratic_trend; V];
            end
            % if trend is included, add to prior mean
            if trend
                V = [V_trend; V];
            end
            % if constant is included, add to prior mean
            if constant
                V = [V_constant; V];
            end
            % convert the vector V into an array
            inv_V_b = b ./ V;
            inv_V = diag(1 ./ V);
            V = diag(V);
            % save as attributes
            self.V = V;
            self.inv_V = inv_V;
            self.inv_V_b = inv_V_b;
        end 
            
        
        function posterior(self)
        
            % creates posterior parameters b_bar and V_bar defined in (3.9.14)

            % unpack
            y = self.y;
            X = self.X;
            XX = self.XX;
            Xy = self.Xy;
            n = self.n;
            inv_sigma_XX = self.inv_sigma_XX;
            inv_sigma_Xy = self.inv_sigma_Xy;          
            inv_V = self.inv_V;
            inv_V_b = self.inv_V_b;
            verbose = self.verbose;
            % V_bar, defined in (3.9.14)
            inv_V_bar = inv_V + inv_sigma_XX;
            V_bar = la.invert_spd_matrix(inv_V_bar);
            % b_bar, defined in (3.9.14)
            b_bar = V_bar * (inv_V_b + inv_sigma_Xy);
            % if verbose, display progress bar
            if verbose
                cu.progress_bar_complete('Model parameters:')
            end
            % save as attributes
            self.inv_V_bar = inv_V_bar;
            self.V_bar = V_bar;
            self.b_bar = b_bar;
        end
        
        
        function parameter_estimates(self)
            
            % point estimates and credibility intervals for model parameters
            % use the multivariate normal distribution defined in (3.9.17)

            % unpack
            V_bar = self.V_bar;
            b_bar = self.b_bar;
            credibility_level = self.credibility_level;
            k = self.k;
            % initiate storage: 4 columns: lower bound, median, upper bound, standard deviation
            estimates_beta = zeros(k,4);
            % critical value of normal distribution for credibility level
            Z = su.normal_icdf((1 + credibility_level) / 2);
            % mean and standard deviation of posterior distribution
            mean = b_bar;
            standard_deviation = sqrt(diag(V_bar));
            % fill estimates
            estimates_beta(:,1) = mean - Z * standard_deviation;
            estimates_beta(:,2) = mean;
            estimates_beta(:,3) = mean + Z * standard_deviation;
            estimates_beta(:,4) = standard_deviation;
            % save as attributes
            self.estimates_beta = estimates_beta;
        end
        
        
        function [negative_log_f_y] = negative_likelihood_simple_V(self, v)
            
            % negative log marginal likelihood for scalar V (common variance)

            % unpack
            k = self.k;
            b = self.b;
            inv_sigma_XX = self.inv_sigma_XX;
            inv_sigma_Xy = self.inv_sigma_Xy;
            % build elements for equation (3.10.10)
            V = v * eye(k);
            inv_V = eye(k) / v;
            inv_V_b = b / v;
            % compute log of marginal likelihood, omitting irrelevant terms
            inv_V_bar = inv_V + inv_sigma_XX;
            V_bar = la.invert_spd_matrix(inv_V_bar);
            b_bar = V_bar * (inv_V_b + inv_sigma_Xy);
            term_1 = -0.5 * la.stable_determinant(V * inv_sigma_XX);
            term_2 = -0.5 * (b' * inv_V_b - b_bar' * inv_V_bar * b_bar);
            % take negative (minimize negative to maximize)
            negative_log_f_y = -(term_1 + term_2);
        end
        
        
        function [negative_log_f_y] = negative_likelihood_full_V(self, v)
            
            % negative log marginal likelihood for vector V (individual variances)

            % unpack
            b = self.b;
            inv_sigma_XX = self.inv_sigma_XX;
            inv_sigma_Xy = self.inv_sigma_Xy;
            % build elements for equation (3.10.10)
            V = diag(v);
            inv_V = diag(1./v);
            inv_V_b = b ./ v;
            % compute log of marginal likelihood, omitting irrelevant terms
            inv_V_bar = inv_V + inv_sigma_XX;
            V_bar = la.invert_spd_matrix(inv_V_bar);
            b_bar = V_bar * (inv_V_b + inv_sigma_Xy);
            term_1 = -0.5 * la.stable_determinant(V * inv_sigma_XX);
            term_2 = -0.5 * (b' * inv_V_b - b_bar' * inv_V_bar * b_bar);
            % take negative (minimize negative to maximize)
            negative_log_f_y = -(term_1 + term_2);
        end 
        
        
    end
    
    
end
    
    
