classdef IndependentBayesianRegression < handle & LinearRegression & BayesianRegression
    

    % Independent Bayesian linear regression, developed in section 9.4
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
    % alpha : float, default = 1e-4
    %     prior shape, defined in (3.9.21)
    %
    % delta : float, default = 1e-4
    %     prior scale, defined in (3.9.21)
    %
    % iterations : int, default = 2000
    %     post burn-in iterations for MCMC algorithm
    %
    % burn : int, default = 1000
    %     burn-in iterations for MCMC algorithm
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
    % alpha : float
    %     prior shape, defined in (3.9.21)
    %
    % delta : float
    %     prior scale, defined in (3.9.21)
    %
    % iterations : int
    %     post burn-in iterations for MCMC algorithm
    %
    % burn : int
    %     burn-in iterations for MCMC algorithm    
    % 
    % credibility_level : float
    %     credibility level (between 0 and 1)
    % 
    % verbose : bool
    %     if true, displays a progress bar during MCMC algorithms
    %
    % y : vector of size (n,1)
    %     explained variables, defined in (3.9.3)
    % 
    % X : matrix of size (n,k)
    %     explanatory variables, defined in (3.9.3)
    %
    % n : int
    %     number of observations, defined in (3.9.1)
    % 
    % k : int
    %     dimension of beta, defined in (3.9.1)
    %
    % alpha_bar : float
    %     posterior scale, defined in (3.9.35)
    %
    % mcmc_beta : matrix of size (k, iterations)
    %     storage of mcmc values for beta
    %
    % mcmc_sigma : vector of size (iterations, 1)
    %     storage of mcmc values for sigma
    %
    % beta_estimates : matrix of size (k,4)
    %     estimates for beta
    %     column 1: point estimate; column 2: interval lower bound; 
    %     column 3: interval upper bound; column 4: standard deviation
    %
    % sigma_estimates : float
    %     posterior estimate for sigma
    %
    % X_hat : matrix of size (m,k)
    %     predictors for the model 
    %
    % m : int
    %     number of predicted observations, defined in (3.10.1)   
    %
    % mcmc_forecasts : matrix of size (m, iterations)
    %     storage of mcmc values for forecasts
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
        alpha
        delta
        iterations
        burn
        credibility_level
        verbose
        alpha_bar
        mcmc_beta
        mcmc_sigma
        beta_estimates
        sigma_estimates
        X_hat
        m
        mcmc_forecasts
        forecast_estimates
        forecast_evaluation_criteria   
        m_y
    end        


    %---------------------------------------------------
    % Methods
    %---------------------------------------------------        
    
    
    methods (Access = public)
        
        
        function self = IndependentBayesianRegression(endogenous, ...
                exogenous, varargin)    
            
            % constructor for the IndependentBayesianRegression class
            
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
            default_alpha = 1e-4;
            default_delta = 1e-4;
            default_iterations = 2000;
            default_burn = 1000;
            default_credibility_level = 0.95;
            default_verbose = false;
            addRequired(parser,'endogenous');
            addRequired(parser,'exogenous');
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
            addParameter(parser, 'alpha', default_alpha);
            addParameter(parser, 'delta', default_delta);
            addParameter(parser, 'iterations', default_iterations);
            addParameter(parser, 'burn', default_burn);
            addParameter(parser,'credibility_level', default_credibility_level);
            addParameter(parser,'verbose', default_verbose);
            parse(parser, endogenous, exogenous, varargin{:});
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
            self.alpha = parser.Results.alpha;
            self.delta = parser.Results.delta;
            self.iterations = parser.Results.iterations;
            self.burn = parser.Results.burn;
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
            % run MCMC algorithm for regression parameters
            self.parameter_mcmc();
            % obtain posterior estimates for regression parameters
            self.parameter_estimates();
        end
        
        
        function [forecast_estimates] = forecast(self, X_hat, credibility_level)
            
            % [estimates_forecasts] = forecast(X_hat, credibility_level)
            % predictions for the linear regression model using algorithm 10.1
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

            % run mcmc algorithm for predictive density
            [mcmc_forecasts, m] = self.forecast_mcmc(X_hat);
            % obtain posterior estimates
            forecast_estimates = self.make_forecast_estimates(mcmc_forecasts, credibility_level);
            % save as attributes
            self.X_hat = X_hat;
            self.m = m;
            self.mcmc_forecasts = mcmc_forecasts;
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
            mcmc_forecasts = self.mcmc_forecasts;
            forecast_estimates = self.forecast_estimates;
            m = self.m;
            iterations = self.iterations;
            % obtain regular forecast evaluation criteria
            y_hat = forecast_estimates(:,1);
            standard_evaluation_criteria = ru.forecast_evaluation_criteria(y_hat, y);
            % obtain Bayesian forecast evaluation criteria
            bayesian_evaluation_criteria = self.bayesian_forecast_evaluation_criteria(y, mcmc_forecasts, iterations, m);
            % merge structures
            forecast_evaluation_criteria = iu.concatenate_structures(standard_evaluation_criteria, bayesian_evaluation_criteria);
            % save as attributes
            self.forecast_evaluation_criteria = forecast_evaluation_criteria;
        end


        function [m_y] = marginal_likelihood(self)

            % marginal_likelihood()
            % log10 marginal likelihood, defined in (3.10.27)
            %
            % parameters:
            % none
            %
            % returns:
            % m_y: float
            %     log10 marginal likelihood            

            % unpack
            y = self.y;
            X = self.X;
            XX = self.XX;
            Xy = self.Xy;
            b = self.b;
            V = self.V;
            inv_V = self.inv_V;
            inv_V_b = self.inv_V_b;
            alpha = self.alpha;
            delta = self.delta;
            alpha_bar = self.alpha_bar;
            beta_estimates = self.beta_estimates;;
            mcmc_sigma = self.mcmc_sigma;
            n = self.n;
            iterations = self.iterations;
            % generate high density values
            beta = beta_estimates(:,1);
            res = y - X * beta;
            delta_bar = delta + res' * res;
            % generate the vector of summation terms
            summation = zeros(iterations,1);
            VXX = V * XX;
            for i = 1:iterations
                inv_sigma = 1 / mcmc_sigma(i);
                inv_V_bar = inv_V + inv_sigma * XX;
                V_bar = la.invert_spd_matrix(inv_V_bar);
                b_bar = V_bar * (inv_V_b + inv_sigma * Xy);
                summation(i) = 0.5 * (la.stable_determinant(inv_sigma * VXX) ...
                               - (beta - b_bar)' * inv_V_bar * (beta - b_bar));
            end
            % compute the marginal likelihood, part by part
            term_1 = -0.5 * n * log(pi) ...
                - 0.5 * (beta - b)' * inv_V * (beta - b) + log(iterations);
            term_2 = - mu.log_sum_exp(summation);
            term_3 = 0.5 * (alpha * log(delta) - alpha_bar * log(delta_bar));
            term_4 = log(gamma(alpha_bar / 2)) - log(gamma(alpha / 2));
            log_f_y = term_1 + term_2 + term_3 + term_4;
            m_y = log_f_y / log(10);
            self.m_y = m_y;
        end
    end
    
    
    methods (Access = protected, Hidden = true)

        
        function posterior(self)
        
            % creates constant posterior element alpha_bar defined in (3.9.35)

            % unpack
            alpha = self.alpha;
            n = self.n;
            % set value
            alpha_bar = alpha + n;
            % save as attribute
            self.alpha_bar = alpha_bar;
        end
        
        
        function parameter_mcmc(self)
            
            % posterior distribution for parameters from algorithm 9.1
            
            % unpack
            X = self.X;
            y = self.y;
            XX = self.XX;
            Xy = self.Xy;
            k = self.k;
            inv_V = self.inv_V;
            inv_V_b = self.inv_V_b;
            alpha_bar = self.alpha_bar;
            delta = self.delta;
            iterations = self.iterations;
            burn = self.burn;
            verbose = self.verbose;
            % preallocate storage space
            mcmc_beta = zeros(k, iterations);
            mcmc_sigma = zeros(iterations, 1);
            total_iterations = iterations + burn;
            % set initial values
            beta = zeros(k, 1);
            inv_sigma = 1;
            % run algorithm 9.1 (Gibbs sampling for the model parameters)
            for i = 1:total_iterations
                % draw beta from its conditional posterior
                beta = self.draw_beta(inv_V, inv_V_b, inv_sigma, XX, Xy);
                % draw sigma from its conditional posterior
                [sigma, inv_sigma] = self.draw_sigma(y, X, beta, delta, alpha_bar);
                % if burn-in sample is over, record value
                if i > burn
                    mcmc_beta(:,i-burn) = beta;
                    mcmc_sigma(i-burn) = sigma;
                end
                % if verbose, display progress bar
                if verbose
                    cu.progress_bar(i, total_iterations, 'Model parameters:');
                end
            end
            % save as attributes
            self.mcmc_beta = mcmc_beta;
            self.mcmc_sigma = mcmc_sigma;
        end
            
            
        function [beta] = draw_beta(self, inv_V, inv_V_b, inv_sigma, XX, Xy)
            
            % draw beta from its conditional posterior defined in (3.9.33)
            
            % posterior V_bar
            inv_V_bar = inv_V + inv_sigma * XX;
            % posterior b_bar
            b_bar_temp = inv_V_b + inv_sigma * Xy;
            % efficient sampling of beta (algorithm 9.4)
            beta = rn.efficient_multivariate_normal(b_bar_temp, inv_V_bar);
        end
        
        
        function [sigma, inv_sigma] = draw_sigma(self, y, X, beta, delta, alpha_bar)
            
            % draw sigma from its conditional posterior defined in (3.9.35)

            % compute residuals
            residuals = y - X * beta;
            % compute delta_bar
            delta_bar = delta + residuals' * residuals;
            % sample sigma
            sigma = rn.inverse_gamma(alpha_bar / 2, delta_bar / 2);
            inv_sigma = 1 / sigma;
        end
        
        
        function parameter_estimates(self)
            
            % point estimates and credibility intervals for model parameters
            % uses quantiles of the empirical posterior distribution

            % unpack
            mcmc_beta = self.mcmc_beta;
            mcmc_sigma = self.mcmc_sigma;
            credibility_level = self.credibility_level;
            k = self.k;
            % initiate storage: 4 columns: lower bound, median, upper bound, standard deviation
            beta_estimates = zeros(k,4);
            % fill estimates
            beta_estimates(:,1) = quantile(mcmc_beta, 0.5, 2);
            beta_estimates(:,2) = quantile(mcmc_beta, (1-credibility_level)/2, 2);
            beta_estimates(:,3) = quantile(mcmc_beta, (1+credibility_level)/2, 2);
            beta_estimates(:,4) = std(mcmc_beta, 0, 2);
            % get point estimate for sigma
            sigma_estimates = quantile(mcmc_sigma, 0.5);
            % save as attributes
            self.beta_estimates = beta_estimates;
            self.sigma_estimates = sigma_estimates; 
        end
        

        function [mcmc_forecasts, m] = forecast_mcmc(self, X_hat)
        
            % posterior predictive distribution from algorithm 10.1

            % unpack
            mcmc_beta = self.mcmc_beta;
            mcmc_sigma = self.mcmc_sigma;
            iterations = self.iterations;
            verbose = self.verbose;     
            n = self.n;
            constant = self.constant;
            trend = self.trend;
            quadratic_trend = self.quadratic_trend;            
            % add constant and trends if included
            m = size(X_hat, 1);
            X_hat = ru.add_intercept_and_trends(X_hat, constant, trend, quadratic_trend, n);
            % initiate storage, loop over simulations and simulate predictions
            mcmc_forecasts = zeros(m, iterations);        
            for i = 1:iterations
                beta = mcmc_beta(:,i);
                sigma = mcmc_sigma(i);
                y_hat = X_hat * beta + sqrt(sigma) * randn(m,1);
                mcmc_forecasts(:,i) = y_hat;
                if verbose
                    cu.progress_bar(i, iterations, 'Predictions:');
                end
            end
        end
        
        
        function [forecast_estimates] = make_forecast_estimates(self, mcmc_forecasts, credibility_level)
            
            % point estimates and credibility intervals for predictions

            m = size(mcmc_forecasts,1);
            % initiate estimate storage; 3 columns: lower bound, median, upper bound
            forecast_estimates = zeros(m, 3);
            % fill estimates
            forecast_estimates(:,1) = quantile(mcmc_forecasts, 0.5, 2);
            forecast_estimates(:,2) = quantile(mcmc_forecasts, (1-credibility_level)/2, 2);
            forecast_estimates(:,3) = quantile(mcmc_forecasts, (1+credibility_level)/2, 2);
        end


        function [bayesian_forecast_evaluation_criteria] = bayesian_forecast_evaluation_criteria(self, y, mcmc_forecasts, iterations, m)
            
            % Bayesian forecast evaluation criteria from equations from equations (3.10.13) and (3.10.15)
    
            log_score = zeros(m,1);
            crps = zeros(m,1);
            for i = 1:m
                % get actual, prediction mean, prediction variance
                y_i = y(i);
                forecasts = mcmc_forecasts(i,:);
                mu_i = mean(forecasts);
                sigma_i = var(forecasts);
                % get log score from equation (3.10.14)
                [log_pdf, ~] = su.normal_pdf(y_i, mu_i, sigma_i);
                log_score(i) = - log_pdf;
                % get CRPS from equation (3.10.17)
                term_1 = sum(abs(forecasts - y_i));
                term_2 = 0;
                for j = 1:iterations
                    term_2 = term_2 + sum(abs(forecasts(j) - forecasts));
                end
                crps(i) = term_1 / iterations - term_2 / (2 * iterations^2);
            end
            log_score = mean(log_score);
            crps = mean(crps);  
            bayesian_forecast_evaluation_criteria = struct;
            bayesian_forecast_evaluation_criteria.log_score = log_score;
            bayesian_forecast_evaluation_criteria.crps = crps;
        end


    end

    
end
