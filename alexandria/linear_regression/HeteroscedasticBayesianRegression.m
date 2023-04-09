classdef HeteroscedasticBayesianRegression < handle & LinearRegression

    
    % Heteroscedastic Bayesian linear regression, developed in section 9.5
    % 
    % Parameters:
    % -----------
    % endogenous : vector of size (n_obs,1)
    %     endogenous variable, defined in (3.9.1)
    % 
    % exogenous : matrix of size (n_obs,n_regressors)
    %     exogenous variables, defined in (3.9.1)
    %
    % heteroscedastic : matrix of size (n_obs,h), default = exogenous
    %     heteroscedasticity variables, defined in (3.9.37)
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
    % g : float or vector of size (h,1), default = 0
    %     prior mean, defined in (3.9.43)
    %
    % Q : float or vector of size (h,1), default = 100
    %     prior variance, defined in (3.9.43)    
    %
    % tau : float, default = 0.001
    %      variance of the random walk shock, defined in (3.9.50)
    %
    % iterations : int, default = 2000
    %     post burn-in iterations for MCMC algorithm
    %
    % burn : int, default = 1000
    %     burn-in iterations for MCMC algorithm
    % 
    % thinning : bool, default = false
    %     if true, thinning is applied to posterior draws from MCMC algorithm
    % 
    % thinning_frequency : int, default = 10
    %     if thinning is true, retains only one out of so many draws from MCMC algorithm
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
    % heteroscedastic : matrix of size (n_obs,h), default = exogenous
    %     heteroscedasticity variables, defined in (3.9.37)
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
    % g : vector of size (h,1)
    %     prior mean, defined in (3.9.43)
    %
    % Q : matrix of size (h,h)
    %     prior variance, defined in (3.9.43)
    %
    % tau : float
    %      variance of the random walk shock, defined in (3.9.50)    
    %
    % iterations : int
    %     post burn-in iterations for MCMC algorithm
    %
    % burn : int
    %     burn-in iterations for MCMC algorithm
    %
    % thinning : bool
    %     if true, thinning is applied to posterior draws from MCMC algorithm
    % 
    % thinning_frequency : int
    %     if thinning is true, retains only one out of so many draws from MCMC algorithm
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
    % Z : matrix of size (n,h)
    %     heteroscedasticity variables, defined in (3.9.39)
    %
    % n : int
    %     number of observations, defined in (3.9.1)
    % 
    % k : int
    %     dimension of beta, defined in (3.9.1)
    % 
    % h : int
    %     dimension of gamma, defined in (3.9.37)    
    %
    % alpha_bar : float
    %     posterior scale, defined in (3.9.35)
    %
    % mcmc_beta : matrix of size (k,iterations)
    %     storage of mcmc values for beta
    %
    % mcmc_sigma : vector of size (iterations,1)
    %     storage of mcmc values for sigma
    %
    % mcmc_gamma : matrix of size (h,iterations)
    %     storage of mcmc values for gamma
    %
    % estimates_beta : matrix of size (k,4)
    %     posterior estimates for beta
    %     column 1: interval lower bound; column 2: median; 
    %     column 3: interval upper bound; column 4: standard deviation
    %
    % estimates_sigma : float
    %     posterior estimate for sigma
    %
    % estimates_gamma : matrix of size (h,3)
    %     posterior estimates for gamma
    %     column 1: interval lower bound; column 2: median; 
    %     column 3: interval upper bound
    %
    % X_hat : matrix of size (m,k)
    %     predictors for the model 
    %
    % Z_hat : matrix of size (m,h)
    %     heteroscedasticity predictors for the model
    %
    % m : int
    %     number of predicted observations, defined in (3.10.1)   
    %
    % mcmc_forecasts : matrix of size (m, iterations)
    %     storage of mcmc values for forecasts
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
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess= protected)
        endogenous
        exogenous
        heteroscedastic
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
        alpha
        delta
        g
        Q
        tau
        iterations
        burn
        thinning
        thinning_frequency
        credibility_level
        verbose
        Z
        h
        alpha_bar
        mcmc_beta
        mcmc_sigma
        mcmc_gamma
        estimates_beta
        estimates_sigma
        estimates_gamma
        X_hat
        Z_hat
        m
        mcmc_forecasts
        estimates_forecasts
        estimates_fit
        estimates_residuals
        insample_evaluation
        forecast_evaluation_criteria
        m_y
    end      
    
    
    properties (GetAccess = private, SetAccess = private)
        inv_V
        inv_V_b
        inv_Q
    end        
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------     
    
    
    methods (Access = public)
        
        
        function self = HeteroscedasticBayesianRegression(endogenous, ...
                exogenous, varargin)    
            
            % constructor for the HeteroscedasticBayesianRegression class
            
            % allow for optional arguments
            parser = inputParser;
            default_heteroscedastic = nan;
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
            default_g = 0;
            default_Q = 100;
            default_tau = 0.001;
            default_iterations = 2000;
            default_burn = 1000;
            default_thinning = false;
            default_thinning_frequency = 10;
            default_credibility_level = 0.95;
            default_verbose = false;
            addRequired(parser, 'endogenous');
            addRequired(parser, 'exogenous');
            addParameter(parser, 'heteroscedastic', default_heteroscedastic);
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
            addParameter(parser, 'g', default_g);
            addParameter(parser, 'Q', default_Q);
            addParameter(parser, 'tau', default_tau);
            addParameter(parser, 'iterations', default_iterations);
            addParameter(parser, 'burn', default_burn);
            addParameter(parser, 'thinning', default_thinning);
            addParameter(parser, 'thinning_frequency', default_thinning_frequency);
            addParameter(parser, 'credibility_level', default_credibility_level);
            addParameter(parser, 'verbose', default_verbose);
            parse(parser, endogenous, exogenous, varargin{:});
            self.endogenous = endogenous;
            self.exogenous = exogenous;
            if isnan(parser.Results.heteroscedastic)
                self.heteroscedastic = exogenous;
            else
                self.heteroscedastic = parser.Results.heteroscedastic;
            end
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
            self.g = parser.Results.g;
            self.Q = parser.Results.Q;
            self.tau = parser.Results.tau;
            self.iterations = parser.Results.iterations;
            self.burn = parser.Results.burn;
            self.thinning = parser.Results.thinning;
            self.thinning_frequency = parser.Results.thinning_frequency;
            self.credibility_level = parser.Results.credibility_level;
            self.verbose = parser.Results.verbose;
            % make regressors
            self.make_regressors();            
        end        
        
        
        function estimate(self)
            
            % estimate()
            % generates posterior estimates for linear regression model parameters beta, sigma and gamma 
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
        
        
        function [estimates_forecasts] = forecast(self, X_hat, credibility_level, varargin)
            
            % [estimates_forecasts] = forecast(X_hat, credibility_level, varargin)
            % predictions for the linear regression model using algorithm 10.2
            %
            % parameters:
            % X_hat : matrix of shape (m,k)
            %     array of predictors
            % credibility_level : float
            %     credibility level for predictions (between 0 and 1)
            % Z_hat : matrix of shape (m,k), optional
            %     array of heteroscedasticity predictors
            %
            % returns:
            % estimates_forecasts : matrix of size (m,3)
            %     posterior estimates for predictions
            %     column 1: interval lower bound; column 2: median; 
            %     column 3: interval upper bound
            
            % by default, set heteroscedastic regressors to exogenous regressors
            parser = inputParser;
            default_Z_hat = nan;
            addRequired(parser, 'X_hat');
            addRequired(parser, 'credibility_level');
            addParameter(parser, 'Z_hat', default_Z_hat);
            parse(parser, X_hat, credibility_level, varargin{:});
            if isnan(parser.Results.Z_hat)
                Z_hat = X_hat;
            else
                Z_hat = parser.Results.Z_hat;
            end       
            % run mcmc algorithm for predictive density
            [mcmc_forecasts, m] = self.forecast_mcmc(X_hat, Z_hat);
            % obtain posterior estimates
            estimates_forecasts = self.forecast_estimates(mcmc_forecasts, credibility_level);
            % save as attributes
            self.X_hat = X_hat;
            self.Z_hat = Z_hat;
            self.m = m;
            self.mcmc_forecasts = mcmc_forecasts;
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
            res = estimates_residuals;
            ssr = res' * res;
            tss = (y - mean(y))' * (y - mean(y));
            r2 = 1 - ssr / tss;
            adj_r2 = 1 - (1 - r2) * (n - 1) / (n - k);
            insample_evaluation = struct;            
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
            mcmc_forecasts = self.mcmc_forecasts;
            estimates_forecasts = self.estimates_forecasts;
            m = self.m;
            iterations = self.iterations;
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
            % log10 marginal likelihood, defined in (3.10.29)
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
            Z = self.Z;
            n = self.n;
            k = self.k;
            h = self.h;
            b = self.b;
            inv_V = self.inv_V;            
            V = self.V;
            g = self.g;
            Q = self.Q;            
            inv_Q = self.inv_Q;
            alpha = self.alpha;
            delta = self.delta;
            mcmc_beta = self.mcmc_beta;
            mcmc_sigma = self.mcmc_sigma;
            mcmc_gamma = self.mcmc_gamma;            
            iterations = self.iterations;
            % generate theta_hat and Sigma_hat
            mcmc_theta = [mcmc_beta; mcmc_sigma; mcmc_gamma];
            theta_hat = mean(mcmc_theta, 2);
            Sigma_hat = cov(mcmc_theta');
            inv_Sigma_hat = la.invert_spd_matrix(Sigma_hat);          
            % generate parameters for truncation of the Chi2
            omega = 0.5;
            bound = su.chi2_icdf(omega, k + h + 1);
            % compute the log of first row of (3.10.29)
            J = iterations;
            term_1 = - log(omega * J);
            term_2 = (n - 1) / 2 * log(2 * pi);
            term_3 = -0.5 * log(det(Sigma_hat));
            term_4 = 0.5 * log(det(V));
            term_5 = 0.5 * log(det(Q));
            term_6 = log(gamma(alpha / 2));
            term_7 = - alpha / 2 * log(delta / 2);
            row_1 = - (term_1 + term_2 + term_3 + term_4 + term_5 + term_6 + term_7);
            % for second row of (3.10.29), compute the log of each term in summation
            summation = zeros(iterations, 1);
            for i=1:iterations
                theta = mcmc_theta(:,i);
                beta = mcmc_beta(:,i);
                sigma = mcmc_sigma(i);
                gama = mcmc_gamma(:,i);
                quadratic_form = (theta - theta_hat)' * inv_Sigma_hat * (theta - theta_hat);
                if quadratic_form > bound
                    summation(i) = -1000;
                else
                    W = exp(Z * gama);
                    term_1 = 0.5 * sum(log(W));
                    term_2 = ((alpha + n) / 2 + 1) * log(sigma);
                    residuals = y - X * beta;
                    inv_sigma_W = 1 ./ (sigma * W);
                    term_3 = 0.5 * (residuals .* inv_sigma_W)' * residuals;
                    term_4 = 0.5 * (beta - b)' * inv_V * (beta - b);
                    term_5 = 0.5 * delta / sigma;
                    term_6 = 0.5 * (gama - g)' * inv_Q * (gama - g);
                    term_7 = - 0.5 * quadratic_form;
                    summation(i) = term_1 + term_2 + term_3 + term_4 + ...
                        term_5 + term_6 + term_7;
                end
            end
            % turn sum of the logs into log of the sum
            row_2 = - mu.log_sum_exp(summation);
            % sum the two rows and convert to log10
            log_f_y = row_1 + row_2;
            m_y = log_f_y / log(10);
            self.m_y = m_y;
        end
    end
                        
    
    methods (Access = protected, Hidden = true) 
        

        function make_regressors(self)
            
            % generates regressors Z defined in (3.9.39)

            % run superclass function
            make_regressors@LinearRegression(self);
            % unpack
            heteroscedastic = self.heteroscedastic;
            % define Z
            Z = heteroscedastic;
            % get dimensions
            h = size(Z,2);
            % save as attributes
            self.Z = Z;
            self.h = h;
        end        
        
        
        function prior(self)
            
            % creates prior elements b, V, g and Q defined in (3.9.10) and (3.9.43)

            % generate b
            self.generate_b();
            % generate V
            self.generate_V();
            % generate g
            self.generate_g();
            % generate Q
            self.generate_Q();
        end
        
        
        function generate_b(self)
        
            % creates prior element b
        
            % unpack
            exogenous = self.exogenous;            
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
            exogenous = self.exogenous;
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
        
        
        function generate_g(self)
        
            % creates prior element g
        
            % unpack
            g = self.g;
            % if g is a scalar, turn it into a vector replicating the value
            if isscalar(g)
                h = self.h;
                g = g * ones(h, 1);
            end
            % save as attribute
            self.g = g;
        end
        
        
        function generate_Q(self)
            
            % creates prior element Q

            % unpack
            Q = self.Q;
            % if Q is a scalar, turn it into diagonal array replicating value
            if isscalar(Q)
                h = self.h;
                inv_Q = eye(h) / Q;
                Q = eye(h) * Q;
            % if Q is a vector, turn it directly into diagonal array
            else
                inv_Q = diag(1 ./ Q);
                Q = diag(Q);
            end
            % save as attributes
            self.Q = Q;
            self.inv_Q = inv_Q;
        end
        
        
        function posterior(self)
        
            % creates constant posterior element alpha_bar defined in (3.9.48)

            % unpack
            alpha = self.alpha;
            n = self.n;
            % set value
            alpha_bar = alpha + n;
            % save as attribute
            self.alpha_bar = alpha_bar;
        end
        
        
        function parameter_mcmc(self)
            
            % posterior distribution for parameters from algorithm 9.2
            
            % unpack
            y = self.y;
            X = self.X;
            Z = self.Z;
            k = self.k;
            h = self.h;
            inv_V = self.inv_V;
            inv_V_b = self.inv_V_b;
            alpha_bar = self.alpha_bar;
            delta = self.delta;
            iterations = self.iterations;
            burn = self.burn;
            thinning = self.thinning;
            thinning_frequency = self.thinning_frequency;
            verbose = self.verbose;
            % if thinning, multiply iterations accordingly
            if thinning
                iterations = iterations * thinning_frequency;
                burn = burn * thinning_frequency;
            end
            % preallocate storage space
            mcmc_beta = zeros(k, iterations);
            mcmc_sigma = zeros(1, iterations);
            mcmc_gamma = zeros(h, iterations);
            total_iterations = iterations + burn;
            % set initial values
            beta = zeros(k, 1);
            inv_sigma = 1;
            gamma = zeros(h, 1);
            inv_W = 1 ./ exp(Z * gamma);
            acceptance_rate = 0;
            % run algorithm 9.2 (Gibbs sampling for the model parameters)
            for i = 1:total_iterations
                % draw beta from its conditional posterior
                [beta, residuals] = self.draw_beta(inv_V, inv_V_b, inv_W, inv_sigma, X, y);
                % draw sigma from its conditional posterior
                [sigma, inv_sigma] = self.draw_sigma(inv_W, delta, alpha_bar, residuals);
                % draw gamma, using a Metropolis_Hastings step
                [gamma, inv_W, accept] = self.draw_gamma(gamma, residuals, inv_sigma, Z);
                % if burn-in sample is over, record value
                if i > burn
                    acceptance_rate = acceptance_rate + accept;
                    mcmc_beta(:,i-burn) = beta;
                    mcmc_sigma(:,i-burn) = sigma;
                    mcmc_gamma(:,i-burn) = gamma;
                end
                % if verbose, display progress bar
                if verbose
                    cu.progress_bar(i, total_iterations, 'Model parameters:');
                end
            end
            % calculate aceptance rate
            acceptance_rate = acceptance_rate / iterations;
            if verbose
                cu.print_message(['Acceptance rate on Metropolis-Hastings algorithm: ' ...
                    num2str(round(100 * acceptance_rate, 2)) '%.']);
            end
            % trim if thinning is applied
            if thinning
                mcmc_beta = mcmc_beta(:,thinning_frequency:thinning_frequency:end);
                mcmc_sigma = mcmc_sigma(thinning_frequency:thinning_frequency:end);
                mcmc_gamma = mcmc_gamma(:,thinning_frequency:thinning_frequency:end);
            end
            % save as attributes
            self.mcmc_beta = mcmc_beta;
            self.mcmc_sigma = mcmc_sigma;
            self.mcmc_gamma = mcmc_gamma;
        end

        
        function [beta, residuals] = draw_beta(self, inv_V, inv_V_b, inv_W, inv_sigma, X, y)
            
            % draw beta from its conditional posterior defined in (3.9.45)
            
            % posterior parameters for beta, defined in (3.9.46)
            XW = X' .* inv_W';
            inv_V_bar = inv_V + inv_sigma * XW * X;
            % posterior b_bar
            b_bar_temp = inv_V_b + inv_sigma * XW * y;
            % efficient sampling of beta (algorithm 9.4)
            beta = rng.efficient_multivariate_normal(b_bar_temp, inv_V_bar);
            % compute residuals
            residuals = y - X * beta;
        end
        
        
        function [sigma, inv_sigma] = draw_sigma(self, inv_W, delta, alpha_bar, residuals)
            
            % draw sigma from its conditional posterior defined in (3.9.35)

            % compute delta_bar
            delta_bar = delta + residuals' .* inv_W' * residuals;
            % sample sigma
            sigma = rng.inverse_gamma(alpha_bar / 2, delta_bar / 2);
            inv_sigma = 1 / sigma;
        end

        
        function [gamma, inv_W, accept] = draw_gamma(self, gamma, residuals, inv_sigma, Z)
            
            % sample gamma, using a Metropolis_Hastings step

            % run the Metropolis-Hastings algorithm
            [gamma, accept] = self.metropolis_hastings(gamma, residuals, inv_sigma);
            inv_W = 1 ./ exp(Z * gamma);
        end
        
        
        function [gamma, accept] = metropolis_hastings(self, gamma, residuals, inv_sigma)
            
            % Metropolis-hastings step, using (3.9.50) and (3.9.51)

            % unpack
            n = self.n;
            h = self.h;
            Z = self.Z;
            tau = self.tau;
            g = self.g;
            inv_Q = self.inv_Q;
            % get candidate value from (3.9.50)
            gamma_tilde = gamma + sqrt(tau) * randn(h, 1);
            % compute acceptance probability from (3.9.51)
            term_1 = ones(1,n) * Z * (gamma_tilde - gamma);
            term_2 = (residuals .* (exp(-Z * gamma_tilde) - exp(-Z * gamma)))' ...
                     * residuals * inv_sigma;
            term_3 = (gamma_tilde - g)' * inv_Q * (gamma_tilde - g);
            term_4 = - (gamma - g)' * inv_Q * (gamma - g);
            prob = min(1, exp(-0.5 * (term_1 + term_2 + term_3 + term_4)));
            % get uniform random number to decide on acceptance or rejection
            u = rand;
            if u <= prob
                gamma = gamma_tilde;
                accept = 1;
            else
                accept = 0;
            end
        end
        
        
        function parameter_estimates(self)
            
            % point estimates and credibility intervals for model parameters
            % uses quantiles of the empirical posterior distribution

            % unpack
            mcmc_beta = self.mcmc_beta;
            mcmc_sigma = self.mcmc_sigma;
            mcmc_gamma = self.mcmc_gamma;
            credibility_level = self.credibility_level;
            k = self.k;
            h = self.h;
            % initiate storage: 4 columns: lower bound, median, upper bound, standard deviation
            estimates_beta = zeros(k,4);
            % fill estimates
            estimates_beta(:,1) = quantile(mcmc_beta, (1-credibility_level)/2, 2);
            estimates_beta(:,2) = quantile(mcmc_beta, 0.5, 2);
            estimates_beta(:,3) = quantile(mcmc_beta, (1+credibility_level)/2, 2);
            estimates_beta(:,4) = std(mcmc_beta, 0, 2);
            % get median for sigma
            estimates_sigma = quantile(mcmc_sigma, 0.5);
            % estimates for gamma
            estimates_gamma = zeros(h,4);
            estimates_gamma(:,1) = quantile(mcmc_gamma, (1-credibility_level)/2, 2);
            estimates_gamma(:,2) = quantile(mcmc_gamma, 0.5, 2);
            estimates_gamma(:,3) = quantile(mcmc_gamma, (1+credibility_level)/2, 2);
            estimates_gamma(:,4) = std(mcmc_gamma, 0, 2);
            % save as attributes
            self.estimates_beta = estimates_beta;
            self.estimates_sigma = estimates_sigma;
            self.estimates_gamma = estimates_gamma;
        end
        
        
        function [mcmc_forecasts, m] = forecast_mcmc(self, X_hat, Z_hat)
        
            % posterior predictive distribution from algorithm 10.2

            % unpack
            mcmc_beta = self.mcmc_beta;
            mcmc_sigma = self.mcmc_sigma;
            mcmc_gamma = self.mcmc_gamma;
            iterations = self.iterations;
            verbose = self.verbose;
            % add constant and trends if included
            m = size(X_hat, 1);
            X_hat = self.add_intercept_and_trends(X_hat, false);
            % initiate storage, loop over simulations and simulate predictions
            mcmc_forecasts = zeros(m, iterations);        
            for i = 1:iterations
                beta = mcmc_beta(:,i);
                sigma = mcmc_sigma(i);
                gamma = mcmc_gamma(:,i);
                W = exp(Z_hat * gamma);
                y_hat = X_hat * beta + sqrt(sigma * W) .* randn(m,1);
                mcmc_forecasts(:,i) = y_hat;
                if verbose
                    cu.progress_bar(i, iterations, 'Predictions:');
                end
            end
        end
            
            
        function [estimates_forecasts] = forecast_estimates(self, mcmc_forecasts, credibility_level)
            
            % point estimates and credibility intervals for predictions

            m = size(mcmc_forecasts,1);
            % initiate estimate storage; 3 columns: lower bound, median, upper bound
            estimates_forecasts = zeros(m, 3);
            % fill estimates
            estimates_forecasts(:,1) = quantile(mcmc_forecasts, (1-credibility_level)/2, 2);
            estimates_forecasts(:,2) = quantile(mcmc_forecasts, 0.5, 2);
            estimates_forecasts(:,3) = quantile(mcmc_forecasts, (1+credibility_level)/2, 2);
        end
    end    
end
