classdef AutocorrelatedBayesianRegression < handle & LinearRegression
    
    
    % Autocorrelated Bayesian linear regression, developed in section 9.6
    % 
    % Parameters:
    % -----------
    % endogenous : vector of size (n,1)
    %     endogenous or explained variable
    % 
    % exogenous : matrix of size (n,k)
    %     exogenous or explanatory variables
    %
    % q : int, default = 1
    %     order of autocorrelation (number of residual lags)
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
    % p : float or vector of size (q,1), default = 0
    %     prior mean, defined in (3.9.62)
    %     if float, value is duplicated to vector of size (q,1)
    %
    % H : float or vector of size (q,1), default = 100
    %     prior variance, defined in (3.9.62)  
    %     if float, value is duplicated to vector of size (q,1)    
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
    % endogenous : vector of size (n,1)
    %     endogenous or explained variable
    % 
    % exogenous : matrix of size (n,k)
    %     exogenous or explanatory variables, defined in (3.9.3)
    %
    % q : int
    %     order of autocorrelation (number of residual lags)
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
    % p : vector of size (q,1)
    %     prior mean, defined in (3.9.62)
    %
    % H : matrix of size (q,q)
    %     prior variance, defined in (3.9.62)     
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
    % T : int
    %     number of observations, defined in (3.9.52)
    % 
    % k : int
    %     dimension of beta, defined in (3.9.1) 
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
    % mcmc_phi : matrix of size (q,iterations)
    %     storage of mcmc values for phi
    %
    % estimates_beta : matrix of size (k,4)
    %     posterior estimates for beta
    %     column 1: interval lower bound; column 2: median; 
    %     column 3: interval upper bound; column 4: standard deviation
    %
    % estimates_sigma : float
    %     posterior estimate for sigma
    %
    % estimates_phi : matrix of size (q,3)
    %     posterior estimates for phi
    %     column 1: interval lower bound; column 2: median; 
    %     column 3: interval upper bound
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
        q
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
        p
        H
        iterations
        burn
        credibility_level
        verbose
        T
        alpha_bar
        mcmc_beta
        mcmc_sigma
        mcmc_phi
        estimates_beta
        estimates_sigma
        estimates_phi
        X_hat
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
        inv_H
        inv_H_p
    end      
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------     
    
    
    methods (Access = public)
        
        
        function self = AutocorrelatedBayesianRegression(endogenous, ...
                exogenous, varargin)  
            
            % constructor for the AutocorrelatedBayesianRegression class

            % allow for optional arguments
            parser = inputParser;
            default_q = 1;            
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
            default_p = 0;
            default_H = 100;
            default_iterations = 2000;
            default_burn = 1000;
            default_credibility_level = 0.95;
            default_verbose = false;
            addRequired(parser, 'endogenous');
            addRequired(parser, 'exogenous');
            addParameter(parser, 'q', default_q);
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
            addParameter(parser, 'p', default_p);
            addParameter(parser, 'H', default_H);       
            addParameter(parser, 'iterations', default_iterations);
            addParameter(parser, 'burn', default_burn);
            addParameter(parser, 'credibility_level', default_credibility_level);
            addParameter(parser, 'verbose', default_verbose);        
            parse(parser, endogenous, exogenous, varargin{:});
            self.endogenous = endogenous;
            self.exogenous = exogenous;
            self.q = parser.Results.q;            
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
            self.p = parser.Results.p;
            self.H = parser.Results.H;
            self.iterations = parser.Results.iterations;
            self.burn = parser.Results.burn;
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
            % predictions for the linear regression model using algorithm 10.3
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
            estimates_forecasts = self.forecast_estimates(mcmc_forecasts, credibility_level);
            % save as attributes
            self.X_hat = X_hat;
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
            T = self.T;            
            % get fit and residual estimates
            estimates_fit = X * beta;
            estimates_residuals = y - X * beta;
            % estimate in-sample prediction criteria from equation (3.10.8)
            res = estimates_residuals;
            ssr = res' * res;
            tss = (y - mean(y))' * (y - mean(y));
            r2 = 1 - ssr / tss;
            adj_r2 = 1 - (1 - r2) * (T - 1) / (T - k);
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
            % log10 marginal likelihood, defined in (3.10.32)
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
            T = self.T;
            k = self.k;
            q = self.q;
            b = self.b;
            inv_V = self.inv_V;            
            V = self.V;
            p = self.p;
            H = self.H;            
            inv_H = self.inv_H;   
            alpha = self.alpha;
            delta = self.delta;
            mcmc_beta = self.mcmc_beta;
            mcmc_sigma = self.mcmc_sigma;
            mcmc_phi = self.mcmc_phi;            
            iterations = self.iterations;
            % generate theta_hat and Sigma_hat
            mcmc_theta = [mcmc_beta; mcmc_sigma; mcmc_phi];
            theta_hat = mean(mcmc_theta, 2);
            Sigma_hat = cov(mcmc_theta');
            inv_Sigma_hat = la.invert_spd_matrix(Sigma_hat);          
            % generate parameters for truncation of the Chi2
            omega = 0.5;
            bound = su.chi2_icdf(omega, k + q + 1);
            % compute the log of first row of (3.10.32)
            J = iterations;
            term_1 = - log(omega * J);
            term_2 = (T - 1) / 2 * log(2 * pi);
            term_3 = -0.5 * log(det(Sigma_hat));            
            term_4 = 0.5 * log(det(V));
            term_5 = 0.5 * log(det(H));
            term_6 = log(gamma(alpha / 2));
            term_7 = - alpha / 2 * log(delta / 2);
            row_1 = - (term_1 + term_2 + term_3 + term_4 + term_5 + term_6 + term_7);
            % for second row of (3.10.32), compute the log of each term in summation
            summation = zeros(iterations, 1);
            for i=1:iterations
                theta = mcmc_theta(:,i);
                beta = mcmc_beta(:,i);
                sigma = mcmc_sigma(i);
                phi = mcmc_phi(:,i);
                quadratic_form = (theta - theta_hat)' * inv_Sigma_hat * (theta - theta_hat);
                if quadratic_form > bound
                    summation(i) = -1000;
                else
                    term_1 = ((alpha + T) / 2 + 1) * log(sigma);
                    residuals = y - X * beta;
                    [epsilon, E] = la.lag_matrix(residuals, q);
                    inv_sigma = 1 / sigma;
                    u = epsilon - E * phi;                    
                    term_2 = 0.5 * u' * u * inv_sigma;
                    term_3 = 0.5 * (beta - b)' * inv_V * (beta - b);
                    term_4 = 0.5 * delta * inv_sigma;
                    term_5 = 0.5 * (phi - p)' * inv_H * (phi - p);
                    term_6 = - 0.5 * quadratic_form;
                    summation(i) = term_1 + term_2 + term_3 + term_4 + term_5 + term_6 ;
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
            
            % generates T defined in (3.9.52)

            % run superclass function
            make_regressors@LinearRegression(self);
            % unpack
            n = self.n;
            q = self.q;
            % get dimensions
            T = n - q;
            % save as attributes
            self.T = T;
        end           
        
        
        function prior(self)
            
            % creates prior elements b, V, p and Z defined in (3.9.10) and (3.9.62)

            % generate b
            self.generate_b();
            % generate V
            self.generate_V();
            % generate p
            self.generate_p();
            % generate H
            self.generate_H();
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
        
        
        function generate_p(self)
        
            % creates prior element p
        
            % unpack
            p = self.p;
            % if p is a scalar, turn it into a vector with p as first value, 0 afterwards
            if isscalar(p)
                q = self.q;
                p = [p; zeros(q-1, 1)];
            end
            % save as attribute
            self.p = p;
        end       
        
        
        function generate_H(self)
            
            % creates prior element H

            % unpack
            p = self.p;
            H = self.H;
            % if H is a scalar, turn it into diagonal array replicating value
            if isscalar(H)
                q = self.q;
                inv_H_p = p / H;
                inv_H = eye(q) / H;
                H = eye(q) * H;
            % if H is a vector, turn it directly into diagonal array
            else
                inv_H_p = p ./ H;
                inv_H = diag(1 ./ H);
                H = diag(H);
            end
            % save as attributes
            self.H = H;
            self.inv_H = inv_H;
            self.inv_H_p = inv_H_p;
        end       
        
        
        function posterior(self)
        
            % creates constant posterior element alpha_bar defined in (3.9.67)

            % unpack
            alpha = self.alpha;
            T = self.T;
            % set value
            alpha_bar = alpha + T;
            % save as attribute
            self.alpha_bar = alpha_bar;
        end        
        
        
        function parameter_mcmc(self)
            
            % posterior distribution for parameters from algorithm 9.3
            
            % unpack
            y = self.y;
            X = self.X;
            k = self.k;
            q = self.q;
            inv_V = self.inv_V;
            inv_V_b = self.inv_V_b;
            alpha_bar = self.alpha_bar;
            delta = self.delta;
            inv_H = self.inv_H;
            inv_H_p = self.inv_H_p;
            iterations = self.iterations;
            burn = self.burn;
            verbose = self.verbose;
            % preallocate storage space
            mcmc_beta = zeros(k, iterations);
            mcmc_sigma = zeros(1, iterations);
            mcmc_phi = zeros(q, iterations);
            total_iterations = iterations + burn;
            % set initial values
            beta = zeros(k, 1);
            inv_sigma = 1;
            phi = zeros(q, 1);
            X_star = la.lag_polynomial(X, phi);
            y_star = la.lag_polynomial(y, phi);
            % run algorithm 9.3 (Gibbs sampling for the model parameters)
            for i = 1:total_iterations
                % draw beta from its conditional posterior
                [beta, res, res_star] = self.draw_beta(inv_V, inv_V_b, inv_sigma, X, y, X_star, y_star);
                % draw sigma from its conditional posterior
                [sigma, inv_sigma] = self.draw_sigma(delta, alpha_bar, res_star);                
                % draw phi from its conditional posterior
                [phi, X_star, y_star] = self.draw_phi(inv_H, inv_H_p, inv_sigma, res, q, X, y);
                % if burn-in sample is over, record value
                if i > burn
                    mcmc_beta(:,i-burn) = beta;
                    mcmc_sigma(:,i-burn) = sigma;
                    mcmc_phi(:,i-burn) = phi;
                end
                % if verbose, display progress bar
                if verbose
                    cu.progress_bar(i, total_iterations, 'Model parameters:');
                end
            end
            % save as attributes
            self.mcmc_beta = mcmc_beta;
            self.mcmc_sigma = mcmc_sigma;
            self.mcmc_phi = mcmc_phi;
        end        
        
        
        function [beta, res, res_star] = draw_beta(self, inv_V, inv_V_b, inv_sigma,X, y, X_star, y_star)
            
            % draw beta from its conditional posterior defined in (3.9.45)
            
            % posterior parameters for beta, defined in (3.9.65)
            inv_V_bar = inv_V + inv_sigma * X_star' * X_star;
            % posterior b_bar
            b_bar_temp = inv_V_b + inv_sigma * X_star' * y_star;
            % efficient sampling of beta (algorithm 9.4)
            beta = rng.efficient_multivariate_normal(b_bar_temp, inv_V_bar);
            % compute residuals, as defined in (3.9.53) and (3.9.57)
            res = y - X * beta;
            res_star = y_star - X_star * beta;
        end       
        
        
        function [sigma, inv_sigma] = draw_sigma(self, delta, alpha_bar, res_star)
            
            % draw sigma from its conditional posterior defined in (3.9.66)

            % compute delta_bar, defined in (3.9.67)
            delta_bar = delta + res_star' * res_star;
            % sample sigma
            sigma = rng.inverse_gamma(alpha_bar / 2, delta_bar / 2);
            inv_sigma = 1 / sigma;
        end        
        
        
        function [phi, X_star, y_star] = draw_phi(self, inv_H, inv_H_p, inv_sigma, res, q, X, y)
            
            % draw phi from its conditional posterior defined in (3.9.69)

            % compute epsilon and E, defined in (3.9.60)
            [epsilon, E] = la.lag_matrix(res, q);
            % posterior parameters for phi, defined in (3.9.70)
            inv_H_bar = inv_H + inv_sigma * E' * E;
            % posterior p_bar
            p_bar_temp = inv_H_p + inv_sigma * E' * epsilon;
            % efficient sampling of phi (algorithm 9.4)
            phi = rng.efficient_multivariate_normal(p_bar_temp, inv_H_bar);
            % update X_star and Y_star
            X_star = la.lag_polynomial(X, phi);
            y_star = la.lag_polynomial(y, phi);             
        end         
        
        
        function parameter_estimates(self)
            
            % point estimates and credibility intervals for model parameters
            % uses quantiles of the empirical posterior distribution

            % unpack
            mcmc_beta = self.mcmc_beta;
            mcmc_sigma = self.mcmc_sigma;
            mcmc_phi = self.mcmc_phi;
            credibility_level = self.credibility_level;
            k = self.k;
            q = self.q;
            % initiate storage: 4 columns: lower bound, median, upper bound, standard deviation
            estimates_beta = zeros(k,4);
            % fill estimates
            estimates_beta(:,1) = quantile(mcmc_beta, (1-credibility_level)/2, 2);
            estimates_beta(:,2) = quantile(mcmc_beta, 0.5, 2);
            estimates_beta(:,3) = quantile(mcmc_beta, (1+credibility_level)/2, 2);
            estimates_beta(:,4) = std(mcmc_beta, 0, 2);
            % get point estimate for sigma
            estimates_sigma = quantile(mcmc_sigma, 0.5);
            % get point estimate for phi
            estimates_phi = zeros(q,4);
            estimates_phi(:,1) = quantile(mcmc_phi, (1-credibility_level)/2, 2);
            estimates_phi(:,2) = quantile(mcmc_phi, 0.5, 2);
            estimates_phi(:,3) = quantile(mcmc_phi, (1+credibility_level)/2, 2);
            estimates_phi(:,4) = std(mcmc_phi, 0, 2);
            % save as attributes
            self.estimates_beta = estimates_beta;
            self.estimates_sigma = estimates_sigma;
            self.estimates_phi = estimates_phi;
        end
        
        
        function [mcmc_forecasts, m] = forecast_mcmc(self, X_hat)
        
            % posterior predictive distribution from algorithm 10.3

            % unpack
            X = self.X;
            y = self.y;
            q = self.q;
            mcmc_beta = self.mcmc_beta;
            mcmc_sigma = self.mcmc_sigma;
            mcmc_phi = self.mcmc_phi;
            iterations = self.iterations;
            verbose = self.verbose;
            % add constant if necessary
            m = size(X_hat, 1);
            X_hat = self.add_intercept_and_trends(X_hat, false);
            % initiate storage, loop over simulations and simulate predictions
            mcmc_forecasts = zeros(m, iterations);        
            for i = 1:iterations
                beta = mcmc_beta(:,i);
                sigma = mcmc_sigma(i);
                phi = mcmc_phi(:,i);
                % get in-sample residuals
                residuals = y - X * beta;                
                % set e_t, defined in (3.9.53), for first out-of-sample period
                e_t = flip(residuals(end+1-q:end));
                % loop over periods and build epsilon_t, defined in (3.9.53)
                epsilon = zeros(m,1);
                for t=1:m
                    % generate epsilon_t
                    u_t = sqrt(sigma) * randn;
                    epsilon_t = e_t' * phi + u_t;
                    % update e_t for next period
                    e_t = [epsilon_t; e_t(1:end-1)];
                    % record value
                    epsilon(t) = epsilon_t;
                end
                % form prediction
                y_hat = X_hat * beta + epsilon;
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