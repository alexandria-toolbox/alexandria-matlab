classdef HierarchicalBayesianRegression < handle & LinearRegression & BayesianRegression
    

    % Hierarchical Bayesian linear regression, developed in section 9.3
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
    % hyperparameter_optimization : bool
    %     if true, performs hyperparameter optimization by marginal likelihood
    %    
    % optimization_type : int, 1 or 2
    %     if 1, simple optimization (scalar v); if 2, full optimization (vector V)  
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
    % b_bar : vector of size (k,1)
    %     posterior mean, defined in (3.9.14)
    %
    % V_bar : matrix of size (k,k)
    %     posterior variance, defined in (3.9.14)  
    %
    % alpha_bar : float
    %     posterior shape, defined in (3.9.24)
    %
    % delta_bar : float
    %     posterior scale, defined in (3.9.24)    
    % 
    % location : vector of size (k,1)
    %     location for the student posterior of beta, defined in (3.9.28)
    %
    % scale : matrix of size (k,k)
    %     scale for the student posterior of beta, defined in (3.9.28)
    %
    % df : float
    %     degrees of freedom for the student posterior of beta, defined in (3.9.28)
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
        credibility_level
        verbose
        b_bar;        
        V_bar;
        alpha_bar
        delta_bar
        hyperparameter_optimization
        optimization_type
        location
        scale
        df
        beta_estimates
        sigma_estimates
        X_hat
        m
        forecast_estimates
        forecast_evaluation_criteria        
        m_y
    end    
    
    
    properties (GetAccess = private, SetAccess = private)  
        sigma
        inv_V_bar
        forecast_location
        forecast_scale
    end    
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------    
    
    
    methods (Access = public)    
    
    
        function self = HierarchicalBayesianRegression(endogenous, ...
                exogenous, varargin)    
            
            % constructor for the HierarchicalBayesianRegression class
            
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
            default_hyperparameter_optimization = false;
            default_optimization_type = 1;
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
            addParameter(parser, 'hyperparameter_optimization', default_hyperparameter_optimization); 
            addParameter(parser, 'optimization_type', default_optimization_type);             
            addParameter(parser,'credibility_level', default_credibility_level);
            addParameter(parser,'verbose', default_verbose);
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
            self.alpha = parser.Results.alpha;
            self.delta = parser.Results.delta;
            self.hyperparameter_optimization = parser.Results.hyperparameter_optimization;  
            self.optimization_type = parser.Results.optimization_type;  
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

            % optimize hyperparameters, if applicable
            self.optimize_hyperparameters();
            % define prior values
            self.prior();
            % define posterior values
            self.posterior();
            % obtain posterior estimates for regression parameters
            self.parameter_estimates();
        end
        
        
        function [forecast_estimates] = forecast(self, X_hat, credibility_level)
            
            % [forecast_estimates] = forecast(self, X_hat, credibility_level)
            % predictions for the linear regression model using (3.10.6)
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
            X = self.X;
            verbose = self.verbose;
            b_bar = self.b_bar;
            V_bar = self.V_bar;
            alpha_bar = self.alpha_bar;
            delta_bar = self.delta_bar;
            n = self.n;
            k = self.k;
            constant = self.constant;
            trend = self.trend;
            quadratic_trend = self.quadratic_trend;
            % add constant and trends if included
            X_hat = ru.add_intercept_and_trends(X_hat, constant, trend, quadratic_trend, n);
            % obtain prediction location, scale and degrees of freedom from (3.10.6)
            m = size(X_hat, 1);
            location = X_hat * b_bar;
            scale = (delta_bar / alpha_bar) * (eye(m) + X_hat * V_bar * X_hat');
            df = alpha_bar;
            % if verbose, display progress bar
            if verbose
                cu.progress_bar_complete('Predictions:')
            end
            % initiate estimate storage; 3 columns: lower bound, median, upper bound
            forecast_estimates = zeros(m, 3);
            % critical value of Student distribution for credibility level
            Z = su.student_icdf((1 + credibility_level) / 2, df);
            % scale in textbook is square of scale for Python: take square root
            sqrt_scale = sqrt(diag(scale));
            % fill estimates
            forecast_estimates(:,1) = location;
            forecast_estimates(:,2) = location - Z * sqrt_scale;
            forecast_estimates(:,3) = location + Z * sqrt_scale;
            % save as attributes
            self.X_hat = X_hat;
            self.m = m;
            self.forecast_location = location;
            self.forecast_scale = scale;
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
            forecast_location = self.forecast_location;
            forecast_scale = self.forecast_scale;
            nu_i = self.alpha_bar;
            forecast_estimates = self.forecast_estimates;
            m = self.m;
            % obtain regular forecast evaluation criteria
            y_hat = forecast_estimates(:,1);
            standard_evaluation_criteria = ru.forecast_evaluation_criteria(y_hat, y);      
            % obtain Bayesian forecast evaluation criteria
            bayesian_evaluation_criteria = self.bayesian_forecast_evaluation_criteria(y, forecast_location, forecast_scale, nu_i, m);
            % merge structures
            forecast_evaluation_criteria = iu.concatenate_structures(standard_evaluation_criteria, bayesian_evaluation_criteria);
            % save as attributes
            self.forecast_evaluation_criteria = forecast_evaluation_criteria;
        end

        
        function [m_y] = marginal_likelihood(self)
            
            % marginal_likelihood()
            % log10 marginal likelihood, defined in (3.10.24)
            %
            % parameters:
            % none
            %
            % returns:
            % m_y: float
            %     log10 marginal likelihood

            % unpack
            n = self.n;
            alpha = self.alpha;
            delta = self.delta;
            V = self.V;
            alpha_bar = self.alpha_bar;
            delta_bar = self.delta_bar;
            XX = self.XX;
            % evaluate the log marginal likelihood from equation (3.10.24)
            term_1 = -(n / 2) * log(pi);
            term_2 = -0.5 * la.stable_determinant(V * XX);
            term_3 = (alpha / 2) * log(delta) - (alpha_bar / 2) * log(delta_bar);
            term_4 = log(gamma(alpha_bar / 2)) - log(gamma(alpha / 2));
            log_f_y = term_1 + term_2 + term_3 + term_4;
            % convert to log10
            m_y = log_f_y / log(10);
            % save as attributes
            self.m_y = m_y;
        end
        

    end
        

    methods (Access = protected, Hidden = true)


        function optimize_hyperparameters(self, type)
            
            % optimize hyperparameter V with marginal likelihood, either scalar v or vector V

            if self.hyperparameter_optimization
                % unpack
                optimization_type = self.optimization_type;
                constant = self.constant;
                trend = self.trend;
                quadratic_trend = self.quadratic_trend;            
                verbose = self.verbose;
                k = self.k;
                % estimate prior elements to get proper estimate of b
                self.prior();
                % optimize: in the simplest case,  V = vI so only scalar v is optimized
                if optimization_type == 1
                    % initial value of optimizer
                    x_0 = ones(2,1);
                    % bounds for parameter values
                    lower_bound = eps * ones(2, 1);
                    upper_bound = 1000 * ones(2, 1);
                    % optimize
                    [solution, fval, exitflag] = minimize(@self.negative_likelihood_simple_V, ...
                                       [x_0], [], [], [], [], [lower_bound], [upper_bound]);
                    % generate V again and save as attribute
                    self.V_constant = solution(1,1);
                    self.V_trend = solution(1,1);
                    self.V_quadratic_trend = solution(1,1);
                    self.V_exogenous = solution(1,1);                
                    self.delta = solution(2,1);
                % in the second case, V is diagonal from vector v, so vector v is optimized
                elseif optimization_type == 2
                    % initial value of optimizer
                    x_0 = ones(k + 1, 1);
                    % bounds for parameter values
                    lower_bound = eps * ones(k + 1, 1);
                    upper_bound = 1000 * ones(k + 1, 1);
                    % optimize
                    [solution, fval, exitflag] = minimize(@self.negative_likelihood_full_V, ...
                                       [x_0], [], [], [], [], [lower_bound], [upper_bound]);
                    % save as attribute
                    if constant
                        self.V_constant = solution(1,1);
                    end
                    if trend
                        self.V_trend = solution(1+constant,1);
                    end
                    if quadratic_trend
                        self.V_quadratic_trend = solution(1+constant+trend,1);
                    end
                    self.V_exogenous = solution(1+constant+trend+quadratic_trend:end-1);
                    self.delta = solution(end,1);
                end
                % if verbose, display progress bar and success/failure of optimization
                if verbose
                    cu.progress_bar_complete('Hyperparameter optimization:');
                    cu.optimization_completion(exitflag == 1);
                end
            end
        end

        
        function posterior(self)
        
            % creates posterior parameters b_bar and V_bar defined in (3.9.14)

            % unpack
            XX = self.XX;
            Xy = self.Xy;
            yy = self.yy;
            n = self.n;            
            b = self.b;
            inv_V = self.inv_V;
            inv_V_b = self.inv_V_b;
            alpha = self.alpha;
            delta = self.delta;
            verbose = self.verbose;
            % V_bar, defined in (3.9.24)
            inv_V_bar = inv_V + XX;
            V_bar = la.invert_spd_matrix(inv_V_bar);
            % b_bar, defined in (3.9.24)
            b_bar = V_bar * (inv_V_b + Xy);
            % alpha_bar, defined in (3.9.24)
            alpha_bar = alpha + n;
            % delta_bar, defined in (3.9.24)
            delta_bar = delta + yy + b' * inv_V_b - b_bar' * inv_V_bar * b_bar;
            % posterior location, defined in (3.9.28)
            location = b_bar;
            % posterior scale, defined in (3.9.28)
            scale = (delta_bar / alpha_bar) * V_bar;
            % posterior degrees of freedom, defined in (3.9.29)
            df = alpha_bar;
            % if verbose, display progress bar
            if verbose
                cu.progress_bar_complete('Model parameters:')
            end
            % save as attributes
            self.inv_V_bar = inv_V_bar;
            self.V_bar = V_bar;
            self.b_bar = b_bar;
            self.alpha_bar = alpha_bar;
            self.delta_bar = delta_bar;
            self.location = location;
            self.scale = scale;
            self.df = df; 
        end        
        

        function parameter_estimates(self)
            
            % point estimates and credibility intervals for model parameters
            % use the multivariate normal distribution defined in (3.9.28)

            % unpack
            alpha_bar = self.alpha_bar;
            delta_bar = self.delta_bar;
            location = self.location;
            scale = self.scale;
            df = self.df;
            credibility_level = self.credibility_level;
            k = self.k;
            % initiate storage: 4 columns: lower bound, median, upper bound, standard deviation
            beta_estimates = zeros(k,4);
            % critical value of Student distribution for credibility level
            Z = su.student_icdf((1 + credibility_level) / 2, df);
            % scale in textbook is square of scale for Python: take square root
            sqrt_scale = sqrt(diag(scale));
            % fill estimates
            beta_estimates(:,1) = location;
            beta_estimates(:,2) = location - Z * sqrt_scale;
            beta_estimates(:,3) = location + Z * sqrt_scale;
            beta_estimates(:,4) = sqrt(df / (df - 2)) * sqrt_scale;
            % get point estimate for sigma (approximate median from Table d.17)
            shape = alpha_bar / 2;
            scale = delta_bar / 2;
            sigma_estimates = (scale * (3 * shape + 0.2)) / (shape * (3 * shape - 0.8));
            % save as attributes
            self.beta_estimates = beta_estimates;
            self.sigma_estimates = sigma_estimates; 
        end
        
        
        function [negative_log_f_y] = negative_likelihood_simple_V(self, x)
            
            % negative log marginal likelihood for scalar V and delta

            % unpack
            v = x(1,1);
            delta = x(2,1);
            k = self.k;
            n = self.n;
            b = self.b;
            alpha = self.alpha;
            XX = self.XX;
            yy = self.yy;
            Xy = self.Xy;
            % build elements for equation (3.10.14)
            V = v * eye(k);
            inv_V = eye(k) / v;
            inv_V_b = b / v;
            inv_V_bar = inv_V + XX;
            V_bar = la.invert_spd_matrix(inv_V_bar);
            b_bar = V_bar * (inv_V_b + Xy);
            alpha_bar= alpha + n;
            delta_bar = delta + yy + b' * inv_V_b - b_bar' * inv_V_bar * b_bar;
            % compute log of marginal likelihood, omitting irrelevant terms
            term_1 = -0.5 * la.stable_determinant(V * XX);
            term_2 = (alpha / 2) * log(delta) - (alpha_bar / 2) * log(delta_bar);
            % take negative (minimize negative to maximize)
            negative_log_f_y = -(term_1 + term_2);
        end
        
        
        function [negative_log_f_y] = negative_likelihood_full_V(self, x)
            
            % negative log marginal likelihood for vector V and delta

            % unpack
            v = x(1:end-1,1);
            delta = x(end,1);
            n = self.n;
            b = self.b;
            alpha = self.alpha;
            XX = self.XX;
            yy = self.yy;
            Xy = self.Xy;
            % build elements for equation (3.10.14)
            V = diag(v);
            inv_V = diag(1./v);
            inv_V_b = b ./ v;
            inv_V_bar = inv_V + XX;
            V_bar = la.invert_spd_matrix(inv_V_bar);
            b_bar = V_bar * (inv_V_b + Xy);
            alpha_bar= alpha + n;
            delta_bar = delta + yy + b' * inv_V_b - b_bar' * inv_V_bar * b_bar;
            % compute log of marginal likelihood, omitting irrelevant terms
            term_1 = -0.5 * la.stable_determinant(V * XX);
            term_2 = (alpha / 2) * log(delta) - (alpha_bar / 2) * log(delta_bar);
            % take negative (minimize negative to maximize)
            negative_log_f_y = -(term_1 + term_2);
        end   


        function [bayesian_forecast_evaluation_criteria] = bayesian_forecast_evaluation_criteria(self, y, forecast_location, forecast_scale, nu_i, m)
            
            % Bayesian forecast evaluation criteria from equations from equations (3.10.13) and (3.10.16)
    
            log_score = zeros(m,1);
            crps = zeros(m,1);
            for i = 1:m
                % get actual, prediction mean, prediction variance
                y_i = y(i);
                mu_i = forecast_location(i);
                sigma_i = forecast_scale(i,i);
                % get log score from equation (3.10.13)
                [log_pdf, ~] = su.student_pdf(y_i, mu_i, sigma_i, nu_i);
                log_score(i) = - log_pdf;
                % get CRPS from equation (3.10.16)
                s_i = sqrt(sigma_i);
                y_tld = (y_i - mu_i) / s_i;
                [~, pdf] = su.student_pdf(y_tld, 0, 1, nu_i);
                cdf = su.student_cdf(y_tld, 0, 1, nu_i);
                term_1 = y_tld * (2 * cdf - 1);
                term_2 = 2 * pdf * (nu_i + y_tld^2) / (nu_i - 1);
                term_3 = - (2 * sqrt(nu_i) * beta(0.5, nu_i - 0.5)) / ((nu_i - 1) * beta(0.5, nu_i / 2)^2);
                crps(i) = s_i * (term_1 + term_2 + term_3);
            end
            log_score = mean(log_score);
            crps = mean(crps);   
            bayesian_forecast_evaluation_criteria = struct;
            bayesian_forecast_evaluation_criteria.log_score = log_score;
            bayesian_forecast_evaluation_criteria.crps = crps;
        end


    end

    
end
    
        
