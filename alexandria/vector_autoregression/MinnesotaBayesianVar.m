classdef MinnesotaBayesianVar < handle & VectorAutoRegression & BayesianVar


    % Litterman Minnesota vector autoregression, developed in section 11.2
    % 
    % Parameters:
    % -----------
    % endogenous : matrix of size (n_obs,n)
    %     endogenous variables, defined in (4.11.1)
    % 
    % exogenous : matrix of size (n_obs,m), default = []
    %     exogenous variables, defined in (4.11.1)
    %
    % structural_identification : int, default = 2
    %     structural identification scheme, as defined in section 13.2
    %     1 = none, 2 = Cholesky, 3 = triangular, 4 = restrictions
    % 
    % restriction_table : matrix
    %     numerical matrix of restrictions for structural identification
    %
    % lags : int, default = 4
    %     number of lags, defined in (4.11.1)
    %
    % constant : bool, default = true
    %     if true, an intercept is included in the VAR model exogenous
    %
    % trend : bool, default = false
    %     if true, a linear trend is included in the VAR model exogenous
    %
    % quadratic_trend : bool, default = false
    %     if true, a quadratic trend is included in the VAR model exogenous
    %
    % ar_coefficients : float or matrix of size (n_endo,1), default = 0.95
    %     prior mean delta for AR coefficients, defined in (4.11.16)
    % 
    % pi1 : float, default = 0.1
    %     overall tightness hyperparameter, defined in (4.11.17)
    % 
    % pi2 : float, default = 0.5
    %     cross-variable shrinkage hyperparameter, defined in (4.11.18)
    % 
    % pi3 : float, default = 1
    %     lag decay hyperparameter, defined in (4.11.17)    
    % 
    % pi4 : float, default = 100
    %     exogenous slackness hyperparameter, defined in (4.11.19)    
    % 
    % pi5 : float, default = 1
    %     sums-of-coefficients hyperparameter, defined in (4.12.6)     
    % 
    % pi6 : float, default = 0.1
    %     initial observation hyperparameter, defined in (4.12.10)   
    % 
    % pi7 : float, default = 0.1
    %     long-run hyperparameter, defined in (4.12.16)      
    %
    % constrained_coefficients : bool, default = false
    %     if true, applies constrained coefficients, as defined in section 12.1
    %
    % constrained_coefficients_table : matrix, default = []
    %     numerical matrix of constrained coefficients
    %
    % sums_of_coefficients : bool, default = false
    %     if true, applies sums-of-coefficients, as defined in section 12.2
    %
    % dummy_initial_observation : bool, default = false
    %     if true, applies dummy initial observation, as defined in section 12.2    
    %
    % long_run_prior : bool, default = false
    %     if true, applies long-run prior, as defined in section 12.2     
    %
    % long_run_table : matrix, default = []
    %     numerical matrix of long-run prior
    %
    % hyperparameter_optimization : bool, default = false
    %     if true, applies hyperparameter optimization by marginal likelihood 
    %
    % stationary_prior : bool, default = false
    %     if true, applies stationary prior, as defined in section 12.4       
    %
    % credibility_level : float, default = 0.95
    %     VAR model credibility level (between 0 and 1)
    %
    % iterations : int, default = 2000
    %     number of Gibbs sampler replications   
    % 
    % verbose : bool, default = false
    %     if true, displays a progress bar 
    % 
    % 
    % Properties
    % ----------
    % endogenous : matrix of size (n_obs,n)
    %     endogenous variables, defined in (4.11.1)
    % 
    % exogenous : matrix of size (n_obs,m)
    %     exogenous variables, defined in (4.11.1)
    %
    % structural_identification : int
    %     structural identification scheme, as defined in section 13.2
    %     1 = none, 2 = Cholesky, 3 = triangular, 4 = restrictions
    % 
    % restriction_table : matrix
    %     numerical matrix of structural identification restrictions
    % 
    % lags : int
    %     number of lags, defined in (4.11.1)
    %
    % constant : bool
    %     if true, an intercept is included in the VAR model exogenous
    %
    % trend : bool
    %     if true, a linear trend is included in the VAR model exogenous
    %
    % quadratic_trend : bool
    %     if true, a quadratic trend is included in the VAR model exogenous
    %
    % ar_coefficients : float or matrix of size (n_endo,1)
    %     prior mean delta for AR coefficients, defined in (4.11.16)
    % 
    % pi1 : float
    %     overall tightness hyperparameter, defined in (4.11.17)
    % 
    % pi2 : float
    %     cross-variable shrinkage hyperparameter, defined in (4.11.18)
    % 
    % pi3 : float
    %     lag decay hyperparameter, defined in (4.11.17)    
    % 
    % pi4 : float
    %     exogenous slackness hyperparameter, defined in (4.11.19)    
    % 
    % pi5 : float
    %     sums-of-coefficients hyperparameter, defined in (4.12.6)     
    % 
    % pi6 : float
    %     initial observation hyperparameter, defined in (4.12.10)   
    % 
    % pi7 : float
    %     long-run hyperparameter, defined in (4.12.16)  
    %
    % constrained_coefficients : bool
    %     if true, applies constrained coefficients, as defined in section 12.1
    %
    % constrained_coefficients_table : matrix
    %     numerical matrix of constrained coefficients
    %
    % sums_of_coefficients : bool
    %     if true, applies sums-of-coefficients, as defined in section 12.2
    %
    % dummy_initial_observation : bool
    %     if true, applies dummy initial observation, as defined in section 12.2    
    %
    % long_run_prior : bool
    %     if true, applies long-run prior, as defined in section 12.2    
    % 
    % J : matrix of size (n,n)
    %     matrix of long-run prior coefficients, defined in (4.12.15)
    %
    % hyperparameter_optimization : bool
    %     if true, applies hyperparameter optimization by marginal likelihood 
    %
    % stationary_prior : bool
    %     if true, applies stationary prior, as defined in section 12.4       
    %
    % credibility_level : float
    %     VAR model credibility level (between 0 and 1)
    %
    % iterations : int
    %     number of Gibbs sampler replications   
    % 
    % verbose : bool
    %     if true, displays a progress bar      
    % 
    % Sigma : matrix of size (n,n)
    %     variance-covariance matrix of VAR residuals, defined in (4.11.1)
    % 
    % b_bar : matrix of size (q,1)
    %     posterior mean of VAR coefficients, defined in (4.11.15)
    % 
    % V_bar : matrix of size (q,q)
    %     posterior mean of VAR coefficients, defined in (4.11.15)    
    %
    % beta_estimates : matrix of size (k,n,3)
    %     estimates of VAR coefficients
    %     page 1: median, page 2: st dev,  page 3: lower bound, page 4: upper bound
    %
    % Sigma_estimates : matrix of size (n,n)
    %     estimates of variance-covariance matrix of VAR residuals
    % 
    % mcmc_beta : matrix of size (k,n,iterations)
    %     MCMC values of VAR coefficients   
    % 
    % mcmc_Sigma : matrix of size (n,n,iterations)
    %     MCMC values of residual variance-covariance matrix
    % 
    % m_y : float
    %     log 10 marginal likelihood, defined in (4.12.21) 
    % 
    % b : matrix of size (q,1)
    %     prior mean of VAR coefficients, defined in (4.11.16)
    % 
    % V : matrix of size (q,q)
    %     prior mean of VAR coefficients, defined in (4.11.20)           
    %
    % mcmc_H :  matrix of size (n,n,iterations)
    %     MCMC values of structural identification matrix, defined in (4.13.5)
    %
    % mcmc_Gamma : matrix of size (iterations,n)
    %     MCMC values of structural shock variance matrix, defined in definition 13.1
    % 
    % Y : matrix of size (T,n)
    %     matrix of in-sample endogenous variables, defined in (4.11.3)
    % 
    % Z : matrix of size (T,m)
    %     matrix of in-sample endogenous variables, defined in (4.11.3)
    % 
    % X : matrix of size (T,k)
    %     matrix of exogenous and lagged regressors, defined in (4.11.3)
    % 
    % n : int
    %     number of endogenous variables, defined in (4.11.1)
    % 
    % m : int
    %     number of exogenous variables, defined in (4.11.1)
    % 
    % p : int
    %     number of lags, defined in (4.11.1)
    % 
    % T : int
    %     number of sample periods, defined in (4.11.1)
    % 
    % k : int
    %     number of VAR coefficients in each equation, defined in (4.11.1)
    % 
    % q : int
    %     total number of VAR coefficients, defined in (4.11.1)
    %
    % delta : matrix of size (n,1)
    %     prior mean delta for AR coefficients, defined in (4.11.16)
    %
    % s : matrix of size (n,1)
    %     individual AR models residual variance, defined in (4.11.18)
    %
    % Y_sum : matrix of size (n,n)
    %     sums-of-coefficients Y matrix, defined in (4.12.6)
    %
    % X_sum : matrix of size (n,k)
    %     sums-of-coefficients X matrix, defined in (4.12.6)
    %
    % Y_obs : matrix of size (1,n)
    %     dummy initial observation Y matrix, defined in (4.12.10)
    %
    % X_obs : matrix of size (1,k)
    %     dummy initial observation X matrix, defined in (4.12.10)
    %
    % Y_lrp : matrix of size (1,n)
    %     long run prior Y matrix, defined in (4.12.16)
    %
    % X_lrp : matrix of size (1,k)
    %     long run prior X matrix, defined in (4.12.16)
    %
    % Y_d : matrix of size (T_d,n)
    %     full Y matrix combining sample data and dummy observations, defined in (4.11.62)
    %
    % X_d : matrix of size (T_d,k)
    %     full X matrix combining sample data and dummy observations, defined in (4.11.62)
    %
    % T_d : int
    %     total number of observations combining sample data and dummy observations, defined in (4.11.62)
    %
    % steady_state_estimates : matrix of size (T,n,3)
    %     estimates of steady-state, defined in (4.12.30)
    %
    % fitted_estimates : matrix of size (T,n,3)
    %     estimates of in-sample fit, defined in (4.11.2)
    %
    % residual_estimates : matrix of size (T,n,3)
    %     estimates of in-sample residuals, defined in (4.11.2)
    %
    % structural_shocks_estimates : matrix of size (T,n,3)
    %     estimates of in-sample structural shocks, defined in definition 13.1
    %
    % insample_evaluation : struct
    %     in-sample evaluation criteria, defined in (4.13.15)-(4.13.17)
    %
    % mcmc_structural_shocks : matrix of size (T,n,iterations)
    %     MCMC values of structural shocks
    %
    % mcmc_forecasts : matrix of size (f_periods,n,iterations)
    %     MCMC values of forecasts
    %
    % forecast_estimates : matrix of size (f_periods,n,3)
    %     forecast estimates, defined in (4.13.12) and (4.13.13)
    %     page 1: median, page 2: lower bound, page 3: upper bound
    %
    % forecast_evaluation_criteria : struct
    %     forecast evaluation criteria, defined in (4.13.18)-(4.13.21)
    %
    % mcmc_irf : matrix of size (n,n,irf_periods,iterations)
    %     MCMC values of impulse response function, defined in section 13.1
    %
    % mcmc_irf_exo : matrix of size (n,m,irf_periods,iterations)
    %     MCMC values of exogenous impulse response function
    %
    % mcmc_structural_irf : matrix of size (n,n,irf_periods,iterations)
    %     MCMC values of structural impulse response function, defined in section 13.2
    %
    % irf_estimates : matrix of size (n,n,irf_periods,3)
    %     posterior estimates of impulse response function, defined in section 13.1 - 13.2
    %     page 1: median, page 2: lower bound, page 3: upper bound    
    %
    % exo_irf_estimates : matrix of size (n,m,irf_periods,3)
    %     posterior estimates of exogenous impulse response function, if any exogenous variable
    %     page 1: median, page 2: lower bound, page 3: upper bound
    %
    % mcmc_fevd : matrix of size (n,n,fevd_periods,iterations)
    %     MCMC values of forecast error variance decompositions, defined in section 13.4
    %
    % fevd_estimates : matrix of size (n,n,fevd_periods,3)
    %     posterior estimates of forecast error variance decomposition, defined in section 13.4
    %     page 1: median, page 2: lower bound, page 3: upper bound 
    %
    % mcmc_hd : matrix of size (n,n,T,iterations)
    %     MCMC values of historical decompositions, defined in section 13.5
    %
    % hd_estimates : matrix of size (n,n,T,3)
    %     posterior estimates of historical decomposition, defined in section 13.5
    %     page 1: median, page 2: lower bound, page 3: upper bound 
    % 
    % mcmc_conditional_forecasts : matrix of size (f_periods,n,iterations)
    %     MCMC values of conditional forecasts, defined in section 14.1
    %
    % conditional_forecast_estimates : matrix of size (f_periods,n,3)
    %     posterior estimates of conditional forecast, defined in section 14.1
    %     page 1: median, page 2: lower bound, page 3: upper bound
    %
    % mcmc_structural_conditional_forecasts : matrix of size (f_periods,n,iterations)
    %     MCMC values of structural conditional forecasts, defined in section 14.2
    %
    % structural_conditional_forecast_estimates : matrix of size (f_periods,n,3)
    %     structural conditional forecast estimates, defined in section 14.2
    %     page 1: median, page 2: lower bound, page 3: upper bound
    %
    % H_estimates : matrix of size (n,n)
    %     posterior estimates of structural matrix, defined in section 13.2
    %
    % Gamma_estimates : matrix of size (1,n)
    %     estimates of structural shock variance matrix, defined in section 13.2
    %
    %
    % Methods
    % ----------
    % estimate
    % marginal_likelihood
    % insample_fit
    % forecast
    % forecast_evaluation
    % impulse_response_function
    % forecast_error_variance_decomposition
    % historical_decomposition
    % conditional_forecast
    % structural_conditional_forecast


    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
     
    
    properties (GetAccess = public, SetAccess= protected)
        endogenous
        exogenous
        structural_identification
        restriction_table
        lags
        constant
        trend
        quadratic_trend
        ar_coefficients
        pi1
        pi2
        pi3
        pi4
        pi5
        pi6
        pi7
        constrained_coefficients
        constrained_coefficients_table
        sums_of_coefficients
        dummy_initial_observation
        long_run_prior
        J
        hyperparameter_optimization
        stationary_prior
        credibility_level
        iterations
        verbose
        Sigma
        b_bar
        V_bar
        beta_estimates
        Sigma_estimates
        mcmc_beta
        mcmc_Sigma
        m_y
    end
    
    
    properties (GetAccess = public, SetAccess = {?BayesianVar})
        b
        V
        mcmc_H
        mcmc_Gamma
    end
    
    
    properties (GetAccess = {?BayesianVar}, SetAccess = private)
        mcmc_chol_Sigma
    end
    
    
    properties (GetAccess = private, SetAccess = private)
        inv_Sigma
        inv_V
        inv_V_b
        inv_V_bar
    end
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods (Access = public)
        
        
        function self = MinnesotaBayesianVar(endogenous, varargin)
            
            % constructor for the MinnesotaBayesianVar class
            
            % allow for optional arguments
            parser = inputParser;
            default_exogenous = [];
            default_structural_identification = 2;
            default_restriction_table = [];
            default_lags = 4;
            default_constant = true;
            default_trend = false;
            default_quadratic_trend = false;
            default_ar_coefficients = 0.95;
            default_pi1 = 0.1;
            default_pi2 = 0.5;
            default_pi3 = 1;
            default_pi4 = 100;
            default_pi5 = 1;
            default_pi6 = 0.1;
            default_pi7 = 0.1;
            default_constrained_coefficients = false;
            default_constrained_coefficients_table = [];
            default_sums_of_coefficients = false;
            default_dummy_initial_observation = false;
            default_long_run_prior = false;
            default_long_run_table = [];
            default_hyperparameter_optimization = false;
            default_stationary_prior = false;
            default_credibility_level = 0.95;
            default_iterations = 2000;
            default_verbose = false;
            addRequired(parser, 'endogenous');
            addParameter(parser, 'exogenous', default_exogenous);
            addParameter(parser, 'structural_identification', default_structural_identification);
            addParameter(parser, 'restriction_table', default_restriction_table);
            addParameter(parser, 'lags', default_lags);
            addParameter(parser, 'constant', default_constant);
            addParameter(parser, 'trend', default_trend);
            addParameter(parser, 'quadratic_trend', default_quadratic_trend);  
            addParameter(parser, 'ar_coefficients', default_ar_coefficients);  
            addParameter(parser, 'pi1', default_pi1);   
            addParameter(parser, 'pi2', default_pi2);  
            addParameter(parser, 'pi3', default_pi3);  
            addParameter(parser, 'pi4', default_pi4);
            addParameter(parser, 'pi5', default_pi5);
            addParameter(parser, 'pi6', default_pi6);
            addParameter(parser, 'pi7', default_pi7);
            addParameter(parser, 'constrained_coefficients', default_constrained_coefficients);
            addParameter(parser, 'constrained_coefficients_table', default_constrained_coefficients_table);
            addParameter(parser, 'sums_of_coefficients', default_sums_of_coefficients);
            addParameter(parser, 'dummy_initial_observation', default_dummy_initial_observation);
            addParameter(parser, 'long_run_prior', default_long_run_prior);
            addParameter(parser, 'long_run_table', default_long_run_table);
            addParameter(parser, 'hyperparameter_optimization', default_hyperparameter_optimization);
            addParameter(parser, 'stationary_prior', default_stationary_prior);
            addParameter(parser, 'credibility_level', default_credibility_level);
            addParameter(parser, 'iterations', default_iterations);
            addParameter(parser, 'verbose', default_verbose);
            parse(parser, endogenous, varargin{:});
            self.endogenous = endogenous;
            self.exogenous = parser.Results.exogenous;
            self.structural_identification = parser.Results.structural_identification;
            self.restriction_table = parser.Results.restriction_table;
            self.lags = parser.Results.lags;
            self.constant = parser.Results.constant;
            self.trend = parser.Results.trend;
            self.quadratic_trend = parser.Results.quadratic_trend;
            self.ar_coefficients = parser.Results.ar_coefficients;
            self.pi1 = parser.Results.pi1;
            self.pi2 = parser.Results.pi2;
            self.pi3 = parser.Results.pi3;
            self.pi4 = parser.Results.pi4;
            self.pi5 = parser.Results.pi5;
            self.pi6 = parser.Results.pi6;
            self.pi7 = parser.Results.pi7;
            self.constrained_coefficients = parser.Results.constrained_coefficients;
            self.constrained_coefficients_table = parser.Results.constrained_coefficients_table;
            self.sums_of_coefficients = parser.Results.sums_of_coefficients;
            self.dummy_initial_observation = parser.Results.dummy_initial_observation;
            self.long_run_prior = parser.Results.long_run_prior;
            self.J = parser.Results.long_run_table;
            self.hyperparameter_optimization = parser.Results.hyperparameter_optimization;
            self.stationary_prior = parser.Results.stationary_prior;
            self.credibility_level = parser.Results.credibility_level;
            self.iterations = parser.Results.iterations;
            self.verbose = parser.Results.verbose;
            % make regressors
            self.make_regressors();
            % make delta
            self.make_delta();
            % make individual residual variance
            self.individual_ar_variances();
        end
        
        
        function estimate(self)
            
            % estimate()
            % generates posterior estimates for Bayesian VAR model parameters beta and Sigma
            %
            % parameters:
            % none
            %
            % returns:
            % none

            % obtain Sigma from OLS estimates
            self.get_Sigma();
            % generate dummy extensions, if applicable
            self.dummy_extensions();
            % optimize hyperparameters, if applicable
            self.optimize_hyperparameters(); 
            % define prior values
            self.prior();
            % apply constrained coefficients, if applicable
            self.make_constrained_coefficients();
            % define posterior values
            self.posterior();
            % obtain posterior estimates for regression parameters
            self.parameter_estimates();
            % run MCMC algorithm (Gibbs sampling) for VAR parameters
            self.parameter_mcmc();
            % estimate structural identification
            self.make_structural_identification();
        end  
        
        
        function [m_y] = marginal_likelihood(self)
            
            % marginal_likelihood()
            % log10 marginal likelihood, defined in (4.12.21)
            %
            % parameters:
            % none
            %
            % returns:
            % m_y: float
            %     log10 marginal likelihood
            
            % unpack
            n = self.n;
            T = self.T;
            Sigma = self.Sigma;
            inv_Sigma = self.inv_Sigma;
            V = self.V;
            XX = self.XX;
            YY = self.YY;
            b = self.b;
            inv_V_b = self.inv_V_b;
            b_bar = self.b_bar;
            inv_V_bar = self.inv_V_bar;    
            % evaluate the log marginal likelihood from equation (4.12.21)
            term_1 = - (n*T/2) * log(2 * pi);
            term_2 = - T / 2 * log(det(Sigma));
            term_3 = - 0.5 * la.stable_determinant(V .* kron(inv_Sigma,XX));
            term_4 = - 0.5 * (b' * inv_V_b - b_bar' * inv_V_bar * b_bar + trace(YY * inv_Sigma));
            log_f_y = term_1 + term_2 + term_3 + term_4;      
            % convert to log10
            m_y = log_f_y / log(10);
            % save as attributes
            self.m_y = m_y;         
        end

 
    end
    
    
    methods (Access = protected, Hidden = true)
        
        
        function get_Sigma(self)
            
            % generates Sigma defined in (4.11.9)
            
            var = MaximumLikelihoodVar(self.endogenous, 'lags', self.lags);
            var.estimate();
            Sigma = var.Sigma;
            inv_Sigma = la.invert_spd_matrix(Sigma);
            self.Sigma = Sigma;
            self.inv_Sigma = inv_Sigma;
        end
        
        
        function optimize_hyperparameters(self)
        
            % optimize_hyperparameters delta, pi1, pi2, pi3 and pi4

            if self.hyperparameter_optimization
                % initial value of optimizer
                x_0 = [0.1 0.5 1 100 repmat(0.9, [1,self.n])]';
                % bounds for parameter values
                lower_bound = [0.01 0.01 1 10 zeros(1,self.n)]';
                upper_bound = [1 1 5 1000 ones(1,self.n)]';
                % optimize
                [solution, fval, exitflag] = minimize(@self.negative_likelihood, ...
                               [x_0], [], [], [], [], [lower_bound], [upper_bound]);
                % update hyperparameters
                self.pi1 = solution(1);
                self.pi2 = solution(2);
                self.pi3 = solution(3);
                self.pi4 = solution(4);
                self.delta = solution(5:end);
                % if verbose, display progress bar and success/failure of optimization
                if self.verbose
                    cu.progress_bar_complete('hyperparameter optimization:');
                    cu.optimization_completion(exitflag == 1);
                end
            end
        end
        
        
        function [negative_log_f_y] = negative_likelihood(self, x)
        
            % negative log marginal likelihood for Minnesota prior

            % unpack
            pi1 = x(1);
            pi2 = x(2);
            pi3 = x(3);
            pi4 = x(4);
            delta = x(5:end);
            % recover prior and posterior elements
            b = vu.make_b(delta, self.n, self.m, self.p);
            V = vu.make_V(self.s, pi1, pi2, pi3, pi4, self.n, self.m, self.p);       
            [inv_V inv_V_b] = vu.make_V_b_inverse(b, V);
            [b_bar V_bar inv_V_bar] = vu.minnesota_posterior(inv_V, ...
                                      inv_V_b, self.XX, self.XY, self.inv_Sigma);
            % compute log of marginal likelihood, omitting irrelevant terms
            term_1 = la.stable_determinant(V .* kron(self.inv_Sigma, self.XX));
            term_2 = (b ./ V)' * b - b_bar' * inv_V_bar * b_bar;
            negative_log_f_y = term_1 + term_2;
        end
        
        
        function prior(self)

            % creates prior elements b and V defined in (4.11.16)-(4.11.20)

            % define b and V
            b = vu.make_b(self.delta, self.n, self.m, self.p);
            V = vu.make_V(self.s, self.pi1, self.pi2, self.pi3, self.pi4, self.n, self.m, self.p);
            self.b = b;
            self.V = V;
        end
        
        
        function posterior(self)
            
            % creates posterior elements b_bar and V_bar defined in (4.11.15)
            
            % unpack
            b = self.b;
            V = self.V;
            XX = self.XX_d;
            XY = self.XY_d;
            inv_Sigma = self.inv_Sigma;
            % generate preliminary posterior elements
            [inv_V inv_V_b] = vu.make_V_b_inverse(b, V);
            % define posterior elements
            [b_bar V_bar inv_V_bar] = vu.minnesota_posterior(inv_V, inv_V_b, XX, XY, inv_Sigma);
            self.b_bar = b_bar;
            self.V_bar = V_bar;
            self.inv_V = inv_V;
            self.inv_V_b = inv_V_b;
            self.inv_V_bar = inv_V_bar;            
        end

        
        function parameter_estimates(self)
            
            % point estimates and credibility intervals for model parameters
            % use the multivariate normal distribution defined in (4.11.14)

            % unpack
            k = self.k;
            n = self.n;
            % mean and standard deviation of posterior distribution
            mean = reshape(self.b_bar, [k n]);
            standard_deviation = reshape(diag(self.V_bar).^0.5, [k,n]);
            % critical value of normal distribution for credibility level
            credibility_level = self.credibility_level;
            Z = su.normal_icdf((1 + credibility_level) / 2);
            % initiate storage: 4 columns: lower bound, median, upper bound, standard deviation
            beta_estimates = zeros(k,n,4);
            % fill estimates
            beta_estimates(:,:,1) = mean;
            beta_estimates(:,:,2) = mean - Z * standard_deviation;
            beta_estimates(:,:,3) = mean + Z * standard_deviation;
            beta_estimates(:,:,4) = standard_deviation;          
            % save as attributes
            self.beta_estimates = beta_estimates;
            self.Sigma_estimates = self.Sigma;
        end        
        
        
        function parameter_mcmc(self)

            % Gibbs sampler for VAR parameters beta and Sigma

            % unpack
            b_bar = self.b_bar;
            V_bar = self.V_bar;
            q = self.q;
            k = self.k;
            n = self.n;
            chol_Sigma = la.cholesky_nspd(self.Sigma);
            chol_V_bar = la.cholesky_nspd(V_bar);
            iterations = self.iterations;
            % if stationary prior is not activated, draw all values at once
            if ~self.stationary_prior
                mcmc_beta = b_bar + chol_V_bar * randn(q,iterations);
                % if verbose, display progress bar
                if self.verbose
                    cu.progress_bar_complete('Model parameters:');
                end
            % if stationary prior is activated, get values one by one and check for stationarity
            elseif self.stationary_prior
                m = self.m;
                p = self.lags;
                iteration = 1;
                mcmc_beta = zeros(q,iterations);
                while iteration <= iterations
                    beta = b_bar + chol_V_bar * randn(q,1);
                    B = reshape(beta, [k n]);
                    stationary = vu.check_stationarity(B, n, m, p);
                    if stationary
                        mcmc_beta(:,iteration) = beta;
                        % if verbose, display progress bar
                        if self.verbose
                            cu.progress_bar(iteration, self.iterations, 'Model parameters:');
                        end
                        iteration = iteration + 1;
                    end   
                end
            end
            mcmc_beta = reshape(mcmc_beta,[k n iterations]);
            self.mcmc_beta = mcmc_beta;
            self.mcmc_Sigma = repmat(self.Sigma,[1 1 iterations]); 
            self.mcmc_chol_Sigma = repmat(chol_Sigma,[1 1 iterations]); 
        end
        

    end
    

end

    