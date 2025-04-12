classdef DummyObservationBayesianVar < handle & VectorAutoRegression & BayesianVar
    

    % Dummy observation vector autoregression, developed in section 11.5
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
    % Y_dum : matrix of size (n,n)
    %     dummy observation prior Y matrix, defined in (4.12.6)
    %
    % X_dum : matrix of size (n,k)
    %     dummy observation prior X matrix, defined in (4.12.6)
    %
    % B_hat : matrix of size (k,n)
    %     posterior mean of VAR coefficients, defined in (4.11.50)
    % 
    % W_hat : matrix of size (k,k)
    %     posterior mean of VAR coefficients, defined in (4.11.50) 
    % 
    % alpha_hat : float
    %     posterior degrees of freedom, defined in (4.11.50)
    % 
    % S_hat : matrix of size (n,n)
    %     posterior scale matrix, defined in (4.11.50) 
    %     
    % alpha_tilde : float
    %     posterior degrees of freedom, defined in (4.11.53)
    % 
    % S_tilde : matrix of size (n,n)
    %     posterior scale matrix, defined in (4.11.53) 
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
    % insample_fit
    % forecast
    % forecast_evaluation
    % impulse_response_function
    % forecast_error_variance_decomposition
    % historical_decomposition
    % conditional_forecast


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
        pi3
        pi4
        pi5
        pi6
        pi7
        sums_of_coefficients
        dummy_initial_observation
        long_run_prior
        J
        stationary_prior
        credibility_level
        iterations
        verbose
        Y_dum
        X_dum
        B_hat
        W_hat
        alpha_hat
        S_hat
        alpha_tilde
        S_tilde
        beta_estimates
        Sigma_estimates
        mcmc_beta
        mcmc_Sigma
    end
    
    
    properties (GetAccess = public, SetAccess = {?BayesianVar})
        mcmc_H
        mcmc_Gamma
    end
    
    
    properties (GetAccess = {?BayesianVar}, SetAccess = private)
        mcmc_chol_Sigma
    end
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods (Access = public)
        
        
        function self = DummyObservationBayesianVar(endogenous, varargin)
            
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
            default_pi3 = 1;
            default_pi4 = 100;
            default_pi5 = 1;
            default_pi6 = 0.1;
            default_pi7 = 0.1;
            default_sums_of_coefficients = false;
            default_dummy_initial_observation = false;
            default_long_run_prior = false;
            default_long_run_table = [];
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
            addParameter(parser, 'pi3', default_pi3);  
            addParameter(parser, 'pi4', default_pi4);
            addParameter(parser, 'pi5', default_pi5);
            addParameter(parser, 'pi6', default_pi6);
            addParameter(parser, 'pi7', default_pi7);
            addParameter(parser, 'sums_of_coefficients', default_sums_of_coefficients);
            addParameter(parser, 'dummy_initial_observation', default_dummy_initial_observation);
            addParameter(parser, 'long_run_prior', default_long_run_prior);
            addParameter(parser, 'long_run_table', default_long_run_table);
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
            self.pi3 = parser.Results.pi3;
            self.pi4 = parser.Results.pi4;
            self.pi5 = parser.Results.pi5;
            self.pi6 = parser.Results.pi6;
            self.pi7 = parser.Results.pi7;
            self.sums_of_coefficients = parser.Results.sums_of_coefficients;
            self.dummy_initial_observation = parser.Results.dummy_initial_observation;
            self.long_run_prior = parser.Results.long_run_prior;
            self.J = parser.Results.long_run_table;
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

            % generate dummy extensions, if applicable
            self.dummy_extensions();
            % define prior values
            self.prior();
            % define posterior values
            self.posterior();
            % obtain posterior estimates for regression parameters
            self.parameter_estimates();
            % run MCMC algorithm (Gibbs sampling) for VAR parameters
            self.parameter_mcmc();
            % estimate structural identification
            self.make_structural_identification();
        end  
        

    end
    
    
    methods (Access = protected, Hidden = true)
        
        
        function prior(self)

            % creates prior elements Y_dum and X_dum defined in (4.11.54)

            % matrices H, J and K
            H = diag(self.delta .* self.s .^ 0.5);
            J = diag((1:self.lags) .^ self.pi3);
            K = diag(self.s .^ 0.5);
            % create the different blocks and concatenate
            block_1 = zeros(self.m,self.n);
            block_2 = H / self.pi1;
            block_3 = zeros(self.n*(self.lags-1),self.n);
            block_4 = K;
            Y_dum = [block_1;block_2;block_3;block_4];
            block_5 = eye(self.m) / (self.pi1 * self.pi4);
            block_6 = zeros(self.m,self.n*self.lags);
            block_7 = zeros(self.n*self.p,self.m);
            block_8 = kron(J,K) / self.pi1;
            block_9 = zeros(self.n,self.m);
            block_10 = zeros(self.n,self.n*self.lags);
            X_dum = [block_5 block_6; block_7 block_8; block_9 block_10];
            self.Y_dum = Y_dum;
            self.X_dum = X_dum;            
            % update data with dummy observations, following (4.11.62)
            self.Y_d = [Y_dum;self.Y_d];
            self.X_d = [X_dum;self.X_d];
            self.T_d = self.T_d + size(Y_dum, 1);
            self.XX_d = self.X_d' * self.X_d;
            self.XY_d = self.X_d' * self.Y_d;
            self.YY_d = self.Y_d' * self.Y_d;             
        end
        
        
        function posterior(self)
            
            % creates posterior elements B_hat, W_hat, alpha_hat, S_hat, alpha_tilde and S_tilde
            % defined in (4.11.50) and (4.11.53)
            
            % unpack
            Y = self.Y_d;
            X = self.X_d;
            XX = self.XX_d;
            XY = self.XY_d;
            T = self.T_d;
            n = self.n;
            k = self.k;
            % define posterior elements
            inv_XX = la.invert_spd_matrix(XX);
            W_hat = inv_XX;
            B_hat = inv_XX * XY;
            alpha_hat = T - k + 2;
            E = Y - X * B_hat;
            S_hat = E' * E;
            alpha_tilde = T - n - k + 3;
            S_tilde = S_hat / alpha_tilde;         
            % save as attributes
            self.B_hat = B_hat;
            self.W_hat = W_hat;
            self.alpha_hat = alpha_hat;
            self.S_hat = S_hat;
            self.alpha_tilde = alpha_tilde;
            self.S_tilde = S_tilde;
        end

        
        function parameter_estimates(self)
            
            % point estimates and credibility intervals for model parameters
            % use the inverse Wishart defined in (4.11.51) and matrix Student defined in (4.11.53)

            % mean and standard deviation of posterior distribution, using matrix Student properties
            mean = self.B_hat;
            scale = reshape(diag(kron(self.S_tilde, self.W_hat)), [self.k self.n]);
            standard_deviation = (self.alpha_tilde / (self.alpha_tilde - 2) * scale).^ 0.5;
            % critical value of Student distribution for credibility level
            credibility_level = self.credibility_level;
            Z = su.student_icdf((1 + credibility_level) / 2, self.alpha_tilde);            
            % initiate storage: 4 columns: lower bound, median, upper bound, standard deviation
            beta_estimates = zeros(self.k,self.n,4);            
            % fill estimates
            beta_estimates(:,:,1) = mean;
            beta_estimates(:,:,2) = mean - Z * scale.^ 0.5;
            beta_estimates(:,:,3) = mean + Z * scale.^ 0.5;
            beta_estimates(:,:,4) = standard_deviation;   
            % posterior estimates for Sigma, using mean of inverse Wishart
            Sigma_estimates = self.S_hat / (self.alpha_hat - self.n - 1);        
            % save as attributes
            self.beta_estimates = beta_estimates;
            self.Sigma_estimates = Sigma_estimates;
        end        
        
        
        function parameter_mcmc(self)

            % Gibbs sampler for VAR parameters beta and Sigma

            % unpack
            B_hat = self.B_hat;
            W_hat = self.W_hat;
            alpha_hat = self.alpha_hat;
            S_hat = self.S_hat;
            iterations = self.iterations;
            n = self.n;
            m = self.m;
            p = self.p;
            k = self.k;
            % initiate storage
            mcmc_beta = zeros(k,n,iterations);
            mcmc_Sigma = zeros(n,n,iterations);
            mcmc_chol_Sigma = zeros(n,n,iterations);
            % compute constant Cholesky factors
            chol_W_hat = la.cholesky_nspd(W_hat);
            chol_S_hat = la.cholesky_nspd(S_hat);
            % loop over iterations, checking for stationarity if stationary prior is activated
            iteration = 1;
            while iteration <= iterations
                Sigma = rn.inverse_wishart(alpha_hat, chol_S_hat, true);
                chol_Sigma = la.cholesky_nspd(Sigma);
                % sample B from matrix normal distribution defined in (4.12.34)
                B = rn.matrix_normal(B_hat, chol_W_hat, Sigma, true, false);
                if self.stationary_prior
                    stationary = vu.check_stationarity(B, n, m, p);
                else
                    stationary = true;
                end
                if stationary
                    mcmc_beta(:,:,iteration) = B;
                    mcmc_Sigma(:,:,iteration) = Sigma;
                    mcmc_chol_Sigma(:,:,iteration) = chol_Sigma;
                    % if verbose, display progress bar
                    if self.verbose
                        cu.progress_bar(iteration, iterations, 'Model parameters:');
                    end
                    iteration = iteration + 1;
                end
            end
            self.mcmc_beta = mcmc_beta;
            self.mcmc_Sigma = mcmc_Sigma;
            self.mcmc_chol_Sigma = mcmc_chol_Sigma;
        end
        

    end
    

end
    
    