classdef VectorAutoregressiveMovingAverage < handle


    % Vector Autoregressive Moving Average, developed in chapter 16
    % 
    % Parameters:
    % -----------
    % endogenous : matrix of size (n_obs,n)
    %     endogenous variables, defined in (5.16.1)
    % 
    % exogenous : matrix of size (n_obs,m), default = []
    %     exogenous variables, defined in (5.16.1)
    %
    % structural_identification : int, default = 2
    %     structural identification scheme, as defined in section 16.5
    %     1 = none, 2 = Cholesky, 3 = triangular, 4 = restrictions
    % 
    % restriction_table : matrix
    %     numerical matrix of restrictions for structural identification
    %
    % lags : int, default = 4
    %     number of lags, defined in (5.16.1)
    %
    % residual_lags : int, default = 1
    %     number of lags in past residuals, defined in (5.16.1)
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
    %     prior mean delta for AR coefficients, defined in (5.16.9)
    % 
    % pi1 : float, default = 0.1
    %     overall tightness hyperparameter, defined in (5.16.9)
    % 
    % pi2 : float, default = 0.5
    %     cross-variable shrinkage hyperparameter, defined in (5.16.9)
    % 
    % pi3 : float, default = 1
    %     lag decay hyperparameter, defined in (5.16.9)    
    % 
    % pi4 : float, default = 100
    %     exogenous slackness hyperparameter, defined in (5.16.9)      
    % 
    % lambda1 : float, default = 0.1
    %     overall tightness hyperparameter, defined in (5.16.10)
    % 
    % lambda2 : float, default = 0.5
    %     cross-variable shrinkage hyperparameter, defined in (5.16.10)
    % 
    % lambda3 : float, default = 1
    %     lag decay hyperparameter, defined in (5.16.10) 
    %
    % credibility_level : float, default = 0.95
    %     VAR model credibility level (between 0 and 1)
    %
    % iterations : int, default = 2000
    %     number of Gibbs sampler replications   
    %
    % burnin : int, default = 1000
    %     number of Gibbs sampler burn-in replications  
    % 
    % verbose : bool, default = false
    %     if true, displays a progress bar 
    % 
    % 
    % Properties
    % ----------
    % endogenous : matrix of size (n_obs,n)
    %     endogenous variables, defined in (5.16.1)
    % 
    % exogenous : matrix of size (n_obs,m), default = []
    %     exogenous variables, defined in (5.16.1)
    %
    % structural_identification : int, default = 2
    %     structural identification scheme, as defined in section 16.5
    %     1 = none, 2 = Cholesky, 3 = triangular, 4 = restrictions
    % 
    % restriction_table : matrix
    %     numerical matrix of restrictions for structural identification
    %
    % lags : int, default = 4
    %     number of lags, defined in (5.16.1)
    %
    % residual_lags : int, default = 1
    %     number of residual lags, defined in (5.16.1)
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
    %     prior mean delta for AR coefficients, defined in (5.16.9)
    % 
    % pi1 : float, default = 0.1
    %     overall tightness hyperparameter, defined in (5.16.9)
    % 
    % pi2 : float, default = 0.5
    %     cross-variable shrinkage hyperparameter, defined in (5.16.9)
    % 
    % pi3 : float, default = 1
    %     lag decay hyperparameter, defined in (5.16.9)    
    % 
    % pi4 : float, default = 100
    %     exogenous slackness hyperparameter, defined in (5.16.9)    
    % 
    % lambda1 : float, default = 0.1
    %     overall tightness hyperparameter, defined in (5.16.10)
    % 
    % lambda2 : float, default = 0.5
    %     cross-variable shrinkage hyperparameter, defined in (5.16.10)
    % 
    % lambda3 : float, default = 1
    %     lag decay hyperparameter, defined in (5.16.10) 
    %
    % credibility_level : float, default = 0.95
    %     VAR model credibility level (between 0 and 1)
    %
    % iterations : int, default = 2000
    %     number of Gibbs sampler replications   
    %
    % burnin : int, default = 1000
    %     number of Gibbs sampler burn-in replications  
    % 
    % verbose : bool, default = false
    %     if true, displays a progress bar 
    %
    % Y : matrix of size (T,n)
    %     matrix of in-sample endogenous variables, defined in (5.16.3)
    % 
    % X : matrix of size (T,k1)
    %     matrix of exogenous and lagged regressors, defined in (5.16.3)
    % 
    % n : int
    %     number of endogenous variables, defined in (5.16.1)
    % 
    % m : int
    %     number of exogenous variables, defined in (5.16.1)
    % 
    % p : int
    %     number of lags, defined in (5.16.1)
    % 
    % q : int
    %     number of residual lags, defined in (5.16.1)
    %     
    % T : int
    %     number of sample periods, defined in (5.16.1)
    % 
    % k1 : int
    %     number of autoregressive coefficients in each equation, defined in (5.16.1)
    % 
    % k2 : int
    %     number of lagged residual coefficients in each equation, defined in (5.16.1)
    % 
    % k : int
    %     total number of coefficients in each equation, defined in (5.16.1)
    % 
    % r1 : int
    %     total number of autoregressive coefficients in the model, defined in (5.16.1)
    % 
    % r2 : int
    %     totla number of lagged residual coefficients in the model, defined in (5.16.1)
    % 
    % r : int
    %     total number of coefficients in the model, defined in (5.16.1)
    % 
    % delta : matrix of size (n,1)
    %     prior mean delta for AR coefficients, defined in (5.16.9)
    %
    % s : matrix of size (n,1)
    %     prior scale matrix, defined in (5.16.11) 
    %
    % b : matrix of size (r1,1)
    %     prior mean of VAR coefficients, defined in (5.16.9)
    % 
    % V : matrix of size (q,q)
    %     prior variance of autoregressive coefficients, defined in (5.16.9)           
    %
    % W : matrix of size (q,q)
    %     prior variance of lagged residual coefficients, defined in (5.16.10)           
    %
    % alpha : float
    %     prior degrees of freedom, defined in (5.16.11)
    % 
    % S : matrix of size (n,1)
    %     prior scale matrix, defined in (5.16.11) 
    %
    % alpha_bar : float
    %     posterior degrees of freedom, defined in (5.16.23)
    %
    % mcmc_beta : matrix of size (k1,n,iterations)
    %     MCMC values of autoregressive coefficients   
    % 
    % mcmc_kappa : matrix of size (k2,n,iterations)
    %     MCMC values of autoregressive coefficients   
    %    
    % mcmc_Z : matrix of size (T,k2,iterations)
    %     MCMC values of lagged residuals  
    %       
    % mcmc_E : matrix of size (T,n,iterations)
    %     MCMC values of current residuals  
    %          
    % mcmc_Sigma : matrix of size (n,n,iterations)
    %     MCMC values of residual variance-covariance matrix
    %  
    % mcmc_chol_Sigma : matrix of size (n,n,iterations)
    %     MCMC values of residual variance-covariance matrix (Cholesky factor)
    % 
    % mcmc_inv_Sigma : matrix of size (n,n,iterations)
    %     MCMC values of residual variance-covariance matrix (inverse)
    % 
    % beta_estimates : matrix of size (k1,n,3)
    %     estimates of VAR coefficients
    %     page 1: median, page 2: st dev,  page 3: lower bound, page 4: upper bound
    %
    % kappa_estimates : matrix of size (k2,n,3)
    %     estimates of lagged residual coefficients
    %     page 1: median, page 2: st dev,  page 3: lower bound, page 4: upper bound
    %  
    % E_estimates : matrix of size (T,n)
    %     estimates of current residuals
    %
    % Z_estimates : matrix of size (T,k2)
    %     estimates of lagged residuals
    %    
    % Sigma_estimates : matrix of size (n,n)
    %     estimates of variance-covariance matrix of VAR residuals
    % 
    % mcmc_H :  matrix of size (n,n,iterations)
    %     MCMC values of structural identification matrix, defined in (4.13.5)
    %
    % mcmc_Gamma : matrix of size (iterations,n)
    %     MCMC values of structural shock variance matrix, defined in definition 13.1
    % 
    % H_estimates : matrix of size (n,n)
    %     posterior estimates of structural matrix, defined in section 13.2
    %
    % Gamma_estimates : matrix of size (1,n)
    %     estimates of structural shock variance matrix, defined in section 13.2
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
    % mcmc_structural_shocks : matrix of size (T,n,iterations)
    %     MCMC values of structural shocks
    %
    % structural_shocks_estimates : matrix of size (T,n,3)
    %     estimates of in-sample structural shocks, defined in definition 13.1
    %
    % insample_evaluation : struct
    %     in-sample evaluation criteria, defined in (4.13.15)-(4.13.17)
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
     
    
    properties (GetAccess = public, SetAccess = protected)
        endogenous
        exogenous
        structural_identification
        restriction_table
        lags
        residual_lags
        constant
        trend
        quadratic_trend
        ar_coefficients
        pi1
        pi2
        pi3
        pi4    
        lambda1
        lambda2
        lambda3
        credibility_level
        iterations
        burnin
        verbose
        Y
        X
        n
        m
        p
        q
        T
        k1
        k2
        k
        r1
        r2
        r
        delta
        s
        b
        V
        W
        alpha
        S
        alpha_bar
        mcmc_beta
        mcmc_kappa
        mcmc_Z
        mcmc_E           
        mcmc_Sigma
        mcmc_chol_Sigma
        mcmc_inv_Sigma
        beta_estimates
        kappa_estimates
        E_estimates
        Z_estimates
        Sigma_estimates
        mcmc_H
        mcmc_Gamma
        H_estimates
        Gamma_estimates 
        steady_state_estimates
        fitted_estimates
        residual_estimates
        mcmc_structural_shocks
        structural_shock_estimates
        insample_evaluation
        mcmc_forecast
        forecast_estimates
        forecast_evaluation_criteria
        mcmc_irf
        mcmc_irf_exo
        mcmc_structural_irf
        irf_estimates
        exo_irf_estimates
        mcmc_fevd
        fevd_estimates
        mcmc_hd
        hd_estimates
        mcmc_conditional_forecast
        conditional_forecast_estimates
    end

    
    properties (GetAccess = private, SetAccess = private)
        XX
        inv_V
        inv_V_b
        inv_W
        F
        M
        mcmc_inv_H
        svar_index
    end

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods (Access = public)
        
        
        function self = VectorAutoregressiveMovingAverage(endogenous, varargin)
            
            % constructor for the VectorAutoregressiveMovingAverage class
            
            % allow for optional arguments
            parser = inputParser;
            default_exogenous = [];
            default_structural_identification = 2;
            default_restriction_table = [];
            default_lags = 4;
            default_residual_lags = 1;
            default_constant = true;
            default_trend = false;
            default_quadratic_trend = false;
            default_ar_coefficients = 0.95;
            default_pi1 = 0.1;
            default_pi2 = 0.5;
            default_pi3 = 1;
            default_pi4 = 100;
            default_lambda1 = 0.1;
            default_lambda2 = 0.5;
            default_lambda3 = 1;
            default_credibility_level = 0.95;
            default_iterations = 3000;
            default_burnin = 1000;
            default_verbose = false;
            addRequired(parser, 'endogenous');
            addParameter(parser, 'exogenous', default_exogenous);
            addParameter(parser, 'structural_identification', default_structural_identification);
            addParameter(parser, 'restriction_table', default_restriction_table);
            addParameter(parser, 'lags', default_lags);
            addParameter(parser, 'residual_lags', default_residual_lags);
            addParameter(parser, 'constant', default_constant);
            addParameter(parser, 'trend', default_trend);
            addParameter(parser, 'quadratic_trend', default_quadratic_trend);  
            addParameter(parser, 'ar_coefficients', default_ar_coefficients);
            addParameter(parser, 'pi1', default_pi1);   
            addParameter(parser, 'pi2', default_pi2);  
            addParameter(parser, 'pi3', default_pi3);  
            addParameter(parser, 'pi4', default_pi4);
            addParameter(parser, 'lambda1', default_lambda1);   
            addParameter(parser, 'lambda2', default_lambda2);  
            addParameter(parser, 'lambda3', default_lambda3);             
            addParameter(parser, 'credibility_level', default_credibility_level);
            addParameter(parser, 'iterations', default_iterations);
            addParameter(parser, 'burnin', default_burnin);
            addParameter(parser, 'verbose', default_verbose);
            parse(parser, endogenous, varargin{:});
            self.endogenous = endogenous;
            self.exogenous = parser.Results.exogenous;
            self.structural_identification = parser.Results.structural_identification;
            self.restriction_table = parser.Results.restriction_table;
            self.lags = parser.Results.lags;
            self.residual_lags = parser.Results.residual_lags;
            self.constant = parser.Results.constant;
            self.trend = parser.Results.trend;
            self.quadratic_trend = parser.Results.quadratic_trend;
            self.ar_coefficients = parser.Results.ar_coefficients;
            self.pi1 = parser.Results.pi1;
            self.pi2 = parser.Results.pi2;
            self.pi3 = parser.Results.pi3;
            self.pi4 = parser.Results.pi4;
            self.lambda1 = parser.Results.lambda1;
            self.lambda2 = parser.Results.lambda2;
            self.lambda3 = parser.Results.lambda3;            
            self.credibility_level = parser.Results.credibility_level;
            self.iterations = parser.Results.iterations;
            self.burnin = parser.Results.burnin;
            self.verbose = parser.Results.verbose;
            % make regressors
            self.make_regressors();
        end


        function estimate(self)
        
            % estimate()
            % generates posterior estimates for Bayesian VARMA model parameters
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
            % run MCMC algorithm for VARMA parameters
            self.parameter_mcmc();
            % obtain posterior estimates
            self.parameter_estimates();
            % estimate structural identification
            self.make_structural_identification();
        end


        function insample_fit(self)
            
            % insample_fit(self)
            % generates in-sample fit and residuals along with evaluation criteria
            %
            % parameters:
            % none
            %
            % returns:
            % none
            
            % compute steady-state
            self.steady_state();
            % compute fitted and residuals
            self.fitted_and_residual();
            % compute in-sample criteria
            self.insample_criteria();
        end   


        function [forecast_estimates] = forecast(self, h, credibility_level, varargin)
            
            % [forecast_estimates] = forecast(self, h, credibility_level, varargin)
            % estimates forecasts for the Bayesian VAR model, using algorithm 13.4
            % 
            % parameters:
            % h : int
            %     number of forecast periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % varargin : either zero or one additional argument, matrix of dimension (h, n_exo)
            %     no argument unless the model includes exogenous other than constant, trend and quadratic trend
            %     if argument, n_exo is the number of additional exogenous variables
            % 
            % returns:
            % forecast_estimates : matrix of dimension (h,n,3)
            %     page 1: median; page 2: interval lower bound; page 3: interval upper bound            
            
            % optional argument
            if isempty(varargin)
                Z_p = [];
            else
                Z_p = varargin{1};
            end  
            % get forecast
            self.make_forecast(h, Z_p);
            % obtain posterior estimates
            self.forecast_posterior_estimates(credibility_level);
            forecast_estimates = self.forecast_estimates;        
        end 


        function [forecast_evaluation_criteria] = forecast_evaluation(self, Y)
            
            % [forecast_evaluation_criteria] = forecast_evaluation(self, Y)
            % forecast evaluation criteria for the Bayesian VARMA model, as defined in (4.13.18)-(4.13.22)
            % 
            % parameters:
            % Y : matrix of size(h,n)
            %     matrix of realised values for forecast evaluation, h being the number of forecast periods
            % 
            % returns:
            % forecast_evaluation_criteria : structure
            %     structure with criteria name as keys and corresponding number as value       
            
            % unpack
            Y_hat = self.forecast_estimates(:,:,1);
            mcmc_forecast = self.mcmc_forecast;
            % obtain regular forecast evaluation criteria from equations (4.13.18) and (4.13.19)
            standard_evaluation_criteria = vu.forecast_evaluation_criteria(Y_hat, Y);
            % obtain Bayesian forecast evaluation criteria from equations (4.13.21) and (4.13.22)
            bayesian_evaluation_criteria = vu.bayesian_forecast_evaluation_criteria(mcmc_forecast, Y);
            % merge structures
            forecast_evaluation_criteria = iu.concatenate_structures(standard_evaluation_criteria, bayesian_evaluation_criteria);
            % save as attributes
            self.forecast_evaluation_criteria = forecast_evaluation_criteria;     
        end   


        function [irf_estimates exo_irf_estimates] = impulse_response_function(self, h, credibility_level)

            % [irf_estimates exo_irf_estimates] = impulse_response_function(h, credibility_level)
            % impulse response functions, as defined in (4.13.1)-(4.13.9)
            % 
            % parameters:
            % h : int
            %     number of IRF periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % 
            % returns:
            % irf_estimates : matrix of size (n,n,h,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound
            % exo_irf_estimates : matrix of size (n,n,h,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound

            % get regular impulse response funtion
            self.make_impulse_response_function(h);
            % get exogenous impuse response function
            self.make_exogenous_impulse_response_function(h);
            % get structural impulse response function
            self.make_structural_impulse_response_function(h);
            % obtain posterior estimates
            self.irf_posterior_estimates(credibility_level);
            irf_estimates = self.irf_estimates;
            exo_irf_estimates = self.exo_irf_estimates;
        end


        function [fevd_estimates] = forecast_error_variance_decomposition(self, h, credibility_level)

            % [fevd_estimates] = forecast_error_variance_decomposition(self, h, credibility_level)
            % forecast error variance decomposition, as defined in (4.13.31)
            % 
            % parameters:
            % h : int
            %     number of FEVD periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % 
            % returns:
            % fevd_estimates : matrix of size (n,n,h,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound

            % get forecast error variance decomposition
            self.make_forecast_error_variance_decomposition(h);
            % obtain posterior estimates
            self.fevd_posterior_estimates(credibility_level);
            fevd_estimates = self.fevd_estimates;
        end


        function [hd_estimates] = historical_decomposition(self, credibility_level)

            % [hd_estimates] = historical_decomposition(self, credibility_level)
            % historical decomposition, as defined in (4.13.34)-(4.13-36)
            % 
            % parameters:
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % 
            % returns:
            % hd_estimates : matrix of size (n,n,T,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound

            % get historical decomposition
            self.make_historical_decomposition();
            % obtain posterior estimates
            self.hd_posterior_estimates(credibility_level);
            hd_estimates = self.hd_estimates;
        end 


        function [conditional_forecast_estimates] = conditional_forecast(self, h, credibility_level, conditions, shocks, conditional_forecast_type, varargin)

            % conditional_forecast(h, credibility_level, conditions, shocks, conditional_forecast_type, varargin)
            % estimates conditional forecasts for the Bayesian VAR model, using algorithms 14.1 and 14.2
            % 
            % parameters:
            % h : int
            %     number of forecast periods           
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % conditions: matrix of size (n_conditions,4)
            %     table defining conditions (column 1: variable, column 2: period, column 3: mean, column 4: variance)   
            % shocks: empty matrix or matrix of size (n,1)
            %     vector defining shocks generating the conditions; should be empty if conditional_forecast_type = 1  
            % conditional_forecast_type : int
            %     conditional forecast type (1 = agnostic, 2 = structural)
            % varargin : either zero or one additional argument, matrix of dimension (h, n_exo)
            %     no argument unless the model includes exogenous other than constant, trend and quadratic trend
            %     if argument, n_exo is the number of additional exogenous variables            
            % 
            % returns:
            % conditional_forecast_estimates : matrix of size (h,n,3)
            %     page 1: median; page 2: interval lower bound; page 3: interval upper bound

            % optional argument
            if isempty(varargin)
                Z_p = [];
            else
                Z_p = varargin{1};
            end 
            % if conditional forecast type is agnostic
            if conditional_forecast_type == 1
                % get conditional forecasts
                self.make_conditional_forecast(h, conditions, Z_p);
            % if instead conditional forecast type is structural
            elseif conditional_forecast_type == 2
                % establish type of shocks
                shock_type = self.check_shock_type(h, conditions, shocks);
                % get structural conditional forecasts
                self.make_structural_conditional_forecast(h, conditions, shocks, Z_p, shock_type);
            end
            % obtain posterior estimates
            self.conditional_forecast_posterior_estimates(credibility_level);
            conditional_forecast_estimates = self.conditional_forecast_estimates;
        end  


    end


    methods (Access = protected, Hidden = true)


        function make_regressors(self)
        
            % generates regressors defined in (5.16.3), along with other dimension elements
                  
            % define regressor matrices
            [Y X] = self.make_regressor_matrices();
            % define dimensions
            [n m p q T k1 k2 k r1 r2 r] = self.generate_dimensions();
            % define estimation terms
            XX = X' * X;
            % save as attributes      
            self.Y = Y;
            self.X = X;
            self.n = n;
            self.m = m;
            self.p = p;
            self.q = q;       
            self.T = T;
            self.k1 = k1;
            self.k2 = k2;
            self.k = k;
            self.r1 = r1;
            self.r2 = r2;
            self.r = r;     
            self.XX = XX;
        end


        function [Y X] = make_regressor_matrices(self)
            
            Y = self.endogenous(self.lags+1:end,:);      
            periods = size(self.endogenous,1) - self.lags;
            X_1 = vu.generate_intercept_and_trends(self.constant, self.trend, self.quadratic_trend, periods, 0);
            X_2 = vu.generate_exogenous_regressors(self.exogenous, self.lags);
            X_3 = vu.generate_lagged_endogenous(self.endogenous, self.lags);
            X = [X_1 X_2 X_3];
        end


        function [n m p q T k1 k2 k r1 r2 r] = generate_dimensions(self)
        
            T = size(self.endogenous,1) - self.lags;
            n = size(self.endogenous,2);
            m = double(self.constant) + double(self.trend) + double(self.quadratic_trend);  
            if ~isempty(self.exogenous)
                m = m + size(self.exogenous,2);
            end
            p = self.lags;
            q = self.residual_lags;
            k1 = m + n * p;
            k2 = n * q;
            k = k1 + k2;
            r1 = n * k1;
            r2 = n * k2;
            r = r1 + r2;
        end


        function prior(self)
            
            % creates prior elements b, V, W, alpha, S and F defined in (5.16.9), (5.16.10), (5.16.11) and (5.16.13)
            
            delta = self.make_delta();
            s = self.individual_ar_variances();
            b = vu.make_b(delta, self.n, self.m, self.p);
            V = vu.make_V(s, self.pi1, self.pi2, self.pi3, self.pi4, self.n, self.m, self.p);
            W = vu.make_V(s, self.lambda1, self.lambda2, self.lambda3, 1, self.n, 0, self.q);
            alpha = vu.make_alpha(self.n);
            S = vu.make_S(s);
            F = diag(ones(self.n*self.q,1),-self.n);
            M = diag([1 zeros(1,self.q)]);
            self.delta = delta;
            self.s = s;
            self.b = b;
            self.V = V;
            self.W = W;
            self.alpha = alpha;
            self.S = S;
            self.F = F;
            self.M = M;
        end


        function [delta] = make_delta(self)
            
            if isscalar(self.ar_coefficients)
                ar_coefficients = repmat(self.ar_coefficients, [self.n 1]);
            else
                ar_coefficients = self.ar_coefficients;
            end
            delta = ar_coefficients;
        end
        
        
        function [s] = individual_ar_variances(self)
            
            s = zeros(self.n,1);
            for i=1:self.n
                ar = MaximumLikelihoodVar(self.endogenous(:,i), 'lags', self.lags);
                ar.estimate();
                s(i) = ar.Sigma;
            end
        end


        function posterior(self)
            
            % creates posterior elements defined in (5.16.17), (5.16.20) and (5.16.23)
            
            [inv_V inv_V_b] = vu.make_V_b_inverse(self.b, self.V);
            inv_W = diag(1 ./ self.W);
            alpha_bar = self.alpha + self.T;
            self.inv_V = inv_V;
            self.inv_V_b = inv_V_b;
            self.inv_W = inv_W;
            self.alpha_bar = alpha_bar;
        end


        function parameter_mcmc(self)
            
            % Gibbs sampler for VAR parameters beta, kappa, Sigma and Z, following algorithm 16.1
            
            % unpack
            Y = self.Y;
            X = self.X;
            inv_V = self.inv_V;
            inv_V_b = self.inv_V_b;
            inv_W = self.inv_W;
            F = self.F;
            M = self.M;
            XX = self.XX;
            alpha_bar = self.alpha_bar;
            S = self.S;
            n = self.n;
            q = self.q;
            k1 = self.k1;
            k2 = self.k2;
            T = self.T;
            iterations = self.iterations;
            burnin = self.burnin;
            verbose = self.verbose;
            
            % preallocate storage space
            mcmc_beta = zeros(k1,n,iterations);
            mcmc_kappa = zeros(k2,n,iterations);
            mcmc_Z = zeros(T,k2,iterations);
            mcmc_E = zeros(T,n,iterations);
            mcmc_Sigma = zeros(n,n,iterations);
            mcmc_chol_Sigma = zeros(n,n,iterations);
            mcmc_inv_Sigma = zeros(n,n,iterations);
    
            % set initial values for parameters
            inv_Sigma = diag(1 ./ S);
            Sigma = diag(S);
            S = diag(S);
            Z = zeros(T,k2);
            K = zeros(k2,n);
            
            % iterate over iterations
            iteration = 1;
            while iteration <= (burnin + iterations)
            
                % step 2: sample beta
                beta = self.draw_beta(Y, X, inv_V, inv_V_b, inv_Sigma, XX, Z, K);
                B = reshape(beta,[k1 n]);

                % step 3: sample kappa
                kappa = self.draw_kappa(Y, X, Z, B, inv_W, inv_Sigma);
                K = reshape(kappa,[k2 n]);

                % step 4: sample Z
                [E Z] = self.draw_Z(Y, X, B, K, F, M, Sigma, n, T, q, k2);

                % step 5: sample Sigma
                [Sigma inv_Sigma chol_Sigma] = self.draw_Sigma(E, S, alpha_bar);

                % save if burn is exceeded
                if iteration > burnin
                    mcmc_beta(:,:,iteration-burnin) = B;
                    mcmc_kappa(:,:,iteration-burnin) = K;
                    mcmc_Z(:,:,iteration-burnin) = Z;
                    mcmc_E(:,:,iteration-burnin) = E;
                    mcmc_Sigma(:,:,iteration-burnin) = Sigma;
                    mcmc_chol_Sigma(:,:,iteration-burnin) = chol_Sigma;
                    mcmc_inv_Sigma(:,:,iteration-burnin) = inv_Sigma;
                end                
                if verbose
                    cu.progress_bar(iteration, iterations+burnin, 'Model parameters:');
                end
                iteration = iteration + 1;
            end

            % save as attributes
            self.mcmc_beta = mcmc_beta;
            self.mcmc_kappa = mcmc_kappa;
            self.mcmc_Z = mcmc_Z;
            self.mcmc_E = mcmc_E;           
            self.mcmc_Sigma = mcmc_Sigma;
            self.mcmc_chol_Sigma = mcmc_chol_Sigma;
            self.mcmc_inv_Sigma = mcmc_inv_Sigma;
        end


        function [beta] = draw_beta(self, Y, X, inv_V, inv_V_b, inv_Sigma, XX, Z, K)

            % draw beta from its conditional posterior defined in (5.16.17)

            % posterior V_bar
            inv_V_bar = inv_V + kron(inv_Sigma, XX);
            % posterior b_bar
            b_bar_temp = inv_V_b + la.vec(X' * (Y - Z * K) * inv_Sigma);
            % efficient sampling of beta (algorithm 9.4)
            beta = rn.efficient_multivariate_normal(b_bar_temp, inv_V_bar);
        end  


        function [kappa] = draw_kappa(self, Y, X, Z, B, inv_W, inv_Sigma)
            
            % draw kappa from its conditional posterior defined in (5.16.20)
            
            % posterior W_bar
            inv_W_bar = inv_W + kron(inv_Sigma, Z' * Z);
            % posterior g_bar
            g_bar_temp = la.vec(Z' * (Y - X * B) * inv_Sigma);  
            % efficient sampling of beta (algorithm 9.4)
            kappa = rn.efficient_multivariate_normal(g_bar_temp, inv_W_bar);
        end

        function [E Z] = draw_Z(self, Y, X, B, K, F, M, Sigma, n, T, q, k2)
        
            % draw Z from its conditional posterior, using (5.16.12)-(5.16.13)
            
            % parameters for state-space representation
            [G Y_GX J Omega] = self.state_space_representation(Y, X, B, K, M, Sigma, n);
            % initial values for algorithm
            [z_00 Upsilon_00] = self.kalman_filter_initial_values(Sigma, k2, n, q);
            % Carter-Kohn algorithm: forward pass
            [Z_tt Z_tt1 Ups_tt Ups_tt1] = self.forward_pass(Y_GX, J, F, Omega, z_00, Upsilon_00, T, n, k2);
            % Carter-Kohn algorithm: backward pass
            [Z] = self.backward_pass(Z_tt, Z_tt1, Ups_tt, Ups_tt1, F, T, n, k2);
            E = Z(:,1:n);
            Z = Z(:,n+1:end);
        end


        function [G Y_GX J Omega] = state_space_representation(self, Y, X, B, K, M, Sigma, n)
            
            G = B';
            Y_GX = Y - X * B;
            J = [eye(n) K']; 
            Omega = kron(M, Sigma);
        end


        function [z_00 Upsilon_00] = kalman_filter_initial_values(self, Sigma, k2, n, q)
            
            z_00 = zeros(k2+n,1);
            Upsilon_00 = kron(eye(q+1), Sigma);
        end

    
        function [Z_tt Z_tt1 Ups_tt Ups_tt1] = forward_pass(self, Y_GX, J, F, Omega, z_00, Upsilon_00, T, n, k2)
            
            [Z_tt Z_tt1 Ups_tt Ups_tt1] = ss.varma_forward_pass(Y_GX, J, F, Omega, z_00, Upsilon_00, T, n, k2+n);
        end


        function [Z] = backward_pass(self, Z_tt, Z_tt1, Ups_tt, Ups_tt1, F, T, n, k2)
            
            [Z] = ss.varma_backward_pass(Z_tt, Z_tt1, Ups_tt, Ups_tt1, F, T, k2+n);
        end


        function [Sigma inv_Sigma chol_Sigma] = draw_Sigma(self, E, S, alpha_bar)

            % draw Sigma from its conditional posterior defined in (5.16.23)

            % posterior S_bar
            S_bar = S + E' * E;
            % sample sigma
            Sigma = rn.inverse_wishart(alpha_bar, S_bar);
            % obtain related elements
            inv_Sigma = la.invert_spd_matrix(Sigma);
            chol_Sigma = la.cholesky_nspd(Sigma);        
        end  


        function parameter_estimates(self)
            
            % point estimates and credibility intervals for model parameters
            % use empirical quantiles from MCMC algorithm

            % posterior estimates for beta
            beta_estimates = zeros(self.k1,self.n,4);
            beta_estimates(:,:,1:3) = vu.posterior_estimates(self.mcmc_beta, self.credibility_level);
            beta_estimates(:,:,4) = std(self.mcmc_beta,0,3);
            % posterior estimates for K
            kappa_estimates = zeros(self.k2,self.n,4);
            kappa_estimates(:,:,1:3) = vu.posterior_estimates(self.mcmc_kappa, self.credibility_level);
            kappa_estimates(:,:,4) = std(self.mcmc_kappa,0,3);
            % posterior estimates for E, Z and Sigma
            E_estimates = quantile(self.mcmc_E,0.5,3);
            Z_estimates = quantile(self.mcmc_Z,0.5,3);      
            Sigma_estimates = quantile(self.mcmc_Sigma,0.5,3);
            self.beta_estimates = beta_estimates;
            self.kappa_estimates = kappa_estimates;
            self.E_estimates = E_estimates;
            self.Z_estimates = Z_estimates;
            self.Sigma_estimates = Sigma_estimates;
        end


        function make_structural_identification(self)

            % structural identification estimates
            
            if self.structural_identification == 2
                self.svar_by_choleski_factorization();
            elseif self.structural_identification == 3
                self.svar_by_triangular_factorization();
            elseif self.structural_identification == 4
                self.svar_by_restrictions();
            end
            if self.structural_identification ~= 1
                self.svar_estimates();  
            end
        end


        function svar_by_choleski_factorization(self)
            
            self.mcmc_H = self.mcmc_chol_Sigma;
            self.mcmc_Gamma = ones(self.iterations,self.n);
            mcmc_inv_H = zeros(self.n,self.n,self.iterations);
            for j=1:self.iterations
                mcmc_inv_H(:,:,j) = la.invert_lower_triangular_matrix(self.mcmc_H(:,:,j));
                if self.verbose
                    cu.progress_bar(j, self.iterations, 'Structural identification:');
                end
            end
            self.mcmc_inv_H = mcmc_inv_H;
            self.svar_index = (1:self.iterations)';
        end


        function svar_by_triangular_factorization(self)
            
            mcmc_H = zeros(self.n,self.n,self.iterations);
            mcmc_inv_H = zeros(self.n,self.n,self.iterations);
            mcmc_Gamma = zeros(self.iterations,self.n);
            for j=1:self.iterations
                [H Gamma] = la.triangular_factorization(self.mcmc_chol_Sigma(:,:,j), true);
                mcmc_H(:,:,j) = H;
                mcmc_inv_H(:,:,j) = la.invert_lower_triangular_matrix(H);
                mcmc_Gamma(j,:) = Gamma;
                if self.verbose
                    cu.progress_bar(j, self.iterations, 'Structural identification:');
                end
            end
            self.mcmc_H = mcmc_H;
            self.mcmc_Gamma = mcmc_Gamma;
            self.mcmc_inv_H = mcmc_inv_H;
            self.svar_index = (1:self.iterations)';
        end


        function svar_by_restrictions(self)
            
            % initiate MCMC elements
            svar_index = zeros(self.iterations,1);
            mcmc_H = zeros(self.n,self.n,self.iterations);
            mcmc_inv_H = zeros(self.n,self.n,self.iterations);
            % create matrices of restriction checks
            [restriction_matrices max_irf_period] = vu.make_restriction_matrices(self.restriction_table, self.p);
            % make preliminary orthogonalised impulse response functions, if relevant
            [mcmc_irf] = vvu.make_varma_restriction_irf(self.mcmc_beta, self.mcmc_kappa, self.mcmc_chol_Sigma, ...
                                                 self.iterations, self.n, self.p, self.q, max_irf_period);  
            % make preliminary structural shocks, if relevant
            [mcmc_shocks] = vvu.make_varma_restriction_shocks(self.mcmc_E, self.mcmc_chol_Sigma, self.T, self.n, self.iterations, restriction_matrices);  
            % loop over iterations, until desired number of total iterations is obtained
            i = 1;
            while i <= self.iterations
                % select a random index in number of iterations
                j = rn.discrete_uniform(1, self.iterations);
                % make a random rotation matrix Q: if no zero restrictions, draw from uniform orthogonal distribution
                if isempty(restriction_matrices{1,1})
                    Q = rn.uniform_orthogonal(self.n);
                % if there are zero restrictions, use the zero uniform orthogonal distribution
                else
                    Q = rn.zero_uniform_orthogonal(self.n, restriction_matrices{1,1}, mcmc_irf(:,:,:,j));
                end
                % check restrictions: IRF, sign
                irf_sign_index = restriction_matrices{2,1};
                if ~isempty(irf_sign_index)
                    irf_sign_coefficients = restriction_matrices{2,2};
                    restriction_satisfied = vu.check_irf_sign(irf_sign_index, irf_sign_coefficients, mcmc_irf(:,:,:,j), Q);
                    if ~restriction_satisfied
                        continue
                    end
                end
                % check restrictions: IRF, magnitude
                irf_magnitude_index = restriction_matrices{3,1};
                if ~isempty(irf_magnitude_index)
                    irf_magnitude_coefficients = restriction_matrices{3,2};
                    restriction_satisfied = vu.check_irf_magnitude(irf_magnitude_index, irf_magnitude_coefficients, mcmc_irf(:,:,:,j), Q);
                    if ~restriction_satisfied
                        continue
                    end
                end
                % check restrictions: structural shocks, sign
                shock_sign_index = restriction_matrices{4,1};
                if ~isempty(shock_sign_index)
                    shock_sign_coefficients = restriction_matrices{4,2};
                    restriction_satisfied = vu.check_shock_sign(shock_sign_index, shock_sign_coefficients, mcmc_shocks(:,:,j), Q);
                    if ~restriction_satisfied
                        continue
                    end
                end
                % check restrictions: structural shocks, magnitude
                shock_magnitude_index = restriction_matrices{5,1};
                if ~isempty(shock_magnitude_index)
                    shock_magnitude_coefficients = restriction_matrices{5,2};
                    restriction_satisfied = vu.check_shock_magnitude(shock_magnitude_index, shock_magnitude_coefficients, mcmc_shocks(:,:,j), Q);
                    if ~restriction_satisfied
                        continue                 
                    end
                end
                % historical decomposition values if any of sign or magnitude restrictions apply
                history_sign_index = restriction_matrices{6,1};
                history_magnitude_index = restriction_matrices{7,1};
                if ~isempty(history_sign_index) || ~isempty(history_magnitude_index)
                    [irf shocks] = vu.make_restriction_irf_and_shocks(mcmc_irf(:,:,:,j), mcmc_shocks(:,:,j), Q, self.n);             
                end
                % check restrictions: historical decomposition, sign
                if ~isempty(history_sign_index)
                    history_sign_coefficients = restriction_matrices{6,2};
                    restriction_satisfied = vu.check_history_sign(history_sign_index, history_sign_coefficients, irf, shocks);
                    if ~restriction_satisfied 
                        continue                
                    end
                end
                % check restrictions: historical decomposition, magnitude
                if ~isempty(history_magnitude_index)
                    history_magnitude_coefficients = restriction_matrices{7,2};
                    restriction_satisfied = vu.check_history_magnitude(history_magnitude_index, history_magnitude_coefficients, irf, shocks);
                    if ~restriction_satisfied
                        continue    
                    end
                end
                % if all restriction passed, keep the draw and record
                H = self.mcmc_chol_Sigma(:,:,j) * Q;
                inv_H = Q' * la.invert_lower_triangular_matrix(self.mcmc_chol_Sigma(:,:,j));
                mcmc_H(:,:,i) = H;
                mcmc_inv_H(:,:,i) = inv_H;
                svar_index(i) = j;                
                if self.verbose
                    cu.progress_bar(i, self.iterations, 'Structural identification:');
                end
                i = i + 1;
            self.mcmc_H = mcmc_H;
            self.mcmc_inv_H = mcmc_inv_H;
            self.mcmc_Gamma = ones(self.iterations,self.n);
            self.svar_index = svar_index; 
            end
        end

        
        function svar_estimates(self)
            
            H_estimates = quantile(self.mcmc_H,0.5,3);
            Gamma_estimates = quantile(self.mcmc_Gamma,0.5,1);
            self.H_estimates = H_estimates;
            self.Gamma_estimates = Gamma_estimates;        
        end  


        function steady_state(self)

            ss = zeros(self.T,self.n,self.iterations);
            for j=1:self.iterations
                ss(:,:,j) = vvu.varma_steady_state(self.X, self.mcmc_beta(:,:,j), self.n, self.m, self.p, self.T);
                if self.verbose
                    cu.progress_bar(j, self.iterations, 'Steady-state:');
                end
            end
            ss_estimates = vu.posterior_estimates(ss, self.credibility_level);
            self.steady_state_estimates = ss_estimates;
        end


        function fitted_and_residual(self)
        
            fitted = zeros(self.T,self.n,self.iterations);
            residual = zeros(self.T,self.n,self.iterations);
            for j=1:self.iterations
                [residual(:,:,j) fitted(:,:,j)] = vvu.varma_fit_and_residuals(self.X, self.mcmc_beta(:,:,j), ...
                                                  self.mcmc_Z(:,:,j), self.mcmc_kappa(:,:,j), self.mcmc_E(:,:,j));
                if self.verbose
                    cu.progress_bar(j, self.iterations, 'Fitted and residual:');
                end
            end
            fitted_estimates = vu.posterior_estimates(fitted, self.credibility_level);
            residual_estimates = vu.posterior_estimates(residual, self.credibility_level); 
            self.fitted_estimates = fitted_estimates;
            self.residual_estimates = residual_estimates;
            if self.structural_identification ~= 1
                structural_shocks = zeros(self.T,self.n,self.iterations);
                for j=1:self.iterations
                    index = self.svar_index(j);
                    structural_shocks(:,:,j) = vu.structural_shocks(residual(:,:,index), self.mcmc_inv_H(:,:,j));
                    if self.verbose
                        cu.progress_bar(j, self.iterations, 'Structural shocks:');
                    end
                end
                structural_shock_estimates = vu.posterior_estimates(structural_shocks, self.credibility_level);
                self.mcmc_structural_shocks = structural_shocks;
                self.structural_shock_estimates = structural_shock_estimates;
            end
        end


        function insample_criteria(self)
        
            insample_evaluation = vu.insample_evaluation_criteria(self.Y, ...
                                  self.residual_estimates(:,:,1), self.T, self.k);
            if self.verbose
                cu.progress_bar_complete('In-sample evaluation criteria:');
            end
            self.insample_evaluation = insample_evaluation;
        end


        function make_forecast(self, h, Z_p)         
        
            [Z_p Y] = vu.make_forecast_regressors(Z_p, self.Y, h, self.p, self.T,...
                       self.exogenous, self.constant, self.trend, self.quadratic_trend);
            mcmc_forecast = zeros(h,self.n,self.iterations);
            for j=1:self.iterations
                mcmc_forecast(:,:,j) = vvu.varma_forecast(self.mcmc_beta(:,:,j), self.mcmc_kappa(:,:,j), ...
                             self.mcmc_chol_Sigma(:,:,j), h, Z_p, Y, self.mcmc_E(end-self.q+1:end,:,j), self.n);
                if self.verbose
                    cu.progress_bar(j, self.iterations, 'Forecasts:');
                end
            end
            self.mcmc_forecast = mcmc_forecast;
        end


        function forecast_posterior_estimates(self, credibility_level)           
            
            % obtain posterior estimates
            mcmc_forecast = self.mcmc_forecast;
            forecast_estimates = vu.posterior_estimates(mcmc_forecast, credibility_level);
            self.forecast_estimates = forecast_estimates;
        end


        function make_impulse_response_function(self, h)        

            mcmc_irf = zeros(self.n, self.n, h, self.iterations);
            for i=1:self.iterations
                % get regular impulse response function for VARMA
                mcmc_irf(:,:,:,i) = vvu.varma_impulse_response_function(self.mcmc_beta(:,:,i), ...
                                    self.mcmc_kappa(:,:,i), self.n, self.p, self.q, h);
                if self.verbose   
                    cu.progress_bar(i, self.iterations, 'Impulse response function:');
                end
            end
            self.mcmc_irf = mcmc_irf;
        end
        

        function make_exogenous_impulse_response_function(self, h)

            if ~isempty(self.exogenous)
                r = size(self.exogenous,2);
                mcmc_irf_exo = zeros(self.n, r, h, self.iterations);
                for i=1:self.iterations
                    % get exogenous IRFs: same as VAR since no shocks are involved in exogenous IRF
                    mcmc_irf_exo(:,:,:,i) = vu.exogenous_impulse_response_function(self.mcmc_beta(:,:,i), self.n, self.m, r, self.p, h);
                    if self.verbose
                        cu.progress_bar(i, self.iterations, 'Exogenous impulse response function:');
                    end
                end
            else
                mcmc_irf_exo = [];
            end
            self.mcmc_irf_exo = mcmc_irf_exo;
        end

        
        function make_structural_impulse_response_function(self, h)

            if self.structural_identification == 1
                if self.verbose
                    cu.progress_bar_complete('Structural impulse response function:');
                end
                self.mcmc_structural_irf = [];                
            else
                mcmc_structural_irf = self.mcmc_irf;
                for i=1:self.iterations
                    index = self.svar_index(i);
                    mcmc_structural_irf(:,:,:,i) = vu.structural_impulse_response_function(...
                                                   self.mcmc_irf(:,:,:,index),self.mcmc_H(:,:,i), self.n);
                    if self.verbose   
                        cu.progress_bar(i, self.iterations, 'Structural impulse response function:');
                    end
                end
                self.mcmc_structural_irf = mcmc_structural_irf;
            end
        end
        
        
        function irf_posterior_estimates(self, credibility_level)

            if self.structural_identification == 1
                mcmc_irf = self.mcmc_irf;
            else
                mcmc_irf = self.mcmc_structural_irf;
            end
            irf_estimates = vu.posterior_estimates_3d(mcmc_irf, credibility_level);
            if ~isempty(self.exogenous)
                exo_irf_estimates = vu.posterior_estimates_3d(self.mcmc_irf_exo, credibility_level);
            else
                exo_irf_estimates = [];
            end
            self.irf_estimates = irf_estimates;
            self.exo_irf_estimates = exo_irf_estimates;
        end


        function make_forecast_error_variance_decomposition(self, h)

            if self.structural_identification == 1
                if self.verbose
                    cu.progress_bar_complete('Forecast error variance decomposition:');
                end
                self.mcmc_fevd = [];                
            else
                if self.structural_identification == 3
                    mcmc_Gamma = num2cell(self.mcmc_Gamma,2);
                else
                    mcmc_Gamma = cell(self.iterations,1);
                end
                mcmc_fevd = zeros(self.n, self.n, h, self.iterations);                
                has_irf = ~isempty(self.mcmc_structural_irf) && size(self.mcmc_structural_irf, 3) >= h;
                for i=1:self.iterations
                    if has_irf
                        structural_irf = self.mcmc_structural_irf(:,:,1:h,i);
                    else
                        index = self.svar_index(i);
                        irf = vvu.varma_impulse_response_function(self.mcmc_beta(:,:,i), ...
                                    self.mcmc_kappa(:,:,i), self.n, self.p, self.q, h);
                        structural_irf = vu.structural_impulse_response_function(irf, self.mcmc_H(:,:,i), self.n);
                    end
                    mcmc_fevd(:,:,:,i) = vu.forecast_error_variance_decomposition(structural_irf, mcmc_Gamma{i}, self.n, h); 
                    if self.verbose   
                        cu.progress_bar(i, self.iterations, 'Forecast error variance decomposition:');
                    end
                end
                self.mcmc_fevd = mcmc_fevd;
            end
        end
        
        
        function fevd_posterior_estimates(self, credibility_level)

            if self.structural_identification == 1
                self.fevd_estimates = [];
            else
                mcmc_fevd = self.mcmc_fevd;
                fevd_estimates = vu.posterior_estimates_3d(mcmc_fevd, credibility_level);
                normalized_fevd_estimates = vu.normalize_fevd_estimates(fevd_estimates);
                self.fevd_estimates = normalized_fevd_estimates;
            end
        end


        function make_historical_decomposition(self)

            if self.structural_identification == 1
                if self.verbose
                    cu.progress_bar_complete('Historical decomposition:');
                end
                self.mcmc_hd = [];     
            else
                mcmc_hd = zeros(self.n, self.n, self.T, self.iterations); 
                has_irf = ~isempty(self.mcmc_structural_irf) && size(self.mcmc_structural_irf, 3) >= self.T;
                has_structural_shocks = ~isempty(self.mcmc_structural_shocks);
                for i=1:self.iterations
                    index = self.svar_index(i);
                    if has_irf
                        structural_irf = self.mcmc_structural_irf(:,:,1:self.T,i);
                    else
                        irf = vvu.varma_impulse_response_function(self.mcmc_beta(:,:,index), ...
                                    self.mcmc_kappa(:,:,index), self.n, self.p, self.q, self.T);
                        structural_irf = vu.structural_impulse_response_function(irf, self.mcmc_H(:,:,i), self.n);
                    end
                    if has_structural_shocks
                        structural_shocks = self.mcmc_structural_shocks(:,:,i); 
                    else
                        [E ~] = vvu.varma_fit_and_residuals(self.X, self.mcmc_beta(:,:,i), ...
                                self.mcmc_Z(:,:,i), self.mcmc_kappa(:,:,i), self.mcmc_E(:,:,index));
                        structural_shocks = vu.structural_shocks(E, self.mcmc_inv_H(:,:,i));
                    end
                    mcmc_hd(:,:,:,i) = vu.historical_decomposition(structural_irf, structural_shocks, self.n, self.T);          
                    if self.verbose
                        cu.progress_bar(i, self.iterations, 'Historical decomposition:');
                    end
                end
                self.mcmc_hd = mcmc_hd;
            end
        end        
        
        
        function hd_posterior_estimates(self, credibility_level)

            if self.structural_identification == 1
                self.hd_estimates = [];
            else
                mcmc_hd = self.mcmc_hd;
                hd_estimates = vu.posterior_estimates_3d(mcmc_hd, credibility_level);
                self.hd_estimates = hd_estimates;
            end
        end  


        function make_conditional_forecast(self, h, conditions, Z_p)        

            % make regressors Z_p and Y
            [Z_p Y] = vu.make_forecast_regressors(Z_p, self.Y, h, self.p, self.T,...
                      self.exogenous, self.constant, self.trend, self.quadratic_trend);
            % make conditional forecast regressors y_bar, Z and omega
            [y_bar Q omega] = vvu.varma_conditional_forecast_regressors_1(conditions, h, self.n, self.p, self.q);
            % initiate storage and loop over iterations
            mcmc_conditional_forecast = zeros(h,self.n,self.iterations);
            for i=1:self.iterations
                % recover iteration-specific regressors
                [mu F K gamma_00 Upsilon_00] = vvu.varma_conditional_forecast_regressors_2(Y, self.mcmc_E(end-self.q+1:end,:,i), ...
                                               self.mcmc_beta(:,:,i), self.mcmc_kappa(:,:,i), self.mcmc_Sigma(:,:,i), ...
                                               conditions, Z_p, self.n, self.m, self.p, self.q, h);
                % run Carter Kohn algorithm to obtain conditional forecasts
                bss = BayesianStateSpaceSampler(y_bar, Q, omega, mu, F, K, 'z_00', gamma_00, ...
                      'Upsilon_00', Upsilon_00, 'kalman_type', 'conditional_forecast');
                bss.carter_kohn_algorithm();
                mcmc_conditional_forecast(:,:,i) = bss.Z(:,1:self.n);
                if self.verbose
                    cu.progress_bar(i, self.iterations, 'Conditional forecasts:');
                end
            self.mcmc_conditional_forecast = mcmc_conditional_forecast;
            end
        end


        function conditional_forecast_posterior_estimates(self, credibility_level)

            if isempty(self.mcmc_conditional_forecast)
                self.conditional_forecast_estimates = [];
            else
                mcmc_conditional_forecast = self.mcmc_conditional_forecast;
                conditional_forecast_estimates = vu.posterior_estimates(mcmc_conditional_forecast, credibility_level);
                self.conditional_forecast_estimates = conditional_forecast_estimates; 
            end
        end   


        function [shock_type] = check_shock_type(self, h, conditions, shocks)

            % check for structural identification
            if self.structural_identification == 1 
                if self.verbose
                    cu.progress_bar_complete('Conditional forecasts:');
                end
                shock_type = 'none';
            else
                % identify shocks
                if sum(shocks) == self.n
                    shock_type = 'all_shocks';
                else
                    shock_type = 'shock-specific';  
                end
            end
        end


        function make_structural_conditional_forecast(self, h, conditions, shocks, Z_p, shock_type)        

            % if there is an issue, return empty mcmc matrix
            if isequal(shock_type, 'none')
                self.mcmc_conditional_forecast = [];
                if self.verbose
                    cu.progress_bar_complete('Conditional forecasts:');
                end
            % if condition type is well defined, proceed
            else
                % make regressors Z_p and Y
                [Z_p Y] = vu.make_forecast_regressors(Z_p, self.Y, h, self.p, self.T,...
                          self.exogenous, self.constant, self.trend, self.quadratic_trend);
                has_irf = ~isempty(self.mcmc_structural_irf) && size(self.mcmc_structural_irf, 3) >= h;
                % make conditional forecast regressors R, y_bar and omega
                [R y_bar omega] = vu.conditional_forecast_regressors_3(conditions, h, self.n);
                if isequal(shock_type, 'shock-specific')
                    [P non_generating] = vu.conditional_forecast_regressors_5(shocks, h, self.n);
                end
                % initiate storage and loop over iterations
                mcmc_conditional_forecast = zeros(h,self.n,self.iterations);
                for i=1:self.iterations 
                    index = self.svar_index(i);
                    % make predictions, absent shocks
                    f = vvu.varma_linear_forecast(self.mcmc_beta(:,:,index), self.mcmc_kappa(:,:,index), ...
                        h, Z_p, Y, self.mcmc_E(end-self.q+1:end,:,index), self.n);
                    % recover structural IRF or estimate them
                    if has_irf
                        structural_irf = self.mcmc_structural_irf(:,:,1:h,i);
                    else
                        irf = vvu.varma_impulse_response_function(self.mcmc_beta(:,:,index), ...
                              self.mcmc_kappa(:,:,index), self.n, self.p, self.q, h);
                        structural_irf = vu.structural_impulse_response_function(irf, self.mcmc_H(:,:,i), self.n);
                    end
                    % recover iteration-specific regressors
                    M = vu.conditional_forecast_regressors_4(structural_irf, self.n, h);
                    % get posterior mean and variance, depending on condition type    
                    if isequal(shock_type, 'all_shocks')
                        [mu_hat Omega_hat] = vu.conditional_forecast_posterior(y_bar, f, M, R, self.mcmc_Gamma(i,:)', omega, self.n, h);
                    elseif isequal(shock_type, 'shock-specific')
                        Gamma_nd = vu.conditional_forecast_regressors_6(self.mcmc_Gamma(i,:)', non_generating, h);                      
                        [mu_hat Omega_hat] = vu.shock_specific_conditional_forecast_posterior(...
                                             y_bar, f, M, R, P, self.mcmc_Gamma(i,:)', Gamma_nd, omega, self.n, h);  
                    end
                    % sample values
                    mcmc_conditional_forecast(:,:,i) = reshape(rn.multivariate_normal(mu_hat, Omega_hat)',[self.n,h])';
                    if self.verbose
                        cu.progress_bar(i, self.iterations, 'Conditional forecasts:');
                    end
                end
                self.mcmc_conditional_forecast = mcmc_conditional_forecast;
            end
        end  


    end


end