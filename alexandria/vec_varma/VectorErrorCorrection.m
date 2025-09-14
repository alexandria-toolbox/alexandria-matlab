classdef VectorErrorCorrection < handle & BayesianVar


    % Vector Error Correction, developed in chapter 15
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
    % max_cointegration_rank : int, default = 1
    %     maximum number of cointegration relations
    %
    % prior_type : int, default = 1
    %     prior for VEC model
    %     1 = uninformative, 2 = horseshoe, 3 = selection
    %
    % error_correction_type : int, default = 1
    %     error correction type for VEC model
    %     1 = general, 2 = reduced-rank
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
    % pi1 : float, default = 0.1
    %     overall tightness hyperparameter, defined in (5.15.18)
    % 
    % pi2 : float, default = 0.5
    %     cross-variable shrinkage hyperparameter, defined in (5.15.18)
    % 
    % pi3 : float, default = 1
    %     lag decay hyperparameter, defined in (5.15.18)    
    % 
    % pi4 : float, default = 100
    %     exogenous slackness hyperparameter, defined in (5.15.18)  
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
    % max_cointegration_rank : int, default = 1
    %     maximum number of cointegration relations
    %
    % prior_type : int, default = 1
    %     prior for VEC model
    %     1 = uninformative, 2 = horseshoe, 3 = selection
    %
    % error_correction_type : int, default = 1
    %     error correction type for VEC model
    %     1 = general, 2 = reduced-rank
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
    % pi1 : float, default = 0.1
    %     overall tightness hyperparameter, defined in (5.15.18)
    % 
    % pi2 : float, default = 0.5
    %     cross-variable shrinkage hyperparameter, defined in (5.15.18)
    % 
    % pi3 : float, default = 1
    %     lag decay hyperparameter, defined in (5.15.18)    
    % 
    % pi4 : float, default = 100
    %     exogenous slackness hyperparameter, defined in (5.15.18)  
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
    %     matrix of in-sample endogenous variables, defined in (4.11.3)
    % 
    % X : matrix of size (T,k)
    %     matrix of exogenous and lagged regressors, defined in (4.11.3)
    % 
    % DY : matrix of size (T,n)
    %     matrix of differenced endogenous variables, defined in (5.15.8)
    % 
    % Y_1 : matrix of size (T,n)
    %     matrix of lagged endogenous variables, defined in (5.15.8)
    % 
    % Z : matrix of size (T,k)
    %     matrix of exogenous and differenced lagged regressors, defined in (5.15.8)
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
    % r : int
    %     maximum number of cointegration relations
    %
    % Q : matrix of size (q,1)
    %     prior variance of autoregressive coefficients, defined in (5.15.18)           
    %
    % alpha : float
    %     prior degrees of freedom, defined in (5.15.19)
    % 
    % S : matrix of size (n,n)
    %     prior scale matrix, defined in (5.15.19)
    % 
    % mcmc_Xi : matrix of size (n,n,iterations)
    %     MCMC values of error correction coefficients   
    % 
    % mcmc_Phi : matrix of size (k,n,iterations)
    %     MCMC values of autoregressive coefficients   
    %     
    % mcmc_Sigma : matrix of size (n,n,iterations)
    %     MCMC values of residual variance-covariance matrix
    %    
    % mcmc_chol_Sigma : matrix of size (n,n,iterations)
    %     MCMC values of residual variance-covariance matrix (Cholesky factor)
    % 
    % mcmc_K : matrix of size (n,r,iterations)
    %     MCMC values of cointegration relations  
    % 
    % mcmc_Lamda : matrix of size (n,r,iterations)
    %     MCMC values of cointegration loadings  
    % 
    % mcmc_beta : matrix of size (k,n,iterations)
    %     MCMC values of VAR coefficients   
    % 
    % mcmc_nu : matrix of size (iterations,1)
    %     MCMC values of nu coefficients, defined in (5.15.40)
    % 
    % mcmc_tau_2 : matrix of size (iterations,1)
    %     MCMC values of tau_2 coefficients, defined in (5.15.39)
    % 
    % mcmc_eta : matrix of size (iterations,1)
    %     MCMC values of eta coefficients, defined in (5.15.42)
    % 
    % mcmc_psi_2 : matrix of size (iterations,1)
    %     MCMC values of psi_2 coefficients, defined in (5.15.39)
    % 
    % mcmc_omega : matrix of size (iterations,1)
    %     MCMC values of nu coefficients, defined in (5.15.63)
    % 
    % mcmc_zeta_2 : matrix of size (iterations,1)
    %     MCMC values of zeta_2 coefficients, defined in (5.15.58)
    % 
    % mcmc_delta : matrix of size (iterations,1)
    %     MCMC values of delta coefficients, defined in (5.15.84)
    % 
    % mcmc_chi : matrix of size (iterations,1)
    %     MCMC values of chi coefficients, defined in (5.15.90)
    % 
    % mcmc_gamma : matrix of size (iterations,1)
    %     MCMC values of gamma coefficients, defined in (5.15.92)
    %     
    % Xi_estimates : matrix of size (n,n,4)
    %     estimates of error correction coefficients
    %     page 1: median, page 2: st dev,  page 3: lower bound, page 4: upper bound
    %
    % Phi_estimates : matrix of size (k,n,3)
    %     estimates of autoregressive coefficients
    %     page 1: median, page 2: st dev,  page 3: lower bound, page 4: upper bound
    %
    % beta_estimates : matrix of size (k,n,3)
    %     estimates of VAR coefficients
    %     page 1: median, page 2: st dev,  page 3: lower bound, page 4: upper bound
    %
    % Sigma_estimates : matrix of size (n,n)
    %     estimates of variance-covariance matrix of VAR residuals
    % 
    % tau_2_estimates : float or matrix of size (r,1)
    %     estimates of tightness coefficients tau_2
    % 
    % psi_2_estimates : float or matrix of size (n,1)
    %     estimates of tightness coefficients psi_2
    % 
    % zeta_2_estimates : float or matrix of size (n,1)
    %     estimates of tightness coefficients zeta_2
    % 
    % K_estimates : matrix of size (n,r,4)
    %     estimates of cointegration relations K
    %     page 1: median, page 2: st dev,  page 3: lower bound, page 4: upper bound
    %
    % Lamda_estimates : matrix of size (n,r,4)
    %     estimates of loadings matrix Lamda
    %     page 1: median, page 2: st dev,  page 3: lower bound, page 4: upper bound    
    %
    % delta_estimates : float or matrix of size (n,n)
    %     estimates of binary coefficients delta (mean of all iterations)
    % 
    % chi_estimates : float or matrix of size (n,n)
    %     estimates of binary coefficients chi (mean of all iterations)
    % 
    % gamma_estimates : float or matrix of size (n,n)
    %     estimates of binary coefficients gamma (mean of all iterations)
    % 
    % mcmc_H :  matrix of size (n,n,iterations)
    %     MCMC values of structural identification matrix, defined in (4.13.5)
    %
    % mcmc_Gamma : matrix of size (iterations,n)
    %     MCMC values of structural shock variance matrix, defined in definition 13.1
    % 
    % s : matrix of size (n,1)
    %     prior scale matrix, defined in (5.15.19) 
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
     
    
    properties (GetAccess = public, SetAccess = protected)
        endogenous
        exogenous
        structural_identification
        restriction_table
        lags
        max_cointegration_rank
        prior_type
        error_correction_type
        constant
        trend
        quadratic_trend
        pi1
        pi2
        pi3
        pi4
        credibility_level
        iterations
        burnin
        verbose
        Y
        X
        DY
        Y_1
        Z
        n
        m
        p
        T
        k
        q
        r
        %s
        Q
        alpha
        S
        mcmc_Xi
        mcmc_Phi
        mcmc_Sigma
        mcmc_chol_Sigma
        mcmc_K
        mcmc_Lamda
        mcmc_beta      
        mcmc_nu
        mcmc_tau_2
        mcmc_eta
        mcmc_psi_2
        mcmc_omega
        mcmc_zeta_2
        mcmc_delta
        mcmc_chi
        mcmc_gamma
        Xi_estimates
        Phi_estimates
        beta_estimates
        Sigma_estimates
        tau_2_estimates
        psi_2_estimates
        zeta_2_estimates
        K_estimates
        Lamda_estimates
        delta_estimates
        chi_estimates
        gamma_estimates
    end


    properties (GetAccess = public, SetAccess = {?BayesianVar})
        b
        V
        mcmc_H
        mcmc_Gamma
    end
    
    
    properties (GetAccess = private, SetAccess = private)
        Z_vec
        p_vec
        T_vec
        k_vec
        q_vec
        ZZ
        ZDY
        ZY
        YY
        YDY
        YZ
        inv_Q
        alpha_bar
        mcmc_inv_Sigma
    end

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods (Access = public)
        
        
        function self = VectorErrorCorrection(endogenous, varargin)
            
            % constructor for the VectorErrorCorrection class
            
            % allow for optional arguments
            parser = inputParser;
            default_exogenous = [];
            default_structural_identification = 2;
            default_restriction_table = [];
            default_lags = 4;
            default_max_cointegration_rank = 1;
            default_prior_type = 1;
            default_error_correction_type = 1;
            default_constant = true;
            default_trend = false;
            default_quadratic_trend = false;
            default_pi1 = 0.1;
            default_pi2 = 0.5;
            default_pi3 = 1;
            default_pi4 = 100;
            default_credibility_level = 0.95;
            default_iterations = 3000;
            default_burnin = 1000;
            default_verbose = false;
            addRequired(parser, 'endogenous');
            addParameter(parser, 'exogenous', default_exogenous);
            addParameter(parser, 'structural_identification', default_structural_identification);
            addParameter(parser, 'restriction_table', default_restriction_table);
            addParameter(parser, 'lags', default_lags);
            addParameter(parser, 'max_cointegration_rank', default_max_cointegration_rank);
            addParameter(parser, 'prior_type', default_prior_type);
            addParameter(parser, 'error_correction_type', default_error_correction_type);
            addParameter(parser, 'constant', default_constant);
            addParameter(parser, 'trend', default_trend);
            addParameter(parser, 'quadratic_trend', default_quadratic_trend);  
            addParameter(parser, 'pi1', default_pi1);   
            addParameter(parser, 'pi2', default_pi2);  
            addParameter(parser, 'pi3', default_pi3);  
            addParameter(parser, 'pi4', default_pi4);
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
            self.max_cointegration_rank = parser.Results.max_cointegration_rank;
            self.prior_type = parser.Results.prior_type;
            self.error_correction_type = parser.Results.error_correction_type;
            self.constant = parser.Results.constant;
            self.trend = parser.Results.trend;
            self.quadratic_trend = parser.Results.quadratic_trend;
            self.pi1 = parser.Results.pi1;
            self.pi2 = parser.Results.pi2;
            self.pi3 = parser.Results.pi3;
            self.pi4 = parser.Results.pi4;
            self.credibility_level = parser.Results.credibility_level;
            self.iterations = parser.Results.iterations;
            self.burnin = parser.Results.burnin;
            self.verbose = parser.Results.verbose;
            % make regressors
            self.make_regressors();
        end

    
        function estimate(self)
        
            % estimate()
            % generates posterior estimates for Bayesian VEC model parameters
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
            % run MCMC algorithm (Gibbs sampling) for VAR parameters
            self.parameter_mcmc(); 
            % obtain posterior estimates for vector error correction parameters
            self.parameter_estimates();
            % estimate structural identification
            self.make_structural_identification();
        end

 
    end
    
    
    methods (Access = protected, Hidden = true)
        
        
        function make_regressors(self)
            
            % generates VEC regressors along with other dimension elements

            % first generate VAR regressors, for later use once VEC is converted to VAR
            [Y Z X] = vvu.make_var_regressors(self.endogenous, self.exogenous, self.lags+1, ...
                                           self.constant, self.trend, self.quadratic_trend);
            % then generate VEC regressors, following (5.15.8)
            [DY Y_1 Z_vec] = vvu.make_vec_regressors(self.endogenous, self.exogenous, self.lags, ...
                                           self.constant, self.trend, self.quadratic_trend);
            % define dimensions
            [n m p_vec p T_vec T k_vec k q_vec q r] = vvu.generate_dimensions(Y, DY, self.exogenous, self.lags, ...
                              self.max_cointegration_rank, self.constant, self.trend, self.quadratic_trend);
            % get individual ar variances
            [s] = vvu.individual_ar_variances(n, self.endogenous, self.lags);
            % define estimation terms
            ZZ = Z_vec' * Z_vec;
            ZDY = Z_vec' * DY;
            ZY = Z_vec' * Y_1;
            YY = Y_1' * Y_1;
            YDY = Y_1' * DY;
            YZ = Y_1' * Z_vec;                        
            % save as attributes      
            self.Y = Y;
            self.X = X;
            self.DY = DY;
            self.Y_1 = Y_1;
            self.Z_vec = Z_vec;
            self.Z = Z;
            self.n = n;
            self.m = m;
            self.p_vec = p_vec;
            self.p = p;
            self.T_vec = T_vec;
            self.T = T;
            self.k_vec = k_vec;
            self.k = k;
            self.q_vec = q_vec;
            self.q = q;
            self.r = r;
            self.s = s;
            self.ZZ = ZZ;
            self.ZDY = ZDY;
            self.ZY = ZY;
            self.YY = YY;
            self.YDY = YDY;
            self.YZ = YZ;
        end


        function prior(self)
            
            % creates prior elements Q, alpha and S defined in (5.15.18) and (5.15.19)
            
            Q = vu.make_V(self.s, self.pi1, self.pi2, self.pi3, self.pi4, self.n, self.m, self.lags);
            alpha = vu.make_alpha(self.n);
            S = diag(vu.make_S(self.s));
            self.Q = Q;
            self.alpha = alpha;
            self.S = S;
        end

    
        function posterior(self)
            
            % creates posterior elements
    
            inv_Q = diag(1 ./ self.Q);
            alpha_bar = self.alpha + self.T;
            self.inv_Q = inv_Q;
            self.alpha_bar = alpha_bar;
        end

    
        function parameter_mcmc(self)
            
            % Gibbs sampler for VEC parameters, depending on prior and error correction type
            
            % if prior is uninformative and error correction is general
            if self.prior_type == 1 && self.error_correction_type == 1
                self.mcmc_uninformative_general();
            % else, if prior is uninformative and error correction is reduced-rank
            elseif self.prior_type == 1 && self.error_correction_type == 2
                self.mcmc_uninformative_reduced_rank();
            % else, if prior is horseshoe and error correction is reduced-rank
            elseif self.prior_type == 2 && self.error_correction_type == 1
                self.mcmc_horseshoe_general();
            % else, if prior is horseshoe and error correction is reduced-rank
            elseif self.prior_type == 2 && self.error_correction_type == 2
                self.mcmc_horseshoe_reduced_rank();  
            % else, if prior is selection and error correction is general
            elseif self.prior_type == 3 && self.error_correction_type == 1
                self.mcmc_selection_general();
            % else, if prior is selection and error correction is reduced-rank
            elseif self.prior_type == 3 && self.error_correction_type == 2
                self.mcmc_selection_reduced_rank();
            end
        end


        function mcmc_uninformative_general(self)
            
            % Gibbs sampler for VEC with uninformative prior and general approach: algorithm 15.1
            
            % unpack
            DY = self.DY;
            Y_1 = self.Y_1;
            Z = self.Z_vec;
            ZZ = self.ZZ;
            ZDY = self.ZDY;
            ZY = self.ZY;
            YY = self.YY;
            YDY = self.YDY;
            YZ = self.YZ;    
            inv_Q = self.inv_Q;
            alpha_bar = self.alpha_bar;
            S = self.S;
            n = self.n;
            m = self.m;
            p = self.p_vec;
            k = self.k_vec;
            iterations = self.iterations;
            burnin = self.burnin;
            verbose = self.verbose;
    
            % other prior values
            inv_U = diag(ones(n*n,1)/10);

            % preallocate storage space
            mcmc_Xi = zeros(n,n,iterations);
            mcmc_Phi = zeros(k,n,iterations);
            mcmc_Sigma = zeros(n,n,iterations);
            mcmc_chol_Sigma = zeros(n,n,iterations);
            mcmc_inv_Sigma = zeros(n,n,iterations);
            mcmc_beta = zeros(k+n,n,iterations);

            % step 1: set initial values for MCMC algorithm
            Phi = zeros(k,n);
            inv_Sigma = diag(1 ./ diag(S));

            % iterate over iterations
            iteration = 1;
            while iteration <= (burnin + iterations)

                % step 2: sample xi
                [xi Xi_T Xi] = self.draw_xi(inv_U, inv_Sigma, YY, YDY, YZ, Phi, n);

                % step 3: sample phi
                [phi Phi] = self.draw_phi(inv_Q, inv_Sigma, ZZ, ZDY, ZY, Xi_T, k, n);

                % step 4: sample Sigma
                [Sigma inv_Sigma chol_Sigma] = self.draw_Sigma(alpha_bar, S, DY, Y_1, Xi_T, Z, Phi);

                % convert into VAR model
                B = vvu.vec_to_var(Xi_T, Phi, n, m, p, k);

                % save if burn is exceeded, and display progress bar
                if iteration > burnin
                    mcmc_Xi(:,:,iteration-burnin) = Xi;
                    mcmc_Phi(:,:,iteration-burnin) = Phi;       
                    mcmc_Sigma(:,:,iteration-burnin) = Sigma;
                    mcmc_chol_Sigma(:,:,iteration-burnin) = chol_Sigma;
                    mcmc_inv_Sigma(:,:,iteration-burnin) = inv_Sigma;
                    mcmc_beta(:,:,iteration-burnin) = B;
                end
                if verbose
                    cu.progress_bar(iteration, iterations+burnin, 'Model parameters:')
                end
                iteration = iteration + 1;
            end

            % save as attributes
            self.mcmc_Xi = mcmc_Xi;
            self.mcmc_Phi = mcmc_Phi;
            self.mcmc_Sigma = mcmc_Sigma;
            self.mcmc_chol_Sigma = mcmc_chol_Sigma;
            self.mcmc_inv_Sigma = mcmc_inv_Sigma;
            self.mcmc_beta = mcmc_beta;
        end


        function mcmc_uninformative_reduced_rank(self)
            
            % Gibbs sampler for VEC with uninformative prior and reduced-rank approach: algorithm 15.2
            
            % unpack
            DY = self.DY;
            Y_1 = self.Y_1;
            Z = self.Z_vec;
            ZZ = self.ZZ;
            ZDY = self.ZDY;
            ZY = self.ZY;
            YY = self.YY;
            YDY = self.YDY;
            YZ = self.YZ;    
            inv_Q = self.inv_Q;
            alpha_bar = self.alpha_bar;
            S = self.S;
            n = self.n;
            m = self.m;
            p = self.p_vec;
            k = self.k_vec;
            r = self.r;
            iterations = self.iterations;
            burnin = self.burnin;
            verbose = self.verbose;
    
            % other prior values
            inv_R = diag(ones(n*r,1)/10);
            inv_P = diag(ones(n*r,1)/10);

            % preallocate storage space
            mcmc_K = zeros(n,r,iterations);
            mcmc_Lamda = zeros(n,r,iterations);
            mcmc_Phi = zeros(k,n,iterations);
            mcmc_Xi = zeros(n,n,iterations);
            mcmc_Sigma = zeros(n,n,iterations);
            mcmc_chol_Sigma = zeros(n,n,iterations);
            mcmc_inv_Sigma = zeros(n,n,iterations);
            mcmc_beta = zeros(k+n,n,iterations);

            % step 1: set initial values for MCMC algorithm
            YYZ = eye(n);
            Lamda = eye(n,r);
            Phi = zeros(k,n);
            inv_Sigma = diag(1 ./ diag(S));

            % iterate over iterations
            iteration = 1;
            while iteration <= (burnin + iterations)

                % step 2: sample kappa
                [kappa K] = self.draw_kappa(inv_R, Lamda, inv_Sigma, YY, YYZ, n, r);
    
                % step 3: sample lambda
                [lamda Lamda_T Lamda Xi_T Xi] = self.draw_lamda(inv_P, inv_Sigma, K, YY, YYZ, n, r);

                % step 4: sample phi
                [phi Phi] = self.draw_phi(inv_Q, inv_Sigma, ZZ, ZDY, ZY, Xi_T, k, n);

                % step 5: sample Sigma
                [Sigma inv_Sigma chol_Sigma] = self.draw_Sigma(alpha_bar, S, DY, Y_1, Xi_T, Z, Phi);
                YYZ = (YDY - YZ * Phi) * inv_Sigma;

                % convert into VAR model
                B = vvu.vec_to_var(Xi_T, Phi, n, m, p, k);

                % save if burn is exceeded, and display progress bar
                if iteration > burnin
                    mcmc_K(:,:,iteration-burnin) = K;
                    mcmc_Lamda(:,:,iteration-burnin) = Lamda;
                    mcmc_Xi(:,:,iteration-burnin) = Xi;
                    mcmc_Phi(:,:,iteration-burnin) = Phi;       
                    mcmc_Sigma(:,:,iteration-burnin) = Sigma;
                    mcmc_chol_Sigma(:,:,iteration-burnin) = chol_Sigma;
                    mcmc_inv_Sigma(:,:,iteration-burnin) = inv_Sigma;
                    mcmc_beta(:,:,iteration-burnin) = B;
                end
                if verbose
                    cu.progress_bar(iteration, iterations+burnin, 'Model parameters:')
                end
                iteration = iteration + 1;
            end

            % save as attributes
            self.mcmc_K = mcmc_K;
            self.mcmc_Lamda = mcmc_Lamda;
            self.mcmc_Xi = mcmc_Xi;
            self.mcmc_Phi = mcmc_Phi;
            self.mcmc_Sigma = mcmc_Sigma;
            self.mcmc_chol_Sigma = mcmc_chol_Sigma;
            self.mcmc_inv_Sigma = mcmc_inv_Sigma;
            self.mcmc_beta = mcmc_beta;
        end


        function mcmc_horseshoe_general(self)
            
            % Gibbs sampler for VEC with horseshoe prior and general approach: algorithm 15.3
            
            % unpack
            DY = self.DY;
            Y_1 = self.Y_1;
            Z = self.Z_vec;
            ZZ = self.ZZ;
            ZDY = self.ZDY;
            ZY = self.ZY;
            YY = self.YY;
            YDY = self.YDY;
            YZ = self.YZ;    
            inv_Q = self.inv_Q;
            alpha_bar = self.alpha_bar;
            S = self.S;
            n = self.n;
            m = self.m;
            p = self.p_vec;
            k = self.k_vec;
            iterations = self.iterations;
            burnin = self.burnin;
            verbose = self.verbose;
    
            % other prior values
            inv_U = diag(ones(n*n,1)/10);
            a_nu = 1;
            a_tau = (n*n+1) / 2;
            a_eta = 1;
            a_psi = 1;

            % preallocate storage space
            mcmc_nu = zeros(iterations,1);
            mcmc_tau_2 = zeros(iterations,1);
            mcmc_eta = zeros(n,n,iterations);
            mcmc_psi_2 = zeros(n,n,iterations);
            mcmc_Xi = zeros(n,n,iterations);
            mcmc_Phi = zeros(k,n,iterations);
            mcmc_Sigma = zeros(n,n,iterations);
            mcmc_chol_Sigma = zeros(n,n,iterations);
            mcmc_inv_Sigma = zeros(n,n,iterations);
            mcmc_beta = zeros(k+n,n,iterations);

            % step 1: set initial values for MCMC algorithm
            tau_2 = 1;
            psi_2 = ones(n,n);
            Xi = ones(n,n);
            Phi = zeros(k,n);
            inv_Sigma = diag(1 ./ diag(S));

            % iterate over iterations
            iteration = 1;
            while iteration <= (burnin + iterations)

                % step 2: sample nu
                [nu] = self.draw_nu(a_nu, tau_2);
                
                % step 3: sample tau_2
                [tau_2] = self.draw_tau_2(a_tau, Xi, psi_2, nu);
                
                % step 4: sample eta
                [eta] = self.draw_eta(a_eta, psi_2, n);
                
                % step 5: sample psi_2
                [psi_2] = self.draw_psi_2(a_psi, Xi, tau_2, eta, n);    
                inv_U = diag(1 ./ la.vec(tau_2 * psi_2'));

                % step 6: sample xi
                [xi Xi_T Xi] = self.draw_xi(inv_U, inv_Sigma, YY, YDY, YZ, Phi, n);

                % step 7: sample phi
                [phi Phi] = self.draw_phi(inv_Q, inv_Sigma, ZZ, ZDY, ZY, Xi_T, k, n);

                % step 8: sample Sigma
                [Sigma inv_Sigma chol_Sigma] = self.draw_Sigma(alpha_bar, S, DY, Y_1, Xi_T, Z, Phi);

                % convert into VAR model
                B = vvu.vec_to_var(Xi_T, Phi, n, m, p, k);

                % save if burn is exceeded, and display progress bar
                if iteration > burnin
                    mcmc_nu(iteration-burnin) = nu;
                    mcmc_tau_2(iteration-burnin) = tau_2;
                    mcmc_eta(:,:,iteration-burnin) = eta;
                    mcmc_psi_2(:,:,iteration-burnin) = psi_2;
                    mcmc_Xi(:,:,iteration-burnin) = Xi;
                    mcmc_Phi(:,:,iteration-burnin) = Phi;       
                    mcmc_Sigma(:,:,iteration-burnin) = Sigma;
                    mcmc_chol_Sigma(:,:,iteration-burnin) = chol_Sigma;
                    mcmc_inv_Sigma(:,:,iteration-burnin) = inv_Sigma;
                    mcmc_beta(:,:,iteration-burnin) = B;
                end
                if verbose
                    cu.progress_bar(iteration, iterations+burnin, 'Model parameters:')
                end
                iteration = iteration + 1;
            end

            % save as attributes
            self.mcmc_nu = mcmc_nu;
            self.mcmc_tau_2 = mcmc_tau_2;
            self.mcmc_eta = mcmc_eta;
            self.mcmc_psi_2 = mcmc_psi_2;
            self.mcmc_Xi = mcmc_Xi;
            self.mcmc_Phi = mcmc_Phi;
            self.mcmc_Sigma = mcmc_Sigma;
            self.mcmc_chol_Sigma = mcmc_chol_Sigma;
            self.mcmc_inv_Sigma = mcmc_inv_Sigma;
            self.mcmc_beta = mcmc_beta;
        end


        function mcmc_horseshoe_reduced_rank(self)
            
            % Gibbs sampler for VEC with horseshoe prior and reduced-rank approach: algorithm 15.4
            
            % unpack
            DY = self.DY;
            Y_1 = self.Y_1;
            Z = self.Z_vec;
            ZZ = self.ZZ;
            ZDY = self.ZDY;
            ZY = self.ZY;
            YY = self.YY;
            YDY = self.YDY;
            YZ = self.YZ;    
            inv_Q = self.inv_Q;
            alpha_bar = self.alpha_bar;
            S = self.S;
            n = self.n;
            m = self.m;
            p = self.p_vec;
            k = self.k_vec;
            r = self.r;
            iterations = self.iterations;
            burnin = self.burnin;
            verbose = self.verbose;
    
            % other prior values
            a_nu = 1;
            a_tau = n + 1 / 2;
            a_eta = 1;
            a_psi = (r + 1) / 2;
            a_omega = 1;
            a_zeta = (r + 1) / 2;

            % preallocate storage space
            mcmc_nu = zeros(r,iterations);
            mcmc_tau_2 = zeros(r,iterations);        
            mcmc_eta = zeros(n,iterations);
            mcmc_psi_2 = zeros(n,iterations);
            mcmc_omega = zeros(n,iterations);
            mcmc_zeta_2 = zeros(n,iterations);
            mcmc_K = zeros(n,r,iterations);
            mcmc_Lamda = zeros(n,r,iterations);
            mcmc_Phi = zeros(k,n,iterations);
            mcmc_Xi = zeros(n,n,iterations);
            mcmc_Sigma = zeros(n,n,iterations);
            mcmc_chol_Sigma = zeros(n,n,iterations);
            mcmc_inv_Sigma = zeros(n,n,iterations);
            mcmc_beta = zeros(k+n,n,iterations);

            % step 1: set initial values for MCMC algorithm
            tau_2 = ones(r,1);
            psi_2 = ones(n,1);
            zeta_2 = ones(n,1);
            K = eye(n,r);
            Lamda = eye(n,r);
            inv_Sigma = diag(1 ./ diag(S));        
            YYZ = eye(n);

            % iterate over iterations
            iteration = 1;
            while iteration <= (burnin + iterations)

                % step 2: sample nu
                [nu] = self.draw_nu_j(a_nu, tau_2, r);
    
                % step 3: sample tau_2
                [tau_2] = self.draw_tau_2_j(a_tau, K, Lamda, psi_2, zeta_2, nu, r);
    
                % step 4: sample eta
                [eta] = self.draw_eta_i(a_eta, psi_2, n);
    
                % step 5: sample psi_2
                [psi_2] = self.draw_psi_2_i(a_psi, K, tau_2, eta, n);       
                inv_R = diag(1 ./ la.vec(psi_2 * tau_2'));
                
                % step 6: sample omega
                [omega] = self.draw_omega_i(a_omega, zeta_2, n);
                
                % step 7: sample zeta_2
                [zeta_2] = self.draw_zeta_2_i(a_zeta, Lamda, tau_2, omega, n);
                inv_P = diag(1 ./ la.vec(tau_2 * zeta_2'));

                % step 8: sample kappa
                [kappa K] = self.draw_kappa(inv_R, Lamda, inv_Sigma, YY, YYZ, n, r);
    
                % step 9: sample lambda
                [lamda Lamda_T Lamda Xi_T Xi] = self.draw_lamda(inv_P, inv_Sigma, K, YY, YYZ, n, r);

                % step 10: sample phi
                [phi Phi] = self.draw_phi(inv_Q, inv_Sigma, ZZ, ZDY, ZY, Xi_T, k, n);

                % step 11: sample Sigma
                [Sigma inv_Sigma chol_Sigma] = self.draw_Sigma(alpha_bar, S, DY, Y_1, Xi_T, Z, Phi);
                YYZ = (YDY - YZ * Phi) * inv_Sigma;

                % convert into VAR model
                B = vvu.vec_to_var(Xi_T, Phi, n, m, p, k);

                % save if burn is exceeded, and display progress bar
                if iteration > burnin
                    mcmc_nu(:,iteration-burnin) = nu;
                    mcmc_tau_2(:,iteration-burnin) = tau_2; 
                    mcmc_eta(:,iteration-burnin) = eta;
                    mcmc_psi_2(:,iteration-burnin) = psi_2;
                    mcmc_omega(:,iteration-burnin) = omega;
                    mcmc_zeta_2(:,iteration-burnin) = zeta_2;
                    mcmc_K(:,:,iteration-burnin) = K;
                    mcmc_Lamda(:,:,iteration-burnin) = Lamda;
                    mcmc_Xi(:,:,iteration-burnin) = Xi;
                    mcmc_Phi(:,:,iteration-burnin) = Phi;       
                    mcmc_Sigma(:,:,iteration-burnin) = Sigma;
                    mcmc_chol_Sigma(:,:,iteration-burnin) = chol_Sigma;
                    mcmc_inv_Sigma(:,:,iteration-burnin) = inv_Sigma;
                    mcmc_beta(:,:,iteration-burnin) = B;
                end
                if verbose
                    cu.progress_bar(iteration, iterations+burnin, 'Model parameters:')
                end
                iteration = iteration + 1;
            end

            % save as attributes
            self.mcmc_nu = mcmc_nu;
            self.mcmc_tau_2 = mcmc_tau_2;
            self.mcmc_eta = mcmc_eta;
            self.mcmc_psi_2 = mcmc_psi_2;
            self.mcmc_omega = mcmc_omega;
            self.mcmc_zeta_2 = mcmc_zeta_2;
            self.mcmc_K = mcmc_K;
            self.mcmc_Lamda = mcmc_Lamda;
            self.mcmc_Xi = mcmc_Xi;
            self.mcmc_Phi = mcmc_Phi;
            self.mcmc_Sigma = mcmc_Sigma;
            self.mcmc_chol_Sigma = mcmc_chol_Sigma;
            self.mcmc_inv_Sigma = mcmc_inv_Sigma;
            self.mcmc_beta = mcmc_beta;
        end

    
        function mcmc_selection_general(self)
            
            % Gibbs sampler for VEC with selection prior and general approach: algorithm 15.5
            
            % unpack
            DY = self.DY;
            Y_1 = self.Y_1;
            Z = self.Z_vec;
            ZZ = self.ZZ;
            ZDY = self.ZDY;
            ZY = self.ZY;
            YY = self.YY;
            YDY = self.YDY;
            YZ = self.YZ;   
            inv_Q = self.inv_Q;
            alpha_bar = self.alpha_bar;
            S = self.S;
            n = self.n;
            m = self.m;
            p = self.p_vec;
            k = self.k_vec;
            iterations = self.iterations;
            burnin = self.burnin;
            verbose = self.verbose;
    
            % other prior values
            mu = 0.5;
            tau_1_2 = 10;
            tau_1 = tau_1_2^0.5;
            tau_2_2 = 0.01;
            tau_2 = tau_2_2^0.5;
            temp_1 = mu / tau_1;
            temp_2 = - 0.5 / tau_1_2;
            temp_3 = (1 - mu) / tau_2;
            temp_4 = - 0.5 / tau_2_2;
        
            % preallocate storage space
            mcmc_delta = zeros(n*n,iterations);
            mcmc_Xi = zeros(n,n,iterations);
            mcmc_Phi = zeros(k,n,iterations);
            mcmc_Sigma = zeros(n,n,iterations);
            mcmc_chol_Sigma = zeros(n,n,iterations);
            mcmc_inv_Sigma = zeros(n,n,iterations);
            mcmc_beta = zeros(k+n,n,iterations);
    
            % step 1: set initial values for MCMC algorithm
            xi = ones(n*n,1);
            Phi = zeros(k,n);
            inv_Sigma = diag(1 ./ diag(S));
            
            % iterate over iterations
            iteration = 1;
            while iteration <= (burnin + iterations)
                
                % step 2: sample delta
                [delta] = self.draw_delta(xi, temp_1, temp_2, temp_3, temp_4, n);
                inv_U = diag(1 ./ (delta * tau_1_2 + (1 - delta) * tau_2_2));
    
                % step 3: sample xi
                [xi Xi_T Xi] = self.draw_xi(inv_U, inv_Sigma, YY, YDY, YZ, Phi, n);
    
                % step 4: sample phi
                [phi Phi] = self.draw_phi(inv_Q, inv_Sigma, ZZ, ZDY, ZY, Xi_T, k, n);
                
                % step 5: sample Sigma
                [Sigma inv_Sigma chol_Sigma] = self.draw_Sigma(alpha_bar, S, DY, Y_1, Xi_T, Z, Phi);
                
                % convert into VAR model
                B = vvu.vec_to_var(Xi_T, Phi, n, m, p, k);
    
                % save if burn is exceeded, and display progress bar
                if iteration > burnin
                    mcmc_delta(:,iteration-burnin) = delta;
                    mcmc_Xi(:,:,iteration-burnin) = Xi;
                    mcmc_Phi(:,:,iteration-burnin) = Phi;           
                    mcmc_Sigma(:,:,iteration-burnin) = Sigma;
                    mcmc_chol_Sigma(:,:,iteration-burnin) = chol_Sigma;
                    mcmc_inv_Sigma(:,:,iteration-burnin) = inv_Sigma;
                    mcmc_beta(:,:,iteration-burnin) = B;
                end
                if verbose
                    cu.progress_bar(iteration, iterations+burnin, 'Model parameters:');
                end
                iteration = iteration + 1;
            end
            % save as attributes
            self.mcmc_delta = mcmc_delta;
            self.mcmc_Xi = mcmc_Xi;
            self.mcmc_Phi = mcmc_Phi;
            self.mcmc_Sigma = mcmc_Sigma;
            self.mcmc_chol_Sigma = mcmc_chol_Sigma;
            self.mcmc_inv_Sigma = mcmc_inv_Sigma;
            self.mcmc_beta = mcmc_beta;
        end
    
    
        function mcmc_selection_reduced_rank(self)
            
            % Gibbs sampler for VEC with selection prior and reduced-rank approach: algorithm 15.6
            
            % unpack
            DY = self.DY;
            Y_1 = self.Y_1;
            Z = self.Z_vec;
            ZZ = self.ZZ;
            ZDY = self.ZDY;
            ZY = self.ZY;
            YY = self.YY;
            YDY = self.YDY;
            YZ = self.YZ;      
            inv_Q = self.inv_Q;
            alpha_bar = self.alpha_bar;
            S = self.S;
            n = self.n;
            m = self.m;
            p = self.p_vec;
            k = self.k_vec;
            r = self.r;
            iterations = self.iterations;
            burnin = self.burnin;
            verbose = self.verbose;
    
            % other prior values
            mu = 0.5;
            tau_1_2 = 10;
            tau_1 = tau_1_2^0.5;
            tau_2_2 = 0.01;
            tau_2 = tau_2_2^0.5;
            temp_1 = mu / tau_1;
            temp_2 = - 0.5 / tau_1_2;
            temp_3 = (1 - mu) / tau_2;
            temp_4 = - 0.5 / tau_2_2;
            
            % preallocate storage space
            mcmc_chi = zeros(n*r,iterations);
            mcmc_gamma = zeros(n*r,iterations);
            mcmc_K = zeros(n,r,iterations);
            mcmc_Lamda = zeros(n,r,iterations);
            mcmc_Phi = zeros(k,n,iterations);
            mcmc_Xi = zeros(n,n,iterations);
            mcmc_Sigma = zeros(n,n,iterations);
            mcmc_chol_Sigma = zeros(n,n,iterations);
            mcmc_inv_Sigma = zeros(n,n,iterations);
            mcmc_beta = zeros(k+n,n,iterations);
    
            % step 1: set initial values for MCMC algorithm
            kappa = ones(n*r,1);
            lamda = ones(n*r,1);
            Phi = zeros(k,n);
            inv_Sigma = diag(1 ./ diag(S));
            YYZ = eye(n);
            Lamda = eye(n,r);
            Phi = zeros(k,n);
            inv_Sigma = diag(1 ./ diag(S));
            
            % iterate over iterations
            iteration = 1;
            while iteration <= (burnin + iterations)
                
                % step 2: sample chi
                [chi] = self.draw_chi(kappa, temp_1, temp_2, temp_3, temp_4, n, r);
                inv_R = diag(1 ./ (chi * tau_1_2 + (1 - chi) * tau_2_2));           
                
                % step 3: sample kappa
                [kappa K] = self.draw_kappa(inv_R, Lamda, inv_Sigma, YY, YYZ, n, r);
    
                % step 4: sample gamma
                [gama] = self.draw_gamma(lamda, temp_1, temp_2, temp_3, temp_4, n, r);
                inv_P = diag(1 ./ (gama * tau_1_2 + (1 - gama) * tau_2_2));
    
                % step 5: sample lambda
                [lamda Lamda_T Lamda Xi_T Xi] = self.draw_lamda(inv_P, inv_Sigma, K, YY, YYZ, n, r);
                
                % step 6: sample phi
                [phi Phi] = self.draw_phi(inv_Q, inv_Sigma, ZZ, ZDY, ZY, Xi_T, k, n);
                
                % step 7: sample Sigma
                [Sigma inv_Sigma chol_Sigma] = self.draw_Sigma(alpha_bar, S, DY, Y_1, Xi_T, Z, Phi);
                YYZ = (YDY - YZ * Phi) * inv_Sigma;
                
                % convert into VAR model
                B = vvu.vec_to_var(Xi_T, Phi, n, m, p, k);
    
                % save if burn is exceeded, and display progress bar
                if iteration > burnin
                    mcmc_chi(:,iteration-burnin) = chi;
                    mcmc_gamma(:,iteration-burnin) = gama;
                    mcmc_K(:,:,iteration-burnin) = K;
                    mcmc_Lamda(:,:,iteration-burnin) = Lamda;
                    mcmc_Xi(:,:,iteration-burnin) = Xi;
                    mcmc_Phi(:,:,iteration-burnin) = Phi;           
                    mcmc_Sigma(:,:,iteration-burnin) = Sigma;
                    mcmc_chol_Sigma(:,:,iteration-burnin) = chol_Sigma;
                    mcmc_inv_Sigma(:,:,iteration-burnin) = inv_Sigma;
                    mcmc_beta(:,:,iteration-burnin) = B;
                end
                if verbose
                    cu.progress_bar(iteration, iterations+burnin, 'Model parameters:');
                end
                iteration = iteration + 1;
            end
                
            % save as attributes
            self.mcmc_chi = mcmc_chi;
            self.mcmc_gamma = mcmc_gamma;
            self.mcmc_K = mcmc_K;
            self.mcmc_Lamda = mcmc_Lamda;
            self.mcmc_Xi = mcmc_Xi;
            self.mcmc_Phi = mcmc_Phi;
            self.mcmc_Sigma = mcmc_Sigma;
            self.mcmc_chol_Sigma = mcmc_chol_Sigma;
            self.mcmc_inv_Sigma = mcmc_inv_Sigma;
            self.mcmc_beta = mcmc_beta;
        end
    

        function [xi Xi_T Xi] = draw_xi(self, inv_U, inv_Sigma, YY, YDY, YZ, Phi, n)
            
            % draw xi from its conditional posterior defined in (5.15.29)
            
            inv_U_bar = inv_U + kron(inv_Sigma, YY);
            d_bar_temp = la.vec((YDY - YZ * Phi) * inv_Sigma);
            xi = rn.efficient_multivariate_normal(d_bar_temp, inv_U_bar);
            Xi_T = reshape(xi,[n n]);
            Xi = Xi_T';
        end

    
        function [phi Phi] = draw_phi(self, inv_Q, inv_Sigma, ZZ, ZDY, ZY, Xi_T, k, n)
            
            % draw phi from its conditional posterior defined in (5.15.22)
            
            inv_Q_bar = inv_Q + kron(inv_Sigma, ZZ);
            f_bar_temp = la.vec((ZDY - ZY * Xi_T) * inv_Sigma);
            phi = rn.efficient_multivariate_normal(f_bar_temp, inv_Q_bar);
            Phi = reshape(phi,[k n]);
        end


        function [Sigma inv_Sigma chol_Sigma] = draw_Sigma(self, alpha_bar, S, DY, Y_1, Xi_T, Z, Phi)
            
            % draw Sigma from its conditional posterior defined in (5.15.25)
            
            residuals = DY - Y_1 * Xi_T - Z * Phi;
            S_bar = S + residuals' * residuals;
            Sigma = rn.inverse_wishart(alpha_bar, S_bar);
            inv_Sigma = la.invert_spd_matrix(Sigma);
            chol_Sigma = la.cholesky_nspd(Sigma);
        end

    
        function [kappa K] = draw_kappa(self, inv_R, Lamda, inv_Sigma, YY, YYZ, n, r)
            
            % draw kappa from its conditional posterior defined in (5.15.35)     
    
            inv_R_bar = inv_R + kron(Lamda' * inv_Sigma * Lamda, YY);
            g_bar_temp = la.vec(YYZ * Lamda);
            kappa = rn.efficient_multivariate_normal(g_bar_temp, inv_R_bar);
            K = reshape(kappa,[n r]);
        end
            
            
        function [lamda Lamda_T Lamda Xi_T Xi] = draw_lamda(self, inv_P, inv_Sigma, K, YY, YYZ, n, r)
            
            % draw lamda from its conditional posterior defined in (5.15.38) 
    
            inv_P_bar = inv_P + kron(inv_Sigma, K' * YY * K);
            h_bar_temp = la.vec(K' * YYZ);
            lamda = rn.efficient_multivariate_normal(h_bar_temp, inv_P_bar);
            Lamda_T = reshape(lamda,[r n]);
            Lamda = Lamda_T';
            Xi_T = K * Lamda_T;
            Xi = Xi_T';
        end   


        function [nu] = draw_nu(self, a_nu, tau_2)
            
            % draw nu from its conditional posterior defined in (5.15.50)     
    
            b_nu = 1 / tau_2 + 1;
            nu = rn.inverse_gamma(a_nu, b_nu);
        end


        function [tau_2] = draw_tau_2(self, a_tau, Xi, psi_2, nu)
    
            % draw tau_2 from its conditional posterior defined in (5.15.47)
    
            b_tau = sum((Xi .* Xi) ./ psi_2, "all") / 2 + 1 / nu;
            tau_2 = rn.inverse_gamma(a_tau, b_tau);
        end


        function [eta] = draw_eta(self, a_eta, psi_2, n)
            
            % draw eta from its conditional posterior defined in (5.15.56)      
    
            b_eta = 1 ./ psi_2 + 1;
            eta = zeros(n,n);
            for i=1:n
                for j=1:n
                    eta(i,j) = rn.inverse_gamma(a_eta, b_eta(i,j));
                end
            end
        end


        function [psi_2] = draw_psi_2(self, a_psi, Xi, tau_2, eta, n)
        
            % draw psi_2 from its conditional posterior defined in (5.15.53)   
    
            b_psi = (Xi .* Xi) / (2 * tau_2) + 1 ./ eta;
            psi_2 = zeros(n,n);
            for i=1:n
                for j=1:n
                    psi_2(i,j) = rn.inverse_gamma(a_psi, b_psi(i,j));
                end
            end
        end


        function [nu] = draw_nu_j(self, a_nu, tau_2, r)
            
            % draw nu_j from its conditional posterior defined in (5.15.71)     
            
            nu = zeros(r,1);
            for j=1:r
                b_nu = 1 / tau_2(j) + 1;
                nu(j) = rn.inverse_gamma(a_nu, b_nu);
            end
        end

            
        function [tau_2] = draw_tau_2_j(self, a_tau, K, Lamda, psi_2, zeta_2, nu, r)
            
            % draw tau_2_j from its conditional posterior defined in (5.15.68)
    
            tau_2 = zeros(r,1);
            for j=1:r
                b_tau = sum((K(:,j) .* K(:,j)) ./ psi_2 + (Lamda(:,j) .* Lamda(:,j)) ./ zeta_2) / 2 + 1 / nu(j);
                tau_2(j) = rn.inverse_gamma(a_tau, b_tau);
            end
        end


        function [eta] = draw_eta_i(self, a_eta, psi_2, n)
    
            % draw eta_i from its conditional posterior defined in (5.15.77)      
    
            eta = zeros(n,1);
            for i=1:n
                b_eta = 1 / psi_2(i) + 1;
                eta(i) = rn.inverse_gamma(a_eta, b_eta);
            end
        end


        function [psi_2] = draw_psi_2_i(self, a_psi, K, tau_2, eta, n)
            
            % draw psi_2_i from its conditional posterior defined in (5.15.74)
    
            psi_2 = zeros(n,1);
            for i=1:n
                b_psi = sum((K(i,:) .* K(i,:)) ./ tau_2') / 2 + 1 / eta(i);
                psi_2(i) = rn.inverse_gamma(a_psi, b_psi);
            end
        end


        function [omega] = draw_omega_i(self, a_omega, zeta_2, n)
    
            % draw omega_i from its conditional posterior defined in (5.15.83)      
    
            omega = zeros(n,1);
            for i=1:n
                b_omega = 1 / zeta_2(i) + 1;
                omega(i) = rn.inverse_gamma(a_omega, b_omega);
            end
        end

    
        function [zeta_2] = draw_zeta_2_i(self, a_zeta, Lamda, tau_2, omega, n)
            
            % draw zeta_2_i from its conditional posterior defined in (5.15.80)
    
            zeta_2 = zeros(n,1);
            for i=1:n
                b_zeta = sum((Lamda(i,:) .* Lamda(i,:)) ./ tau_2') / 2 + 1 / omega(i);
                zeta_2(i) = rn.inverse_gamma(a_zeta, b_zeta);
            end
        end


        function [delta] = draw_delta(self, xi, temp_1, temp_2, temp_3, temp_4, n)
            
            % draw delta from its conditional posterior defined in (5.15.89)
            
            xi_2 = xi .* xi;
            a_1 = temp_1 * exp(temp_2 * xi_2);
            a_2 = temp_3 * exp(temp_4 * xi_2);
            mu_xi = a_1 ./ (a_1 + a_2);
            delta = rn.multivariate_bernoulli(mu_xi, n*n);
        end

    
        function [chi] = draw_chi(self, kappa, temp_1, temp_2, temp_3, temp_4, n, r)
            
            % draw chi from its conditional posterior defined in (5.15.97)
            
            kappa_2 = kappa .* kappa;
            a_1 = temp_1 * exp(temp_2 * kappa_2);
            a_2 = temp_3 * exp(temp_4 * kappa_2);
            mu_kappa = a_1 ./ (a_1 + a_2);
            chi = rn.multivariate_bernoulli(mu_kappa, n*r);
        end
    
    
        function [gama] = draw_gamma(self, lamda, temp_1, temp_2, temp_3, temp_4, n, r)
            
            % draw gamma from its conditional posterior defined in (5.15.100)
            
            lamda_2 = lamda .* lamda;
            a_1 = temp_1 * exp(temp_2 * lamda_2);
            a_2 = temp_3 * exp(temp_4 * lamda_2);
            mu_lamda = a_1 ./ (a_1 + a_2);
            gama = rn.multivariate_bernoulli(mu_lamda, n*r);     
        end


        function parameter_estimates(self)
            
            % point estimates and credibility intervals for model parameters
            % use empirical quantiles from MCMC algorithm
            
            % common parameters
            mcmc_Xi = self.mcmc_Xi;
            mcmc_Phi = self.mcmc_Phi;      
            mcmc_Sigma = self.mcmc_Sigma;
            mcmc_beta = self.mcmc_beta;
            credibility_level = self.credibility_level;
            k_vec = self.k_vec;
            k = self.k;
            n = self.n;
            Xi_estimates = zeros(n,n,4);
            Phi_estimates = zeros(k_vec,n,4);
            beta_estimates = zeros(k,n,4);
            Xi_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_Xi, credibility_level);
            Xi_estimates(:,:,4) = std(mcmc_Xi,0,3);
            Phi_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_Phi, credibility_level);
            Phi_estimates(:,:,4) = std(mcmc_Phi,0,3);
            beta_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_beta, credibility_level);
            beta_estimates(:,:,4) = std(mcmc_beta,0,3);
            Sigma_estimates = quantile(mcmc_Sigma,0.5,3);
            self.Xi_estimates = Xi_estimates;
            self.Phi_estimates = Phi_estimates;
            self.beta_estimates = beta_estimates;
            self.Sigma_estimates = Sigma_estimates;
            % model-specific parameters
            if self.prior_type == 1 && self.error_correction_type == 2
                r = self.r;
                mcmc_K = self.mcmc_K;
                K_estimates = zeros(n,r,4);
                K_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_K, credibility_level);
                K_estimates(:,:,4) = std(mcmc_K,0,3);
                mcmc_Lamda = self.mcmc_Lamda;
                Lamda_estimates = zeros(n,r,4);
                Lamda_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_Lamda, credibility_level);
                Lamda_estimates(:,:,4) = std(mcmc_Lamda,0,3);
                self.K_estimates = K_estimates;
                self.Lamda_estimates = Lamda_estimates;
            elseif self.prior_type == 2 && self.error_correction_type == 1
                mcmc_tau_2 = self.mcmc_tau_2;
                tau_2_estimates = quantile(mcmc_tau_2,0.5);
                mcmc_psi_2 = self.mcmc_psi_2;
                psi_2_estimates = zeros(n,n,4);
                psi_2_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_psi_2, credibility_level);
                psi_2_estimates(:,:,4) = std(mcmc_psi_2,0,3);
                self.tau_2_estimates = tau_2_estimates;
                self.psi_2_estimates = psi_2_estimates;
            elseif self.prior_type == 2 && self.error_correction_type == 2
                r = self.r;
                mcmc_tau_2 = self.mcmc_tau_2;
                tau_2_estimates = quantile(mcmc_tau_2,0.5,2);
                mcmc_psi_2 = self.mcmc_psi_2;
                psi_2_estimates = quantile(mcmc_psi_2,0.5,2);
                mcmc_zeta_2 = self.mcmc_zeta_2;
                zeta_2_estimates = quantile(mcmc_zeta_2,0.5,2);
                mcmc_K = self.mcmc_K;
                K_estimates = zeros(n,r,4);
                K_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_K, credibility_level);
                K_estimates(:,:,4) = std(mcmc_K,0,3);
                mcmc_Lamda = self.mcmc_Lamda;
                Lamda_estimates = zeros(n,r,4);
                Lamda_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_Lamda, credibility_level);
                Lamda_estimates(:,:,4) = std(mcmc_Lamda,0,3);
                self.tau_2_estimates = tau_2_estimates;
                self.psi_2_estimates = psi_2_estimates;
                self.zeta_2_estimates = zeta_2_estimates;
                self.K_estimates = K_estimates;
                self.Lamda_estimates = Lamda_estimates;
            elseif self.prior_type == 3 && self.error_correction_type == 1 
                mcmc_delta = self.mcmc_delta;
                delta_estimates = reshape(mean(mcmc_delta,2),[n n])';  
                self.delta_estimates = delta_estimates;
            elseif self.prior_type == 3 && self.error_correction_type == 2
                r = self.r;
                mcmc_chi = self.mcmc_chi;
                chi_estimates = reshape(mean(mcmc_chi,2),[n r]);              
                mcmc_K = self.mcmc_K;
                K_estimates = zeros(n,r,4);
                K_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_K, credibility_level);
                K_estimates(:,:,4) = std(mcmc_K,0,3);
                mcmc_gamma = self.mcmc_gamma;
                gamma_estimates = reshape(mean(mcmc_gamma,2),[r n])';            
                mcmc_Lamda = self.mcmc_Lamda;
                Lamda_estimates = zeros(n,r,4);
                Lamda_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_Lamda, credibility_level);
                Lamda_estimates(:,:,4) = std(mcmc_Lamda,0,3);
                self.chi_estimates = chi_estimates;
                self.K_estimates = K_estimates;
                self.gamma_estimates = gamma_estimates;
                self.Lamda_estimates = Lamda_estimates ;   
            end
        end


    end
    

end
