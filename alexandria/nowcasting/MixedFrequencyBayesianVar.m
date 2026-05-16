classdef MixedFrequencyBayesianVar < handle & BayesianVar


    % Mixed frequency Bayesian VAR, developed in chapter 17
    % 
    % Parameters:
    % -----------
    % 
    % endogenous : matrix of size (n_obs,n)
    %     endogenous variables, defined in (6.17.1)
    % 
    % exogenous : matrix of size (n_obs,m), default = []
    %     exogenous variables, defined in (6.17.1)
    % 
    % decomposition : bool
    %     if true, applies frequency decomposition as developed in section 17.3
    % 
    % decomposition_table : ndarray
    %     numerical matrix of frequency decomposition       
    % 
    % structural_identification : int, default = 2
    %     structural identification scheme, as defined in section 13.2
    %     1 = none, 2 = Cholesky, 3 = triangular, 4 = restrictions
    % 
    % restriction_table : ndarray
    %     numerical matrix of restrictions for structural identification
    % 
    % lags : int, default = 4
    %     number of lags, defined in (6.17.1)
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
    %     prior mean delta for AR coefficients, defined in (6.17.8)
    % 
    % pi1 : float, default = 0.1
    %     overall tightness hyperparameter, defined in (6.17.8)
    % 
    % pi2 : float, default = 0.5
    %     cross-variable shrinkage hyperparameter, defined in (6.17.8)
    % 
    % pi3 : float, default = 1
    %     lag decay hyperparameter, defined in (6.17.8)    
    % 
    % pi4 : float, default = 100
    %     exogenous slackness hyperparameter, defined in (6.17.8)             
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
    % Attributes
    % ----------
    % endogenous : matrix of size (n_obs,n)
    %     endogenous variables, defined in (6.17.1)
    % 
    % exogenous : matrix of size (n_obs,m), default = []
    %     exogenous variables, defined in (6.17.1)
    % 
    % decomposition : bool
    %     if true, applies frequency decomposition as developed in section 17.3
    % 
    % decomposition_table : ndarray
    %     numerical matrix of frequency decomposition       
    % 
    % structural_identification : int, default = 2
    %     structural identification scheme, as defined in section 13.2
    %     1 = none, 2 = Cholesky, 3 = triangular, 4 = restrictions
    % 
    % restriction_table : ndarray
    %     numerical matrix of restrictions for structural identification
    % 
    % lags : int, default = 4
    %     number of lags, defined in (6.17.1)
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
    %     prior mean delta for AR coefficients, defined in (6.17.8)
    % 
    % pi1 : float, default = 0.1
    %     overall tightness hyperparameter, defined in (6.17.8)
    % 
    % pi2 : float, default = 0.5
    %     cross-variable shrinkage hyperparameter, defined in (6.17.8)
    % 
    % pi3 : float, default = 1
    %     lag decay hyperparameter, defined in (6.17.8)    
    % 
    % pi4 : float, default = 100
    %     exogenous slackness hyperparameter, defined in (6.17.8)             
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
    % Z : matrix of size (T,m)
    %     exogenous variables, defined in (6.17.3)    
    % 
    % Z_ : matrix of size (T+p,m)
    %     exogenous variables with initial conditions  
    % 
    % n : int
    %     number of endogenous variables, defined in (6.17.1)
    % 
    % m : int
    %     number of exogenous variables, defined in (6.17.1)
    % 
    % p : int
    %     number of lags, defined in (6.17.1)
    % 
    % T : int
    %     number of sample periods, defined in (6.17.1)
    % 
    % k : int
    %     number of VAR coefficients in each equation, defined in (6.17.1)
    % 
    % q : int
    %     total number of VAR coefficients, defined in (6.17.1)    
    % 
    % yo_ : list of len (T+p)
    %     list of observed endogenous variables, defined in (6.17.16)    
    % 
    % L_ : list of len (T+p)
    %     list of selection matrices, defined in (6.17.21)      
    % 
    % b : matrix of size (q,1)
    %     prior mean of VAR coefficients, defined in (6.17.8)
    % 
    % V : matrix of size (q,q)
    %     prior mean of VAR coefficients, defined in (6.17.8)       
    % 
    % alpha : float
    %     prior degrees of freedom, defined in (6.17.9)
    % 
    % S : matrix of size (n,n)
    %     prior scale matrix, defined in (6.17.9) 
    % 
    % alpha_bar : float
    %     posterior degrees of freedom, defined in (6.17.15)      
    % 
    % mcmc_beta : matrix of size (k,n,iterations)
    %     MCMC values of VAR coefficients   
    % 
    % mcmc_Sigma : matrix of size (n,n,iterations)
    %     MCMC values of residual variance-covariance matrix     
    % 
    % mcmc_W : matrix of size (T,n,iterations)
    %     MCMC values of latent endogenous variables    
    % 
    % beta_estimates : matrix of size (k,n,4)
    %     estimates of VAR coefficients
    %     page 1: median, page 2: st dev,  page 3: lower bound, page 4: upper bound
    % 
    % Sigma_estimates : matrix of size (n,n)
    %     estimates of variance-covariance matrix of VAR residuals    
    % 
    % W_estimates : matrix of size (T,n,3)
    %     estimates of VAR coefficients
    %     page 1: median, page 2: lower bound, page 3: upper bound
    % 
    % Y : matrix of size (T,n)
    %     matrix of in-sample endogenous variables, obtained from W
    % 
    % X : matrix of size (T,k)
    %     matrix of exogenous and lagged regressors, defined in (6.17.3)   
    % 
    % mcmc_H :  matrix of size (n,n,iterations)
    %     MCMC values of structural identification matrix, defined in (4.13.5)
    % 
    % mcmc_Gamma : matrix of size (iterations,n)
    %     MCMC values of structural shock variance matrix, defined in definition 13.1    
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
    % insample_evaluation : dict
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
    % forecast_evaluation_criteria : dict
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
        decomposition
        decomposition_table
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
        credibility_level
        iterations
        burnin
        verbose
        Z
        Z_
        n
        m
        p
        T
        k
        q
        yo_
        L_
        b
        V
        alpha
        S
        alpha_bar
        mcmc_beta
        mcmc_Sigma
        mcmc_W
        beta_estimates
        Sigma_estimates
        W_estimates
        Y
        X
    end
    

    properties (GetAccess = public, SetAccess = {?BayesianVar})
        mcmc_H
        mcmc_Gamma
    end


    properties (GetAccess = private, SetAccess = private)
        T_
        r
        inv_V
        inv_V_b
        mcmc_X
        mcmc_inv_Sigma
    end
    

    properties (GetAccess = {?BayesianVar}, SetAccess = private)
        mcmc_chol_Sigma
    end


    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods (Access = public)


        function self = MixedFrequencyBayesianVar(endogenous, varargin)
            
            % constructor for the MixedFrequencyBayesianVar class
            
            % allow for optional arguments
            parser = inputParser;
            default_exogenous = [];
            default_decomposition = false;
            default_decomposition_table = [];
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
            default_credibility_level = 0.95;
            default_iterations = 2000;
            default_burnin = 1000;
            default_verbose = false;
            addRequired(parser, 'endogenous');
            addParameter(parser, 'exogenous', default_exogenous);
            addParameter(parser, 'decomposition', default_decomposition);
            addParameter(parser, 'decomposition_table', default_decomposition_table);            
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
            addParameter(parser, 'credibility_level', default_credibility_level);
            addParameter(parser, 'iterations', default_iterations);
            addParameter(parser, 'burnin', default_burnin);            
            addParameter(parser, 'verbose', default_verbose);
            parse(parser, endogenous, varargin{:});
            self.endogenous = endogenous;
            self.exogenous = parser.Results.exogenous;
            self.decomposition = parser.Results.decomposition;
            self.decomposition_table = parser.Results.decomposition_table;
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
            self.credibility_level = parser.Results.credibility_level;
            self.iterations = parser.Results.iterations;
            self.burnin = parser.Results.burnin;
            self.verbose = parser.Results.verbose;
            % make regressors
            self.make_regressors();
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
    
            % define prior values
            self.prior();
            % define posterior values
            self.posterior();
            % run MCMC algorithm (Gibbs sampling) for VAR parameters
            self.parameter_mcmc();   
            % obtain posterior estimates for regression parameters
            self.parameter_estimates();
            % estimate structural identification
            self.make_structural_identification();
        end


    end


    methods (Access = protected, Hidden = true)


        function make_regressors(self)
            
            % creates regressors and dimensions

            % make exogenous regressors
            [self.Z self.Z_] = self.make_exogenous_regressors();
            % generate dimensions
            [self.n self.m self.p self.T self.k self.q] = self.generate_dimensions();
            % generate delta
            [self.delta] = self.make_delta();
            % generate s
            [self.s] = self.individual_ar_variances();
            % generate state-space regressors
            [self.yo_ self.L_ self.T_ self.r] = self.make_state_space_regressors();           
        end  


        function [Z, Z_] = make_exogenous_regressors(self)
    
            % creates exogenous regressors

            periods = size(self.endogenous,1);
            X_1 = vu.generate_intercept_and_trends(self.constant, self.trend, self.quadratic_trend, periods, 0);   
            X_2 = vu.generate_exogenous_regressors(self.exogenous, 0);
            Z_ = [X_1 X_2];
            Z = Z_(self.lags+1:end,:);
        end


        function [n m p T k q] = generate_dimensions(self)
    
            % creates VAR dimension  
    
            T = size(self.endogenous,1) - self.lags;
            n = size(self.endogenous,2);
            p = self.lags;
            m = self.constant + self.trend + self.quadratic_trend;
            if ~isempty(self.exogenous)
                m = m + size(self.exogenous,2);
            end
            k = m + n * p;
            q = n * k;
        end


        function [delta] = make_delta(self)

            % creates delta hyperparameter
            
            if isscalar(self.ar_coefficients)
                ar_coefficients = repmat(self.ar_coefficients, [self.n 1]);
            else
                ar_coefficients = self.ar_coefficients;
            end
            delta = ar_coefficients;
        end


        function [s] = individual_ar_variances(self)

            % creates individual AR variances
            
            s = zeros(self.n,1);
            for i=1:self.n
                endogenous = self.endogenous(:,i);
                endogenous = endogenous(~isnan(endogenous));
                ar = MaximumLikelihoodVar(endogenous, 'lags', self.lags);
                ar.estimate();
                s(i) = ar.Sigma;
            end
            self.s = s;
        end


        function [yo_ L_ T_ r] = make_state_space_regressors(self)
            
            % creates period-specific parameters for state-space formulation

            % identify maximum lag for state-space sampler
            if self.decomposition
                r = max(self.p, max(self.decomposition_table));
            else
                r = self.p;
            end
            % create full selection matrix
            L = zeros(self.n,self.n*r);
            eye_matrix = eye(self.n);
            zero_matrix = zeros(self.n,self.n);
            if self.decomposition
                decomposition_table = self.decomposition_table;
            else
                decomposition_table = ones(self.n,1);
            end
            for i=1:self.n
                periods = decomposition_table(i);
                temp = [repmat(eye_matrix,[1 periods]) repmat(zero_matrix,[1 r-periods])];
                L(i,:) = temp(i,:);
            end
            % initiate observation and selection matrices
            T_ = size(self.endogenous,1);
            yo_ = cell(T_,1);
            L_ = cell(T_,1);
            % loop over periods to obtain period-specific matrices
            for t=1:T_
                endogenous_t = self.endogenous(t,:);
                non_nan_entries = find(~isnan(endogenous_t));
                yo_{t} = endogenous_t(non_nan_entries)';
                L_{t} = L(non_nan_entries,:);
            end
        end


        function prior(self)
            
            % creates prior elements b and V
            
            self.b = vu.make_b(self.delta, self.n, self.m, self.p);
            self.V = vu.make_V(self.s, self.pi1, self.pi2, self.pi3, self.pi4, self.n, self.m, self.p);
            self.alpha = vu.make_alpha(self.n);
            self.S = vu.make_S(self.s);
        end


        function posterior(self)
            
            % creates posterior elements
            
            % generate preliminary posterior elements
            [inv_V inv_V_b] = vu.make_V_b_inverse(self.b, self.V);
            % generate posterior alpha_bar
            alpha_bar = self.alpha + self.T;
            self.alpha_bar = alpha_bar;
            self.inv_V = inv_V;
            self.inv_V_b = inv_V_b;
        end


        function parameter_mcmc(self)
            
            % Gibbs sampler for VAR parameters beta, Sigma and W, following algorithm 17.1
    
            % unpack
            endogenous = self.endogenous;
            Z = self.Z;
            Z_ = self.Z_;
            inv_V = self.inv_V;
            inv_V_b = self.inv_V_b;
            alpha_bar = self.alpha_bar;
            S = self.S;
            n = self.n;
            m = self.m;
            p = self.p;
            k = self.k;
            T = self.T;
            yo_ = self.yo_;
            L_ = self.L_;
            T_ = self.T_;
            r = self.r;
            iterations = self.iterations;
            burnin = self.burnin;
            verbose = self.verbose;
            
            % preallocate storage space
            mcmc_beta = zeros(k,n,iterations);
            mcmc_Sigma = zeros(n,n,iterations);
            mcmc_W = zeros(T,n,iterations);
            mcmc_X = zeros(T,k,iterations);
            mcmc_chol_Sigma = zeros(n,n,iterations);
            mcmc_inv_Sigma = zeros(n,n,iterations);
            % set initial values
            inv_Sigma = diag(1./S);
            S = diag(S);
            W = randn(T,n);
            X = randn(T,k);

            % iterate over iterations
            iteration = 1;
            while iteration <= (burnin + iterations)

                % step 2: sample beta
                beta = self.draw_beta(inv_V, inv_V_b, inv_Sigma, W, X);
                B = reshape(beta,[k n]);

                % step 3: sample Sigma
                [Sigma inv_Sigma chol_Sigma] = self.draw_Sigma(W, X, B, S, alpha_bar);

                % step 4: sample gamma
                [Gamma W X] = self.draw_gamma(Z, Z_, B, Sigma, yo_, L_, n, m, p, T, T_);

                % save if burn is exceeded
                if iteration > burnin

                    % save parameter values
                    mcmc_beta(:,:,iteration-burnin) = B;
                    mcmc_Sigma(:,:,iteration-burnin) = Sigma;
                    mcmc_W(:,:,iteration-burnin) = W;
                    mcmc_X(:,:,iteration-burnin) = X;
                    mcmc_chol_Sigma(:,:,iteration-burnin) = chol_Sigma;
                    mcmc_inv_Sigma(:,:,iteration-burnin) = inv_Sigma;
                end

                % if verbose, display progress bar
                if verbose
                    cu.progress_bar(iteration, iterations+burnin, 'Model parameters:'); 
                end

                % update iterations    
                iteration = iteration + 1;                
            end

            % save as attributes
            self.mcmc_beta = mcmc_beta;
            self.mcmc_Sigma = mcmc_Sigma;
            self.mcmc_W = mcmc_W;
            self.mcmc_X = mcmc_X;
            self.mcmc_chol_Sigma = mcmc_chol_Sigma;
            self.mcmc_inv_Sigma = mcmc_inv_Sigma;
        end
            
    
        function [beta] = draw_beta(self, inv_V, inv_V_b, inv_Sigma, W, X)
            
            % draw beta from its conditional posterior defined in (6.17.11)
            
            % posterior V_bar
            inv_V_bar = inv_V + kron(inv_Sigma, X' * X);
            % posterior b_bar
            b_bar_temp = inv_V_b + la.vec(X' * W * inv_Sigma);
            % efficient sampling of beta (algorithm 9.4)
            beta = rn.efficient_multivariate_normal(b_bar_temp, inv_V_bar);
        end

    
        function [Sigma inv_Sigma chol_Sigma] = draw_Sigma(self, W, X, B, S, alpha_bar)
            
            % draw Sigma from its conditional posterior defined in (6.17.14)
    
            % compute residuals
            residuals = W - X * B;
            % compute S_bar
            S_bar = S + residuals' * residuals;
            % sample sigma
            Sigma = rn.inverse_wishart(alpha_bar, S_bar);
            inv_Sigma = la.invert_spd_matrix(Sigma);
            chol_Sigma = la.cholesky_nspd(Sigma);
        end


        function [Gamma W X] = draw_gamma(self, Z, Z_, B, Sigma, yo_, L_, n, m, p, T, T_)
            
            % draw Gamma from its conditional posterior, as defined in (6.17.20)-(6.17.21)
            
            % compute state-space parameters
            [F mu_ Upsilon] = nwu.make_mfbvar_state_regressors(Z_, B, Sigma, n, m, p, T_);
            % get initial values for algorithm
            [gamma_00 Upsilon_00] = nwu.mfbvar_kalman_initial_values(Sigma, n, p);
            % run forward pass
            [Gamma_tt Gamma_tt1 Ups_tt Ups_tt1] = ss.mfbvar_forward_pass(yo_, L_, ...
                                   F, mu_, n, p, Upsilon, T_, gamma_00, Upsilon_00);
            % run backward pass
            [Gamma] = ss.mfbvar_backward_pass(Gamma_tt, Gamma_tt1, Ups_tt, Ups_tt1, F, T_, n, p);
            % update state regressors
            [W X] = nwu.update_mfbvar_state_regressors(Z, Gamma, n, p);
        end


        function parameter_estimates(self)
            
            % point estimates and credibility intervals for model parameters
            % use empirical quantiles from MCMC algorithm

            % unpack
            mcmc_beta = self.mcmc_beta;
            mcmc_Sigma = self.mcmc_Sigma;
            mcmc_W = self.mcmc_W;
            mcmc_X = self.mcmc_X;
            credibility_level = self.credibility_level;
            k = self.k;
            n = self.n;
            % initiate storage: 4 columns: lower bound, median, upper bound, standard deviation
            beta_estimates = zeros(k,n,4);
            % fill estimates
            beta_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_beta, credibility_level);
            beta_estimates(:,:,4) = std(mcmc_beta,0,3);  
            Sigma_estimates = quantile(mcmc_Sigma,0.5,3);
            W_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_W, credibility_level);
            X = quantile(mcmc_X,0.5,3);
            Y = W_estimates(:,:,1);
            self.beta_estimates = beta_estimates;
            self.Sigma_estimates = Sigma_estimates;
            self.W_estimates = W_estimates;
            self.Y = Y;
            self.X = X;           
        end 


    end


end
