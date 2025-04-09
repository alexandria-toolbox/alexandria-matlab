classdef BayesianProxySvar < handle & VectorAutoRegression & BayesianVar
    

    % Bayesian proxy-SVAR, developed in section 14.5
    % 
    % Parameters:
    % -----------
    % endogenous : matrix of size (n_obs,n_endogenous)
    %     endogenous variables, defined in (4.11.1)
    % 
    % proxys : matrix of size (n_obs,n_proxys)
    %     proxy variables, defined in (4.14.42)
    %     
    % exogenous : matrix of size (n_obs,n_exogenous), default = []
    %     exogenous variables, defined in (4.11.1)
    %
    % structural_identification : int, default = 1
    %     structural identification scheme, additional to proxy-SVAR
    %     1 = none, 4 = restrictions
    % 
    % restriction_table : matrix
    %     numerical matrix of restrictions for structural identification
    %
    % lamda : float, default = 0.2
    %     relevance parameter, defined in (4.14.54)
    %  
    % proxy_prior : int, default = 1
    %     prior scheme for normal-generalized-normal prior
    %     1 = uninformative, 2 = Minnesota
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
    % endogenous : matrix of size (n_obs,n_endogenous)
    %     endogenous variables, defined in (4.11.1)
    % 
    % proxys : matrix of size (n_obs,n_proxys)
    %     endogenous variables, defined in (4.14.42)
    % 
    % exogenous : matrix of size (n_obs,n_exogenous)
    %     exogenous variables, defined in (4.11.1)
    %
    % structural_identification : int
    %     structural identification scheme, additional to proxy-SVAR
    %     1 = none, 4 = restrictions
    % 
    % restriction_table : matrix
    %     numerical matrix of structural identification restrictions
    % 
    % lamda : float, default = 0.2
    %     relevance parameter, defined in (4.14.54)
    %  
    % proxy_prior : int, default = 1
    %     prior scheme for normal-generalized-normal prior
    %     1 = uninformative, 2 = Minnesota
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
    % pi4 : float
    %     exogenous slackness hyperparameter, defined in (4.11.19)    
    %
    % credibility_level : float
    %     VAR model credibility level (between 0 and 1)
    %
    % iterations : int
    %     number of Gibbs sampler replications   
    % 
    % burnin : int
    %     number of Gibbs sampler burn-in replications   
    %     
    % verbose : bool
    %     if true, displays a progress bar      
    % 
    % R : matrix of size (T,h)
    %     in-sample matrix of proxy variables, defined in (4.14.42) 
    % 
    % Y_bar : matrix of size (T,n+h)
    %     in-sample matrix of endogenous and proxy regressors, defined in (4.14.49) 
    % 
    % X_bar : matrix of size (T,m+(n+h)*p)
    %     in-sample matrix of full regressors, defined in (4.14.49) 
    %     
    % h : int
    %     number of proxy variables
    %
    % n_bar : int
    %     number of endogenous and proxy variables (n+h)
    %
    % k_bar : int
    %     total number of proxy-SVAR regressors (m+(n+h)*p)
    %    
    % alpha : float
    %     prior degrees of freedom, defined in (4.14.50)
    % 
    % inv_W : matrix of size (k_bar,k_bar)
    %     prior variance of VAR coefficients, defined in (4.14.50)           
    %
    % B : matrix of size (k_bar,n_bar)
    %     prior mean of VAR coefficients, defined in (4.14.50)           
    %
    % S : matrix of size (n_bar,n_bar)
    %     prior scale matrix, defined in (4.14.50) 
    % 
    % alpha_bar : float
    %     posterior degrees of freedom, defined in (4.11.33)
    %
    % W_bar : matrix of size (k_bar,k_bar)
    %     posterior variance of VAR coefficients, defined in (4.14.50)           
    %
    % B_bar : matrix of size (k_bar,n_bar)
    %     posterior mean of VAR coefficients, defined in (4.14.50)           
    %
    % S_bar : matrix of size (n_bar,n_bar)
    %     posterior scale matrix, defined in (4.14.50) 
    % 
    % mcmc_beta : matrix of size (k,n,iterations)
    %     MCMC values of VAR coefficients   
    % 
    % mcmc_Sigma : matrix of size (n,n,iterations)
    %     MCMC values of residual variance-covariance matrix
    %  
    % mcmc_V : matrix of size (h,h,iterations)
    %     MCMC values of covariance matrix V, defined in (4.14.46)
    %  
    % mcmc_min_eigenvalue : matrix of size (iterations,1)
    %     MCMC values of minimum eigenvalue, defined in (4.14.54)
    %  
    % beta_estimates : matrix of size (k,n,3)
    %     estimates of VAR coefficients
    %     page 1: median, page 2: st dev,  page 3: lower bound, page 4: upper bound
    %
    % Sigma_estimates : matrix of size (n,n)
    %     estimates of variance-covariance matrix of VAR residuals
    % 
    % V_estimates : matrix of size (h,h)
    %     posterior estimates of covariance matrix V, defined in (4.14.46)
    % 
    % min_eigenvalue_estimates : float
    %     posterior estimate of minimum eigenvalue, defined in (4.14.54)
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
        proxys
        exogenous
        structural_identification
        restriction_table
        lamda
        proxy_prior
        lags
        constant
        trend
        quadratic_trend
        ar_coefficients
        pi1
        pi3
        pi4        
        credibility_level
        iterations
        burnin
        verbose
        R
        Y_bar
        X_bar
        h
        n_bar
        k_bar
        alpha
        inv_W
        B
        S
        alpha_bar
        W_bar
        B_bar
        S_bar
        mcmc_beta
        mcmc_Sigma
        mcmc_V
        mcmc_min_eigenvalue
        beta_estimates
        Sigma_estimates
        V_estimates
        min_eigenvalue_estimates
    end
    
    
    properties (GetAccess = public, SetAccess = {?BayesianVar})
        mcmc_H
        mcmc_Gamma
    end
    
    
    properties (GetAccess = {?BayesianVar}, SetAccess = private)
        mcmc_chol_Sigma
    end    
    
    
    properties (GetAccess = private, SetAccess = private)
        delta_bar
        s_bar
        inv_W_bar
        U
        V
        K
        P
        F
        FU
        C
        x_dim
    end
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods (Access = public)
        
        
        function self = BayesianProxySvar(endogenous, proxys, varargin)
            
            % constructor for the BayesianProxySvar class
            
            % allow for optional arguments
            parser = inputParser;
            default_exogenous = [];
            default_structural_identification = 1;
            default_restriction_table = [];
            default_lamda = 0.2;
            default_proxy_prior = 1;
            default_lags = 4;
            default_constant = true;
            default_trend = false;
            default_quadratic_trend = false;
            default_ar_coefficients = 0.95;
            default_pi1 = 0.1;
            default_pi3 = 1;
            default_pi4 = 100;            
            default_credibility_level = 0.95;
            default_iterations = 2000;
            default_burnin = 1000;
            default_verbose = false;            
            addRequired(parser, 'endogenous');
            addRequired(parser, 'proxys');            
            addParameter(parser, 'exogenous', default_exogenous);
            addParameter(parser, 'structural_identification', default_structural_identification);
            addParameter(parser, 'restriction_table', default_restriction_table);
            addParameter(parser, 'lamda', default_lamda);
            addParameter(parser, 'proxy_prior', default_proxy_prior);
            addParameter(parser, 'lags', default_lags);
            addParameter(parser, 'constant', default_constant);
            addParameter(parser, 'trend', default_trend);
            addParameter(parser, 'quadratic_trend', default_quadratic_trend);  
            addParameter(parser, 'ar_coefficients', default_ar_coefficients);  
            addParameter(parser, 'pi1', default_pi1);   
            addParameter(parser, 'pi3', default_pi3);  
            addParameter(parser, 'pi4', default_pi4);            
            addParameter(parser, 'credibility_level', default_credibility_level);
            addParameter(parser, 'iterations', default_iterations);
            addParameter(parser, 'burnin', default_burnin);
            addParameter(parser, 'verbose', default_verbose);            
            parse(parser, endogenous, proxys, varargin{:});
            self.endogenous = endogenous;
            self.proxys = proxys;
            self.exogenous = parser.Results.exogenous;
            self.structural_identification = parser.Results.structural_identification;
            self.restriction_table = parser.Results.restriction_table;
            self.lamda = parser.Results.lamda;
            self.proxy_prior = parser.Results.proxy_prior;
            self.lags = parser.Results.lags;
            self.constant = parser.Results.constant;
            self.trend = parser.Results.trend;
            self.quadratic_trend = parser.Results.quadratic_trend;
            self.ar_coefficients = parser.Results.ar_coefficients;
            self.pi1 = parser.Results.pi1;
            self.pi3 = parser.Results.pi3;
            self.pi4 = parser.Results.pi4;            
            self.credibility_level = parser.Results.credibility_level;
            self.iterations = parser.Results.iterations;
            self.burnin = parser.Results.burnin;
            self.verbose = parser.Results.verbose;
            % make regressors
            self.make_regressors();
            % make proxy SVAR regressors
            self.make_proxy_svar_regressors();
            % make delta
            self.make_delta();
            % make individual residual variance
            self.individual_ar_variances();
            % complete with proxy elements
            self.proxy_delta_and_ar_variances();          
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
            % define posterior parameters
            self.posterior();
            % complete with proxy elements
            self.orthogonal_triangular_block_parameters();
            % run MCMC algorithm (Gibbs sampling) for proxy SVAR parameters
            self.parameter_mcmc();
            % obtain posterior estimates for regression parameters
            self.parameter_estimates();         
        end         
        
        
    end
           
        
    methods (Access = protected, Hidden = true)
        
        
        function make_proxy_svar_regressors(self)
            
            % generates proxy regressors Y_bar and X_bar as defined in (4.14.49)
            
            R = self.proxys(self.lags+1:end,:);
            Y_bar = [self.Y R];
            X_1 = vu.generate_intercept_and_trends(self.constant, self.trend, self.quadratic_trend, self.T, 0);
            X_2 = vu.generate_exogenous_regressors(self.exogenous, self.lags);
            X_3 = vu.generate_lagged_endogenous([self.endogenous,self.proxys], self.lags);
            X_bar = [X_1,X_2,X_3];
            h = size(self.proxys,2);
            n_bar = self.n + h;
            k_bar = size(X_bar,2);
            self.R = R;
            self.Y_bar = Y_bar;
            self.X_bar = X_bar;
            self.h = h;
            self.n_bar = n_bar;
            self.k_bar = k_bar;
        end        

        
        function proxy_delta_and_ar_variances(self)      

            % generates proxy elements for delta and AR variances

            delta_bar = [self.delta; 0.9 * ones(self.h,1)];
            s = zeros(self.h,1);
            for i=1:self.h
                ar = MaximumLikelihoodVar(self.proxys(:,i), 'lags', self.lags);
                ar.estimate();
                s(i) = ar.Sigma;
            end
            s_bar = [self.s; s];
            self.delta_bar = delta_bar;
            self.s_bar = s_bar;
        end
        
        
        function prior(self)
        
            % creates prior elements alpha, inv_W, B and S defined in (4.14.54)

            % if prior is naive, set all elements to uninformative
            if self.proxy_prior == 1
                alpha = self.n_bar;
                inv_W = zeros(self.k_bar,1);
                B = zeros(self.k_bar,self.n_bar);
                S = zeros(self.n_bar,self.n_bar);
            % if prior is Minnesota, implement parameters that replicate the normal Wishart prior
            elseif self.proxy_prior == 2
                alpha = 0;
                inv_W = 1 ./ vu.make_W(self.s_bar, self.pi1, self.pi3, self.pi4, self.n_bar, self.m, self.p);
                B = vu.make_B(self.delta_bar, self.n_bar, self.m, self.p);
                S = vu.make_S(self.s_bar);
            end
            self.alpha = alpha;
            self.inv_W = inv_W;
            self.B = B;
            self.S = S;
        end

        
        function posterior(self)

            % creates posterior elements alpha_bar, W_bar, B_bar and S_bar defined in (4.14.54)   

            alpha_bar = self.alpha + self.T;
            inv_W_bar = diag(self.inv_W) + self.X_bar' * self.X_bar;
            W_bar = la.invert_spd_matrix(inv_W_bar);
            B_bar = W_bar * (diag(self.inv_W) * self.B + self.X_bar' * self.Y_bar);
            S_bar = diag(self.S) + self.Y_bar' * self.Y_bar + self.B' * diag(self.inv_W) * self.B ...
            - B_bar' * inv_W_bar * B_bar;
            self.alpha_bar = alpha_bar;
            self.W_bar = W_bar;
            self.inv_W_bar = inv_W_bar;
            self.B_bar = B_bar;
            self.S_bar = S_bar;
        end

        
        function orthogonal_triangular_block_parameters(self)

            % creates U, V, H, P, Q and F defined in algorithm 14.9, along with associated elements

            U = cell(self.n_bar,1);
            V = cell(self.n_bar,1);
            K = cell(self.n_bar,1);
            P = cell(self.n_bar,1);
            F = cell(self.n_bar,1);
            FU = cell(self.n_bar,1);
            C = cell(self.n_bar,1);
            for j=1:self.n_bar
                U_j = eye(self.n_bar,j);            
                if j <= self.n
                    V_j = blkdiag(eye(self.m),kron(eye(self.p),eye(self.n_bar,self.n)));
                else
                    V_j = eye(self.k_bar);
                end            
                inv_H_j = V_j' * self.inv_W_bar * V_j;
                H_j = la.invert_spd_matrix(inv_H_j);
                K_j = la.cholesky_nspd(H_j);
                P_j = H_j * V_j' * self.inv_W_bar * self.B_bar * U_j;
                Q_j = self.alpha_bar * la.invert_spd_matrix(U_j' * self.S_bar * U_j + ...
                      U_j' * self.B_bar' * self.inv_W_bar * self.B_bar * U_j - P_j' * inv_H_j * P_j);            
                F_j= la.cholesky_nspd(Q_j);
                FU_j = F_j' * U_j';
                C_j = setdiff(1:self.n_bar,j);
                U{j} = U_j;
                V{j} = V_j;
                K{j} = K_j;
                P{j} = P_j;
                F{j} = F_j;
                FU{j} = FU_j;
                C{j} = C_j;
            end
            z = zeros(self.n,1);
            z(1:(self.n-self.h)) = self.h;
            x_dim = ones(self.n,1) + self.n - (1:self.n)' - z;
            self.U = U;
            self.V = V;
            self.K = K;
            self.P = P;
            self.F = F;
            self.FU = FU;
            self.C = C;
            self.x_dim = x_dim;            
        end
        
        
        function parameter_mcmc(self)

            % Gibbs sampler for proxy SVAR parameters H_bar_0 and H_bar_+, following algorithm 14.12

            % step 1: posterior parameters
            U_j = self.U;
            V_j = self.V;
            K_j = self.K;
            P_j = self.P;
            F_j = self.F;
            FU_j = self.FU;
            C_j = self.C;
            x_dim = self.x_dim;
            
            % unpack other parameters
            Y = self.Y;
            X = self.X;
            n = self.n;
            h = self.h;
            p = self.p;
            k = self.k;
            n_bar = self.n_bar;
            k_bar = self.k_bar;
            alpha_bar = self.alpha_bar;
            lamda = self.lamda;
            structural_identification = self.structural_identification;
            restriction_table = self.restriction_table;        
            iterations = self.iterations;
            burnin = self.burnin;
            verbose = self.verbose;     
            
            % preallocate storage space
            mcmc_beta = zeros(k,n,iterations);
            mcmc_Sigma = zeros(n,n,iterations);        
            mcmc_chol_Sigma = zeros(n,n,iterations);
            mcmc_H = zeros(n,n,iterations);
            mcmc_inv_H = zeros(n,n,iterations);
            mcmc_V = zeros(h,h,iterations);
            mcmc_min_eigenvalue = zeros(iterations,1);            
            
            % create matrices of restriction and checks, if applicable
            [restriction_matrices covariance_restriction_matrices max_irf_period max_zero_irf_period no_zero_restrictions ...
            shock_history_restrictions] = self.make_restriction_matrices_and_checks(restriction_table, p);

            % step 2: set initial value for Lambda_0
            Lambda_bar_0 = eye(n_bar);        
            
            % iterate over iterations
            iteration = 1;
            while iteration <= (burnin + iterations) 
                
                % step 3: sample Lambda_bar_0 and Lambda_bar_plus
                [Lambda_bar_0 inv_Lambda_bar_0 Lambda_bar_plus] = self.draw_Lambda_bar_0_and_Lambda_bar_plus(U_j, V_j, ...
                                               K_j, P_j, F_j, FU_j, C_j, Lambda_bar_0, alpha_bar, n_bar, k_bar);      
                
                % for next steps: recover reduced-form parameters
                [B Sigma chol_Sigma] = self.recover_reduced_form_parameters(Lambda_bar_0, ...
                                       inv_Lambda_bar_0, Lambda_bar_plus, n, V_j{1});                
                
                % step 4: sample Q
                [Q] = self.draw_Q(inv_Lambda_bar_0, chol_Sigma, B, n, h, p, x_dim, ...
                      no_zero_restrictions, max_zero_irf_period, restriction_matrices);              

                % step 5: obtain SVAR parameters
                [H_bar_0 H_bar_plus H, inv_H] = self.get_svar_parameters(Lambda_bar_0, inv_Lambda_bar_0, Lambda_bar_plus, Q, n);

                % step 6: check relevance
                [V min_eigenvalue relevance_satisfied] = self.check_relevance(H_bar_0, inv_Lambda_bar_0, Q, lamda, n, h);
                if ~ relevance_satisfied
                    continue  
                end
                  
                % step 7: check restrictions
                if structural_identification == 4
                    [restriction_satisfied] = self.check_restrictions(restriction_matrices, covariance_restriction_matrices, ...
                                              shock_history_restrictions, max_irf_period, Y, X, B, H, inv_H, V, n, h, p);
                    if ~ restriction_satisfied
                        continue                
                    end
                end
                  
                % save if burn is exceeded
                if iteration > burnin
                    mcmc_beta(:,:,iteration-burnin) = B;
                    mcmc_Sigma(:,:,iteration-burnin) = Sigma;
                    mcmc_chol_Sigma(:,:,iteration-burnin) = chol_Sigma;
                    mcmc_H(:,:,iteration-burnin) = H;
                    mcmc_inv_H(:,:,iteration-burnin) = inv_H;
                    mcmc_V(:,:,iteration-burnin) = V;
                    mcmc_min_eigenvalue(iteration-burnin) = min_eigenvalue;
                end
                
                % display progress bar
                if verbose
                    cu.progress_bar(iteration, iterations+burnin, 'Model parameters:') 
                end
                iteration = iteration + 1;      
            end
            
            % save as attributes
            self.mcmc_beta = mcmc_beta;
            self.mcmc_Sigma = mcmc_Sigma;
            self.mcmc_chol_Sigma = mcmc_chol_Sigma;
            self.mcmc_H = mcmc_H;
            self.mcmc_inv_H = mcmc_inv_H;
            self.mcmc_Gamma = ones(iterations,n);
            self.mcmc_V = mcmc_V;
            self.mcmc_min_eigenvalue = mcmc_min_eigenvalue;
            self.svar_index = (1:iterations)';
        end
        

        function [restriction_matrices covariance_restriction_matrices max_irf_period max_zero_irf_period ...
                  no_zero_restrictions shock_history_restrictions] = make_restriction_matrices_and_checks(self, restriction_table, p)

            % elements for later check of restrictions

            if isempty(restriction_table)
                restriction_matrices = [];
                covariance_restriction_matrices = [];
                max_irf_period = 0;
                max_zero_irf_period = 0;
                no_zero_restrictions = true;
                shock_history_restrictions = false;
            else
                [restriction_matrices max_irf_period] = vu.make_restriction_matrices(restriction_table, p);
                [covariance_restriction_matrices] = vu.make_covariance_restriction_matrices(restriction_table);
                [no_zero_restrictions max_zero_irf_period shock_history_restrictions] = self.make_restriction_checks(restriction_table);
            end
        end


        function [no_zero_restrictions max_zero_irf_period shock_history_restrictions] = make_restriction_checks(self, restriction_table)

            % boolean for zero/shock resrictions and max IRF period for zero restrictions

            zero_restriction_table = restriction_table(restriction_table(:,1)==1,:);
            zero_restriction_number = size(zero_restriction_table,1);        
            no_zero_restrictions = (zero_restriction_number == 0);
            if no_zero_restrictions
                max_zero_irf_period = 0;
            else
                max_zero_irf_period = max(zero_restriction_table(:,3));
            end
            shock_restriction_table = restriction_table(restriction_table(:,1)==3,:);
            history_restriction_table = restriction_table(restriction_table(:,1)==4,:);
            shock_history_restriction_table = [shock_restriction_table;history_restriction_table];
            shock_history_restriction_number = size(shock_history_restriction_table,1);    
            shock_history_restrictions = (shock_history_restriction_number ~= 0); 
        end

        
        function [Lambda_bar_0, inv_Lambda_bar_0, Lambda_bar_plus] = draw_Lambda_bar_0_and_Lambda_bar_plus(self, U, V, ...
                 K, P, F, FU, C, Lambda_bar_0, alpha_bar, n_bar, k_bar)

            % generation of proxy SVAR parameters Lambda_bar_0 and Lambda_bar_plus, following algorithm 14.9

            Lambda_bar_plus = zeros(n_bar,k_bar);
            for j=1:n_bar
                % step 3
                temp = Lambda_bar_0(C{j},:);
                kernel = null(temp);
                z = kernel(:,1);                
                % step 4
                w_1 = FU{j} * z;
                w_1 = w_1 / norm(w_1);               
                % step 5
                w = zeros(j,j);
                w(:,1) = w_1;
                for i=2:j
                    w_i = zeros(j,1);
                    c_i = w_1(1:i)' * w_1(1:i);
                    c_i_1 = w_1(1:i-1)' * w_1(1:i-1);
                    w_i(i) = - c_i_1;
                    w_i(1:i-1) = w_1(1:i-1) * w_1(i);
                    w_i = w_i / sqrt(c_i * c_i_1);
                    w(:,i) = w_i;
                end
                % step 6
                beta = zeros(j,1);
                s = randn(alpha_bar+1,1) / sqrt(alpha_bar);
                r = s' * s;
                beta(1) = (1 - 2 * (rand>0.5)) * sqrt(r);              
                % step 7
                beta(2:end) = randn(j-1,1) / sqrt(alpha_bar);
                % step 8
                gamma_0 = F{j} * (w * beta);
                gamma_0 = sign(gamma_0(j)) * gamma_0;            
                % step 9
                gamma_plus = K{j} * randn(size(P{j},1),1) + P{j} * gamma_0;
                % step 10
                lambda_0 = U{j} * gamma_0;
                lambda_plus = V{j} * gamma_plus;
                Lambda_bar_0(j,:) = lambda_0;
                Lambda_bar_plus(j,:) = lambda_plus;
            end
            inv_Lambda_bar_0 = la.invert_lower_triangular_matrix(Lambda_bar_0);
        end     

        
        function [B Sigma chol_Sigma] = recover_reduced_form_parameters(self, Lambda_bar_0, inv_Lambda_bar_0, Lambda_bar_plus, n, V_j)

            % recover reduced-form parameters

            inv_Lambda_0 = inv_Lambda_bar_0(1:n,1:n);
            chol_Sigma = inv_Lambda_0;
            Sigma = chol_Sigma * chol_Sigma';
            Lambda_plus = Lambda_bar_plus(1:n,:) * V_j;
            B = (inv_Lambda_0 * Lambda_plus)';
        end        
        
        
        function [Q] = draw_Q(self, inv_Lambda_bar_0, chol_Sigma, B, n, h, p, x_dim, no_zero_restrictions, max_zero_irf_period, restriction_matrices)

            % Q rotation matrix, with or without additional zero restrictions

            if no_zero_restrictions
                Q = self.draw_regular_Q(inv_Lambda_bar_0, x_dim, n, h);
            else
                irf = vu.impulse_response_function(B, n, p, max_zero_irf_period);
                structural_irf = vu.structural_impulse_response_function(irf, chol_Sigma, n);
                Q = self.draw_zero_restriction_Q(inv_Lambda_bar_0, x_dim, n, h, restriction_matrices{1,1}, structural_irf);
            end
        end       
        
        
        function [Q] = draw_regular_Q(self, inv_Lambda_bar_0, x_dim, n, h)

            % generation of proxy SVAR parameter Q, following algorithm 14.10    

            G = inv_Lambda_bar_0(n+1:end,1:n);
            % initiate Q_1 and iterate
            Q_1 = zeros(n,n); 
            for j=1:n
                % step 1
                x1_j = randn(x_dim(j),1);
                w1_j = x1_j / norm(x1_j);           
                if j <= n-h
                    M1_j = [Q_1(:,1:j-1) G']';
                else
                    M1_j = Q_1(:,1:j-1)';
                end
                K1_j = null(M1_j);
                Q_1(:,j) = K1_j * w1_j;
            end
            % step 3
            Q_2 = rn.uniform_orthogonal(h);
            % step 4
            Q = blkdiag(Q_1, Q_2)';    
        end
        
        
        function [Q] = draw_zero_restriction_Q(self, inv_Lambda_bar_0, x_dim, n, h, restriction_matrix, irf)

            % generation of proxy SVAR parameter Q, following algorithm 14.12     

            zero_restriction_shocks = unique(restriction_matrix(:,2));
            G = inv_Lambda_bar_0(n+1:end,1:n);
            % initiate Q_1 and iterate
            Q_1 = zeros(n,n);   
            for j=1:n            
                % G_j for zero restrictions
                if ismember(j,zero_restriction_shocks)
                    shock_restrictions = restriction_matrix(restriction_matrix(:,2)==j,:);
                    h_j = size(shock_restrictions,1);
                    G_j = zeros(h_j,n);
                    for i=1:h_j
                        G_j(i,:) = irf(shock_restrictions(i,1),: ,shock_restrictions(i,3));
                    end
                else
                    G_j = [];
                    h_j = 0;
                end
                % step 1
                x1_j = randn(x_dim(j)-h_j,1);
                w1_j = x1_j / norm(x1_j);
                % step 2
                if j <= n-h
                    M1_j = [Q_1(:,1:j-1) G' G_j']';
                else
                    M1_j = [Q_1(:,1:j-1) G_j']';
                end
                K1_j = null(M1_j);
                Q_1(:,j) = K1_j * w1_j;                
            end
            % step 3
            Q_2 = rn.uniform_orthogonal(h);
            % step 4
            Q = blkdiag(Q_1, Q_2)';    
        end

        
        function [H_bar_0 H_bar_plus H inv_H] = get_svar_parameters(self, Lambda_bar_0, inv_Lambda_bar_0, Lambda_bar_plus, Q, n)

            % SVAR parameters from orthogonal triangular-block parameterization

            H_bar_0 = Q * Lambda_bar_0;
            H_bar_plus = Q * Lambda_bar_plus;
            inv_Lambda_0 = inv_Lambda_bar_0(1:n,1:n);
            Q_1 = Q(1:n,1:n);
            H = inv_Lambda_0 * Q_1';
            inv_H = H_bar_0(1:n,1:n);
        end   
        
        
        function [V min_eigenvalue relevance_satisfied] = check_relevance(self, H_bar_0, inv_Lambda_bar_0, Q, lamda, n, h)
        
            % check relevance conditions

            % create Gamma_01, inv_Gamma_02 and inv_H_0
            Gamma_01 = H_bar_0(n+1:end,1:n);
            Q_1 = Q(1:n,1:n);
            Q_2 = Q(n+1:end,n+1:end);
            inv_Lambda_0 = inv_Lambda_bar_0(1:n,1:n);
            inv_Lambda_02 = inv_Lambda_bar_0(n+1:end,n+1:end);            
            inv_H_0 = inv_Lambda_0 * Q_1';
            inv_Gamma_02 = inv_Lambda_02 * Q_2';            
            % create V from (4.14.46)
            E_r_xi = - inv_Gamma_02 * Gamma_01 * inv_H_0;
            V = E_r_xi(:,n-h+1:end);
            % relevance matrix P, from (4.14.54)
            VV = V * V';
            P = (inv_Gamma_02 * inv_Gamma_02' + VV) \ VV;
            % minimum eigenvalue
            eigenvalues = eig(P);
            min_eigenvalue = min(eigenvalues);
            if min_eigenvalue > lamda
                relevance_satisfied = true;
            else
                relevance_satisfied = false;
            end
        end
        
        
        function [restriction_satisfied] = check_restrictions(self, restriction_matrices, covariance_restriction_matrices, shock_history_restrictions, max_irf_period, Y, X, B, H, inv_H, V, n, h, p)

            % check of all structural identification restrictions

            % if no restriction, stop restriction check
            if isempty(restriction_matrices)
                restriction_satisfied = true;
                return
            end
            % preliminary IRFs and shocks
            if max_irf_period ~= 0
                irf = vu.impulse_response_function(B, n, p, max_irf_period);
            end
            if shock_history_restrictions     
                [E ~] = vu.fit_and_residuals(Y, X, B);
            end
            % check restrictions: IRF, sign
            irf_sign_index = restriction_matrices{2,1};
            if ~isempty(irf_sign_index)
                irf_sign_coefficients = restriction_matrices{2,2};
                restriction_satisfied = vu.check_irf_sign(irf_sign_index, irf_sign_coefficients, irf, H);
                if ~restriction_satisfied
                    return
                end
            end
            % check restrictions: IRF, magnitude
            irf_magnitude_index = restriction_matrices{3,1};
            if ~isempty(irf_magnitude_index)
                irf_magnitude_coefficients = restriction_matrices{3,2};
                restriction_satisfied = vu.check_irf_magnitude(irf_magnitude_index, irf_magnitude_coefficients, irf, H);
                if ~restriction_satisfied
                    return
                end
            end
            % check restrictions: structural shocks, sign
            shock_sign_index = restriction_matrices{4,1};
            if ~isempty(shock_sign_index)
                shock_sign_coefficients = restriction_matrices{4,2};
                restriction_satisfied = vu.check_shock_sign(shock_sign_index, shock_sign_coefficients, E, inv_H');
                if ~restriction_satisfied
                    return
                end
            end
            % check restrictions: structural shocks, magnitude
            shock_magnitude_index = restriction_matrices{5,1};
            if ~isempty(shock_magnitude_index)
                shock_magnitude_coefficients = restriction_matrices{5,2};
                restriction_satisfied = vu.check_shock_magnitude(shock_magnitude_index, shock_magnitude_coefficients, E, inv_H');
                if ~restriction_satisfied
                    return
                end
            end
            % historical decomposition values if any of sign or magnitude restrictions apply
            history_sign_index = restriction_matrices{6,1};
            history_magnitude_index = restriction_matrices{7,1};
            if ~isempty(history_sign_index) || ~isempty(history_magnitude_index)
                structural_irf = vu.structural_impulse_response_function(irf, H, n);
                structural_shocks = E * inv_H';
            end
            % check restrictions: historical decomposition, sign
            if ~isempty(history_sign_index)
                history_sign_coefficients = restriction_matrices{6,2};
                restriction_satisfied = vu.check_history_sign(history_sign_index, ...
                                        history_sign_coefficients, structural_irf, structural_shocks);
                if ~restriction_satisfied 
                    return
                end
            end
            % check restrictions: historical decomposition, magnitude
            if ~isempty(history_magnitude_index)
                history_magnitude_coefficients = restriction_matrices{7,2};
                restriction_satisfied = vu.check_history_magnitude(history_magnitude_index, ...
                                        history_magnitude_coefficients, structural_irf, structural_shocks);
                if ~restriction_satisfied
                    return
                end
            end
            % check restrictions: covariance, sign
            covariance_sign_index = covariance_restriction_matrices{1,1};
            if ~isempty(covariance_sign_index)
                covariance_sign_coefficients = covariance_restriction_matrices{1,2};
                restriction_satisfied = vu.check_covariance_sign(covariance_sign_index, covariance_sign_coefficients, V, n, h);
                if ~restriction_satisfied
                    return
                end
            end
            % check restrictions: covariance, magnitude
            covariance_magnitude_index = covariance_restriction_matrices{2,1};
            if ~isempty(covariance_magnitude_index)
                covariance_magnitude_coefficients = covariance_restriction_matrices{2,2};
                restriction_satisfied = vu.check_covariance_magnitude(covariance_magnitude_index, covariance_magnitude_coefficients, V, n, h);
                if ~restriction_satisfied
                    return
                end
            end
            restriction_satisfied = true;
        end

        
        function parameter_estimates(self)

            % point estimates and credibility intervals for model parameters
            % use empirical quantiles from MCMC algorithm

            % unpack
            mcmc_beta = self.mcmc_beta;
            mcmc_Sigma = self.mcmc_Sigma;
            mcmc_H = self.mcmc_H;
            mcmc_Gamma = self.mcmc_Gamma;
            mcmc_V = self.mcmc_V;
            mcmc_min_eigenvalue = self.mcmc_min_eigenvalue;
            k = self.k;
            n = self.n;
            credibility_level = self.credibility_level;
            % initiate and fill storage for beta, 4 columns: lower bound, median, upper bound, standard deviation
            beta_estimates = zeros(k,n,4);
            beta_estimates(:,:,1:3) = vu.posterior_estimates(mcmc_beta, credibility_level);
            beta_estimates(:,:,4) = std(mcmc_beta,0,3);
            % other estimates for Sigma, H, Gamma
            Sigma_estimates = quantile(mcmc_Sigma,0.5,3);
            H_estimates = quantile(mcmc_H,0.5,3);
            Gamma_estimates = quantile(mcmc_Gamma,0.5);
            % proxy SVAR specific estimates
            V_estimates = quantile(mcmc_V,0.5,3);
            min_eigenvalue_estimates = quantile(mcmc_min_eigenvalue,0.5);
            % save as attributes
            self.beta_estimates = beta_estimates;
            self.Sigma_estimates = Sigma_estimates;
            self.H_estimates = H_estimates;
            self.Gamma_estimates = Gamma_estimates;
            self.V_estimates = V_estimates;
            self.min_eigenvalue_estimates = min_eigenvalue_estimates; 
        end
        
        
    end
    
    
end
        
        