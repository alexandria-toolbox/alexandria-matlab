classdef BayesianDynamicFactorModel < handle


    % Bayesian dynamic factor model, developed in chapter 18
    % 
    % Parameters:
    % -----------
    % endogenous : matrix of size (n_obs,n)
    %     endogenous variables, defined in (6.18.1)
    % 
    % factors : int, default = 3
    %     number of latent factors, defined in (6.18.1)
    % 
    % loadings_lags : int, default = 2
    %     number of loadings lags, defined in (6.18.2)
    % 
    % factor_lags : int, default = 2
    %     number of factor lags, defined in (6.18.6)
    % 
    % residual_lags : int, default = 1
    %     number of residual lags, defined in (6.18.7)
    % 
    % sigma : float, default = 0.1
    %     variance on residual shock e_it, defined in (6.18.7)        
    % 
    % omega : float, default = 0.1
    %     variance on factor shock xi_t, defined in (6.18.6) 
    % 
    % delta1 : float, default = 0.1
    %     overall tightness hyperparameter on lambda_ij, defined in (6.18.21)
    % 
    % pi1 : float, default = 0.1
    %     overall tightness hyperparameter, defined in (6.18.25)
    % 
    % pi2 : float, default = 0.5
    %     cross-variable shrinkage hyperparameter, defined in (6.18.25)
    % 
    % pi3 : float, default = 1
    %     lag decay hyperparameter, defined in (6.18.25)  
    % 
    % omega_1 : float, default = 0.1
    %     overall tightness hyperparameter on gamma_i, defined in (6.18.27)     
    % 
    % credibility_level : float, default = 0.95
    %     VAR model credibility level (between 0 and 1)
    % 
    % burnin : int, default = 1000
    %     number of Gibbs sampler burn-in replications  
    % 
    % iterations : int, default = 2000
    %     number of Gibbs sampler replications   
    % 
    % verbose : bool, default = False
    %     if True, displays a progress bar 
    % 
    % 
    % Attributes
    % ----------
    % endogenous : matrix of size (n_obs,n)
    %     endogenous variables, defined in (6.18.1)
    % 
    % factors : int, default = 3
    %     number of latent factors, defined in (6.18.1)
    % 
    % loadings_lags : int, default = 2
    %     number of loadings lags, defined in (6.18.2)
    % 
    % factor_lags : int, default = 2
    %     number of factor lags, defined in (6.18.6)
    % 
    % residual_lags : int, default = 1
    %     number of residual lags, defined in (6.18.7)
    % 
    % sigma : float, default = 0.1
    %     variance on residual shock e_it, defined in (6.18.7)        
    % 
    % omega : float, default = 0.1
    %     variance on factor shock xi_t, defined in (6.18.6) 
    % 
    % delta1 : float, default = 0.1
    %     overall tightness hyperparameter on lambda_ij, defined in (6.18.21)
    % 
    % pi1 : float, default = 0.1
    %     overall tightness hyperparameter, defined in (6.18.25)
    % 
    % pi2 : float, default = 0.5
    %     cross-variable shrinkage hyperparameter, defined in (6.18.25)
    % 
    % pi3 : float, default = 1
    %     lag decay hyperparameter, defined in (6.18.25)  
    % 
    % omega_1 : float, default = 0.1
    %     overall tightness hyperparameter on gamma_i, defined in (6.18.27)     
    % 
    % credibility_level : float, default = 0.95
    %     VAR model credibility level (between 0 and 1)
    % 
    % burnin : int, default = 1000
    %     number of Gibbs sampler burn-in replications  
    % 
    % iterations : int, default = 2000
    %     number of Gibbs sampler replications   
    % 
    % verbose : bool, default = False
    %     if True, displays a progress bar         
    % 
    % y : matrix of size (n_obs,n)
    %     standardized endogenous variables
    % 
    % c : matrix of size (n_obs,1)
    %     mean of endogenous variables, prior to standardization   
    % 
    % S : matrix of size (n_obs,1)
    %     standard deviation of endogenous variables, prior to standardization  
    % 
    % T : int
    %     number of sample periods, defined in (6.18.1)  
    % 
    % n : int
    %     number of endogenous variables, defined in (6.18.1)
    % 
    % t_ : cell of dimension (n)
    %     cell of observed sample periods t_i, defined in (6.18.9) 
    % 
    % T_ : cell of dimension (n)
    %     cell of observed sample length T_i, defined in (6.18.9)   
    % 
    % x : cell of dimension (n)
    %     cell of sample observations x_i, defined in (6.18.10)         
    % 
    % x_ : cell of dimension (T)
    %     cell of sample observations, periods by periods
    % 
    % J : cell of dimension (T)
    %     selection matrices J_t, defined in (6.18.3)
    % 
    % d : cell of dimension (T)
    %     number of observations d_t for each sample period, defined in (6.18.3)
    % 
    % m : int
    %     number of latent factors, defined in (6.18.1)
    % 
    % q : int
    %     number of loadings lags, defined in (6.18.2)
    % 
    % p : int
    %     number of factor lags, defined in (6.18.6)
    % 
    % r : int
    %     number of residual lags, defined in (6.18.7)
    % 
    % s : int
    %     max(q,p), defined in (6.18.39)    
    % 
    % k : int
    %     factor VAR coefficients per equation, defined in (6.18.15)
    % 
    % l : int
    %     regression coefficients per loading equation, defined in (6.18.11)
    % 
    % mcmc_lambda : matrix of size (n,l,iterations)
    %     MCMC values of loadings coefficients lambda      
    % 
    % mcmc_beta : matrix of size (k,m,iterations)
    %     MCMC values of VAR coefficients beta         
    % 
    % mcmc_gamma : matrix of size (n,r,iterations)
    %     MCMC values of AR coefficients gamma           
    % 
    % mcmc_f : matrix of size (T,m*(p+1),iterations)
    %     MCMC values of latent factors f            
    % 
    % mcmc_eps : matrix of size (T,n*(r+1),iterations)
    %     MCMC values of latent residuals epsilon
    % 
    % lambda_estimates : matrix of size (n,l,4)
    %     posterior estimates for lambda
    %     page 1: median, page 2: lower bound, page 3: upper bound, page 4: standard deviation
    % 
    % beta_estimates : matrix of size (k,m,4)
    %     posterior estimates for beta
    %     page 1: median, page 2: lower bound, page 3: upper bound, page 4: standard deviation
    % 
    % gamma_estimates : matrix of size (n,r,4)
    %     posterior estimates for gamma
    %     page 1: median, page 2: lower bound, page 3: upper bound, page 4: standard deviation
    % 
    % f_estimates : matrix of size (T,m,4)
    %     posterior estimates for factors
    %     page 1: median, page 2: lower bound, page 3: upper bound, page 4: standard deviation
    %
    % mcmc_forecasts : matrix of size (f_periods,n,iterations)
    %     MCMC values of forecasts
    % 
    % mcmc_f_forecasts : matrix of size (f_periods,m,iterations)
    %     MCMC values of factor forecasts
    % 
    % forecast_estimates : matrix of size (f_periods,n,3)
    %     forecast estimates
    %     page 1: median, page 2: lower bound, page 3: upper bound
    % 
    % f_forecast_estimates : matrix of size (f_periods,m,3)
    %     forecast estimates for latent factors
    %     page 1: median, page 2: lower bound, page 3: upper bound
    % 
    % mcmc_irf : matrix of size (n,m+1,irf_periods,iterations)
    %     MCMC values of impulse response function
    % 
    % irf_estimates : matrix of size (n,m+1,irf_periods,3)
    %     posterior estimates of impulse response function
    %     page 1: median, page 2: lower bound, page 3: upper bound    
    % 
    % mcmc_fevd : matrix of size (n,m+1,fevd_periods,iterations)
    %     MCMC values of forecast error variance decompositions
    % 
    % fevd_estimates : matrix of size (n,m+1,fevd_periods,3)
    %     posterior estimates of forecast error variance decomposition
    %     page 1: median, page 2: lower bound, page 3: upper bound 
    % 
    % mcmc_hd : matrix of size (n,m+1,T,iterations)
    %     MCMC values of historical decompositions
    % 
    % hd_estimates : matrix of size (n,m+1,T,3)
    %     posterior estimates of historical decomposition
    %     page 1: median, page 2: lower bound, page 3: upper bound 
    % 
    % fitted_estimates : matrix of size (T,n,3)
    %     estimates of in-sample fit
    % 
    % residual_estimates : matrix of size (T,n,3)
    %     estimates of in-sample residuals
    % 
    % factor_residual_estimates : matrix of size (T,m,3)
    %     estimates of in-sample residuals     
    % 
    % insample_evaluation : struct
    %     in-sample evaluation criteria  
    %
    % forecast_evaluation_criteria : struct
    %     forecast evaluation criteria
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


    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
     
    
    properties (GetAccess = public, SetAccess= protected)
        endogenous
        factors
        loadings_lags
        factor_lags
        residual_lags
        sigma
        omega
        delta1
        delta2 
        pi1
        pi2
        pi3
        omega1
        credibility_level
        burnin
        iterations
        verbose
        y
        c
        S
        T
        n
        t_
        T_
        x
        x_
        J
        d
        m
        q
        p
        r
        s
        k
        l
        mcmc_lambda
        mcmc_beta
        mcmc_gamma
        mcmc_f
        mcmc_eps
        lambda_estimates
        beta_estimates
        gamma_estimates
        f_estimates
        mcmc_forecast
        mcmc_f_forecast
        forecast_estimates
        f_forecast_estimates
        mcmc_irf
        irf_estimates
        mcmc_fevd
        fevd_estimates
        mcmc_hd
        hd_estimates  
        fitted_estimates
        residual_estimates
        factor_residual_estimates
        insample_evaluation
        forecast_evaluation_criteria
    end
    
    
    properties (GetAccess = private, SetAccess = private)
        inv_U
        inv_U_h
        inv_V
        inv_Q
        mcmc_W
        mcmc_E_
        mcmc_xi
        mcmc_e
        mcmc_Y
    end
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods (Access = public)
        
        
        function self = BayesianDynamicFactorModel(endogenous, varargin)
            
            % constructor for the BayesianDynamicFactorModel class
            
            % allow for optional arguments
            parser = inputParser;
            default_factors = 3;
            default_loadings_lags = 2;
            default_factor_lags = 2;
            default_residual_lags = 1;
            default_sigma = 0.1;
            default_omega = 0.1;
            default_delta1 = 0.1;
            default_pi1 = 0.1;
            default_pi2 = 0.5;
            default_pi3 = 1;
            default_omega1 = 0.1;
            default_credibility_level = 0.95;
            default_burnin = 1000;
            default_iterations = 2000;
            default_verbose = false;
            addRequired(parser, 'endogenous');
            addParameter(parser, 'factors', default_factors);
            addParameter(parser, 'loadings_lags', default_loadings_lags);
            addParameter(parser, 'factor_lags', default_factor_lags);
            addParameter(parser, 'residual_lags', default_residual_lags);
            addParameter(parser, 'sigma', default_sigma); 
            addParameter(parser, 'omega', default_omega); 
            addParameter(parser, 'delta1', default_delta1);
            addParameter(parser, 'pi1', default_pi1);   
            addParameter(parser, 'pi2', default_pi2);  
            addParameter(parser, 'pi3', default_pi3);  
            addParameter(parser, 'omega1', default_omega1);  
            addParameter(parser, 'credibility_level', default_credibility_level);
            addParameter(parser, 'burnin', default_burnin);
            addParameter(parser, 'iterations', default_iterations);
            addParameter(parser, 'verbose', default_verbose);
            parse(parser, endogenous, varargin{:});
            self.endogenous = endogenous;
            self.factors = parser.Results.factors;
            self.loadings_lags = parser.Results.loadings_lags;
            self.factor_lags = parser.Results.factor_lags;
            self.residual_lags = parser.Results.residual_lags;
            self.sigma = parser.Results.sigma;
            self.omega = parser.Results.omega;
            self.delta1 = parser.Results.delta1;
            self.delta2 = 1e-10;
            self.pi1 = parser.Results.pi1;
            self.pi2 = parser.Results.pi2;
            self.pi3 = parser.Results.pi3;
            self.omega1 = parser.Results.omega1;
            self.credibility_level = parser.Results.credibility_level;
            self.burnin = parser.Results.burnin;
            self.iterations = parser.Results.iterations;
            self.verbose = parser.Results.verbose;
            % make regressors
            self.make_regressors();
        end
        
        
        function estimate(self)
            
            % estimate()
            % generates posterior estimates for the Bayesian dynamic factor model
            %
            % parameters:
            % none
            %
            % returns:
            % none


            % define prior values
            self.prior();
            % run MCMC algorithm (Gibbs sampling) for VAR parameters
            self.parameter_mcmc();
            % compute bridge equations
            self.bridge_equations();            
            % obtain posterior estimates for dfm parameters
            self.parameter_estimates();            
        end  
        

        function insample_fit(self)
            
            % insample_fit(self)
            % generates in-sample fit and residuals
            %
            % parameters:
            % none
            %
            % returns:
            % none
            
            % compute fitted and residuals
            self.fitted_and_residual();
            % compute in-sample criteria
            self.insample_criteria();          
        end   


        function [forecast_estimates] = forecast(self, h, credibility_level)
            
            % forecast(h, credibility_level)
            % estimates forecasts for the Bayesian dynamic factor model, using algorithm 18.3
            % 
            % parameters:
            % h : int
            %     number of forecast periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % 
            % returns:
            % forecast_estimates : matrix of size (h,n,3)
            %     page 1: median; page 2: interval lower bound; page 3: interval upper bound
            
            % get forecast
            self.make_forecast(h);
            % obtain posterior estimates
            self.forecast_posterior_estimates(credibility_level);
            forecast_estimates = self.forecast_estimates;
        end


        function [forecast_evaluation_criteria] = forecast_evaluation(self, Y)
            
            % [forecast_evaluation_criteria] = forecast_evaluation(self, Y)
            % forecast evaluation criteria for the Bayesian DFM
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
            % obtain regular forecast evaluation criteria 
            standard_evaluation_criteria = vu.forecast_evaluation_criteria(Y_hat, Y);
            % obtain Bayesian forecast evaluation criteria 
            bayesian_evaluation_criteria = vu.bayesian_forecast_evaluation_criteria(mcmc_forecast, Y);
            % merge structures
            forecast_evaluation_criteria = iu.concatenate_structures(standard_evaluation_criteria, bayesian_evaluation_criteria);
            % save as attributes
            self.forecast_evaluation_criteria = forecast_evaluation_criteria;     
        end 

    
        function [irf_estimates] = impulse_response_function(self, h, credibility_level)

            % impulse_response_function(h, credibility_level)
            % impulse response functions, as defined in (6.18.46)-(6.18.47)
            % 
            % parameters:
            % h : int
            %     number of IRF periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % 
            % returns:
            % irf_estimates : ndarray of shape (n,m+1,h,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound       
            
            % get regular impulse response funtion
            self.make_impulse_response_function(h);
            % obtain posterior estimates
            self.irf_posterior_estimates(credibility_level);  
            irf_estimates = self.irf_estimates;
        end


        function [fevd_estimates] = forecast_error_variance_decomposition(self, h, credibility_level)
            
            % forecast_error_variance_decomposition(self, h, credibility_level)
            % forecast error variance decomposition, as defined in (6.18.48)-(6.18.50)
            % 
            % parameters:
            % h : int
            %     number of FEVD periods
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % 
            % returns:
            % fevd_estimates : ndarray of shape (n,m+1,h,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound
            
            % get forecast error variance decomposition
            self.make_forecast_error_variance_decomposition(h);
            % obtain posterior estimates
            self.fevd_posterior_estimates(credibility_level);
            fevd_estimates = self.fevd_estimates;
        end

    
        function [hd_estimates] = historical_decomposition(self, credibility_level)

            % historical_decomposition(self, credibility_level)
            % historical decomposition, as defined in (6.18.51)-(6.18.53)
            % 
            % parameters:
            % credibility_level: float between 0 and 1
            %     credibility level for forecast credibility bands
            % 
            % returns:
            % hd_estimates : ndarray of shape (n,m+1,T,3)
            %     first 3 dimensions are variable, shock, period; 4th dimension is median, lower bound, upper bound
            
            % get historical decomposition
            self.make_historical_decomposition();
            % obtain posterior estimates
            self.hd_posterior_estimates(credibility_level);
            hd_estimates = self.hd_estimates;
        end


    end
    
    
    methods (Access = protected, Hidden = true)


        function make_regressors(self)
            
            % generates regressors, hyperparameters and dimensions
            
            % verify data integrity
            nwu.check_data_integrity(self.endogenous);
            % standardize data
            [y c S] = nwu.standardize_data(self.endogenous);
            % create regressors
            [T n t_ T_ x x_ J d] = nwu.make_dfm_regressors(y);
            % get dimensions
            [m q p r s k l] = nwu.make_dfm_dimensions(self.factors, ...
                              self.loadings_lags, self.factor_lags, self.residual_lags);
            % save as attributes
            self.y = y;
            self.c = c;
            self.S = S;
            self.T = T;
            self.n = n;
            self.t_ = t_;
            self.T_ = T_;
            self.x = x;
            self.x_ = x_;
            self.J = J;
            self.d = d;
            self.m = m;
            self.q = q;
            self.p = p;
            self.r = r;
            self.s = s;
            self.k = k;
            self.l = l;
        end


        function prior(self)
    
            % create prior elements for the dynamic factor model 
    
            % define prior hyperparameters for lambda
            [self.inv_U self.inv_U_h] = nwu.make_lambda_prior(self.n, self.m, ...
                                        self.q, self.l, self.delta1, self.delta2);
            % define prior hyperparameters for beta
            [self.inv_V] = nwu.make_beta_prior(self.m, self.p, self.pi1, self.pi2, self.pi3); 
            % define prior hyperparameters for gamma
            self.inv_Q = nwu.make_gamma_prior(self.r, self.omega1);
        end

        
        function parameter_mcmc(self)

            % Gibbs sampler for DFM parameters

            % unpack
            n = self.n;
            m = self.m;
            p = self.p;
            q = self.q;
            r = self.r;
            s = self.s;
            k = self.k;
            l = self.l;
            T = self.T;
            sigma = self.sigma;
            omega = self.omega;
            inv_U = self.inv_U;
            inv_U_h = self.inv_U_h;
            inv_V = self.inv_V;
            inv_Q = self.inv_Q;
            t_ = self.t_;
            x = self.x;
            x_ = self.x_;
            J = self.J;
            c = self.c;
            S = self.S;
            burnin = self.burnin;
            iterations = self.iterations;
            verbose = self.verbose;
            % preallocate storage space
            mcmc_lambda = zeros(n,l,iterations);
            mcmc_beta = zeros(k,m,iterations);
            mcmc_gamma = zeros(n,r,iterations);
            mcmc_W = zeros(T,l,iterations);
            mcmc_f = zeros(T,m*(s+1),iterations);
            mcmc_eps = zeros(T,n*(r+1),iterations);
            mcmc_E_ = cell(iterations,1);
            mcmc_xi = zeros(T,m,iterations);
            mcmc_e = zeros(T,n,iterations);
            mcmc_Y = zeros(T,n,iterations);

            % step 1: set initial values
            z = randn(T,m*(s+1)+n*(r+1));
            gamma = zeros(n,r);
            beta = zeros(k,m);
            [W F Z eps_ E_ E xi e] = nwu.make_state_regressors(beta, gamma, z, m, q, p, s, n, r, T);

            % iterate over iterations
            iteration = 1;
            while iteration <= (burnin + iterations)  
                
                % step 2: draw lambda_i
                lamda = self.draw_lambda(x, sigma, inv_U, inv_U_h, W, E_, gamma, t_, n, l);

                % step3: draw beta_j
                beta = self.draw_beta(F, omega, inv_V, Z, m, k);

                % step 4: draw gamma_i
                gamma = self.draw_gamma(eps_, E_, sigma, inv_Q, t_, n, r);
    
                % step 5: draw z
                z = self.draw_z(x_, J, lamda, beta, gamma, sigma, omega, n, m, q, p, s, r, T);
                [W F Z eps_ E_ E xi e] = nwu.make_state_regressors(beta, gamma, z, m, q, p, s, n, r, T);

                % recover nowcasts
                Y = self.nowcasts(W, lamda, E, c, S);

                % save if burn is exceeded
                if iteration > burnin
                    
                    % save parameter values
                    mcmc_lambda(:,:,iteration-burnin) = lamda;
                    mcmc_beta(:,:,iteration-burnin) = beta;
                    mcmc_gamma(:,:,iteration-burnin) = gamma;
                    mcmc_W(:,:,iteration-burnin) = W;
                    mcmc_f(:,:,iteration-burnin) = z(:,1:m*(s+1));
                    mcmc_eps(:,:,iteration-burnin) = z(:,end-n*(r+1)+1:end);
                    mcmc_E_{iteration-burnin} = E_;
                    mcmc_xi(:,:,iteration-burnin) = xi;
                    mcmc_e(:,:,iteration-burnin) = e;
                    mcmc_Y(:,:,iteration-burnin) = Y;
                end

                % if verbose, display progress bar
                if verbose
                    cu.progress_bar(iteration, iterations+burnin, 'Model parameters:'); 
                end

                % update iterations
                iteration = iteration + 1;

            end

            % save as attributes
            self.mcmc_lambda = mcmc_lambda;
            self.mcmc_beta = mcmc_beta;
            self.mcmc_gamma = mcmc_gamma;
            self.mcmc_f = mcmc_f;
            self.mcmc_eps = mcmc_eps;
            self.mcmc_W = mcmc_W;
            self.mcmc_E_ = mcmc_E_;
            self.mcmc_xi = mcmc_xi;
            self.mcmc_e = mcmc_e;
            self.mcmc_Y = mcmc_Y;
        end


        function [lamda] = draw_lambda(self, x, sigma, inv_U, inv_U_h, W, E_, gamma, t_, n, l)
            
            % draw lambda from its conditional posterior defined in (6.18.30)

            lamda = zeros(n,l);
            lamda(1:l,1:l) = eye(l);
            for i=l+1:n
                % get regressors
                x_i = x{i};
                E_i = E_{i}(t_{i},:);
                gamma_i = gamma(i,:)';
                W_i = W(t_{i},:);
                inv_U_i = inv_U{i};
                inv_U_h_i = inv_U_h{i};
                % posterior U_i_bar
                inv_U_i_bar = inv_U_i + W_i' * W_i / sigma;
                % posterior h_i_bar
                h_i_bar_temp = inv_U_h_i + W_i' * (x_i - E_i * gamma_i) / sigma;
                % efficient sampling of lambda_i (algorithm 9.4)
                lambda_i = rn.efficient_multivariate_normal(h_i_bar_temp, inv_U_i_bar);
                lamda(i,:) = lambda_i;
            end
        end
            
        function [beta] = draw_beta(self, F, omega, inv_V, Z, m, k)
    
            % draw beta from its conditional posterior defined in (6.18.33)

            beta = zeros(k,m);
            for j=1:m
                % get regressors
                f_j = F(:,j);
                inv_V_j = inv_V{j};
                % posterior V_j_bar
                inv_V_j_bar = inv_V_j + Z' * Z / omega;
                % posterior b_j_bar
                b_j_bar_temp = Z' * f_j / omega;
                % efficient sampling of beta (algorithm 9.4)
                beta_j = rn.efficient_multivariate_normal(b_j_bar_temp, inv_V_j_bar);
                beta(:,j) = beta_j;
            end
        end

    
        function [gamma] = draw_gamma(self, eps_, E_, sigma, inv_Q, t_, n, r)
            
            % draw gamma from its conditional posterior defined in (6.18.36)

            gamma = zeros(n,r);
            if r > 0
                for i=1:n
                    % get regressors
                    eps_i = eps_{i}(t_{i});
                    E_i = E_{i}(t_{i},:);
                    % posterior Q_i_bar
                    inv_Q_i_bar = inv_Q + E_i' * E_i / sigma; 
                    % posterior g_i_bar
                    g_i_bar_temp = E_i' * eps_i / sigma;
                    % efficient sampling of beta (algorithm 9.4)
                    gamma_i = rn.efficient_multivariate_normal(g_i_bar_temp, inv_Q_i_bar);
                    gamma(i,:) = gamma_i;
                end
            end
        end

    
        function [Z] = draw_z(self, x_, J, lamda, beta, gamma, sigma, omega, n, m, q, p, s, r, T)
            
            % draw z from the state-space representation (6.18.39)
            
            % get parameters for state-space representation
            [Lambda B Upsilon] = nwu.dfm_state_space_representation(lamda, beta, gamma, ...
                                 sigma, omega, n, m, q, p, s, r);
            % get initial values for algorithm
            [z_00 Upsilon_00] = nwu.dfm_kalman_initial_values(sigma, omega, s, r, m, n); 
            % run forward pass
            [Z_tt Z_tt1 Ups_tt Ups_tt1] = ss.dfm_forward_pass(x_, Lambda, B, Upsilon, ...
                                          z_00, Upsilon_00, T, m, n, s, r, J);
            % run backward pass
            [Z] = ss.dfm_backward_pass(Z_tt, Z_tt1, Ups_tt, Ups_tt1, B, T, m, s, n, r);
        end

     
        function [Y] = nowcasts(self, W, lamda, E, c, S)
            
            % draw Y from the dynamic factor model (6.18.2)
            
            % recover all values jointly from compact representation
            Y_hat = W * lamda' + E;
            % rescale the fitted values
            Y = Y_hat .* S + c;
        end


        function bridge_equations(self)
            
            % bridge equations for DFM
            
            % unpack
            n = self.n;
            r = self.r;       
            l = self.l;
            T = self.T;
            sigma = self.sigma;  
            inv_U = self.inv_U;
            inv_U_h = self.inv_U_h;
            inv_Q = self.inv_Q;
            t_ = self.t_;
            x = self.x;
            x_ = self.x_;
            J = self.J;     
            iterations = self.iterations;
            verbose = self.verbose;
            mcmc_W = self.mcmc_W;
            gamma = self.mcmc_gamma(1:l,:,1);
            E_ = self.mcmc_E_{1};
            
            % iterate over iterations
            for iteration=1:iterations  
                
                % get iteration parameters
                W = mcmc_W(:,:,iteration);

                % step 1: update lambda
                [lamda] = self.update_lambda(iteration, x, sigma, inv_U, inv_U_h, W, E_, gamma, t_, n, l);

                % step 2: update epsilon
                [z eps_ E_ E e] = self.update_epsilon(x_, lamda, W, gamma, J, sigma, l, r, T);

                % step 3: update gamma
                [gamma] = self.update_gamma(l, r, eps_, t_, E_, inv_Q, sigma);

                % update MCMC values
                self.update_mcmc_values(lamda, gamma, z, E_, e, n, l, r, iteration);

                % if verbose, display progress bar
                if verbose
                    cu.progress_bar(iteration, iterations, 'Bridge equations:')
                end

                % update iterations
                iteration = iteration + 1;     
            end
        end


        function [lamda] = update_lambda(self, iteration, x, sigma, inv_U, inv_U_h, W, E_, gamma, t_, n, l)
            
            % loadings lambda for the restricted variables
    
            lamda = zeros(l,l);
            for i=1:l            
                % get regressors
                x_i = x{i};
                E_i = E_{i}(t_{i},:);
                gamma_i = gamma(i,:)';
                W_i = W(t_{i},:);
                inv_U_i = inv_U{i};
                inv_U_h_i = inv_U_h{i};
                % posterior U_i_bar
                inv_U_i_bar = inv_U_i + W_i' * W_i / sigma;            
                % posterior h_i_bar    
                h_i_bar_temp = inv_U_h_i + W_i' * (x_i - E_i * gamma_i) / sigma;
                % efficient sampling of lambda_i (algorithm 9.4)
                lambda_i = rn.efficient_multivariate_normal(h_i_bar_temp, inv_U_i_bar);
                % update lambda
                lamda(i,:) = lambda_i;   
            end
        end

    
        function [z eps_ E_ E e] = update_epsilon(self, x_, lamda, W, gamma, J, sigma, l, r, T)
            
            % residuals epsilon for the restricted variables
            
            % get parameters for state-space representation
            [A B Upsilon] = nwu.epsilon_state_space_representation(gamma, sigma, l, r);
            % get initial values for algorithm
            [z_00 Upsilon_00] = nwu.epsilon_kalman_initial_values(sigma, r, l);
            % run forward pass
            [Z_tt Z_tt1 Ups_tt Ups_tt1] = ss.epsilon_forward_pass(x_, lamda, W, A, B, Upsilon, z_00, Upsilon_00, T, l, r, J);
            % run backward pass
            [z] = ss.epsilon_backward_pass(Z_tt, Z_tt1, Ups_tt, Ups_tt1, B, T, l, r);
            % update regressors
            [eps_ E_ E e] = nwu.update_epsilon_regressors(gamma, z, l, r, T);
        end


        function [gamma] = update_gamma(self, l, r, eps_, t_, E_, inv_Q, sigma)
            
            % AR coefficients gamma for the restricted variables

            gamma = zeros(l,r);
            if r > 0
                for i=1:l
                    % get regressors
                    eps_i = eps_{i}(t_{i});
                    E_i = E_{i}(t_{i},:);
                    % posterior Q_i_bar
                    inv_Q_i_bar = inv_Q + E_i' * E_i / sigma;
                    % posterior g_i_bar
                    g_i_bar_temp = E_i' * eps_i / sigma;
                    % efficient sampling of beta (algorithm 9.4)
                    gamma_i = rn.efficient_multivariate_normal(g_i_bar_temp, inv_Q_i_bar);
                    gamma(i,:) = gamma_i;
                end
            end
        end


        function update_mcmc_values(self, lamda, gamma, z, E_, e, n, l, r, iteration)
    
            % update mcmc values of restricted variables

            self.mcmc_lambda(1:l,:,iteration) = lamda;
            self.mcmc_gamma(1:l,:,iteration) = gamma; 
            for i=1:r+1
                self.mcmc_eps(:,(i-1)*n+1:(i-1)*n+l,iteration) = z(:,(i-1)*l+1:i*l);
            end
            for i=1:l
                self.mcmc_E_{iteration}{i} = E_{i};
            end
            self.mcmc_e(:,1:l,iteration) = e;
        end
            

        function parameter_estimates(self)
            
            % point estimates and credibility intervals for model parameters
    
            % unpack
            mcmc_lambda = self.mcmc_lambda;
            mcmc_beta = self.mcmc_beta;
            mcmc_gamma = self.mcmc_gamma;
            mcmc_f = self.mcmc_f;
            mcmc_eps = self.mcmc_eps;
            mcmc_Y = self.mcmc_Y;
            credibility_level = self.credibility_level;
            m = self.m;
            % recover posterior estimates
            lambda_estimates = nwu.posterior_estimates(mcmc_lambda, credibility_level);
            beta_estimates = nwu.posterior_estimates(mcmc_beta, credibility_level);
            gamma_estimates = nwu.posterior_estimates(mcmc_gamma, credibility_level);
            f_estimates = nwu.posterior_estimates(mcmc_f(:,1:m,:), credibility_level);

            % save as attributes
            self.lambda_estimates = lambda_estimates;
            self.beta_estimates = beta_estimates;
            self.gamma_estimates = gamma_estimates;
            self.f_estimates = f_estimates;
        end


        function fitted_and_residual(self)
        
            % point estimates and credibility intervals for fitted and residuals

            mcmc_eps = self.mcmc_eps;
            mcmc_xi = self.mcmc_xi;
            mcmc_Y = self.mcmc_Y;   
            n = self.n;
            m = self.m;
            credibility_level = self.credibility_level;
            eps_estimates = nwu.posterior_estimates(mcmc_eps(:,1:n,:), credibility_level);
            xi_estimates = nwu.posterior_estimates(mcmc_xi(:,1:m,:), credibility_level);
            Y_estimates = nwu.posterior_estimates(mcmc_Y, credibility_level);
            self.fitted_estimates = Y_estimates;
            self.residual_estimates = eps_estimates;
            self.factor_residual_estimates = xi_estimates;
        end


        function insample_criteria(self)
            
            % in-sample fit evaluation criteria
            
            insample_evaluation = vu.insample_evaluation_criteria(self.fitted_estimates(:,:,1), ...
                                  quantile(self.mcmc_e, 0.5, 3), self.T, self.l);      
            if self.verbose
                cu.progress_bar_complete('In-sample evaluation criteria:');
            end
            self.insample_evaluation = insample_evaluation;
        end


        function make_forecast(self, h)  
            
            % forecasts for the dynamic factor model
            
            % initiate storage and loop over iterations
            mcmc_forecast = zeros(h,self.n,self.iterations);
            mcmc_f_forecast = zeros(h,self.m,self.iterations);
            for i=1:self.iterations
                % make MCMC simulation for beta and Sigma
                [mcmc_forecast(:,:,i) mcmc_f_forecast(:,:,i)] = nwu.dfm_forecast(self.mcmc_lambda(:,:,i), ...
                    self.mcmc_beta(:,:,i), self.mcmc_gamma(:,:,i), self.sigma, self.omega, ...
                    self.mcmc_f(end,:,i), self.mcmc_eps(end,:,i), h, self.m, self.n, self.p, self.r, self.l);
                if self.verbose
                    cu.progress_bar(i, self.iterations, 'Forecasts:');
                end
            end
            self.mcmc_forecast = mcmc_forecast;
            self.mcmc_f_forecast = mcmc_f_forecast;
        end


        function forecast_posterior_estimates(self, credibility_level)
            
            % posterior estimates for forecasts
            
            % obtain posterior estimates
            mcmc_forecast = self.mcmc_forecast;
            mcmc_f_forecast = self.mcmc_f_forecast;
            forecast_estimates = nwu.posterior_estimates(mcmc_forecast, credibility_level);
            f_forecast_estimates = nwu.posterior_estimates(mcmc_f_forecast, credibility_level);
            self.forecast_estimates = forecast_estimates(:,:,1:3);
            self.f_forecast_estimates = f_forecast_estimates(:,:,1:3);
        end

    
        function make_impulse_response_function(self, h)

            % impulse response function for the dynamic factor model
            
            mcmc_irf = zeros(self.n, self.m+1, h, self.iterations);
            for i=1:self.iterations
                % get impulse response function
                mcmc_irf(:,:,:,i) = nwu.dfm_impulse_response_function(self.mcmc_lambda(:,:,i), ...
                                    self.mcmc_beta(:,:,i), self.mcmc_gamma(:,:,i), self.n, ...
                                    self.m, self.q, self.p, self.r, h);
                if self.verbose    
                    cu.progress_bar(i, self.iterations, 'Impulse response function:');
                end
            end
            self.mcmc_irf = mcmc_irf;
        end


        function irf_posterior_estimates(self, credibility_level)
            
            % posterior estimates for impulse response function
            
            mcmc_irf = self.mcmc_irf;
            irf_estimates = nwu.posterior_estimates_3d(mcmc_irf, credibility_level);
            self.irf_estimates = irf_estimates;
        end


        function make_forecast_error_variance_decomposition(self, h)
            
            % forecast error variance decomposition for the dynamic factor model    
    
            mcmc_fevd = zeros(self.n, self.m+1, h, self.iterations);
            has_irf = ~isempty(self.mcmc_irf) && size(self.mcmc_irf, 3) >= h;
            for i=1:self.iterations
                % recover structural IRF or estimate them
                if has_irf
                    irf = self.mcmc_irf(:,:,1:h,i);
                else
                    irf = nwu.dfm_impulse_response_function(self.mcmc_lambda(:,:,i), self.mcmc_beta(:,:,i), ...
                          self.mcmc_gamma(:,:,i), self.n, self.m, self.q, self.p, self.r, h);
                end
                % recover fevd
                mcmc_fevd(:,:,:,i) = nwu.dfm_forecast_error_variance_decomposition(irf, ...
                                     self.sigma, self.omega, self.n, self.m, h);
                if self.verbose   
                    cu.progress_bar(i, self.iterations, 'Forecast error variance decomposition:');
                end
            end
            self.mcmc_fevd = mcmc_fevd;
        end


        function fevd_posterior_estimates(self, credibility_level)

            % posterior estimates for forecast error variance decomposition

            mcmc_fevd = self.mcmc_fevd;
            fevd_estimates = vu.posterior_estimates_3d(mcmc_fevd, credibility_level);
            normalized_fevd_estimates = nwu.normalize_fevd_estimates(fevd_estimates);
            self.fevd_estimates = normalized_fevd_estimates;
        end

    
        function make_historical_decomposition(self)
            
            % historical decomposition for the dynamic factor model
            
            % initiate MCMC storage for HD and loop over iterations
            mcmc_hd = zeros(self.n, self.m+1, self.T, self.iterations);
            has_irf = ~isempty(self.mcmc_irf) && size(self.mcmc_irf, 3) >= self.T;
            for i=1:self.iterations
                % recover structural IRF or estimate them
                if has_irf
                    irf = self.mcmc_irf(:,:,1:self.T,i);
                else
                    irf = nwu.dfm_impulse_response_function(self.mcmc_lambda(:,:,i), self.mcmc_beta(:,:,i), ...
                          self.mcmc_gamma(:,:,i), self.n, self.m, self.q, self.p, self.r, self.T);
                end
                xi = self.mcmc_xi(:,:,i);
                e = self.mcmc_e(:,:,i);
                % get historical decomposition
                mcmc_hd(:,:,:,i) = nwu.dfm_historical_decomposition(irf, xi, e, self.n, self.m, self.T); 
                if self.verbose
                    cu.progress_bar(i, self.iterations, 'Historical decomposition:');
                end
            end
            self.mcmc_hd = mcmc_hd;
        end


        function hd_posterior_estimates(self, credibility_level)

            % posterior estimates for historical decomposition

            mcmc_hd = self.mcmc_hd;
            hd_estimates = vu.posterior_estimates_3d(mcmc_hd, credibility_level);
            self.hd_estimates = hd_estimates;
        end


    end
    

end

    