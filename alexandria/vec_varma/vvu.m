classdef vvu
    

    % vvu stands for vec varma utilities
    % A class containing static methods for vec and varma utilities

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)


        function [Y Z X] = make_var_regressors(endogenous, exogenous, lags, ...
                         constant, trend, quadratic_trend)

            % make_var_regressors(endogenous, exogenous, lags, constant, trend, quadratic_trend)
            % generates VAR regressors, defined in (4.11.3)
            % 
            % parameters:
            % endogenous : matrix of size (periods,n_endogenous)
            %     matrix of endogenous regressors        
            % exogenous : matrix of size (periods,n_exogenous)
            %     matrix of exogenous regressors     
            % lags: int
            %     number of lags in VAR model        
            % constant : bool
            %     if true, a constant is added to the VAR model
            % trend : bool
            %     if true, a linear trend is added to the VAR model
            % quadratic_trend : bool
            %     if true, a quadratic trend is added to the VAR model
            %     
            % returns:
            % Y : matrix of size (T,n)
            %     matrix of endogenous variables      
            % Z : matrix of size (T,n)
            %     matrix of exogenous variables 
            % X : matrix of size (T,k)
            %     matrix of regressors

            % get Y matrix
            Y = endogenous(lags+1:end,:);
            % get X matrix
            periods = size(endogenous,1) - lags;
            X_1 = vu.generate_intercept_and_trends(constant, trend, quadratic_trend, periods, 0);
            X_2 = vu.generate_exogenous_regressors(exogenous, lags);
            X_3 = vu.generate_lagged_endogenous(endogenous, lags);
            Z = [X_1,X_2];
            X = [X_1,X_2,X_3];
        end


        function [DY Y_1 Z] = make_vec_regressors(endogenous, exogenous, ... 
                              lags, constant, trend, quadratic_trend)
            
            % make_vec_regressors(endogenous, exogenous, lags, constant, trend, quadratic_trend)
            % generates VEC regressors, defined in (5.15.8)
            %
            % parameters:
            % endogenous : matrix of size (periods,n_endogenous)
            %     matrix of endogenous regressors        
            % exogenous : matrix of size (periods,n_exogenous)
            %     matrix of exogenous regressors     
            % lags: int
            %     number of lags in VAR model        
            % constant : bool
            %     if true, a constant is added to the VAR model
            % trend : bool
            %     if true, a linear trend is added to the VAR model
            % quadratic_trend : bool
            %     if true, a quadratic trend is added to the VAR model
            %     
            % returns:
            % DY : matrix of size (T,n)
            %     matrix of differenced endogenous variables   
            % Y_1 : matrix of size (T,n)
            %     matrix of lagged endogenous variables, one period           
            % Z : matrix of size (T,k)
            %     matrix of endogenous and lagged regressors   
            
            % get DY matrix
            diff_endogenous = endogenous(2:end,:) - endogenous(1:end-1,:);
            DY = diff_endogenous(lags+1:end,:);
            % get Y_1 matrix
            Y_1 = endogenous(lags+1:end-1,:);
            % get Z matrix
            periods = size(endogenous, 1) - lags - 1;
            Z_1 = vu.generate_intercept_and_trends(constant, trend, quadratic_trend, periods, 0);
            Z_2 = vu.generate_exogenous_regressors(exogenous, lags+1);
            Z_3 = vu.generate_lagged_endogenous(diff_endogenous, lags);
            Z = [Z_1,Z_2,Z_3];
        end

        
        function [n m p_vec p T_vec T k_vec k q_vec q r] = generate_dimensions(Y, DY, exogenous,  ...
                 lags, max_cointegration_rank, constant, trend, quadratic_trend)
        
            % generate_dimensions(DY, exogenous, lags, max_cointegration_rank, constant, trend, quadratic_trend)
            % generate VEC dimensions, defined in (5.15.8)
            % 
            % parameters:
            % Y : matrix of size (T,n)
            %     matrix of endogenous variables             
            % DY : matrix of size (T,n)
            %     matrix of differenced endogenous variables              
            % exogenous : matrix of size (periods,n_exogenous)
            %     matrix of exogenous regressors     
            % lags: int
            %     number of lags in VAR model     
            % max_cointegration_rank: int
            %     maximum cointegration rank
            % constant : bool
            %     if true, a constant is added to the VAR model
            % trend : bool
            %     if true, a linear trend is added to the VAR model
            % quadratic_trend : bool
            %     if true, a quadratic trend is added to the VAR model
            %    
            % returns:
            % n : int
            %     number of endogenous variables         
            % m : int
            %     number of exogenous variables  
            % p_vec : int
            %     number of lags of VEC model
            % p : int
            %     number of lags of equivalent VAR model
            % T_vec : int
            %     number of sample periods of VEC model
            % T : int
            %     number of sample periods of equivalent VAR model
            % k_vec : int
            %     number of coefficients in each equation of VEC model
            % k : int
            %     number of coefficients in each equation of equivalent VAR model
            % q_vec : int
            %     total number of coefficients of VEC model
            % q : int
            %     total number of coefficients of equivalent VAR model
            % r : int
            %     maximum cointegration rank
        
            T_vec = size(DY, 1);
            T = size(Y,1);
            n = size(DY, 2);
            p_vec = lags;
            p = lags + 1;
            m = double(constant) + double(trend) + double(quadratic_trend);    
            if ~isempty(exogenous)
                m = m + size(exogenous,2);
            end
            k_vec = m + n * p_vec;
            k = m + n * p;
            q_vec = n * k_vec;
            q = n * k;
            r = max_cointegration_rank;
        end


        function [s] = individual_ar_variances(n, endogenous, lags)
            
            % individual_ar_variances(n, endogenous, lags)
            % generates residual variance for each variable
            % 
            % parameters:  
            % n : int
            %     number of endogenous variables  
            % endogenous : matrix of size (periods,n_endogenous)
            %     matrix of endogenous regressors               
            % lags: int
            %     number of lags in VAR model   
            % 
            % returns:
            % s : matrix of size (n,1)
            %     vector of individual residual variances
            
            diff_endogenous = endogenous(2:end,:) - endogenous(1:end-1,:);
            s = zeros(n,1);
            for i=1:n
                ar = MaximumLikelihoodVar(diff_endogenous(:,i), 'lags', lags);
                ar.estimate();
                s(i) = ar.Sigma;
            end
        end


        function [B] = vec_to_var(Xi_T, Phi, n, m, p, k)
            
            % vec_to_var(Xi_T, Phi, n, m, p, k)
            % converts VEC to VAR using equation (5.15.103)
            % 
            % parameters:  
            % Xi_T : matrix of size (n,n)
            %     transpose of error correction matrix Xi, defined in (5.15.3)
            % Phi : matrix of size (k,n)
            %     matrix of VEC coefficients, defined in (5.15.8)          
            % n : int
            %     number of endogenous variables         
            % m : int
            %     number of exogenous variables         
            % p : int
            %     number of lags
            % k : int
            %     number of coefficients in each VEC equation 
            % 
            % returns:
            % B : matrix of size (k+n,n)
            %     matrix of equivalent VAR coefficients, defined in (4.11.3) and (5.15.103)
            
            % initiate storage
            A_T = zeros(n,n,p+1);
            B = zeros(k+n,n);
            % last VAR lag
            A_T(:,:,end) = - Phi(end-n+1:end,:);
            B(end-n+1:end,:) = A_T(:,:,end);
            % VAR lags from 2 to p-1
            for i=p:-1:2 
                A_T(:,:,i) = - Phi(m+1+(i-2)*n:m+(i-1)*n,:) - sum(A_T(:,:,i+1:end),3);
                B(m+1+(i-1)*n:m+i*n,:) = A_T(:,:,i);
            end
            % first VAR lag
            A_T(:,:,1) = Xi_T + eye(n) - sum(A_T(:,:,2:end),3);
            B(m+1:m+n,:) = A_T(:,:,1);
            % exogenous regressors
            B(1:m,:) = Phi(1:m,:);
        end


        function [steady_state_estimates] = vec_steady_state(Y,n,T)
            
            % vec_steady_state(Y);
            % pseudo steady-state estimates for VEC model
            % 
            % parameters:  
            % Y : matrix of size (T,n)
            %     matrix of endogenous variables 
            % n : int
            %     number of endogenous variables 
            % T : int
            %     number of sample periods
            %    
            % returns:
            % steady_state_estimates : matrix of size (T,n,3)
            %     matrix of steady-state estimates
            
            steady_state_estimates = zeros(T,n,3);
            steady_state_mean = mean(Y,1);
            steady_state_std = std(Y,0,1);
            lower_bound = steady_state_mean - 0.1 * steady_state_std;
            upper_bound = steady_state_mean + 0.1 * steady_state_std;
            steady_state_estimates(:,:,1) = repmat(steady_state_mean,[T 1]);
            steady_state_estimates(:,:,2) = repmat(lower_bound,[T 1]);
            steady_state_estimates(:,:,3) = repmat(upper_bound,[T 1]);
        end


        function [mcmc_irf] = make_varma_restriction_irf(mcmc_beta, mcmc_kappa, mcmc_chol_Sigma, iterations, n, p, q, max_irf_period)

            % make_varma_restriction_irf(mcmc_beta, mcmc_kappa, mcmc_chol_Sigma, iterations, n, p, q, max_irf_period)
            % creates preliminary orthogonalized IRFs for restriction algorithm
            % 
            % parameters:
            % mcmc_beta : matrix of size (k1, n, iterations)
            %     matrix of mcmc values for beta
            % mcmc_kappa : matrix of size (k2, n, iterations)
            %     matrix of mcmc values for K           
            % mcmc_chol_Sigma : matrix of size (n, n, iterations)
            %     matrix of mcmc values for h(Sigma)
            % iterations: int
            %     number of MCMC iterations
            % n : int
            %     number of endogenous variables       
            % p : int
            %     number of lags
            % q : int
            %     number of residual lags            
            % max_irf_period : int
            %     maximum number of periods for which IRFs will have to be computed in later algorithms        
            % 
            % returns:
            % mcmc_irf : matrix of size (n, n, iterations)
            %     matrix of mcmc values for preliminary orthogonalized IRFs 

            if max_irf_period == 0
                mcmc_irf = [];
            else
                mcmc_irf = zeros(n, n, max_irf_period, iterations);
                for i=1:iterations
                    irf = vvu.varma_impulse_response_function(mcmc_beta(:,:,i), mcmc_kappa(:,:,i), n, p, q, max_irf_period);
                    structural_irf = vu.structural_impulse_response_function(irf, mcmc_chol_Sigma(:,:,i), n);            
                    mcmc_irf(:,:,:,i) = structural_irf;
                end
            end
        end


        function [mcmc_shocks] = make_varma_restriction_shocks(mcmc_E, mcmc_chol_Sigma, T, n, iterations, restriction_matrices)

            % make_varma_restriction_shocks(mcmc_E, mcmc_chol_Sigma, T, n, iterations, restriction_matrices)
            % creates preliminary structural shocks for restriction algorithm
            % 
            % parameters:
            % mcmc_E : matrix of size (T, n, iterations)
            %     matrix of mcmc values for residuals
            % mcmc_chol_Sigma : matrix of size (n, n, iterations)
            %     matrix of mcmc values for h(Sigma)   
            % T : int
            %     number of sample periods   
            % n : int
            %     number of endogenous variables        
            % iterations: int
            %     number of MCMC iterations        
            % restriction_matrices : cell array of dimension (7,2)
            %     each cell entry stores matrices of restriction and coefficient values     
            % 
            % returns:
            % mcmc_shocks : matrix of size (T, n, iterations)
            %     matrix of mcmc values for preliminary structural shocks

            if isempty(restriction_matrices{4,1}) && isempty(restriction_matrices{5,1})
                mcmc_shocks = [];
            else
                mcmc_shocks = zeros(T, n, iterations);
                for i=1:iterations
                    E = mcmc_E(:,:,i);
                    Xi = E / mcmc_chol_Sigma(:,:,i)';
                    mcmc_shocks(:,:,i) = Xi;
                end
            end
        end  


        function [steady_state] = varma_steady_state(X, B, n, m, p, T)
            
            % [steady_state] = steady_state(Z, B, n, m, p, T)
            % computes the steady-state of the VAR model, using equation (4.12.30)
            %
            % parameters:
            % X : matrix of size (T,k1)
            %     matrix of exogenous variables
            % B : matrix of size (k1,n)
            %     matrix of VAR coefficients
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables            
            % p : int
            %     number of lags            
            % T : int
            %     number of sample periods
            %
            % returns:
            % steady_state : matrix of size (T,n)
            %     matrix of steady-state values
                      
            if m == 0
                steady_state = zeros(T,n);            
            else
                Z = X(:,1:m);
                C = B(1:m,:);
                A = eye(n);
                for i=1:p
                    A = A - B(1+m+(i-1)*n:m+i*n,:);
                end
                steady_state = Z * C / A;
            end
        end  


        function [E, Y_hat] = varma_fit_and_residuals(X, B, Z, K, E)
            
            % varma_fit_and_residuals(X, B, Z, K, E)
            % generates fit and residuals for a VARMA model, using (5.16.2)
            %
            % parameters:
            % X : matrix of size (T,k1)
            %     matrix of regressors
            % B : matrix of size (k,n)
            %     matrix of VAR coefficients
            % Z : matrix of size (T,k2)
            %     matrix of lagged residuals
            % K : matrix of size (k2,n)
            %     matrix of residual coefficients        
            % E : matrix of size (T,n)
            %     matrix of residuals 
            %
            % returns:
            % Y_hat : matrix of size (T,n)
            %     matrix of fitted endogenous
            % E : matrix of size (T,n)
            %     matrix of VAR residuals
            
            Y_hat = X * B + Z * K;        
        end


        function [Y_p] = varma_forecast(B, K, chol_Sigma, h, Z_p, Y, E, n)

            % varma_forecast(B, K, chol_Sigma, h, Z_p, Y, E, n)
            % products simulated forecasts
            % 
            % parameters:
            % B : matrix of size (k1,n)
            %     matrix of VAR coefficients
            % K : matrix of size (k2,n)
            %     matrix of residual coefficients
            % chol_Sigma : ndarray of shape (n,n)
            %     Cholesky factor of residual variance-covariance matrix Sigma
            % h : int
            %     number of forecast periods
            % Z_p : matrix of shape (h,m)
            %     matrix of exogenous regressors for forecasts
            % Y : matrix of shape (p,n)
            %     matrix of endogenous variables
            % E : matrix of size (q,n)
            %     matrix of residuals             
            % n : int
            %     number of endogenous variables
            % 
            % returns:
            % Y_p : matrix of size (h,n)
            %     matrix of simulated forecast values
        
            % E = (chol_Sigma * randn(n,h))';
            Y_p = zeros(h,n);
            for j=1:h
                % get lagged endogenous regressors
                X = la.vec(fliplr(Y'))';
                % add exogenous regressors, if any
                if ~isempty(Z_p)
                    X = [Z_p(j,:) X];
                end
                % get lagged residual regressors
                Z = la.vec(fliplr(E'))';
                % recover residuals
                e = (chol_Sigma * randn(n,1))';
                % generate forecasts
                y = X * B + Z * K + e;
                % update Y, E and Y_p
                Y = [Y(2:end,:);y];
                E = [E(2:end,:);e];
                Y_p(j,:) = y;
            end   
        end


        function [irf] = varma_impulse_response_function(B, K, n, p, q, h)

            % [irf] = varma_impulse_response_function(B, K, n, p, q, h)
            % generates impulse response function for given VARMA coefficients
            % using equations (5.16.31)-(5.16.33)
            %
            % parameters:
            % B : matrix of size (k1,n)
            %     matrix of VAR coefficients
            % K : matrix of size (k2,n)
            %     matrix of residual coefficients
            % n : int
            %     number of endogenous variables
            % p : int
            %     number of lags
            % q : int
            %     number of residual lags
            % h : int
            %     number of irf periods (in addition to impact period)
            %
            % returns:
            % irf : matrix of size (n,n,h)
            %     matrix of impulse response functions

            B = B(end+1-n*p:end,:);
            Yh = eye(n);
            irf = cat(3, Yh, zeros(n,n,h-1));
            Xh = zeros(n,n*p);
            Zh = eye(n,n*q);
            for i=2:h
                Xh = [Yh Xh(:,1:end-n)];
                Yh = Xh * B + Zh * K;
                Zh = [zeros(n,n) Zh(:,1:end-n)];
                irf(:,:,i) = Yh';
            end
        end


        function [y_bar Q omega] = varma_conditional_forecast_regressors_1(conditions, h, n, p, q)

            % varma_conditional_forecast_regressors_1(conditions, h, n, p, q)
            % first set of elements for conditional forecasts: iteration-invariant 
            % 
            % parameters:
            % conditions : matrix of size (nconditions,4)
            %     matrix of conditions (one row per condition: variable, period, mean, variance)
            % h : int
            %     number of forecast periods         
            % n : int
            %     number of endogenous variables 
            % p : int
            %     number of lags   
            % q : int
            %     number of residual lags              
            % 
            % returns:
            % y_bar : matrix of size (h,n)
            %     matrix of mean values for conditions
            % Q : matrix of size (n,n*(p+q))
            %     selection matrix for conditional forecasts state-space representation        
            % omega : matrix of size (h,n)
            %     matrix of variance values for conditions

            y_bar = zeros(h,n);
            omega = 10e6 * ones(h,n);
            for i=1:size(conditions,1)
                variable = conditions(i,1);
                period = conditions(i,2);
                mean = conditions(i,3);
                variance = max(1e-10,conditions(i,4));
                y_bar(period,variable) = mean;
                omega(period,variable) = variance;
            end
            Q = eye(n,n*(p+q));
        end


        function [mu F K gamma_00 Upsilon_00] = varma_conditional_forecast_regressors_2(Y, E, B, Kappa, Sigma, conditions, Z_p, n, m, p, q, h)

            % varma_conditional_forecast_regressors_2(Y, E, B, Kappa, Sigma, conditions, Z_p, n, m, p, q, h)
            % second set of elements for conditional forecasts: iteration-specific 
            % 
            % parameters:    
            % Y : matrix of size (p,n)
            %     matrix of initial conditions for exogenous   
            % E : matrix of size (q,n)
            %     matrix of initial conditions for residuals      
            % B : matrix of size (k1,n)
            %     matrix of VAR coefficients
            % Kappa : matrix of size (k2,n)
            %     matrix of residual coefficients 
            % Sigma : matrix of size (n,n)
            %     variance-covariance matrix of VAR residuals   
            % conditions : matrix of sizee (nconditions,4)
            %     matrix of conditions (one row per condition: variable, period, mean, variance)        
            % Z_p : matrix of size (h,m)
            %     matrix of exogenous regressors for forecasts
            % n : int
            %     number of endogenous variables 
            % m : int
            %     number of exogenous variables          
            % p : int
            %     number of lags 
            % q : int
            %     number of residual lags         
            % h : int
            %     number of forecast periods 
            % 
            % returns:
            % mu : matrix of size (h,n*(p+q))
            %     matrix of intercepts for state variables
            % F : matrix of size (n*p,n*p)
            %     companion form matrix
            % K : matrix of size (n*p,n*p,h)
            %     variance-covariance matrix for state errors
            % gamma_00 : matrix of size (n*p,1)
            %     initial conditions (mean) for the space vector gamma_hat            
            % Upsilon_00 : matrix of size (n*p,)
            %     initial conditions (variance) for the space vector gamma_hat

            F = vvu.make_varma_companion_form(B, Kappa, n, m, p, q);
            mu = zeros(h,n*(p+q));
            mu(:,1:n) = Z_p * B(1:m,:);
            K = zeros(n*(p+q),n*(p+q),h);
            selection = zeros(p+q,p+q);
            selection(1,1) = 1;
            selection(1,p+1) = 1;
            selection(p+1,1) = 1;
            selection(p+1,p+1) = 1;
            for i=1:h
                temp = Sigma;
                period_conditions = conditions(conditions(:,2) == i,:);
                condition_variables = period_conditions(:,1);
                for j=1:size(condition_variables,1)
                    variable = condition_variables(j);
                    temp(variable,variable) = 100;
                end
                K(:,:,i) = kron(selection, temp);
            end
            gamma_00 = [la.vec(fliplr(Y'));la.vec(fliplr(E'))];
            Upsilon_00 = kron(selection, Sigma) + 1e-10 * eye(n*(p+q));
        end


        function [F] = make_varma_companion_form(B, K, n, m, p, q)

            % make_varma_companion_form(B, K, n, m, p, q)
            % creates companion form matix F as defined in (5.16.35)
            % 
            % parameters:
            % B : matrix of size (k1,n)
            %     matrix of VAR coefficients
            % K : matrix of size (k2,n)
            %     matrix of residual coefficients 
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables        
            % p : int
            %     number of lags
            % q : int
            %     number of residual lags  
            % 
            % returns:
            % F : matrix of size (n*p,n*p)
            %     companion form matrix

            block_1 = [B(m+1:end,:)' K'];
            block_2 = eye(n*(p+q-1),n*(p+q));
            F = [block_1; block_2];
            F(end+1-q*n:end-(q-1)*n,end+1-(q+1)*n:end-q*n) = zeros(n,n);
        end  


        function [Y_p] = varma_linear_forecast(B, K, h, Z_p, Y, E, n)

            % varma_linear_forecast(B, K, h, Z_p, Y, E, n)
            % best linear forecasts for VARMA model, absent shocks
            % 
            % parameters:
            % B : matrix of size (k1,n)
            %     matrix of VAR coefficients
            % K : matrix of size (k2,n)
            %     matrix of residual coefficients 
            % h : int
            %     number of forecast periods
            % Z_p : matrix of shape (h,m)
            %     matrix of exogenous regressors for forecasts
            % Y : matrix of shape (T,n)
            %     matrix of endogenous variables
            % E : matrix of size (q,n)
            %     matrix of initial conditions for residuals             
            % n : int
            %     number of endogenous variables
            % 
            % returns:
            % Y_p : matrix of size (h,n)
            %     matrix of simulated forecast values
        
            Y_p = zeros(h,n);
            for j=1:h
                % get lagged endogenous regressors
                X = la.vec(fliplr(Y'))';
                % add exogenous regressors, if any
                if ~isempty(Z_p)
                    X = [Z_p(j,:) X];
                end
                % get lagged residuals
                Z = la.vec(fliplr(E'))';              
                % generate forecasts
                y = X * B + Z * K;
                % update Y, E and Y_p
                Y = [Y(2:end,:);y];
                E = [E(2:end,:);zeros(1,n)];
                Y_p(j,:) = y;
            end   
        end  


    end
    
    
end
        