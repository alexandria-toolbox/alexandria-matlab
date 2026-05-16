classdef nwu
    

    % nwu stands for nowcasting utilities
    % a class containing static methods for nowcasting utilities

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)

        
        function check_data_integrity(data)

            % check_data_integrity(data)
            % verifies that data is valid, i.e. does not contain all-NaN rows or columns
            % 
            % parameters:
            % data : matrix of size (n_periods,n_exogenous)
            %     matrix of endogenous data to verify
            %    
            % returns:
            % none
            
            % check for rows that are all-NaN
            all_nan_rows = find(all(isnan(data), 2));
            if ~isempty(all_nan_rows)
                row = all_nan_rows(1);
                error(['Data error: row ' num2str(row) ' (and possibly others) of endogenous data is all-NaN']);
            end
            % check for columns that are all-NaN
            all_nan_columns = find(all(isnan(data), 1));
            if ~isempty(all_nan_columns)
                column = all_nan_columns(1);
                error(['Data error: column ' num2str(column) ' (and possibly others) of endogenous data is all-NaN']);
            end
        end


        function [y c S] = standardize_data(endogenous)     
                
            % standardize_data(endogenous)
            % standardize data, making it zero mean and unit standard deviation
            % 
            % parameters:
            % endogenous: matrix of size (n_periods,n_endogenous)
            %     array of endogenous data to process 
            %     
            % returns:
            % y: matrix of size (n_periods,n_endogenous)
            %     array of standardized endogenous data  
            % c: matrix of size (n,)
            %     vector of empirical mean     
            % S: matrix of size (n,)
            %     vector of empirical standard deviation          
            
            c = mean(endogenous,1,'omitmissing');
            S = std(endogenous,1,'omitmissing');
            y = (endogenous - c) ./ S;
        end


        function [T n t_ T_ x x_ J d] = make_dfm_regressors(data)
            
            % make_dfm_regressors(data)
            % get dynamic factor model regressors and dimensions
            % 
            % parameters:
            % data : matrix of size (n_periods,n_exogenous)
            %     array of endogenous data to process 
            % 
            % returns:
            % T: int
            %     total number of sample periods, defined in (6.18.1)
            % n: int
            %     number of endogenous variables, defined in (6.18.1)
            % t_: cell of size (n)
            %     cell of non-NaN periods for each endogenous variable, defined in (6.18.9)
            % T_: cell of size (n)
            %     cell of non-NaN dimension for each endogenous variable, defined in (6.18.9)        
            % x: cell of size (n)
            %     cell of non-NaN data for each endogenous variable, defined in (6.18.9)  
            % x_: cell of size (T)
            %     cell of non-NaN data for each sample period, defined in (6.18.10)         
            % J: cell of size (T)
            %     cell of non-NaN variables for each sample period, defined in (6.18.3)
            % d: cell of size (T)
            %     cell of non-NaN dimension for each sample period, defined in (6.18.3)            

            T = size(data,1);
            n = size(data,2);
            t_ = cell(n,1);
            T_ = cell(n,1);
            x = cell(n,1);
            J = cell(T,1);
            d = cell(T,1);
            x_ = cell(T,1);
            for i=1:n
                data_i = data(:,i);
                non_nan_entries = find(~isnan(data_i));
                t_{i} = non_nan_entries;
                T_{i} = size(non_nan_entries,1);
                x{i} = data_i(non_nan_entries);
            end
            for t=1:T
                data_t = data(t,:);
                non_nan_entries = find(~isnan(data_t));
                J{t} = non_nan_entries;
                d{t} = size(non_nan_entries,2);
                x_{t} = data_t(non_nan_entries);
            end
        end


        function [m q p r s k l] = make_dfm_dimensions(factors, loadings_lags, factor_lags, residual_lags)

            % make_dfm_dimensions(factors, loadings_lags, factor_lags)
            % generate additional dimensions for dynamic factor model
            % 
            % parameters:
            % factors : int
            %     number of fundamental factors in the model
            % loadings_lags : int
            %     number of lags in the loadings equation
            % factor_lags : int
            %     number of lags in the dynamic equation        
            % 
            % returns:
            % m: int
            %     number of fundamental factors in the model, defined in (6.18.1)        
            % q: int
            %     number of lags in the loadings equation, defined in (6.18.2)        
            % p: int
            %     number of lags in the VAR model, defined in (6.18.6)  
            % r: int
            %     number of lags in the residual AR process, defined in (6.18.7) 
            % s: int
            %     maximum number of factor lags, defined in (6.18.38)
            % k: int
            %     total number of VAR regressors, defined in (6.18.15) 
            % l: int
            %     total number of loadings regressors, defined in (6.18.11)      
            
            m = factors;
            q = loadings_lags;
            p = factor_lags;
            r = residual_lags;
            s = max(q,p);
            k = m * p;
            l = m * (q + 1);
        end


        function [inv_U inv_U_h] = make_lambda_prior(n, m, q, l, delta1, delta2)
            
            % make_lambda_prior(n, m, q, l, delta1, delta2)
            % generate prior terms for loadings lambda
            % 
            % parameters:
            % 
            % n: int
            %     number of endogenous variables, defined in (6.18.1)        
            % m: int
            %     number of fundamental factors in the model, defined in (6.18.1)        
            % q: int
            %     number of lags in the loadings equation, defined in (6.18.2)           
            % l: int
            %     total number of loadings regressors, defined in (6.18.11)         
            % 
            % factors : int
            %     number of fundamental factors in the model
            % loadings_lags : int
            %     number of lags in the loadings equation
            % factor_lags : int
            %     number of lags in the dynamic equation        
            % delta1: float
            %     overall tightness hyperparameter, defined in (6.18.21)     
            % delta2: float
            %     hyperparameter for identification coefficients, defined in (6.18.23)
            % 
            % returns:
            % inv_U: cell of size (n)
            %     cell of inverse prior variances, defined in (6.18.30) 
            % inv_U_h: cell of size (n)
            %     cell of inverse prior means, defined in (6.18.30)   
        
            % prior mean g
            h = zeros(n,l);
            % prior variance U
            U = zeros(n,m);
            % variance for Lambda0
            U(:,1:m) = (2 * delta1)^2;
            % variance for lagged values
            for j=1:q
                temp = ones(n,m) * (delta1 / j)^2;
                U = [U temp];
            end
            inv_U = cell(n,1);
            inv_U_h = cell(n,1);
            for i=1:n
                inv_u_i = 1 ./ U(i,:);
                inv_U{i} = diag(inv_u_i);
                h_i = h(i,:);
                inv_U_h{i} = (inv_u_i .* h_i)';
            end
        end


        function [inv_V] = make_beta_prior(m, p, pi1, pi2, pi3)
        
            % make_beta_prior(q, p, pi1, pi2, pi3)
            % generate prior terms for VAR coefficients beta
            % 
            % parameters:
            % m: int
            %     number of fundamental factors in the model, defined in (6.18.1) 
            % p: int
            %     number of lags in the VAR model, defined in (6.18.6)  
            % pi1: float
            %     overall tightness hyperparameter, defined in (6.18.24)     
            % pi2: float
            %     cross-variable shrinkage hyperparameter, defined in (6.18.24)
            % pi3: float
            %     lag decay hyperparameter, defined in (6.18.24)
            % 
            % returns:
            % inv_V: cell of size (m)
            %     cell of inverse prior variance, defined in (6.18.33)     
                    
            V = (pi1 * kron(1./(1:p)'.^pi3, pi2 * ones(m,m) + (1-pi2) * eye(m)))'.^2;
            inv_V = cell(m,1);
            for j=1:m
                inv_v_j = 1 ./ V(j,:);
                inv_V{j} = diag(inv_v_j);
            end
        end


        function [inv_Q] = make_gamma_prior(r, omega1)

            % make_gamma_prior(r, omega1)
            % generate prior terms for residual AR coefficients gamma
            % 
            % parameters:
            % r: int
            %     number of lags in the residual AR process, defined in (6.18.7) 
            % omega1: float
            %     overall tightness hyperparameter, defined in (6.18.26)  
            % 
            % returns:
            % inv_Q: matrix of size (r,r)
            %     inverse prior variance, defined in (6.18.36)   
            
            q = (omega1 ./ (1:r)').^2;
            inv_q = 1 ./ q;
            inv_Q = diag(inv_q);
        end


        function [W F Z eps_ E_ E xi e] = make_state_regressors(beta, gamma, z, m, q, p, s, n, r, T)

            % make_state_regressors(beta, gamma, z, m, q, p, s, n, r, T)
            % make factor and residual regressors for dynamic factor models
            % 
            % parameters:
            % beta: matrix of size (k,m)
            %     matrix of coefficients beta_j, defined in (6.18.16)
            % gamma: matrix of size (n,r)
            %     matrix of residual AR coefficients gamma_i, defined in (6.18.13)            
            % z: matrix of size (T,m(s+1)+n(r+1))
            %     array of state variables, defined in (6.18.39)
            % m: int
            %     number of fundamental factors in the model, defined in (6.18.1) 
            % q: int
            %     number of lags in the loadings equation, defined in (6.18.2)          
            % p: int
            %     number of lags in the VAR model, defined in (6.18.6)  
            % s: int
            %     maximum number of factor lags, defined in (6.18.38)
            % n: int
            %     number of endogenous variables, defined in (6.18.1)  
            % r: int
            %     number of lags in the residual AR process, defined in (6.18.7) 
            % T: int
            %     total number of sample periods, defined in (6.18.1)            
            % 
            % returns:
            % W: matrix of size (T,m(q+1))
            %     factor regressors for loadings, defined in (6.18.11)   
            % F: matrix of size (T,m)
            %     factor regressors for VAR, defined in (6.18.15)
            % Z: matrix of size (T,m(p+1))
            %     lagged factor regressors for VAR, defined in (6.18.15)
            % eps_: cell of size (n)
            %     cell of residuals, defined in (6.18.13)   
            % E_: cell of size (n)
            %     cell of lagged residual regressors, defined in (6.18.13) 
            % E: matrix of size (T,n)
            %     current period residuals, defined in (6.18.13)
            % xi: matrix of size (T,m)
            %     factor VAR residuals, defined in (6.18.15) 
            % e: matrix of size (T,n)
            %     residual AR disturbances, defined in (6.18.13)     

            W = z(:,1:m*(q+1));
            F = z(:,1:m);
            Z = z(:,m+1:m*(p+1));
            xi = F - Z * beta;
            eps = z(:,m*(s+1)+1:end);
            E = eps(:,1:n);
            eps_ = cell(n,1);
            E_ = cell(n,1);
            e = zeros(T,n);
            index = n * (1:r)';
            for i=1:n
                eps_{i} = eps(:,i);
                E_{i} = eps(:,index+i);
                e(:,i) = eps_{i} - E_{i} * gamma(i,:)';
            end
        end

        
        function [Lambda B Upsilon] = dfm_state_space_representation(lamda, beta, gamma, sigma, omega, n, m, q, p, s, r)
            
            % dfm_state_space_representation(lamda, beta, gamma, sigma, omega, n, m, q, p, s, r)
            % create state-space representation matrices
            % 
            % parameters:
            % lamda: matrix of size (n,l)
            %     matrix of coefficients lambda_i, defined in (6.18.10)
            % beta: matrix of size (k,m)
            %     matrix of coefficients beta_j, defined in (6.18.16)
            % gamma: matrix of size (n,r)
            %     matrix of residual AR coefficients gamma_i, defined in (6.18.13)
            % sigma: float
            %     variance on residuals, defined in (6.18.7)
            % omega: float
            %     variance on factors, defined in (6.18.6)             
            % n: int
            %     number of endogenous variables, defined in (6.18.1)          
            % m: int
            %     number of fundamental factors in the model, defined in (6.18.1)        
            % q: int
            %     number of lags in the loadings equation, defined in (6.18.2)          
            % p: int
            %     number of lags in the VAR model, defined in (6.18.6)  
            % s: int
            %     maximum number of factor lags, defined in (6.18.38)        
            % r: int
            %     number of lags in the residual AR process, defined in (6.18.7) 
            % 
            % returns:
            % Lambda: matrix of size (n,m(s+1)+n(r+1))
            %     equivalent of A_t matrix, defined in (6.18.39)
            % B: matrix of size (m(s+1)+n(r+1),m(s+1)+n(r+1))
            %     companion form matrix for VAR dynamics, defined in (6.18.39) 
            % Upsilon: matrix of size (m(s+1)+n(r+1),m(s+1)+n(r+1))
            %     covariance matrix for VAR dynamics, defined in (6.18.39)         

            temp_1 = [lamda zeros(n,m*(s-q))];
            temp_2 = [eye(n) zeros(n,n*r)];
            Lambda = [temp_1 temp_2];
            temp_1 = [beta' zeros(m,(s-p+1)*m) zeros(m,(r+1)*n)];
            temp_2 = [eye(s*m) zeros(s*m,m) zeros(s*m,(r+1)*n)];
            temp_3 = zeros(n,(s+1)*m);
            for j=1:r
                temp_3 = [temp_3 diag(gamma(:,j))];
            end
            temp_3 = [temp_3 zeros(n,n)];
            temp_4 = [zeros(n*r,(s+1)*m) eye(r*n) zeros(n*r,n)];
            B = [temp_1;temp_2;temp_3;temp_4];
            Upsilon = diag([ones(m,1)*omega;zeros(s*m,1);ones(n,1)*sigma;zeros(r*n,1)]);
        end


        function [z_00 Upsilon_00] = dfm_kalman_initial_values( sigma, omega, s, r, m, n)

            % dfm_kalman_initial_values( sigma, omega, s, r, m, n)
            % create initial values for state-space representation, defined in (6.18.39)
            % 
            % parameters:
            % sigma: float
            %     variance on residuals, defined in (6.18.7)
            % omega: float
            %     variance on factors, defined in (6.18.6)              
            % q: int
            %     number of fundamental factors in the model, defined in (6.18.1)
            % s: int
            %     maximum number of factor lags, defined in (6.18.38)        
            % r: int
            %     number of lags in the residual AR process, defined in (6.18.7) 
            % m: int
            %     number of fundamental factors in the model, defined in (6.18.1)  
            % n: int
            %     number of endogenous variables, defined in (6.18.1) 
            % 
            % returns:
            % z_00: matrix of size (m(s+1)+n(r+1),)
            %     vector of mean for initial factor value
            % Upsilon_00: matrix of size (m(s+1)+n(r+1),m(s+1)+n(r+1))
            %     matrix of variance for initial factor value        
        
            z_00 = zeros(m*(s+1)+n*(r+1),1);
            Upsilon_00 = diag([5*ones(m,1)*omega;zeros(s*m,1);ones(n,1)*sigma;zeros(r*n,1)]) ...
                         + 1e-10 * eye(m*(s+1)+n*(r+1));
        end


        function [A B Upsilon] = epsilon_state_space_representation(gamma, sigma, l, r)
            
            % epsilon_state_space_representation(gamma, sigma, l, r)
            % create state-space representation matrices for the residuals
            % 
            % parameters:
            % gamma: matrix of size (l,r)
            %     matrix of residual AR coefficients gamma_i for first l variables, defined in (6.18.13)         
            % sigma: float
            %     variance on residuals, defined in (6.18.7) 
            % l: int
            %     total number of loadings regressors, defined in (6.18.11)         
            % r: int
            %     number of lags in the residual AR process, defined in (6.18.7)      
            % 
            % returns:
            % A: matrix of shape (l,l(r+1))
            %     equivalent of A_t matrix
            % B: matrix of shape (l(r+1),l(r+1))
            %     companion form matrix B for residual dynamics 
            % Upsilon: matrix of shape (l(r+1),l(r+1))
            %     covariance matrix for residual dynamics              
        
            A = [eye(l) zeros(l,l*r)];
            if r == 0
                B = zeros(l,l);
            elseif r > 0
                temp_1 = zeros(l,l*(r+1));
                for j=1:r
                    temp_1(:,(j-1)*l+1:j*l) = diag(gamma(:,j));
                end
                temp_2 = [eye(r*l) zeros(l*r,l)];
                B = [temp_1; temp_2];
            end
            Upsilon = diag([ones(l,1) * sigma; zeros(r*l,1)]);    
        end


        function [z_00 Upsilon_00] = epsilon_kalman_initial_values(sigma, r, l)
            
            % eps_kalman_initial_values(sigma, r, l)
            % create initial values for state-space representation  
            % 
            % parameters:        
            % sigma: float
            %     variance on residuals, defined in (6.18.7) 
            % l: int
            %     total number of loadings regressors, defined in (6.18.11)         
            % r: int
            %     number of lags in the residual AR process, defined in (6.18.7)   
            % 
            % returns:
            % z_00: matrix of size (l(r+1),)
            %     vector of mean for initial factor value
            % Upsilon_00: matrix of size (l(r+1),l(r+1))
            %     matrix of variance for initial factor value               
        
            z_00 = zeros(l*(r+1),1);
            Upsilon_00 = diag([3*ones(l,1)*sigma; zeros(r*l,1)]);
        end


        function [eps_ E_ E e] = update_epsilon_regressors(gamma, z, l, r, T)
            
            % update_epsilon_regressors(gamma, z, l, r, T)
            % update residual regressors  
            % 
            % parameters:        
            % z : matrix of size (T,l(r+1))
            %     matrix of sampled values for the state variables   
            % l: int
            %     total number of loadings regressors, defined in (6.18.11)         
            % r: int
            %     number of lags in the residual AR process, defined in (6.18.7)     
            % T: int
            %     total number of sample periods, defined in (6.18.1)
            % 
            % returns:        
            % eps_: struct of size (n)
            %     structure of residuals, defined in (6.18.13)   
            % E_: struct of size (n)
            %     structure of lagged residual regressors, defined in (6.18.13) 
            % E: matrix of size (T,n)
            %     current period residuals, defined in (6.18.13)     
            % e: matrix of size (T,n)
            %     residual AR disturbances, defined in (6.18.13)               
        
            E = z(:,1:l);
            eps_ = cell(l,1);
            E_ = cell(l,1);
            e = zeros(T,l);
            index = l * (1:r)';
            for i=1:l
                eps_{i} = z(:,i);
                E_{i} = z(:,index+i);
                e(:,i) = eps_{i} - E_{i} * gamma(i,:)';
            end
        end


        function [posterior_estimates] = posterior_estimates(X, credibility_level)
            
            % posterior_estimates(X, credibility_level)
            % median, lower bound and upper bound of credibility interval
            % 
            % parameters:
            % X : matrix of shape (n,m,iterations)
            %     matrix of MCMC draws
            % credibility_level : float between 0 and 1
            %     credibility level for credibility interval
            % 
            % returns:
            % posterior_estimates : matrix of shape (n,m,4)
            %     matrix of posterior estimates
        
            posterior_estimates = zeros(size(X,1),size(X,2),4);
            posterior_estimates(:,:,1) = quantile(X,0.5,3);
            posterior_estimates(:,:,2) = quantile(X,(1-credibility_level)/2,3);
            posterior_estimates(:,:,3) = quantile(X,(1+credibility_level)/2,3);
            posterior_estimates(:,:,4) = std(X,0,3);
        end


        function [X_p F_p] = dfm_forecast(lamda, beta, gamma, sigma, omega, f, eps, h, m, n, p, r, l)
        
            % dfm_forecast(lamda, beta, gamma, sigma, omega, f, eps, h, m, n, p, r, l)
            % forecasts for the dynamic factor model, using algorithm 18.3
            % 
            % parameters:
            % lamda: matrix of size (n,l)
            %     matrix of coefficients lambda_i, defined in (6.18.10)
            % beta: matrix of size (k,m)
            %     matrix of coefficients beta_j, defined in (6.18.16)
            % gamma: matrix of size (n,r)
            %     matrix of residual AR coefficients gamma_i, defined in (6.18.13) 
            % sigma: float
            %     variance on residuals, defined in (6.18.7)
            % omega: float
            %     variance on factors, defined in (6.18.6)  
            % f: matrix of size (l,1)
            %     matrix of sample factors     
            % eps: ndarray of size (n*(1+r),1)
            %     matrix of sample residuals
            % h: int
            %     number of forecast periods        
            % m: int
            %     number of fundamental factors in the model, defined in (6.18.1)        
            % n: int
            %     number of endogenous variables, defined in (6.18.1)          
            % p: int
            %     number of lags in the VAR model, defined in (6.18.6)  
            % r: int
            %     number of lags in the residual AR process, defined in (6.18.7)  
            % l: int
            %     total number of loadings regressors, defined in (6.18.11) 
            % 
            % returns:
            % posterior_estimates : matrix of size (n,m,4)
            %     matrix of posterior estimates
        
            % initiate storage and shocks
            X_p = zeros(h,n);
            F_p = zeros(h,m);
            Xi = sqrt(omega) * randn(h,m);
            e = sqrt(sigma) * randn(h,n);
            for t=1:h
                % predict factors using (6.18.15)
                Z = f(1:p*m);
                F = Z * beta + Xi(t,:);
                f = [F f];
                % predict residuals using (6.18.12)
                E = reshape(eps(1:n*r),[n r]);
                eps_ = sum(E .* gamma,2)' + e(t,:);
                eps = [eps_ eps];
                % predict endogenous using (6.18.10)
                W = f(1:l);
                x = W * lamda' + eps_;
                X_p(t,:) = x;
                F_p(t,:) = F;
            end
        end


        function [irf] = dfm_impulse_response_function(lamda, beta, gamma, n, m, q, p, r, h)
            
            % dfm_impulse_response_function(lamda, beta, gamma, n, m, q, p, r, h)
            % impulse response function for dfm, using algorithm 18.4
            % 
            % parameters:
            % lamda: matrix of size (n,l)
            %     matrix of coefficients lambda_i, defined in (6.18.10)
            % beta: matrix of size (k,m)
            %     matrix of coefficients beta_j, defined in (6.18.16)
            % gamma: matrix of size (n,r)
            %     matrix of residual AR coefficients gamma_i, defined in (6.18.13) 
            % n: int
            %     number of endogenous variables, defined in (6.18.1)            
            % m: int
            %     number of fundamental factors in the model, defined in (6.18.1)  
            % q: int
            %     number of fundamental factors in the model, defined in (6.18.1)        
            % p: int
            %     number of lags in the VAR model, defined in (6.18.6)  
            % r: int
            %     number of lags in the residual AR process, defined in (6.18.7) 
            % h: int
            %     number of irf periods 
            % 
            % returns:
            % irf : matrix of size (n,m+1,h)
            %     matrix of dfm impulse response function
        
            irf = zeros(n,m+1,h);
            [factor_irf] = nwu.factor_impulse_response_function(beta, m, p, h);
            [residual_irf] = nwu.residual_impulse_response_function(gamma, n, r, h);
            Lambda_0 = lamda(:,1:m);
            full_factor_irf = Lambda_0 * factor_irf;
            for i=1:q
                Lambda_i = lamda(:,i*m+1:(i+1)*m);
                factor_irf_i = [zeros(m,i*m) factor_irf(:,1:end-m*i)];
                full_factor_irf = full_factor_irf + Lambda_i * factor_irf_i;
            end
            for i=1:h
                irf(:,1:m,i) = full_factor_irf(:,(i-1)*m+1:i*m);
                irf(:,m+1,i) = residual_irf(:,i);
            end
        end


        function [irf] = factor_impulse_response_function(beta, m, p, h)
            
            % factor_impulse_response_function(beta, m, p, h)
            % impulse response function for factor VAR
            % 
            % parameters:
            % beta: matrix of size (k,m)
            %     matrix of coefficients beta_j, defined in (6.18.16)           
            % m: int
            %     number of fundamental factors in the model, defined in (6.18.1)         
            % p: int
            %     number of lags in the VAR model, defined in (6.18.6)  
            % h: int
            %     number of irf periods 
            % 
            % returns:
            % irf : matrix of size (m,m,h)
            %     matrix of factor impulse response function
            
            Yh = eye(m);
            irf = {Yh};
            Xh = zeros(m,m*p);   
            for i=2:h
                Xh = [Yh Xh(:,1:end-m)];
                Yh = Xh * beta;
                irf = [irf Yh];
            end
            irf = cell2mat(irf);
        end


        function [irf] = residual_impulse_response_function(gamma, n, r, h)
        
            % residual_impulse_response_function(gamma, n, r, h)
            % impulse response function for residual AR models
            % 
            % parameters:
            % gamma: ndarray of shape (n,r)
            %     matrix of residual AR coefficients gamma_i, defined in (6.18.13) 
            % n: int
            %     number of endogenous variables, defined in (6.18.1)            
            % r: int
            %     number of lags in the residual AR process, defined in (6.18.7) 
            % h: int
            %     number of irf periods 
            % 
            % returns:
            % irf : ndarray of shape (n,h)
            %     matrix of residual impulse response function
            
            Yh = ones(n,1);
            irf = [Yh zeros(n,h-1)];
            Xh = zeros(n,r);
            for i=2:h
                Xh = [Yh Xh(:,1:r-1)];
                Yh = sum(Xh .* gamma, 2);
                irf(:,i) = Yh;
            end
        end


        function [posterior_estimates] = posterior_estimates_3d(X, credibility_level)

            % posterior_estimates(X, credibility_level)
            % median, lower bound and upper bounf of credibility interval
            % 
            % parameters:
            % X : matrix of size (n,m,h,iterations)
            %     matrix of MCMC draws
            % credibility_level : float between 0 and 1
            %     credibility level for credibility interval
            % 
            % returns:
            % posterior_estimates : ndarray of shape (n,m,h,3)
            %     matrix of posterior estimates
        
            posterior_estimates = zeros(size(X,1),size(X,2),size(X,3),3);
            posterior_estimates(:,:,:,1) = quantile(X,0.5,4);
            posterior_estimates(:,:,:,2) = quantile(X,(1-credibility_level)/2,4);
            posterior_estimates(:,:,:,3) = quantile(X,(1+credibility_level)/2,4); 
        end  


        function [fevd] = dfm_forecast_error_variance_decomposition(irf, sigma, omega, n, m, h)

            % dfm_forecast_error_variance_decomposition(irf, sigma, omega, n, m, h)
            % products FEVD for the dynamic factor model, using (6.18.50)
            % 
            % parameters:
            % irf : matrix of size (n,m+1,h)
            %     matrix of impulse response functions
            % sigma: float
            %     variance on residuals, defined in (6.18.7)
            % omega: float
            %     variance on factors, defined in (6.18.6)      
            % n : int
            %     number of endogenous variables  
            % m: int
            %     number of fundamental factors in the model, defined in (6.18.1)         
            % h : int
            %     number of forecast periods
            % 
            % returns:
            % fevd : matrix of size (n,m+1,h)
            %     matrix of forecast error variance decomposition

            cum_squared_irf = cumsum(irf.^2, 3);
            variances = repmat([omega*ones(n,m) sigma*ones(n,1)], [1 1 h]);
            cum_squared_irf =  variances .* cum_squared_irf;
            total_variance = repmat(sum(cum_squared_irf,2), [1 m+1 1]);
            fevd = cum_squared_irf ./ total_variance;
        end


        function [normalized_fevd_estimates] = normalize_fevd_estimates(fevd_estimates)

            % [normalized_fevd_estimates] = normalize_fevd_estimates(fevd_estimates)
            % normalizes FEVD estimates so that they sum up to 1
            % 
            % parameters:
            % fevd_estimates : matrix of size (n,n,h,iterations)
            %     matrix of posterior FEVD estimates
            % 
            % returns:
            % normalized_fevd_estimates : matrix of size (n,n,h,iterations)
            %     matrix of normalized posterior FEVD estimates

            point_estimate_contribution = sum(fevd_estimates(:,:,:,1), 2);
            total_contribution = repmat(point_estimate_contribution, [1 size(fevd_estimates,2)]);
            estimates_contribution = repmat(total_contribution, [1 1 1 3]);
            normalized_fevd_estimates = fevd_estimates ./ estimates_contribution;
        end


        function [hd] = dfm_historical_decomposition(irf, xi, e, n, m, T)

            % dfm_historical_decomposition(irf, xi, e, n, m, T)
            % products historical decomposition for the dynamic factor model, using (6.18.51)-(6.18.53)
            % 
            % parameters:
            % irf : matrix of size (n,m+1,h)
            %     matrix of impulse response functions
            % xi : matrix of size (T,m)
            %     matrix of factor shocks    
            % e : matrix of size (T,n)
            %     matrix of residual shocks              
            % n : int
            %     number of endogenous variables  
            % m: int
            %     number of fundamental factors in the model, defined in (6.18.1)         
            % T: int
            %     total number of sample periods, defined in (6.18.1)
            % 
            % returns:
            % hd : matrix of size (n,m+1,T)
            %     matrix of historical decomposition

            reshaped_xi = permute(repmat(flip(xi,1),[1 1 n]), [3 2 1]);
            reshaped_e = permute(reshape(flip(e,1), [T n 1]), [2 3 1]);
            reshaped_shocks = cat(2,reshaped_xi,reshaped_e);
            hd = zeros(n,m+1,T);
            for j=1:T
                hd(:,:,j) = sum(irf(:,:,1:j) .* reshaped_shocks(:,:,end-j+1:end), 3);
            end
        end


        function [insample_evaluation] = dfm_insample_evaluation_criteria(Y, E, T, k)
            
            % [insample_evaluation] = dfm_insample_evaluation_criteria(Y, E, T, k)
            % computes ssr, R2 and adjusted R2 for each VAR equation
            %
            % parameters:
            % Y : matrix of size (T,n)
            %     matrix of endogenous variables
            % E : matrix of size (T,n)
            %     matrix of residuals
            % T : int
            %     number of sample periods
            % k : int
            %     number of coefficients in each VAR equation
            %
            % returns:
            % insample_evaluation : struct
            %     structure of insample evaluation criteria
                      
            ssr = diag(E' * E);
            Z = Y - mean(Y);
            tss = diag(Z' * Z);
            r2 = 1 - ssr ./ tss;
            adj_r2 = 1 - (1 - r2) * (T - 1) / (T - k);
            insample_evaluation = struct;
            insample_evaluation.ssr = ssr;
            insample_evaluation.r2 = r2;
            insample_evaluation.adj_r2 = adj_r2;
        end


        function [F, mu_, Upsilon] = make_mfbvar_state_regressors(Z_, B, Sigma, n, m, p, T_)

            % make_mfbvar_state_regressors(Z_, B, Sigma, n, m, p, T_)
            % make state-space parameters for MF-BVAR model, defined in (6.17.20)-(6.17.22)
            % 
            % parameters:
            % Z_ : matrix of size (T_,m)
            %     matrix of endogenous variables, defined in (6.17.1), including initial conditions
            % B : matrix of size (k,n)
            %     matrix of VAR coefficients, defined in (6.17.2)
            % Sigma : matrix of size (n,n)
            %     variance-covariance matrix of residuals, defined in (6.17.1)
            % n : int
            %     number of endogenous variables, defined in (6.17.1)
            % m : int
            %     number of exogenous variables, defined in (6.17.1)
            % p : int
            %     number of lags in the VAR model, defined in (6.17.1)
            % T_ : int
            %     number of sample periods, including initial conditions
            % 
            % returns:
            % F : matrix of size (n*p,n*p)
            %     companion matrix for VAR coefficients
            % mu_ : matrix of size (n*p,)
            %     exogenous vector of state equation
            % Upsilon : matrix of size (n*p,n*p)
            %     variance covariance matrix of state equation     

            F = diag(ones(n*(p-1),1),-n);
            F(1:n,:) = B(m+1:end,:)';
            C = B(1:m,:);
            mu_ = [(Z_*C) zeros(T_,(p-1)*n)];
            Upsilon = 1e-12 * eye(n*p);
            Upsilon(1:n,1:n) = Sigma;
        end


        function [gamma_00 Upsilon_00] = mfbvar_kalman_initial_values(Sigma, n, p)
        
            % mfbvar_kalman_initial_values(Sigma, n, p)
            % create initial values for state-space representation, defined in (6.17.20)-(6.17.21)
            % 
            % parameters:
            % Sigma : matrix of size (n,n)
            %     variance-covariance matrix of residuals, defined in (6.17.1)
            % n : int
            %     number of endogenous variables, defined in (6.17.1)
            % p : int
            %     number of lags in the VAR model, defined in (6.17.1)
            % 
            % returns:
            % gamma_00 : matrix of size (n*p,)
            %     initial values for state vector
            % Upsilon_00 : matrix of size (n*p,n*p)
            %     initial value for variance covariance matrix of state vector  
            
            gamma_00 = zeros(n*p,1);    
            Upsilon_00 = 1e-12 * eye(n*p);
            Upsilon_00(1:n,1:n) = 10 * Sigma;
        end


        function [W X] = update_mfbvar_state_regressors(Z, Gamma, n, p)
            
            % update_mfbvar_state_regressors(Z, Gamma, n, p)
            % define state regressors W and X from Gamma, as defined in (6.17.3)
            % 
            % parameters:
            % Z : ndarray of shape (T,m)
            %     matrix of endogenous variables, defined in (6.17.1)        
            % Gamma : ndarray of shape (T_,n*r)
            %     full matrix of state vectors gamma, including initial conditions
            % n : int
            %     number of endogenous variables, defined in (6.17.1)
            % p : int
            %     number of lags in the VAR model, defined in (6.17.1)
            % 
            % returns:
            % W : ndarray of shape (T,n)
            %     current period states
            % X : ndarray of shape (T,k)
            %     matrix of lagged and endogenous regressors, defined in (6.17.2)
            
            W = Gamma(p+1:end,1:n);
            X = [Z Gamma(p+1:end,n+1:p*n) Gamma(1:end-p,1:n)];
        end


    end
    
    
end
        