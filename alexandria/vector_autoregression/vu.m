classdef vu
    

    % vu stands for var utilities
    % A class containing static methods for vector autoregression utilities

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
        
        
        function [X] = generate_intercept_and_trends(constant, trend, quadratic_trend, periods, shift)
            
            % [X] = generate_intercept_and_trends(constant, trend, quadratic_trend, periods, shift)
            % generates automated exogenous regressors (constant, trend, quadratic trend)
            %
            % parameters:
            % constant : bool
            %     if true, a constant is added to the VAR model
            % trend : bool
            %     if true, a linear trend is added to the VAR model
            % quadratic_trend : bool
            %     if true, a quadratic trend is added to the VAR model
            % periods : int
            %     number of periods for which regressors are generated
            % shift : int
            %     number of periods to add to regressors (e.g. for predictions)
            %
            % returns:
            % X : matrix of shape (periods,)
            %     matrix of automated exogenous regressors           

            X = [];
            if constant
               constant_column = ones(periods,1);
               X = [X constant_column];
            end
            if trend
               trend_column = shift + (1:periods)';
               X = [X trend_column];
            end
            if quadratic_trend
               quadratic_trend_column = (shift + (1:periods)).^2';
               X = [X quadratic_trend_column];
            end  
        end
        

        function [X] = generate_exogenous_regressors(exogenous, lags)
            
            % [X] = generate_exogenous_regressors(exogenous, lags)
            % generate matrix of in-sample exogenous regressors
            %
            % parameters:
            % exogenous : matrix of size (periods,n_exogenous)
            %     matrix of exogenous regressors
            % lags : int
            %     number of lags in VAR model
            %
            % returns:
            % X : matrix of shape (periods,n_exogenous)
            %     matrix of exogenous regressors              
            
            if isempty(exogenous)
                X = [];
            else
                X = exogenous(lags+1:end,:);
            end
        end
        
        
        function [X] = generate_lagged_endogenous(endogenous, lags)
            
            % [X] = generate_lagged_endogenous(self, endogenous, lags)
            % generate in-sample matrix of lagged endogenous regressors
            %
            % parameters:
            % endogenous : matrix of shape (periods,n_endogenous)
            %     matrix of endogenous regressors
            % lags : int
            %     number of lags in VAR model
            %
            % returns:
            % X : matrix of shape (periods,n_endogenous)
            %     matrix of endogenous regressors 

            X = [];
            for i=1:lags
                X = [X endogenous(lags+1-i:end-i,:)];
            end
        end
        
        
        function [E, Y_hat] = fit_and_residuals(Y, X, B)
            
            % [E, Y_hat] = fit_and_residuals(Y, X, B)
            % fit_and_residuals(Y, X, B)
            % generates fit and residuals for a VAR model, using (4.11.2)
            %
            % parameters:
            % Y : matrix of size (T,n)
            %     matrix of endogenous variables
            % X : matrix of size (T,k)
            %     matrix of VAR regressors
            % B : matrix of size (k,n)
            %     matrix of VAR coefficients
            %
            % returns:
            % Y_hat : matrix of size (T,n)
            %     matrix of fitted endogenous
            % E : matrix of size (T,n)
            %     matrix of VAR residuals
            
            Y_hat = X * B;
            E = Y - Y_hat;         
        end
        
        
        function [insample_evaluation] = insample_evaluation_criteria(Y, E, T, k)
            
            % [insample_evaluation] = insample_evaluation_criteria(Y, E, T, k)
            % computes ssr, R2 and adjusted R2 for each VAR equation
            %
            % parameters:
            % Y : matrix of size (T,n)
            %     matrix of endogenous variables
            % E : matrix of size (T,n)
            %     matrix of VAR residuals
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
        
        
        function [steady_state] = steady_state(Z, B, n, m, p, T)
            
            % [steady_state] = steady_state(Z, B, n, m, p, T)
            % computes the steady-state of the VAR model, using equation (4.12.30)
            %
            % parameters:
            % Z : matrix of size (T,m)
            %     matrix of exogenous variables
            % B : matrix of size (m+n*p,n)
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
                C = B(1:m,:);
                A = eye(n);
                for i=1:p
                    A = A - B(1+m+(i-1)*n:m+i*n,:);
                end
                steady_state = Z * C / A;
            end
        end        

        
        function [Xi] = structural_shocks(E, inv_H)
            
            % [Xi] = structural_shocks(E, inv_H)
            % computes the structural shocks of the VAR model, using equation (4.13.9)
            %
            % parameters:
            % E : matrix of size (T,n)
            %     matrix of VAR residuals
            % inv_H : matrix of size (n,n)
            %     inverse of structural matrix H
            %
            % returns:
            % Xi : matrix of size (n,n)
            %     matrix of structural shocks
                      
            Xi = E * inv_H';
        end         
        
        
        function [Z_p Y] = make_forecast_regressors(Z_p, Y, h, p, T, exogenous, constant, trend, quadratic_trend)       

            % [Z_p, Y] = make_forecast_regressors(Z_p, Y, periods, p, T, constant, trend, quadratic_trend)
            % create regressors for forecast estimation
            %
            % parameters:
            % Z_p : matrix of size (periods,m) or empty matrix
            %     matrix of predicted exogenous            
            % Y : matrix of size (T,n)
            %     matrix of endogenous variables
            % h : int
            %     number of forecast periods
            % p : int
            %     number of lags            
            % T : int
            %     number of sample periods
            % exogenous : ndarray of shape (periods,n_exogenous)
            %     matrix of exogenous regressors
            % constant : bool
            %     if true, a constant is added to the VAR model
            % trend : bool
            %     if true, a linear trend is added to the VAR model
            % quadratic_trend : bool
            %     if true, a quadratic trend is added to the VAR model            
            %
            % returns:
            % Z_p : matrix of size (periods,m)
            %     full matrix of endogenous regressors
            % Y : matrix of size (lags,n)
            %     matrix of endogenous variables for initial conditions
        
            temp = vu.generate_intercept_and_trends(constant, trend, quadratic_trend, h, T);
            % if no exogenous, return empty matrix
            if isempty(exogenous)
                Z_p = [];
            elseif isempty(Z_p)
                Z_p = repmat(exogenous(end,:),[h 1]);
            end
            if ~isempty(Z_p)
                Z_p = [temp Z_p];
            elseif any([constant trend quadratic_trend])
                Z_p = temp;
            else
                Z_p = [];
            end
            Y = Y(end-p+1:end,:);  
        end
        
        
        function [forecast_evaluation_criteria] = forecast_evaluation_criteria(Y_hat, Y)
    
            % forecast_evaluation_criteria(Y_hat, Y)
            % estimates RMSE, MAE, MAPE, Theil-U and bias for forecasts
            % 
            % parameters:
            % Y_p : matrix of shape (periods,n)
            %     matrix of predicted endogenous
            % Y : matrix of shape (periods,n)
            %     matrix of actual endogenous values
            % 
            % returns:
            % forecast_evaluation_criteria : struct
            %     structure storing forecast evaluation criteria
    
            % check dimensions
            if ~isequal(size(Y), size(Y_hat))
                error(['Cannot calculate forecast evaluation criteria. Forecasts and actual values have different dimensions.']);
            end
            h = size(Y_hat,1);
            % calculate forecast error
            err = Y - Y_hat;
            % calculate RMSE, MAE and MAPE from (4.13.18)
            rmse = sqrt(sum(err .* err, 1) / h);
            mae = sum(abs(err) / h);
            mape = 100 * sum(abs(err ./ Y),1) / h;
            % calculate Theil-U and bias from (4.13.19)
            theil_u = sqrt(sum(err .* err,1)) ./ (sqrt(sum(Y .* Y,1)) + sqrt(sum(Y_hat .* Y_hat,1)));
            bias = sum(err,1) ./ sum(abs(err),1);
            % store in structure
            forecast_evaluation_criteria = struct;
            forecast_evaluation_criteria.rmse = rmse;
            forecast_evaluation_criteria.mae = mae;
            forecast_evaluation_criteria.mape = mape;
            forecast_evaluation_criteria.theil_u = theil_u;
            forecast_evaluation_criteria.bias = bias;
        end
        
        
        function [forecast_evaluation_criteria] = bayesian_forecast_evaluation_criteria(mcmc_forecast, Y)

            % bayesian_forecast_evaluation_criteria(mcmc_forecast, Y)
            % estimates log scores and CRPS for forecasts
            % 
            % parameters:
            % mcmc_forecast : matrix of size (periods,n,iterations)
            %     matrix of Gibbs sampler forecast values
            % Y : matrix of size (periods,n)
            %     matrix of actual endogenous values
            % 
            % returns:
            % bayesian forecast_evaluation_criteria : structure
            %     structure storing forecast evaluation criteria

        
            % check dimensions
            if ~isequal(size(Y), [size(mcmc_forecast,1) size(mcmc_forecast,2)])
                error(['Cannot calculate forecast evaluation criteria. Forecasts and actual values have different dimensions.']);
            end
            h = size(mcmc_forecast,1);
            n = size(mcmc_forecast,2);
            % initiate storage
            log_pdf = zeros(h,n);
            crps = zeros(h,n);
            % log scores for individual periods
            mu_hat = mean(mcmc_forecast,3);
            sigma_hat = var(mcmc_forecast,0,3);
            for i=1:h
                for j=1:n
                    [log_pdf(i,j) ~] = su.normal_pdf(Y(i,j), mu_hat(i,j), sigma_hat(i,j));
                end
            end
            log_score = - log_pdf;          
            % log scores for joint periods
            if h == 1
                joint_log_score = log_score;
            else
                joint_log_pdf = zeros(n,1);
                flipped_mcmc_forecast = permute(mcmc_forecast,[1 3 2]);
                for i=1:n
                    Sigma = cov(flipped_mcmc_forecast(:,:,i)');
                    [joint_log_pdf(i,1) ~] = su.multivariate_normal_pdf(Y(:,i), mu_hat(:,i), Sigma);
                end
                joint_log_score = - joint_log_pdf;
            end
            % CRPS for individual periods
            for i=1:h
                for j=1:n
                    crps(i,j) = vu.make_crps(Y(i,j), flipped_mcmc_forecast(i,:,j));
                end
            end
            % CRPS for joint periods
            joint_crps = sum(crps);
            % store in structure
            forecast_evaluation_criteria = struct;
            forecast_evaluation_criteria.log_score = log_score;
            forecast_evaluation_criteria.joint_log_score = joint_log_score;
            forecast_evaluation_criteria.crps = crps;
            forecast_evaluation_criteria.joint_crps = joint_crps;
        end
             
            
        function [crps] = make_crps(y, y_hat)

            % make_crps(y, y_hat)
            % continuous rank probability score for prediction y
            % 
            % parameters:
            % y : float
            %     actual (observed) value for the foecast
            % y_hat : matrix of size (1,iteration)
            %     vector of MCMC simulated values for predictions
            % 
            % returns:
            % crps : float
            %     crps value

            J = size(y_hat,2);
            term_1 = sum(abs(y_hat - y)) / J;
            temp = repmat(y_hat, [J 1]);
            term_2 = sum(sum(abs(temp - temp'))) / (2 * J * J);
            crps = term_1 - term_2;
        end
        
        
        function [Y_sum X_sum] = sums_of_coefficients_extension(sums_of_coefficients, pi5, Y, n, m, p)
            
            % sums_of_coefficients_extension(sums_of_coefficients, pi5, Y, n, m, p)
            % generates dummy extension matrices Y_sum and X_sum as defined in (4.12.6)
            % 
            % parameters:
            % sums_of_coefficients : bool
            %     if true, the sums-of-coefficients extension is added to the model
            % pi5 : float
            %     prior shrinkage for the sums-of-coefficients extension
            % Y : matrix of shape (T,n)
            %     matrix of endogenous variables
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables        
            % p : int
            %     number of lags
            % 
            % returns:
            % Y_sum : matrix of shape (n,n) or empty
            %     Y matrix for sums-of-coefficients extension
            % X_sum : matrix of shape (n,m+n*p) or empty
            %     X matrix for sums-of-coefficients extension        

            if sums_of_coefficients
                Y_sum = diag(mean(Y)) / pi5;
                X_sum = [zeros(n,m) kron(ones(1,p),Y_sum)];
            else
                Y_sum = [];
                X_sum = [];
            end
        end
        
        
        function [Y_obs X_obs] =  dummy_initial_observation_extension(dummy_initial_observation, pi6, Y, X, n, m, p)

            % dummy_initial_observation_extension(dummy_initial_observation, pi6, Y, X, n, m, p)
            % generates dummy extension matrices Y_obs and X_obs as defined in (4.12.10)
            %
            % parameters:
            % dummy_initial_observation : bool
            %     if true, the dummy initial observation extension is added to the model
            % pi6 : float
            %     prior shrinkage for the dummy initial observation extension
            % Y : matrix of shape (T,n)
            %     matrix of endogenous variables
            % X : matrix of shape (T,k)
            %     matrix of VAR regressors        
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables        
            % p : int
            %     number of lags
            %
            % returns:
            % Y_obs : matrix of shape (1,n) or empty
            %     Y matrix for dummy initial observation extension
            % X_obs : matrix of shape (1,m+n*p) or empty
            %     X matrix for dummy initial observation extension        

            if dummy_initial_observation
                Y_obs = mean(Y) / pi6;
                X_obs = [mean(X(:,1:m))/pi6 kron(ones(1,p),Y_obs)];
            else
                Y_obs = [];
                X_obs = [];
            end
        end
        
        
        function [Y_lrp X_lrp] = long_run_prior_extension(long_run_prior, pi7, J, Y, n, m, p)

            % long_run_prior_extension(long_run_prior, pi7, H, Y, n, m, p)
            % generates dummy extension matrices Y_lrp and X_lrp as defined in (4.12.16)
            %
            % parameters:
            % long_run_prior : bool
            %     if True, the long run prior extension is added to the model
            % pi7 : float
            %     prior shrinkage for the long run prior extension
            % J : ndarray of shape (n,n)
            %     matrix of long-run relations          
            % Y : ndarray of shape (T,n)
            %     matrix of endogenous variables
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables        
            % p : int
            %    number of lags
            %
            % returns:
            % Y_lrp : matrix of shape (n,n) or empty
            %     Y matrix for long run prior extension
            % X_lrp : matrix of shape (n,m+n*p) or empty
            %    X matrix for long run prior extension        

            if long_run_prior
                Y_lrp = diag(J * mean(Y)' / pi7) / J';
                X_lrp = [zeros(n,m) kron(ones(1,p),Y_lrp)];
            else
                Y_lrp = [];
                X_lrp = [];
            end
        end

        
        function [b] = make_b(delta, n, m, p)
            
            % [b] = make_b(delta, n, m, p)
            % generates Minnesota prior parameter b as defined in (4.11.16)
            %
            % parameters:
            % delta : matrix of size (n,1)
            %     matrix of prior AR coefficients
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables            
            % p : int
            %     number of lags
            %
            % returns:
            % b : matrix of size (q,1)
            %     prior mean for beta     

            b = la.vec([zeros(m,n);diag(delta);zeros(n*(p-1),n)]);
        end        
        
        
        function [V] = make_V(s, pi1, pi2, pi3, pi4, n, m, p)
            
            % [V] = make_V(s, pi1, pi2, pi3, pi4, n, m, p)
            % generates Minnesota prior parameter V as defined in (4.11.17)-(4.11.20)
            %
            % parameters:
            % s : matrix of size (n,1)
            %     matrix of individual residual variances 
            % pi1 : float
            %     overall tightness           
            % pi2 : float
            %     cross-variable shrinkage      
            % pi3 : float
            %     lag decay         
            % pi4 : float
            %     exogenous slackness           
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables            
            % p : int
            %     number of lags
            %
            % returns:    
            % V : matrix of size (q,1)
            %     prior variance (diagonal term only) for beta

            scale = [repmat(s',[m 1]); repmat(repmat(s',[n 1])./repmat(s,[1 n]),[p 1])];
            shrinkage = (pi1 * [pi4 * ones(m,n); kron(1./(1:p)'.^pi3, ...
                         pi2 * ones(n,n) + (1-pi2) * eye(n))]).^2;
            V = la.vec(scale .* shrinkage);
        end        

        
        function [inv_V inv_V_b] = make_V_b_inverse(b, V)
            
            % [inv_V inv_V_b] = make_V_b_inverse(b, V)
            % generates elements inv_V and inv_B * b used in (4.11.15) and (4.11.43)
            %
            % parameters:
            % V : matrix of size (q,1)
            %     prior variance (diagonal term only) for beta
            % b : matrix of size (q,1)
            %     prior mean for beta             
            %
            % returns:    
            % inv_V : matrix of size (q,q)
            %     inverse prior variance for beta
            % inv_V_b : matrix of size (q,1)
            %     product inv_V * b          

            inv_V = diag(1 ./ V);
            inv_V_b = b ./ V;
        end   
        
        
        function [b_bar V_bar inv_V_bar] = minnesota_posterior(inv_V, inv_V_b, XX, XY, inv_Sigma)

            % minnesota_posterior(V, b, XX, XY, inv_Sigma)
            % generates Minnesota posterior parameters b_bar and V_bar as defined in (4.11.15)
            % 
            % parameters:
            % inv_V : matrix of size (q,q)
            %     inverse prior variance for beta
            % inv_V_b : matrix of size (q,1)
            %     product inv_V * b 
            % XX : matrix of size (k,k)
            %     matrix product X' X     
            % XY : matrix of size (k,n)
            %     matrix product X' Y         
            % inv_Sigma : matrix of size (n,n)
            %     inverse of residual variance-covariance matrix Sigma     
            % 
            % returns:
            % b_bar : matrix of size (q,)
            %     posterior mean for beta  
            % V_bar : matrix of size (q,q)
            %     posterior variance for beta 
            % inv_V_bar : matrix of size (q,q)
            %     inverse posterior variance for beta            

            inv_V_bar = inv_V + kron(inv_Sigma, XX);
            V_bar = la.invert_spd_matrix(inv_V_bar);
            b_bar = V_bar * (inv_V_b + la.vec(XY * inv_Sigma));
        end       
        
        
        function [B] = make_B(delta, n, m, p)
            
            % make_B(delta, n, m, p)
            % generates normal-Wishart prior parameter B as defined in (4.11.33)
            %
            % parameters:
            % delta : matrix of size (n,1)
            %     matrix of prior AR coefficients
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables            
            % p : int
            %     number of lags
            %
            % returns:
            % B : matrix of size (k,n)
            %     prior mean for beta     

            B = [zeros(m,n);diag(delta);zeros(n*(p-1),n)];
        end    

        
        function [W] = make_W(s, pi1, pi3, pi4, n, m, p)
            
            % make_W(s, pi1, pi3, pi4, n, m, p)
            % generates normal-Wishart prior parameter W as defined in (4.11.27)
            %
            % parameters:
            % s : matrix of size (n,1)
            %     matrix of individual residual variances 
            % pi1 : float
            %     overall tightness                
            % pi3 : float
            %     lag decay         
            % pi4 : float
            %     exogenous slackness           
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables            
            % p : int
            %     number of lags
            %
            % returns:    
            % W : matrix of size (k,1)
            %     prior variance (diagonal term only) for beta

            scale = [ones(m,1); repmat(1./s,[p 1])];
            shrinkage = (pi1 * [pi4 * ones(m,1); 1./ kron((1:p)', ones(n,1)).^ pi3]).^2;
            W = scale .* shrinkage;
        end
        
        
        function [alpha] = make_alpha(n)
            
            % make_alpha(n)
            % generates normal-Wishart prior parameter alpha as defined in (4.11.30) 
            %
            % parameters:           
            % n : int
            %     number of endogenous variables
            %
            % returns:
            % alpha : float
            %     prior degrees of freedom for Sigma         

            alpha = n + 2;
        end        
        
        
        function [S] = make_S(s)
            
            % make_S(s)
            % generates normal-Wishart prior parameter S as defined in (4.11.30)
            %
            % parameters:
            % s : matrix of size (n,1)
            %     matrix of individual residual variances    
            %
            % returns:  
            % S : matrix of size (n,1)
            %     prior scale (diagonal term only) for Sigma         

            S = s;
        end    
        

        function [B_bar W_bar alpha_bar S_bar alpha_hat S_hat] = ...
                normal_wishart_posterior(B, W, alpha, S, n, T, XX, XY, YY)
            
            % normal_wishart_posterior(B, W, alpha, S, n, T, XX, XY, YY)
            % generates normal-Wishart posterior parameters as defined in (4.11.33) and (4.11.38)
            %
            % parameters:
            % B : matrix of size (k,n)
            %     prior mean for beta 
            % W : matrix of size (k,1)
            %     prior variance (diagonal term only) for beta            
            % alpha : float
            %     prior degrees of freedom for Sigma  
            % S : matrix of size (n,1)
            %     prior scale (diagonal term only) for Sigma    
            % n : int
            %     number of endogenous variables
            % T : int
            %     number of sample periods
            % XX : matrix of size (k,k)
            %     matrix product X' X     
            % XY : matrix of size (k,n)
            %     matrix product X' Y                
            % YY : matrix of size (n,n)
            %     matrix product Y' Y                
            %
            % returns:  
            % B_bar : matrix of shape (k,n)
            %     posterior mean for beta  
            % W_bar : matrix of shape (k,k)
            %     posterior variance for beta 
            % alpha_bar : float
            %     posterior degrees of freedom for Sigma  
            % S_bar : matrix of shape (n,n)
            %     posterior scale for Sigma  
            % alpha_hat : float
            %     posterior degrees of freedom for B
            % S_hat : matrix of shape (n,n)
            %     posterior scale for B        

            inv_W_bar = diag(1 ./ W) + XX;
            W_bar = la.invert_spd_matrix(inv_W_bar);
            B_bar = W_bar * (B ./ W + XY);
            alpha_bar = alpha + T;
            S_bar = diag(S) + YY + B' * (B ./ W) - B_bar' * inv_W_bar * B_bar;
            alpha_hat = alpha + T - n + 1;
            S_hat = S_bar / alpha_hat;
        end         
        

        function [new_b new_V] = make_constrained_coefficients(B, V, n, m, k, ...
                 lags, constant, trend, quadratic_trend, constrained_coefficients_table)
            
            % make_constrained_coefficients(B, V, n, m, k, ...
            % lags, constant, trend, quadratic_trend, constrained_coefficients_table)
            %
            % parameters:
            % B : matrix of size (k,n)
            %     prior mean for beta (reshaped as k * n)
            % V : matrix of size (k,n)
            %     prior variance for beta (reshaped as k * n)     
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables
            % k : int
            %     number of coefficients in each VAR equation            
            % lags : int
            %     number of lags            
            % constant : bool
            %     if true, a constant is added to the VAR model
            % trend : bool
            %     if true, a linear trend is added to the VAR model
            % quadratic_trend : bool
            %     if true, a quadratic trend is added to the VAR model             
            % constrained_coefficients_table : matrix of size (n_constraints,5)
            %     table defining the constraints on VAR coefficients            
            % 
            % returns:
            % new_b : matrix of size (q,1)
            %     prior mean for beta with constraints applied
            % new_V : matrix of size (q,1)
            %     prior variance for beta with constraints applied
            
            new_B = B;
            new_V = V;
            for j=1:size(constrained_coefficients_table,1)
                variable = constrained_coefficients_table(j,1);
                responding = constrained_coefficients_table(j,2);
                lag = constrained_coefficients_table(j,3);
                mean = constrained_coefficients_table(j,4);
                variance = constrained_coefficients_table(j,5);
                if responding == 0.1
                    new_B(1,variable) = mean;
                    new_V(1,variable) = variance;
                elseif responding == 0.2
                    new_B(1+constant,variable) = mean;
                    new_V(1+constant,variable) = variance; 
                elseif responding == 0.3
                    new_B(1+constant+trend,variable) = mean;
                    new_V(1+constant+trend,variable) = variance;
                elseif responding < 0
                    responding = - responding;
                    new_B(constant+trend+quadratic_trend+responding,variable) = mean;
                    new_V(constant+trend+quadratic_trend+responding,variable) = variance; 
                else
                    row = m + (lag-1) * n + responding;
                    new_B(row,variable) = mean;
                    new_V(row,variable) = variance; 
                end
            end
            new_b = la.vec(new_B);
            new_V = la.vec(new_V);
        end
        

        function [stationary] = check_stationarity(B, n, m, p)

            % check_stationarity(B, n, m, p)
            % check for stability of VAR model as in definition 12.1, using companion form (4.12.27)
            % 
            % parameters:
            % B : matrix of shape (m+n*p,n)
            %     matrix of VAR coefficients
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables        
            % p : int
            %     number of lags
            % 
            % returns:
            % stationary : bool
            %     if true, the VAR model is stationary

            F = vu.make_companion_form(B, n, m, p);
            eigenvalues =  sort(abs(eig(F)),'descend');
            max_eigenvalue = eigenvalues(1);
            if max_eigenvalue < 0.999
                stationary = true;
            else
                stationary = false;
            end
        end
         
        
        function [F] = make_companion_form(B, n, m, p)

            % make_companion_form(B, n, m, p)
            % creates companion form matix F as defined in (4.12.27)-(4.12.28)
            % 
            % parameters:
            % B : matrix of shape (m+n*p,n)
            %     matrix of VAR coefficients
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables        
            % p : int
            %     number of lags
            % 
            % returns:
            % F : matrix of size (n*p,n*p)
            %     companion form matrix

            block_1 = B(m+1:end,:)';
            block_2 = eye(n*(p-1));
            block_3 = zeros(n*(p-1),n);
            F = [block_1; block_2 block_3];
        end        
        

        function [irf] = impulse_response_function(B, n, p, h)

            % [irf] = impulse_response_function(B, n, p, h)
            % generates impulse response function for a given matrix B of VAR coefficients
            % using equations (4.13.2)-(4.13.4)
            %
            % parameters:
            % B : matrix of size (m+n*p,n)
            %     matrix of VAR coefficients
            % n : int
            %     number of endogenous variables
            % p : int
            %     number of lags
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
            for i=2:h
                Xh = [Yh Xh(:,1:end-n)];
                Yh = Xh * B;
                irf(:,:,i) = Yh';
            end
        end


        function [irf] = exogenous_impulse_response_function(B, n, m, r, p, h)

            % [irf] = exogenous_impulse_response_function(B, n, m, r, p, h)
            % generates exogenous impulse response function for a given matrix B of VAR coefficients
            % using equations (4.13.2) with exogenous regressors
            %
            % parameters:
            % B : matrix of size (m+n*p,n)
            %     matrix of VAR coefficients
            % n : int
            %     number of endogenous variables
            % m : int
            %     number of exogenous variables 
            % r : int
            %     number of exogenous variables other than constant, trend and quadratic trends 
            % p : int
            %     number of lags
            % h : int
            %     number of irf periods (in addition to impact period)
            %
            % returns:
            % irf : matrix of size (n,r,h)
            %     matrix of impulse response functions

            C = B(end+1-n*p-r:end-n*p,:);
            B = B(end+1-n*p:end,:);
            Yh = C;
            irf = cat(3, Yh, zeros(r,n,h-1));
            Xh = zeros(r,n*p);
            for i=2:h
                Xh = [Yh Xh(:,1:end-n)];
                Yh = Xh * B;
                irf(:,:,i) = Yh;
            end
            irf = permute(irf,[2 1 3]);
        end        
        

        function [structural_irf] = structural_impulse_response_function(irf, H, n)

            % structural_impulse_response_function(irf, H)
            % generates structural impulse response function using equation (4.13.9)
            % uses a vectorized form of (4.13.9) on all IRF periods to gain efficiency
            % 
            % parameters:
            % irf : matrix of size (n,n,h)
            %     matrix of impulse response functions
            % H : matrix of size (n,n)
            %     structural identification matrix
            % 
            % returns:
            % structural_irf : matrix of size (n,n,h)
            %     matrix of structural impulse response functions

            temp = reshape(permute(irf,[1 3 2]),[],n) * H;
            structural_irf = permute(reshape(temp,n,[],n),[1 3 2]);
        end
       

        function [posterior_estimates] = posterior_estimates(X, credibility_level)

            % posterior_estimates(X, credibility_level)
            % median, lower bound and upper bounf of credibility interval
            % 
            % parameters:
            % X : matrix of size (n,m, iterations)
            %     matrix of MCMC draws
            % credibility_level : float between 0 and 1
            %     credibility level for credibility interval
            % 
            % returns:
            % posterior_estimates : ndarray of shape (n,m,3)
            %     matrix of posterior estimates
        
            posterior_estimates = zeros(size(X,1),size(X,2),3);
            posterior_estimates(:,:,1) = quantile(X,0.5,3);
            posterior_estimates(:,:,2) = quantile(X,(1-credibility_level)/2,3);
            posterior_estimates(:,:,3) = quantile(X,(1+credibility_level)/2,3); 
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
        
        
        function [Y_p] = forecast(B, chol_Sigma, h, Z_p, Y, n)

            % [Y_p] = forecast(B, chol_Sigma, h, Z_p, Y, n)
            % products simulated forecasts
            % 
            % parameters:
            % B : matrix of size (k,n)
            %     matrix of VAR coefficients
            % chol_Sigma : ndarray of shape (n,n)
            %     Cholesky factor of residual variance-covariance matrix Sigma
            % h : int
            %     number of forecast periods
            % Z_p : matrix of shape (h,m)
            %     matrix of exogenous regressors for forecasts
            % Y : matrix of shape (T,n)
            %     matrix of endogenous variables
            % n : int
            %     number of endogenous variables
            % 
            % returns:
            % Y_p : matrix of size (h,n)
            %     matrix of simulated forecast values
        
            E = (chol_Sigma * randn(n,h))';
            Y_p = zeros(h,n);
            for j=1:h
                % get lagged endogenous regressors
                X = la.vec(fliplr(Y'))';
                % add exogenous regressors, if any
                if ~isempty(Z_p)
                    X = [Z_p(j,:) X];
                end
                % recover residuals
                e = E(j,:);
                % generate forecasts
                y = X * B + e;
                % update Y and Y_p
                Y = [Y(2:end,:);y];
                Y_p(j,:) = y;
            end   
        end
        
        
        function [fevd] = forecast_error_variance_decomposition(structural_irf, Gamma, n, h)

            % [fevd] = forecast_error_variance_decomposition(structural_irf, Gamma, n, h)
            % products forecast error variance decomposition from structural IRFs
            % 
            % parameters:
            % structural_irf : matrix of size (n,n,h)
            %     matrix of structural impulse response functions
            % Gamma : empty array or matrix of size (1,n)
            %     structural shock variance (empty array if variance is 1)
            % n : int
            %     number of endogenous variables          
            % h : int
            %     number of forecast periods
            % 
            % returns:
            % fevd : matrix of size (n,n,h)
            %     matrix of forecast error variance decomposition

            cum_squared_irf = cumsum(structural_irf.^2, 3);
            if ~isempty(Gamma)
                reshaped_Gamma = repmat(Gamma, [n 1 h]);
                cum_squared_irf =  reshaped_Gamma .* cum_squared_irf;
            end
            total_variance = repmat(sum(cum_squared_irf,2), [1 n 1]);
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
        
        
        function [hd] = historical_decomposition(structural_irf, structural_shocks, n, T)

            % [hd] = historical_decomposition(structural_irf, structural_shocks, n, T)
            % products historical decomposition from structural shocks and IRFs, using algorithm 13.6
            % 
            % parameters:
            % structural_irf : matrix of size (n,n,T)
            %     matrix of structural impulse response functions
            % structural_shocks : matrix of size (T,n)
            %     matrix of structural shocks        
            % n : int
            %     number of endogenous variables          
            % T : int
            %     number of sample periods  
            %
            % returns:
            % hd : matrix of size (n,n,T)
            %     matrix of historical decomposition   

            reshaped_shocks = flip(permute(repmat(structural_shocks,[1 1 n]), [3 2 1]), 3);
            hd = zeros(n,n,T);
            for j=1:T
                hd(:,:,j) = sum(structural_irf(:,:,1:j) .* reshaped_shocks(:,:,end-j+1:end), 3);
            end
        end


        function [y_bar Q omega gamma_00] = conditional_forecast_regressors_1(conditions, h, Y, n, p)

            % conditional_forecast_regressors_1(conditions, h, Y, n, p)
            % first set of elements for conditional forecasts: iteration-invariant 
            % 
            % parameters:
            % conditions : matrix of size (nconditions,4)
            %     matrix of conditions (one row per condition: variable, period, mean, variance)
            % h : int
            %     number of forecast periods       
            % Y : matrix of size (p,n)
            %     matrix of initial conditions for exogenous   
            % n : int
            %     number of endogenous variables 
            % p : int
            %     number of lags        
            % 
            % returns:
            % y_bar : matrix of size (h,n)
            %     matrix of mean values for conditions
            % Q : matrix of size (n,n*p)
            %     selection matrix for conditional forecasts state-space representation        
            % omega : matrix of size (h,n)
            %     matrix of variance values for conditions
            % gamma_00 : matrix of size (n*p,1)
            %     initial conditions (mean) for the space vector gamma_hat

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
            Q = [eye(n) zeros(n,n*(p-1))];
            gamma_00 = la.vec(fliplr(Y'));
        end
        
        
        function [mu F K Upsilon_00] = conditional_forecast_regressors_2(B, Sigma, conditions, Z_p, n, m, p, h)

            % [mu F K Upsilon_00] = conditional_forecast_regressors_2(B, Sigma, conditions, Z_p, n, m, p, h)
            % second set of elements for conditional forecasts: iteration-specific 
            % 
            % parameters:    
            % B : matrix of size (m+n*p,n)
            %     matrix of VAR coefficients
            % Sigma : matrix of size (n,n)
            %     variance-covariance matrix of VAR residuals 
            % conditions : matrix of size (nconditions,4)
            %     matrix of conditions (one row per condition: variable, period, mean, variance)            
            % Z_p : matrix of size (h,m)
            %     matrix of exogenous regressors for forecasts
            % n : int
            %     number of endogenous variables 
            % m : int
            %     number of exogenous variables          
            % p : int
            %     number of lags 
            % h : int
            %     number of forecast periods 
            % 
            % returns:
            % mu : matrix of size (h,n*p)
            %     matrix of intercepts for state variables
            % F : matrix of size (n*p,n*p)
            %     companion form matrix
            % K : matrix of size (n*p,n*p,h)
            %     variance-covariance matrix for state errors
            % Upsilon_00 : matrix of size (n*p,)
            %     initial conditions (variance) for the space vector gamma_hat

            F = vu.make_companion_form(B, n, m, p);
            mu = zeros(h,n*p);
            mu(:,1:n) = Z_p * B(1:m,:);
            K = zeros(n*p,n*p,h);
            for i=1:h
                K(1:n,1:n,i) = Sigma;
                temp = conditions(conditions(:,2) == i,:);
                condition_variables = temp(:,1);
                for j=1:size(condition_variables,1)
                    variable = condition_variables(j);
                    K(variable,variable,i) = 100;
                end
            end
            Upsilon_00 = 1e-10 * eye(n*p);
            Upsilon_00(1:n,1:n) = Sigma;
        end
        

        function [R y_bar omega] = conditional_forecast_regressors_3(conditions, h, n)

            % conditional_forecast_regressors_3(conditions, h, n)
            % first set of elements for structural conditional forecasts: iteration-invariant 
            % 
            % parameters:
            % conditions : matrix of size (nconditions,4)
            %     matrix of conditions (one row per condition: variable, period, mean, variance)
            % h : int
            %     number of forecast periods 
            % n : int
            %     number of endogenous variables       
            % 
            % returns:
            % R : matrix of size (n_conditions,n*h)
            %     selection matrix for conditional forecasts
            % y_bar : matrix of size (n_conditions,1)
            %     vector of mean values for conditions
            % omega : matrix of size (n_conditions,1)
            %     vector of variance values for conditions

            k = size(conditions,1);
            R = zeros(k,n*h);
            y_bar = zeros(k,1);
            omega = zeros(k,1);
            for i=1:k
                variable = conditions(i,1);
                period = conditions(i,2);
                mean = conditions(i,3);
                variance = max(1e-10,conditions(i,4));
                y_bar(i) = mean;
                omega(i) = variance;
                R(i,n*(period-1)+variable) = 1;
            end
        end
     

        function [M] = conditional_forecast_regressors_4(structural_irf, n, h)

            % conditional_forecast_regressors_4(structural_irf, n, h)
            % second set of elements for structural conditional forecasts, as in (4.14.14)
            %
            % parameters:
            % structural_irf : matrix of size (n,n,h)
            %     matrix of structural impulse response functions
            % n : int
            %     number of endogenous variables         
            % h : int
            %     number of forecast periods 
            % 
            % returns:
            % M : matrix of size (n*h,n*h)
            %     matrix of stacked impulse response function 

            temp = zeros(n,n*h);
            M = zeros(n*h,n*h);
            temp(:,end-n+1:end) = structural_irf(:,:,1);
            M(1:n,1:n) = structural_irf(:,:,1);
            for i=1:h-1
                temp(:,end-(i+1)*n+1:end-i*n) = structural_irf(:,:,i+1);
                M(i*n+1:(i+1)*n,1:(i+1)*n) = temp(:,end-(i+1)*n+1:end);
            end
        end
        
        
        function [P non_generating] = conditional_forecast_regressors_5(shocks, h, n)

            % conditional_forecast_regressors_5(shocks, h, n)
            % first set of elements for shock-specific structural conditional forecasts: iteration-invariant 
            %
            % parameters:
            % shocks : matrix of size (n,1)
            %     vector of generating shocks: 1 is generating, 0 is non-generating
            % h : int
            %     number of forecast periods 
            % n : int
            %     number of endogenous variables       
            %
            % returns:
            % P : matrix of size (m,n*h)
            %     selection matrix for the m non-generating shocks
            % non_generating : matrix of size (n,1)
            %     vector of non-generating shocks: 1 is non-generating, 0 is generating

            non_generating = ones(n,1) - shocks;
            P = diag(repmat(non_generating,[h 1]));
            P = P(any(P,2),:);
        end


        function [Gamma_nd] = conditional_forecast_regressors_6(gamma, non_generating, h)

            % conditional_forecast_regressors_6(gamma, non_generating, h)
            % second set of elements for structural conditional forecasts, as in (4.14.21)
            % 
            % parameters:
            % gamma : ndarray of shape (n,1)
            %     vector of variance values for structural shocks
            % non_generating : ndarray of shape (n,1)
            %     vector of non-generating shocks: 1 is non-generating, 0 is generating        
            % h : int
            %     number of forecast periods      
            %
            % returns:
            % Gamma_nd : ndarray of shape (m,1)
            %     variance vector for the m non-generating shocks

            non_generating_variances = gamma .* non_generating;
            non_generating_variances = non_generating_variances(non_generating_variances ~= 0);
            Gamma_nd = repmat(non_generating_variances, [h 1]);
        end


        function [Y_p] = linear_forecast(B, h, Z_p, Y, n)

            % [Y_p] = linear_forecast(B, h, Z_p, Y, n)
            % best linear forecasts, absent shocks
            % 
            % parameters:
            % B : matrix of size (k,n)
            %     matrix of VAR coefficients
            % h : int
            %     number of forecast periods
            % Z_p : matrix of shape (h,m)
            %     matrix of exogenous regressors for forecasts
            % Y : matrix of shape (T,n)
            %     matrix of endogenous variables
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
                % generate forecasts
                y = X * B;
                % update Y and Y_p
                Y = [Y(2:end,:);y];
                Y_p(j,:) = y;
            end   
        end    
        
        
        function [mu_hat Omega_hat] = conditional_forecast_posterior(y_bar, f, M, R, gamma, omega, n, h)

            % [mu_hat Omega_hat] = conditional_forecast_posterior(y_bar, f, M, R, gamma, omega, n, h)
            % posterior parameters for the structural conditional forecasts, as in (4.14.20)
            %
            % parameters:
            % y_bar : matrix of size (n_conditions,1)
            %     vector of mean values for conditions
            % f : matrix of size (h,n)
            %     matrix of simulated forecast values   
            % M : matrix of size (n*h,n*h)
            %     matrix of stacked impulse response function         
            % R : matrix of size (n_conditions,n*h)
            %     selection matrix for conditional forecasts
            % gamma : matrix of size (n,1)
            %     vector of variance values for structural shocks
            % omega : matrix of size (n_conditions,1)
            %     vector of variance values for conditions
            % n : int
            %     number of endogenous variables        
            % h : int
            %     number of forecast periods
            %
            % returns:
            % mu_hat : matrix of size (n*h,1)
            %     vector of posterior mean values
            % Omega_hat : matrix of size (n*h,n*h)
            %     matrix of posterior variance-covariance values        

            f_T = la.vec(f');
            D = R * M;
            D_star = D' / (D * D');
            mu_hat = f_T + M * D_star * (y_bar - R * f_T);
            temp = eye(n*h) - D_star * D;
            Omega_hat = M * (D_star * (omega .* D_star') + temp * (repmat(gamma,[h 1]) .* temp)) * M';  
        end
        

        function [mu_hat Omega_hat] = shock_specific_conditional_forecast_posterior(y_bar, f, M, R, P, gamma, Gamma_nd, omega, n, h)

            % [mu_hat Omega_hat] = shock_specific_conditional_forecast_posterior(y_bar, f, M, R, P, gamma, Gamma_nd, omega, n, h)
            % posterior parameters for the structural conditional forecasts, as in (4.14.20) and (4.14.26)
            % 
            % parameters:
            % y_bar : matrix of size (n_conditions,1)
            %     vector of mean values for conditions
            % f : matrix of size (h,n)
            %     matrix of simulated forecast values   
            % M : matrix of size (n*h,n*h)
            %     matrix of stacked impulse response function         
            % R : matrix of size (n_conditions,n*h)
            %     selection matrix for conditional forecasts
            % P : matrix of size (m,n*h)
            %     selection matrix for the m non-generating shocks        
            % gamma : matrix of size (n,1)
            %     vector of variance values for structural shocks
            % Gamma_nd : ndarray of shape (m,1)
            %     variance vector for the m non-generating shocks        
            % omega : ndarray of shape (n_conditions,1)
            %     vector of variance values for conditions
            % n : int
            %     number of endogenous variables        
            % h : int
            %     number of forecast periods
            %
            % returns:
            % mu_hat : matrix of size (n*h,1)
            %     vector of posterior mean values
            % Omega_hat : matrix of size (n*h,n*h)
            %     matrix of posterior variance-covariance values        

            Q = P / M;
            Z = [R; Q];
            f_T = la.vec(f');
            g_T = [y_bar; Q * f_T];
            xi = [omega; Gamma_nd];
            D = Z * M;
            D_star = D' / (D * D');
            mu_hat = f_T + M * D_star * (g_T - Z * f_T);
            temp = eye(n*h) - D_star * D;
            Omega_hat = M * (D_star * (xi .* D_star') + temp * (repmat(gamma,[h 1]) .* temp)) * M';    
        end       
        
        
        function [restriction_matrices max_irf_period] = make_restriction_matrices(restriction_table, p)

            % [restriction_matrices max_period] = make_restriction_matrices(restriction_table, p)
            % creates restriction and coefficient matrices to check restriction validity in later algorithms
            % 
            % parameters:
            % restriction_table : matrix of size (n_restrictions, 3+n_endogenous)
            %     matrix of restrictions
            % p : int
            %     number of lags 
            % 
            % returns:
            % restriction_matrices : cell array of dimension (7,2)
            %     each cell entry stores matrices of restriction and coefficient values
            % max_irf_period : int
            %     maximum number of periods for which IRFs will have to be computed in later algorithms     
  
            restriction_matrices = cell(7,2);
            max_irf_period = 0;
            % zero restrictions
            zero_restrictions = restriction_table(restriction_table(:,1)==1,:);
            restriction_number = size(zero_restrictions,1);            
            if restriction_number ~= 0
                indices = zeros(restriction_number,3);
                for i=1:restriction_number
                    shock = find(zero_restrictions(i,4:end));
                    period = zero_restrictions(i,3);
                    indices(i,1) = zero_restrictions(i,2);
                    indices(i,2) = shock;
                    indices(i,3) = period;
                    max_irf_period = max(max_irf_period, period);
                end
                restriction_matrices{1,1} = indices; 
            end
            % restrictions on IRFs: sign restrictions
            restrictions = restriction_table(restriction_table(:,1)==2,:);  
            [rows columns] = find(restrictions(:,4:end));
            occurrences_1 = unique(rows);
            occurrences_2 = sum(rows(:)==unique(rows)');
            sign_occurrences = occurrences_1(occurrences_2==1);
            sign_restrictions = restrictions(sign_occurrences,:);
            restriction_number = size(sign_occurrences,1);
            if restriction_number ~= 0
                indices = zeros(restriction_number,3);
                coefficients = zeros(restriction_number,1);
                for i=1:restriction_number
                    shock = find(sign_restrictions(i,4:end));
                    period = sign_restrictions(i,3);
                    indices(i,1) = sign_restrictions(i,2);
                    indices(i,2) = shock;
                    indices(i,3) = period;
                    coefficients(i,1) = sign_restrictions(i,3+shock);
                    max_irf_period = max(max_irf_period, period);
                end
                restriction_matrices{2,1} = indices; 
                restriction_matrices{2,2} = coefficients; 
            end 
            % restrictions on IRFs: magnitude restrictions
            magnitude_occurrences = occurrences_1(occurrences_2==2);
            magnitude_restrictions = restrictions(magnitude_occurrences,:);
            restriction_number = size(magnitude_occurrences,1);
            if restriction_number ~= 0        
                indices = zeros(restriction_number,4);
                coefficients = zeros(restriction_number,2);
                for i=1:restriction_number
                    shocks = find(magnitude_restrictions(i,4:end));
                    period = magnitude_restrictions(i,3);
                    indices(i,1) = magnitude_restrictions(i,2);
                    indices(i,2:3) = shocks;
                    indices(i,4) = period;                
                    coefficients(i,:) = magnitude_restrictions(i,3+shocks);
                    max_irf_period = max(max_irf_period, period);
                end
                restriction_matrices{3,1} = indices; 
                restriction_matrices{3,2} = coefficients;
            end
            % restrictions on shocks: sign restrictions
            restrictions = restriction_table(restriction_table(:,1)==3,:);  
            [rows columns] = find(restrictions(:,4:end));
            occurrences_1 = unique(rows);
            occurrences_2 = sum(rows(:)==unique(rows)');
            sign_occurrences = occurrences_1(occurrences_2==1);
            sign_restrictions = restrictions(sign_occurrences,:);
            restriction_number = size(sign_occurrences,1);                
            if restriction_number ~= 0        
                indices = zeros(restriction_number,2);
                coefficients = zeros(restriction_number,1);
                for i=1:restriction_number
                    shock = find(sign_restrictions(i,4:end));
                    period = sign_restrictions(i,3) - p;           
                    indices(i,1) = period;
                    indices(i,2) = shock;
                    coefficients(i,1) = sign_restrictions(i,3+shock);
                    max_irf_period = max(max_irf_period, period);
                end
                restriction_matrices{4,1} = indices; 
                restriction_matrices{4,2} = coefficients;  
            end
            % restrictions on shocks: magnitude restrictions
            magnitude_occurrences = occurrences_1(occurrences_2==2);
            magnitude_restrictions = restrictions(magnitude_occurrences,:);            
            restriction_number = size(magnitude_occurrences,1); 
            if restriction_number ~= 0    
                indices = zeros(restriction_number,3);
                coefficients = zeros(restriction_number,2);            
                for i=1:restriction_number
                    shocks = find(magnitude_restrictions(i,4:end));
                    period = magnitude_restrictions(i,3) - p;
                    indices(i,1) = period;
                    indices(i,2:end) = shocks;              
                    coefficients(i,:) = magnitude_restrictions(i,3+shocks);
                    max_irf_period = max(max_irf_period, period);            
                end
                restriction_matrices{5,1} = indices; 
                restriction_matrices{5,2} = coefficients;  
            end            
            % restrictions on historical decomposition: sign restrictions
            restrictions = restriction_table(restriction_table(:,1)==4,:);  
            [rows columns] = find(restrictions(:,4:end));
            occurrences_1 = unique(rows);
            occurrences_2 = sum(rows(:)==unique(rows)');
            sign_occurrences = occurrences_1(occurrences_2==1);
            sign_restrictions = restrictions(sign_occurrences,:);
            restriction_number = size(sign_occurrences,1);                
            if restriction_number ~= 0        
                indices = zeros(restriction_number,3);
                coefficients = zeros(restriction_number,1);
                for i=1:restriction_number
                    shock = find(sign_restrictions(i,4:end));
                    period = sign_restrictions(i,3) - p;           
                    indices(i,1) = sign_restrictions(i,2);
                    indices(i,2) = shock;
                    indices(i,3) = period;
                    coefficients(i,1) = sign_restrictions(i,3+shock);
                    max_irf_period = max(max_irf_period, period);
                end
                restriction_matrices{6,1} = indices; 
                restriction_matrices{6,2} = coefficients;  
            end         
            % restrictions on historical decomposition: magnitude restrictions
            magnitude_occurrences = occurrences_1(occurrences_2==2);
            magnitude_restrictions = restrictions(magnitude_occurrences,:);            
            restriction_number = size(magnitude_occurrences,1); 
            if restriction_number ~= 0    
                indices = zeros(restriction_number,4);
                coefficients = zeros(restriction_number,2);            
                for i=1:restriction_number
                    shocks = find(magnitude_restrictions(i,4:end));
                    period = magnitude_restrictions(i,3) - p;
                    indices(i,1) = magnitude_restrictions(i,2);
                    indices(i,2:3) = shocks;
                    indices(i,4) = period;              
                    coefficients(i,:) = magnitude_restrictions(i,3+shocks);
                    max_irf_period = max(max_irf_period, period);            
                end
                restriction_matrices{7,1} = indices; 
                restriction_matrices{7,2} = coefficients;  
            end                  
        end
        
        
        function [restriction_matrices] = make_covariance_restriction_matrices(restriction_table)   

            % make_covariance_restriction_matrices(restriction_table)
            % creates restriction and coefficient matrices to check restriction validity in later algorithms
            % 
            % parameters:
            % restriction_table : matrix of size (n_restrictions, 3+n_endogenous)
            %     matrix of restrictions
            % 
            % returns:
            % restriction_matrices : cell of length 2
            %     each list entry stores matrices of restriction and coefficient values   

            restriction_matrices = cell(2,2);
            % restrictions on covariance: sign restrictions
            restrictions = restriction_table(restriction_table(:,1)==5,:);
            [rows columns] = find(restrictions(:,4:end));
            occurrences_1 = unique(rows);
            occurrences_2 = sum(rows(:)==unique(rows)');            
            sign_occurrences = occurrences_1(occurrences_2==1);
            sign_restrictions = restrictions(sign_occurrences,:);            
            restriction_number = size(sign_occurrences,1);
            if restriction_number ~= 0
                indices = zeros(restriction_number,2);
                coefficients = zeros(restriction_number,1);  
                for i=1:restriction_number
                    shock = find(sign_restrictions(i,4:end));
                    indices(i,1) = sign_restrictions(i,2);
                    indices(i,2) = shock;
                    coefficients(i,1) = sign_restrictions(i,3+shock);
                end
                restriction_matrices{1,1} = indices; 
                restriction_matrices{1,2} = coefficients;                
            end
            % restrictions on covariance: magnitude restrictions
            magnitude_occurrences = occurrences_1(occurrences_2==2);
            magnitude_restrictions = restrictions(magnitude_occurrences,:);            
            restriction_number = size(magnitude_occurrences,1);
            if restriction_number ~= 0
                indices = zeros(restriction_number,3);
                coefficients = zeros(restriction_number,2);
                for i=1:restriction_number
                    shocks = find(magnitude_restrictions(i,4:end));
                    indices(i,1) = magnitude_restrictions(i,2);
                    indices(i,2:3) = shocks;
                    coefficients(i,:) = magnitude_restrictions(i,3+shocks);
                end
                restriction_matrices{2,1} = indices; 
                restriction_matrices{2,2} = coefficients;
            end   
        end
        
        
        function [mcmc_irf] = make_restriction_irf(mcmc_beta, mcmc_chol_Sigma, iterations, n, p, max_irf_period)

            % [mcmc_irf] = make_restriction_irf(mcmc_beta, mcmc_chol_Sigma, iterations, n, p, max_irf_period)
            % creates preliminary orthogonalized IRFs for restriction algorithm
            % 
            % parameters:
            % mcmc_beta : matrix of size (k, n, iterations)
            %     matrix of mcmc values for beta
            % mcmc_chol_Sigma : matrix of size (n, n, iterations)
            %     matrix of mcmc values for h(Sigma)
            % iterations: int
            %     number of MCMC iterations
            % n : int
            %     number of endogenous variables       
            % p : int
            %     number of lags
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
                    irf = vu.impulse_response_function(mcmc_beta(:,:,i), n, p, max_irf_period);
                    structural_irf = vu.structural_impulse_response_function(irf, mcmc_chol_Sigma(:,:,i), n);            
                    mcmc_irf(:,:,:,i) = structural_irf;
                end
            end
        end
       
        
        function [mcmc_shocks] = make_restriction_shocks(mcmc_beta, mcmc_chol_Sigma, Y, X, T, n, iterations, restriction_matrices)

            % make_restriction_shocks(mcmc_beta, mcmc_chol_Sigma, Y, X, iterations, restriction_matrices)
            % creates preliminary structural shocks for restriction algorithm
            % 
            % parameters:
            % mcmc_beta : matrix of size (k, n, iterations)
            %     matrix of mcmc values for beta
            % mcmc_chol_Sigma : matrix of size (n, n, iterations)
            %     matrix of mcmc values for h(Sigma)
            % Y : matrix of size (T,n)
            %     matrix of endogenous variables
            % X : matrix of size (T,k)
            %     matrix of VAR regressors    
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

            if isempty(restriction_matrices{4,1}) && isempty(restriction_matrices{5,1}) ...
               && isempty(restriction_matrices{6,1}) && isempty(restriction_matrices{7,1}) 
                mcmc_shocks = [];
            else
                mcmc_shocks = zeros(T, n, iterations);
                for i=1:iterations
                    [E ~] = vu.fit_and_residuals(Y, X, mcmc_beta(:,:,i));
                    Xi = E / mcmc_chol_Sigma(:,:,i)';
                    mcmc_shocks(:,:,i) = Xi;
                end
            end
        end  
        

        function [restriction_satisfied] = check_irf_sign(indices, coefficients, irf, Q)

            % [restriction_satisfied] = check_irf_sign(indices, coefficients, irf, Q)
            % checks whether IRF sign restrictions are satisfied
            % 
            % parameters:
            % indices : matrix of size (n_restrictions,3)
            %     matrix of IRF indices (variable, shock, period) to check
            % coefficients : matrix of size (n_restrictions,1)
            %     matrix of IRF coefficients to apply for positivity/negativity restrictions
            % irf : matrix of size (n,n,n_periods)
            %     matrix of preliminary orthogonalized IRFs 
            % Q : matrix of size (n,n)
            %     uniform orthogonal matrix 
            % 
            % returns:
            % restriction_satisfied : bool
            %     true if all restrictions are satisfied, false otherwise

            for i=1:size(indices,1)
                irf_value = irf(indices(i,1),:,indices(i,3)) * Q(:,indices(i,2));
                restriction = irf_value * coefficients(i);
                if restriction < 0
                    restriction_satisfied = false;
                    return 
                end
            end
            restriction_satisfied = true; 
        end
        
        
        function [restriction_satisfied] = check_irf_magnitude(indices, coefficients, irf, Q)

            % [restriction_satisfied] = check_irf_magnitude(indices, coefficients, irf, Q)
            % checks whether IRF magnitude restrictions are satisfied
            % 
            % parameters:
            % indices : matrix of size (n_restrictions,4)
            %     matrix of IRF indices (variable, shock1, shock2, period) to check
            % coefficients : matrix of size (n_restrictions,2)
            %     matrix of IRF coefficients to apply for positivity/negativity restrictions
            % irf : matrix of size (n,n,n_periods)
            %     matrix of preliminary orthogonalized IRFs 
            % Q : matrix of size (n,n)
            %     uniform orthogonal matrix 
            % 
            % returns:
            % restriction_satisfied : bool
            %     true if all restrictions are satisfied, false otherwise

            for i=1:size(indices,1)
                irf_values = abs(irf(indices(i,1),:,indices(i,4)) * Q(:,indices(i,[2 3])));
                restriction = sum(irf_values .* coefficients(i,:));
                if restriction < 0
                    restriction_satisfied = false;
                    return
                end
            end
            restriction_satisfied = true;
        end

        
        function [restriction_satisfied] = check_shock_sign(indices, coefficients, shocks, Q)

            % [restriction_satisfied] = check_shock_sign(indices, coefficients, shocks, Q)
            % checks whether shock sign restrictions are satisfied
            % 
            % parameters:
            % indices : matrix of size (n_restrictions,2)
            %     matrix of shock indices (period, shock) to check
            % coefficients : matrix of size (n_restrictions,1)
            %     matrix of shock coefficients to apply for positivity/negativity restrictions
            % shocks : matrix of size (T,n)
            %     matrix of preliminary structural shocks
            % Q : matrix of size (n,n)
            %     uniform orthogonal matrix 
            % 
            % returns:
            % restriction_satisfied : bool
            %     true if all restrictions are satisfied, false otherwise
 
            for i=1:size(indices,1)
                shock_value = shocks(indices(i,1),:) * Q(:,indices(i,2));
                restriction = shock_value * coefficients(i);
                if restriction < 0
                    restriction_satisfied = false;
                    return
                end
            restriction_satisfied = true;
            end
        end
               
        
        function [restriction_satisfied] = check_shock_magnitude(indices, coefficients, shocks, Q)

            % [restriction_satisfied] = check_shock_magnitude(indices, coefficients, shocks, Q)
            % checks whether shock magnitude restrictions are satisfied
            % 
            % parameters:
            % indices : matrix of size (n_restrictions,3)
            %     matrix of shock indices (period,shock1, shock2) to check
            % coefficients : matrix of size (n_restrictions,2)
            %     matrix of shock coefficients to apply for positivity/negativity restrictions
            % irf : matrix of size (n,n,n_periods)
            %     matrix of preliminary orthogonalized IRFs 
            % Q : matrix of size (n,n)
            %     uniform orthogonal matrix 
            % 
            % returns:
            % restriction_satisfied : bool
            %     true if all restrictions are satisfied, false otherwise

            for i=1:size(indices,1)
                shock_values = abs(shocks(indices(i,1),:) * Q(:,indices(i,[2 3])));
                restriction = sum(shock_values .* coefficients(i,:));
                if restriction < 0
                    restriction_satisfied = false;
                    return
                end
            end
            restriction_satisfied = true;         
        end
        
        
        function [structural_irf structural_shocks] = make_restriction_irf_and_shocks(irf, shocks, Q, n)

            % make_restriction_irf_and_shocks(irf, shocks, Q, n)
            % generates structural IRFs and shocks for a given Q matrix
            % 
            % parameters:
            % irf : matrix of size (n,n,n_periods)
            %     matrix of preliminary orthogonalized IRFs 
            % shocks : matrix of size (T,n)
            %     matrix of preliminary structural shocks      
            % Q : matrix of size (n,n)
            %     uniform orthogonal matrix 
            % n : int
            %     number of endogenous variables 
            % 
            % returns:
            % structural_irf : matrix of size (n,n,n_periods)
            %     matrix of final orthogonalized IRFs 
            % structural_shocks : matrix of size (T,n)
            %     matrix of final orthogonalized IRFs         

            structural_irf = vu.structural_impulse_response_function(irf, Q, n); 
            structural_shocks = shocks * Q;
        end  
        
        
        function [restriction_satisfied] = check_history_sign(indices, coefficients, irf, shocks)

            % check_history_sign(indices, coefficients, irf, shocks)
            % checks whether historical sign restrictions are satisfied
            % 
            % parameters:
            % indices : matrix of size (n_restrictions,3)
            %     matrix of IRF indices (variable, shock, period) to check
            % coefficients : matrix of size (n_restrictions,1)
            %     matrix of historical decomposition coefficients to apply for positivity/negativity restrictions
            % irf : matrix of size (n,n,n_periods)
            %     matrix of orthogonalized IRFs 
            % shocks : matrix of size (T,n)
            %     matrix of structural shocks     
            % 
            % returns:
            % restriction_satisfied : bool
            %     true if all restrictions are satisfied, false otherwise

            for i=1:size(indices,1)
                restriction_irf = la.vec(irf(indices(i,1),indices(i,2),1:indices(i,3)));
                restriction_shocks = flipud(shocks(1:indices(i,3),indices(i,2)));
                hd_value = sum(restriction_irf .* restriction_shocks);
                restriction = hd_value * coefficients(i);
                if restriction < 0
                    restriction_satisfied = false;
                    return
                end
            end
            restriction_satisfied = true;         
        end
        
        
        function [restriction_satisfied] = check_history_magnitude(indices, coefficients, irf, shocks)

            % [restriction_satisfied] = check_history_magnitude(indices, coefficients, irf, shocks)
            % checks whether historical magnitude restrictions are satisfied
            % 
            % parameters:
            % indices : matrix of size (n_restrictions,4)
            %     matrix of IRF indices (variable, shock1, shock2, period) to check
            % coefficients : matrix of size (n_restrictions,2)
            %     matrix of historical decomposition coefficients to apply for positivity/negativity restrictions
            % irf : matrix of size (n,n,n_periods)
            %     matrix of orthogonalized IRFs 
            % shocks : matrix of size (T,n)
            %     matrix of structural shocks     
            % 
            % returns:
            % restriction_satisfied : bool
            %     true if all restrictions are satisfied, false otherwise

            for i=1:size(indices,1)
                restriction_irf_1 = la.vec(irf(indices(i,1),indices(i,2),1:indices(i,4)));
                restriction_shocks_1 = flipud(shocks(1:indices(i,4),indices(i,2)));
                hd_value_1 = sum(restriction_irf_1 .* restriction_shocks_1);
                restriction_1 = abs(hd_value_1) * coefficients(i,1);
                restriction_irf_2 = la.vec(irf(indices(i,1),indices(i,3),1:indices(i,4)));
                restriction_shocks_2 = flipud(shocks(1:indices(i,4),indices(i,3)));
                hd_value_2 = sum(restriction_irf_2 .* restriction_shocks_2);
                restriction_2 = abs(hd_value_2) * coefficients(i,2);
                restriction = restriction_1 + restriction_2;
                if restriction < 0
                    restriction_satisfied = false;
                    return
                end
            end
            restriction_satisfied = true;   
        end
        
        
        function [restriction_satisfied]  = check_covariance_sign(indices, coefficients, V, n, h)

            % check_covariance_sign(indices, coefficients, V, n, h)
            % checks whether covariance sign restrictions are satisfied
            % 
            % parameters:
            % indices : matrix of size (n_restrictions,2)
            %     matrix of covariance indices (period, shock) to check
            % coefficients : matrix of size (n_restrictions,1)
            %     matrix of covariance coefficients to apply for positivity/negativity restrictions
            % V : matrix of size (h,h)
            %     matrix of covariance between proxys and structural shocks
            % n : int
            %     number of endogenous variables         
            % h : int
            %     number of proxy variables         
            % 
            % returns:
            % restriction_satisfied : bool
            %     true if all restrictions are satisfied, false otherwise
            
            E_r_xi = [zeros(h,n-h) V];
            for i=1:size(indices,1)
                covariance_value = E_r_xi(indices(i,1),indices(i,2));
                restriction = (covariance_value * coefficients(i));            
                if restriction < 0
                    restriction_satisfied = false;
                    return
                end
            end
            restriction_satisfied = true;   
        end            
                
                
        function [restriction_satisfied] = check_covariance_magnitude(indices, coefficients, V, n, h)

            % check_covariance_magnitude(indices, coefficients, V, n, h)
            % checks whether covariance magnitude restrictions are satisfied
            % 
            % parameters:
            % indices : matrix of size (n_restrictions,3)
            %     matrix of covariance indices (variable, shock1, shock2) to check
            % coefficients : matrix of size (n_restrictions,2)
            %     matrix of covariance coefficients to apply for magnitude restrictions
            % V : matrix of size (h,h)
            %     matrix of covariance between proxys and structural shocks
            % n : int
            %     number of endogenous variables         
            % h : int
            %     number of proxy variables  
            % 
            % returns:
            % restriction_satisfied : bool
            %     true if all restrictions are satisfied, false otherwise

            E_r_xi = [zeros(h,n-h) V];
            for i=1:size(indices,1)
                covariance_values = abs(E_r_xi(indices(i,1),indices(i,[2 3])));
                restriction = sum(covariance_values .* coefficients(i,:));      
                if restriction < 0
                    restriction_satisfied = false;
                    return
                end
            end
            restriction_satisfied = true;   
        end   


        function [mcmc_B] = ols_var_mcmc_beta(B, Sigma, XX, k, n, q)

            % ols_var_mcmc_beta(B, Sigma, XX, k, n)
            % generates pseudo MCMC draws for beta for an OLS VAR model
            % 
            % parameters:
            % B : matrix of size (k,n)
            %     matrix of VAR coefficients
            % Sigma : matrix of size (n,n)
            %     variance-covariance matrix of VAR residuals
            % XX : matrix of size (k,k)
            %     covariance matrix of regressors X
            % k : int
            %     number of VAR coefficients by equation
            % n : int
            %     number of endogenous variables   
            % q : int
            %     total number of VAR coefficients
            % 
            % returns:
            % mcmc_B : matrix of size (k,n,500)
            %     matrix of pseudo MCMC draws

            Q = kron(Sigma, la.invert_spd_matrix(XX));
            mcmc_beta = la.vec(B) + la.cholesky_nspd(Q) * randn(q, 500);
            mcmc_B = reshape(mcmc_beta,[k n 500]);
        end


        function [secondary_constrained_coefficients_table] = ...
                 rework_constraint_table(constrained_coefficients_table, lags)

            % rework_constraint_table(constrained_coefficients_table, lags)
            % update constraint table by switching '-1' lags to 'all lags' entries
            % 
            % parameters:
            % constrained_coefficients_table : matrix of size (n_constraint,5)
            %     matrix of constraints
            % lags : int
            %     number of lags, either -1 or positive
            % 
            % returns:
            % secondary_constrained_coefficients_table : matrix of size (n_constraint,5)
            %     matrix of updated constraints
            
            secondary_constrained_coefficients_table =  [];
            for i=1:(size(constrained_coefficients_table,1))
                row = constrained_coefficients_table(i,:);
                lag = constrained_coefficients_table(i,3);
                if lag == -1
                    temp = repmat(row,[lags 1]);
                    temp(:,3) = (1:lags);
                    secondary_constrained_coefficients_table = [...
                    secondary_constrained_coefficients_table; temp];
                else
                    secondary_constrained_coefficients_table = [...
                    secondary_constrained_coefficients_table; row];
                end
            end
        end

        
    end
    
    
end
        