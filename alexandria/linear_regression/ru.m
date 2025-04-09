classdef ru
    

    % ru stands for var utilities
    % A class containing static methods for linear regression utilities

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)


        function [X] = add_intercept_and_trends(exogenous, constant, trend, quadratic_trend, shift)

            % add_intercept_and_trends(X, constant, trend, quadratic_trend, shift)
            % add intercept and trends to regressor matrix X
            % 
            % parameters:
            % exogenous : ndarray of shape (periods,n_exogenous)
            %     matrix of regressors  
            % constant : bool
            %     if True, a constant is added to the VAR model
            % trend : bool
            %     if True, a linear trend is added to the VAR model
            % quadratic_trend : bool
            %     if True, a quadratic trend is added to the VAR model
            % shift : int
            %     sample period shift: 0 for in-sample data, sample periods otherwise
            %     
            % returns:
            % X : ndarray of shape (periods,n_exogenous_total)
            %     matrix of regressors    
        
            % sample periods
            X = exogenous;
            periods = size(X,1);
            % quadratic trend
            if quadratic_trend
                quadratic_trend_column = (shift + (1:periods)).^2';
                X = [quadratic_trend_column X];
            end
            % trend
            if trend
                trend_column = shift + (1:periods)';
                X = [trend_column X];
            end
            % intercept
            if constant
                constant_column = ones(periods,1);
                X = [constant_column X];
            end
        end


        function [b] = generate_b(b_exogenous, b_constant, b_trend, b_quadratic_trend, constant, trend, quadratic_trend, n_exogenous)
        
            % generate_b(b_exogenous, b_constant, b_trend, b_quadratic_trend, constant, trend, quadratic_trend, n_exogenous)
            % generates prior mean vector b
            % 
            % parameters:
            % b_exogenous : float or matrix of size (n_exogenous,1)
            %     prior mean for exogenous variables
            % b_constant : float
            %     prior mean for constant       
            % b_trend : float
            %     prior mean for linear trend     
            % b_quadratic_trend : float
            %     prior mean for quadratic trend
            % constant : bool
            %     if True, a constant is added to the VAR model
            % trend : bool
            %     if True, a linear trend is added to the VAR model
            % quadratic_trend : bool
            %     if True, a quadratic trend is added to the VAR model
            % n_exogenous : int
            %     number of exogenous regressors, excluding constant, trend and quadratic trend
            %     
            % returns:
            % b : matrix of size (n_exogenous_total,1)
            %     prior mean for exogenous variables, comprising exogenous regressors, constant, trend and quadratic trend

            % if b_exogenous is a scalar, turn it into a vector replicating the value
            if isscalar(b_exogenous)
                b_exogenous = b_exogenous * ones(n_exogenous, 1);
            end
            b = b_exogenous;
            % if quadratic trend is included, add to prior mean
            if quadratic_trend
                b = [b_quadratic_trend; b];
            end
            % if trend is included, add to prior mean
            if trend
                b = [b_trend; b];
            end
            % if constant is included, add to prior mean
            if constant
                b = [b_constant; b];
            end
        end 


        function [V] = generate_V(V_exogenous, V_constant, V_trend, V_quadratic_trend, constant, trend, quadratic_trend, n_exogenous)
        
            % generate_V(V_exogenous, V_constant, V_trend, V_quadratic_trend, constant, trend, quadratic_trend, n_exogenous)
            % generates prior mean vector b
            % 
            % parameters:
            % b_exogenous : float or matrix of size (n_exogenous,1)
            %     prior mean for exogenous variables
            % b_constant : float
            %     prior mean for constant       
            % b_trend : float
            %     prior mean for linear trend     
            % b_quadratic_trend : float
            %     prior mean for quadratic trend
            % constant : bool
            %     if True, a constant is added to the VAR model
            % trend : bool
            %     if True, a linear trend is added to the VAR model
            % quadratic_trend : bool
            %     if True, a quadratic trend is added to the VAR model
            % n_exogenous : int
            %     number of exogenous regressors, excluding constant, trend and quadratic trend
            %     
            % returns:
            % b : matrix of size (n_exogenous_total,1)
            %     prior mean for exogenous variables, comprising exogenous regressors, constant, trend and quadratic trend
                 
            % if V_exogenous is a scalar, turn it into a vector replicating the value
            if isscalar(V_exogenous)
                V_exogenous = V_exogenous * ones(n_exogenous, 1);
            end
            V = V_exogenous;
            % if quadratic trend is included, add to prior mean
            if quadratic_trend
                V = [V_quadratic_trend; V];
            end
            % if trend is included, add to prior mean
            if trend
                V = [V_trend; V];
            end
            % if constant is included, add to prior mean
            if constant
                V = [V_constant; V];
            end
        end


        function [V inv_V inv_V_b] = generate_b_and_V_arrays(b, V)
            
            % generate_b_and_V_arrays(b, V)
            % generates arrays related to b and V
            % 
            % parameters:
            % b : matrix of size (k,1)
            %     prior mean for exogenous variables
            % V : matrix of size (k,1)
            %     prior variance for exogenous variables
            %     
            % returns:
            % V : matrix of size (k,k)
            %     prior variance matrix for exogenous variables  
            % inv_V : matrix of size (k,k)
            %     inverse of prior variance matrix for exogenous variables
            % inv_V_b : matrix of size (k,)
            %     product inv(V) * b 

            % convert the vector V into an array
            inv_V_b = b ./ V;
            inv_V = diag(1 ./ V);
            V = diag(V);
        end

        
        function [beta_hat sigma_hat] = ols_regression(y, X, XX, Xy, n)
            
            % ols_regression(y, X, XX, Xy, n)
            % maximum likelihood estimates for beta and sigma, from (3.9.7)
            % 
            % parameters:
            % y : matrix of size (n,)
            %     endogenous (explained) variable
            % X : matrix of size (n,k)
            %     exogenous (explanatory) variables    
            % XX : matrix of size (k,k)
            %     regressors variance matrix     
            % Xy : matrix of size (k,)
            %     regressors covariance matrix  
            % n : int
            %     number of sample observations
            %     
            % returns:
            % beta_hat : matrix of size (k,)
            %     sample estimate for beta
            % sigma_hat : float
            %     sample estimate for sigma

            beta_hat = XX \ Xy;
            res = y - X * beta_hat;
            sigma_hat = res' * res / n;        
        end


        function [fitted residual] = fitted_and_residual(y, X, beta)
            
            % fitted_and_residual(y, X, beta)
            % generates in-sample fitted and residuals
            % 
            % parameters:
            % y : matrix of size (n,1)
            %     endogenous (explained) variable
            % X : matrix of size (n,k)
            %     exogenous (explanatory) variables
            % beta : matrix of size (k,1)
            %     regression coefficients
            %     
            % returns:
            % fitted : matrix of size (n,1)
            %     in-sample fitted values
            % residual : matrix of size (n,1)
            %     in-sample residual
            
            fitted = X * beta;
            residual = y - X * beta;
        end
        
        
        function [insample_evaluation] = insample_evaluation_criteria(y, res, n, k)
            
            % insample_evaluation_criteria(y, res, n, k)
            % generates in-sample evaluation criteria
            % 
            % parameters:
            % y : ndarray of shape (n,)
            %     endogenous (explained) variable
            % res : ndarray of shape (n,)
            %     in-sample residual
            % n : int
            %     number of sample observations
            % k : int
            %     dimension of VAR coefficients        
            %     
            % returns:
            % insample_evaluation : dict
            %     in-sample evaluation criteria

            ssr = res' * res;
            tss = (y - mean(y))' * (y - mean(y));
            r2 = 1 - ssr / tss;
            adj_r2 = 1 - (1 - r2) * (n - 1) / (n - k);
            insample_evaluation = struct;            
            insample_evaluation.ssr = ssr;
            insample_evaluation.r2 = r2;
            insample_evaluation.adj_r2 = adj_r2;
        end


        function [insample_evaluation] = ml_insample_evaluation_criteria(y, res, n, k, sigma)
            
            % ml_insample_evaluation_criteria(y, res, n, k, sigma)
            % generates in-sample evaluation criteria for maximum lilekihood regression
            % 
            % parameters:
            % y : ndarray of shape (n,)
            %     endogenous (explained) variable
            % res : ndarray of shape (n,)
            %     in-sample residual
            % n : int
            %     number of sample observations
            % k : int
            %     dimension of VAR coefficients
            % sigma : float
            %     residual variance
            %     
            % returns:
            % insample_evaluation : dict
            %     in-sample evaluation criteria

            ssr = res' * res;
            tss = (y - mean(y))' * (y - mean(y));
            r2 = 1 - ssr / tss;
            adj_r2 = 1 - (1 - r2) * (n - 1) / (n - k);
            aic = 2 * k / n + log(sigma);
            bic = k * log(n) / n + log(sigma);
            insample_evaluation = struct;            
            insample_evaluation.ssr = ssr;
            insample_evaluation.r2 = r2;
            insample_evaluation.adj_r2 = adj_r2;
            insample_evaluation.aic = aic;
            insample_evaluation.bic = bic;
        end

        
        function [forecast_evaluation_criteria] = forecast_evaluation_criteria(y_hat, y)
        
            % forecast_evaluation_criterion(y_hat, y)
            % forecast evaluation criteria from equations (3.10.11) and (3.10.12)
            % 
            % parameters:
            % y_hat : matrix of size (m,1)
            %     array of forecast values for forecast evaluation            
            % y : matrix of size (m,1)
            %     array of realised values for forecast evaluation
            %     
            % returns:
            % forecast_evaluation_criteria : struct
            %     forecast evaluation criteria
            
            err = y - y_hat;
            m = size(y_hat,1);
            rmse = sqrt(err' * err / m);
            mae = sum(abs(err)) / m;
            mape = 100 * sum(abs(err ./ y)) / m;
            theil_u = sqrt(err' * err) / (sqrt(y' * y) + sqrt(y_hat' * y_hat));
            bias = sum(err) / sum(abs(err));
            forecast_evaluation_criteria = struct;  
            forecast_evaluation_criteria.rmse = rmse;
            forecast_evaluation_criteria.mae = mae;
            forecast_evaluation_criteria.mape = mape;
            forecast_evaluation_criteria.theil_u = theil_u;
            forecast_evaluation_criteria.bias = bias;
        end

        
    end
    
    
end   