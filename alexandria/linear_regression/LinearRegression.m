classdef LinearRegression < handle
    
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
        y
        X          
        n
        k
    end    
    
    properties (GetAccess = protected, SetAccess = private)
        n_exogenous
        XX
        Xy
        yy
    end 

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = protected, Hidden = true)    
    
    
        function self = LinearRegression()
        end    
    
    
        function make_regressors(self)
            
            % generates regressors X and y defined in (3.9.2)

            % unpack
            endogenous = self.endogenous;
            exogenous = self.exogenous;
            % define y
            y = endogenous;
            % define X, adding constant and trends if included
            X = self.add_intercept_and_trends(exogenous, true);
            % get dimensions
            n_exogenous = size(exogenous, 2);
            n = size(X,1);
            k = size(X,2);
            % define terms for posterior distribution
            XX = X' * X;
            Xy = X' * y;
            yy = y' * y;
            % save as attributes
            self.y = y;
            self.X = X;
            self.n_exogenous = n_exogenous;
            self.n = n;
            self.k = k;
            self.XX = XX;
            self.Xy = Xy;
            self.yy = yy;
        end      
    
    
        function [X] = add_intercept_and_trends(self, X, in_sample)
            
            % add constant, trend and quadratic trend to regressors if selected
            
            % unpack
            constant = self.constant;
            trend = self.trend;
            quadratic_trend = self.quadratic_trend;
            n_rows = size(X,1);
            % consider quadratic trend
            if quadratic_trend
                % if in_sample, quadratic trend starts at 1
                if in_sample
                    quadratic_trend_column = (1:n_rows).^2';
                % if out of sample, quadratic trend starts at (n+1)^2
                else
                    n = self.n;
                    quadratic_trend_column = (n + (1:n_rows)).^2';
                end
                X = [quadratic_trend_column X];
            end            
            % consider trend
            if trend
                % if in_sample, trend starts at 1
                if in_sample
                    trend_column = (1:n_rows)';
                % if out of sample, trend starts at n+1
                else
                    n = self.n;
                    trend_column = n + (1:n_rows)';
                end
                X = [trend_column X];
            end
            % consider intercept
            if constant
                constant_column = ones(n_rows,1);
                X = [constant_column X];
            end
        end      
    
    
        function [beta_hat, sigma_hat] = ols_regression(self, y, X, XX, Xy, n)
            
            % maximum likelihood estimates for beta and sigma, from (3.9.7)
            
            beta_hat = XX \ Xy;
            res = y - X * beta_hat;
            sigma_hat = res' * res / n;
        end     
    

    end
    
    
end
