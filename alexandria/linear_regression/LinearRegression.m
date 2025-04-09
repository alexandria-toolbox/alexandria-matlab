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

    properties (GetAccess = {?BayesianRegression}, SetAccess = private)
        n_exogenous
    end
    
    properties (GetAccess = protected, SetAccess = private)
        XX
        Xy
        yy
    end 

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)


        function self = LinearRegression()
        end 


    end


    methods (Access = protected, Hidden = true)    
    

        function make_regressors(self)
            
            % generates regressors X and y defined in (3.9.2)

            % define y
            y = self.endogenous;
            % define X, adding constant and trends if included
            X = ru.add_intercept_and_trends(self.exogenous, self.constant, self.trend, self.quadratic_trend, 0);
            % get dimensions
            n_exogenous = size(self.exogenous, 2);
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


        function [beta_hat sigma_hat] = ols_regression(self)
            
            % maximum likelihood estimates for beta and sigma, from (3.9.7)
            
            [beta_hat sigma_hat]  = ru.ols_regression(self.y, self.X, self.XX, self.Xy, self.n);
        end 


    end
    
    
end
