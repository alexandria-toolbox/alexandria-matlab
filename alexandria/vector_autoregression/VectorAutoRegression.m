classdef VectorAutoRegression < handle
    
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
        Y
        Z
        X
        n
        m
        p
        T
        k
        q
    end    
    
    properties (GetAccess = protected, SetAccess = private)
        XX
        XY
        YY
    end 

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = protected, Hidden = true)    
    
    
        function self = VectorAutoregression()
        end    
    
    
        function make_regressors(self)
            
            % generates regressors Y, X as defined in (4.11.3), along with other dimension elements

            % define y            
            Y = self.make_endogenous_matrix();
            % define X
            [Z, X] = self.make_regressor_matrix();
            % define dimensions
            [n, m, p, T, k, q] = self.generate_dimensions();
            % define estimation terms
            XX = X' * X;
            XY = X' * Y;
            YY = Y' * Y;          
            % save as attributes
            self.Y = Y;
            self.Z = Z;
            self.X = X;
            self.n = n;
            self.m = m;
            self.p = p;
            self.T = T;
            self.k = k;
            self.q = q;
            self.XX = XX;
            self.XY = XY;
            self.YY = YY;            
        end      

        
        function [Y] = make_endogenous_matrix(self)
        
            % unpack, recover endogenous after trimming inital conditions
            endogenous = self.endogenous;
            lags = self.lags;
            Y = endogenous(lags+1:end,:);
        end   
        
        
        function [Z, X] = make_regressor_matrix(self)

            % unpack
            endogenous = self.endogenous;
            exogenous = self.exogenous;
            lags = self.lags;
            constant = self.constant;
            trend = self.trend;
            quadratic_trend = self.quadratic_trend;
            periods = size(endogenous,1) - lags;
            % get automated exogenous: constant, trend, quadratic trend
            X_1 = vu.generate_intercept_and_trends(constant, trend, quadratic_trend, periods, 0);
            % recover other exogenous
            X_2 = vu.generate_exogenous_regressors(exogenous, lags);
            % get lagged endogenous
            X_3 = vu.generate_lagged_endogenous(endogenous, lags);
            % concat to obtain final regressor matrix
            Z = [X_1 X_2];
            X = [X_1 X_2 X_3];
        end      

        
        function [n, m, p, T, k, q] = generate_dimensions(self)
        
            endogenous = self.endogenous;
            exogenous = self.exogenous;
            lags = self.lags;
            constant = self.constant;
            trend = self.trend;
            quadratic_trend = self.quadratic_trend;
            T = size(endogenous,1) - lags;
            n = size(endogenous,2);
            p = lags;
            m = double(constant) + double(trend) + double(quadratic_trend);    
            if ~isempty(exogenous)
                m = m + size(exogenous,2);
            end
            k = m + n * p;
            q = n * k;
        end
            

        function [B_hat, Sigma_hat] = ols_var(self, Y, X, XX, XY, T)
            
            % maximum likelihood estimates for B and Sigma, from (4.11.9)
            
            B_hat = XX \ XY;
            E_hat = Y - X * B_hat;
            Sigma_hat = E_hat' * E_hat / T;
        end     

        
    end
    
end
