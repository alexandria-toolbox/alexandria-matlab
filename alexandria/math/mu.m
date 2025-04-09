classdef mu
    

    % mu stands for mathematics utilities
    % A class containing static methods for mathematics utilities (excluding linear algebra)

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
       
        
        function [log_sum_x] = log_sum_exp(log_x)
            
            % [log_sum_x] = log_sum_exp(log_x)
            % log-sum-exp function that computes log(sum(x)) as a function of sum(log(x))
            % using equation (3.10.30)
            %
            % parameters:
            % log_x : matrix of size (n,1)
            %     vector containing log(x) values
            %
            % returns:
            % log_sum_x : float
            %     log of the sum of the x's
            
            log_x_bar = max(log_x);
            diff = log_x - log_x_bar;
            log_sum_x = log_x_bar + log(sum(exp(diff)));
        end
        
        
        function [y] = log_multivariate_gamma(x, n)
            
            % log_multivariate_gamma(x, n)
            % log of multivariate Gamme function
            %
            % parameters:
            % x : float
            %     argument of Gamma function
            % n : int
            %     dimension of function
            %
            % returns:
            % y : float
            %     log value returned by the multivariate gamma function
            
            y = (n * (n-1) / 4) * log(pi);
            for j=1:n
                y = y + gammaln(x + (1 - j) / 2);
            end
        end
     
    end
    
end
        
