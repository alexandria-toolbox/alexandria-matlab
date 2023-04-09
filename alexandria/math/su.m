classdef su
        

    % su stands for statistics utilities
    % A class containing static methods for statistics utilities

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
        
        
        function [log_val, val] = normal_pdf(x, mu, sigma)
            
            % [log_val, val] = normal_pdf(x, mu, sigma)
            % log-pdf and pdf for the normal distribution
            %
            % parameters:
            % x : float
            %     x value at which the pdf must be calculated
            % mu : float
            %     mean parameter
            % sigma : float
            %     variance parameter (positive)
            %
            % returns:
            % log_val : float
            %     log-density of normal distribution at x
            % val : float
            %     density of normal distribution at x
            
            log_val = - 0.5 * (log(2 * pi) + log(sigma) + (x - mu)^2 / sigma);
            val = exp(log_val);
        end
        
        
        function [val] = normal_cdf(x, mu, sigma)
            
            % [val] = normal_cdf(x, mu, sigma)
            % cumulative distribution function for the normal distribution
            %
            % parameters:
            % x : float
            %     x value at which the cdf must be calculated
            % mu : float
            %     mean parameter
            % sigma : float
            %     variance parameter (positive)
            %
            % returns:
            % val : float
            %     cdf of normal distribution at x
            
            val = 0.5 * (1 + erf((x - mu) / sqrt(2 * sigma)));
        end
        

        function [val] = normal_icdf(p)
    
            % [val] = normal_icdf(p)
            % inverse cdf for the normal distribution, open source equivalent of norminv
            %
            % parameters:
            % p : float
            %     probability of the cdf (between 0 and 1)
            %
            % returns:
            % val : float
            %     inverse cdf value
            %
            % references:
            %     Abramowitz and Stegun (1964), "Handbook of Mathematical Functions"
            %     sections 7.1 and 26.2
            
            val = - sqrt(2) * erfcinv(2 * p);
        end
        
        
        function [log_val, val] = student_pdf(x, mu, sigma, v)
            
            % [log_val, val] = student_pdf(x, mu, sigma, nu)
            % log-pdf and pdf for the student distribution
            %
            % parameters:
            % x : float
            %     x value at which the pdf must be calculated
            % mu : float
            %     location parameter
            % sigma : float
            %     scale parameter (positive)
            % v : float
            %     degrees of freedom (positive)
            %
            % returns:
            % log_val : float
            %     log-density of student distribution at x
            % val : float
            %     density of student distribution at x
            
            term_1 = gammaln((v + 1) / 2) - gammaln(v / 2);
            term_2 = - 0.5 * (log(v) + log(pi) + log(sigma));
            term_3 = - 0.5 * (v + 1) * log(1 + (x - mu)^2 / (v * sigma));
            log_val = term_1 +  term_2 + term_3;
            val = exp(log_val);
        end
        
        
        function [val] = student_cdf(x, mu, sigma, v)
            
            % [val] = student_cdf(x, mu, sigma, v)
            % cumulative distribution function for the student distribution
            %
            % parameters:
            % x : float
            %     x value at which the cdf must be calculated
            % mu : float
            %     location parameter
            % sigma : float
            %     scale parameter (positive)
            % v : float
            %     degrees of freedom (positive)
            %
            % returns:
            % val : float
            %     cdf value
            %
            % references:
            %     Abramowitz and Stegun (1964), "Handbook of Mathematical Functions"
            %     formula 26.5.2 and 26.7.1
            
            % standardize
            x = (x - mu) / sqrt(sigma);
            x_2 = x^2;
            % for small df, use lower tail of incomplete Beta function
            if v <= x_2
                val = betainc(v / (v + x_2), v / 2, 0.5, 'lower') / 2;
            % for large df, revert and use upper tail
            else
                val = betainc(x_2 / (v + x_2), 0.5, v / 2, 'upper') / 2;
            end
            % if x is positive, cdf is 1 - val
            if x > 0
                val = 1 - val;
            end
        end
            
        
        function [val] = student_icdf(p, v)
    
            % [val] = student_icdf(p, v)
            % inverse cdf for the Student distribution, open source equivalent of tinv
            %
            % parameters:
            % p : float
            %     probability of the cdf (between 0 and 1)
            % v : float
            %     degrees of freedom (positive)
            %
            % returns:
            % val : float
            %     inverse cdf value
            %
            % references:
            % for small integer df:
            %     William Shaw: "New methods for simulating the Student T-distribution
            %     - Direct use of the inverse cumulative distribution function"
            % for inverse cdf with inverse regularized beta function:
            %     https://www.wolframalpha.com/input/?i=inverse+cdf+student%27s+t+distribution
            % for large degrees of freedom:
            %     Abramowitz and Stegun (1964), "Handbook of Mathematical Functions"
            %     formula 26.7.5

            % if p = 0, value is - inf
            if p == 0
                val = -Inf;
            % if p = 1, value is inf
            elseif p == 1
                val = Inf;
            % if p = 0.5, value is 0
            elseif p == 0.5
                val = 0;
            % otherwise, no trivial value, 
            % however, some low df have simple formulas
            % for df = 1
            elseif v == 1
                val = tan(pi * (p - 0.5));
            % for df = 2
            elseif v == 2
                val = (2 * p - 1) / sqrt(2 * p * (1 - p));
            % for df = 4
            elseif v == 4
                alph = 4 * p * (1 - p);
                q = cos(acos(sqrt(alph)) / 3) / sqrt(alph);
                val = sign(p - 0.5) * 2 * sqrt(q - 1);
            % otherwise, for small df, use the general inverse cdf formula 
            % relying on the inverse regularized beta function
            elseif v < 1000
                if p < 0.5
                    val = - sqrt(v) * sqrt(1 / betaincinv(2 * p, v / 2, 0.5) - 1);
                else
                    val = sqrt(v) * sqrt(1 / betaincinv(2 * (1 - p), v / 2, 0.5) - 1);
                end
            % otherwise, go for the large df approximation:
            else
            x = norminv(p);
            val = x + (x + x^3) / (4 * v) + ...
                   (5 * x^5 + 16 * x^3 + 3 * x) / (96 * v^2) + ...
                   (3 * x^7 + 19 * x^5 + 17 * x^3 - 15 * x) / (384 * v^3) +...
                   (79 * x^9 + 776 * x^7 + 1482 * x^5 - 1920 * x^3 - 945 * x) / (92160 * v^4);
            end
        end                    
            
            
        function [val] = gamma_icdf(p, a, b)
    
            % [val] = gamma_icdf(p, a)
            % inverse cdf for the standard gamma distribution, open source equivalent of gaminv
            %
            % parameters:
            % p : float
            %     probability of the cdf (between 0 and 1)
            % a : float
            %     shape parameter (positive)
            % b : float
            %     scale parameter (positive)
            %
            % returns:
            % val : float
            %     the inverse cdf value
            %
            % references:
            %     Abramowitz and Stegun (1964), "Handbook of Mathematical Functions"
            %     section 26.1
            
            val = gammaincinv(p,a) * b;            
        end
        
        
        function [val] = chi2_icdf(p, nu)
    
            % [val] = chi2_icdf(p, nu)
            % inverse cdf for the chi2 distribution, open source equivalent of chi2inv
            %
            % parameters:
            % p : float
            %     probability of the cdf (between 0 and 1)
            % nu : float
            %     degrees of freedom (positive)
            %
            % returns:
            % val : float
            %     the inverse cdf value
            %
            % references:
            %     Abramowitz and Stegun (1964), "Handbook of Mathematical Functions"
            %     section 26.4
            
            val = su.gamma_icdf(p, nu / 2, 2);           
        end
        
        
    end
    
    
end
            
            
            
            
            
            
            
            
            
            
            
