classdef rng
    

    % rng stands for random number generators
    % A class containing static methods for random number generation

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
        
        function [x] = normal(mu, sigma)
    
            % [x] = normal(mu, sigma)
            % random number generator for the normal distribution
            % based on algorithm d.11
            %
            % parameters:
            % mu : float
            %     mean
            % sigma : float
            %     variance (positive)
            %
            % returns:
            % x : float 
            %     pseudo-random number from the normal distribution

            x = mu + sqrt(sigma) * randn;
        end
        
        
        function [x] = multivariate_normal(mu, Sigma)
    
            % [x] = multivariate_normal(mu, Sigma)
            % random number generator for the multivariate normal distribution
            % based on algorithm d.13
            %
            % parameters:
            % mu : matrix of shape (n,1)
            %     mean
            % Sigma : matrix of shape (n,n)
            %     variance-covariance (symmetric positive definite)
            %
            % returns:
            % x : matrix of shape (n,1)
            %     pseudo-random number from the multivariate normal distribution

            x = mu + la.cholesky_nspd(Sigma) * randn(size(Sigma,1), 1); 
        end

        
        function [x] = efficient_multivariate_normal(m, inv_Sigma)
            
            % [x] = efficient_multivariate_normal(m, inv_Sigma)
            % efficient random number generator for the multivariate normal distribution, using algorithm 9.4
            %
            % parameters:
            % m: matrix of size (n,1)
            %     partial term for the distribution mean             
            % inv_Sigma: matrix of size (n,n)
            %     inverse of variance matrix of the distribution
            %
            % returns:
            % x: matrix of size (n,1)
            %     pseudo-random vector from the multivariate normal distribution

            G = la.cholesky_nspd(inv_Sigma);
            zeta = randn(size(G,1), 1);
            temp = G \ m + zeta;
            x = G' \ temp;
        end          
        
        
        function [X] = matrix_normal(M, Sigma, Omega)
            % [X] = matrix_normal(M, Sigma, Omega)
            % random number generator for the matrix normal distribution
            % based on algorithm d.15
            %
            % parameters:
            % M : matrix of shape(n,m)
            %     location
            % Sigma : matrix of shape (n,n)
            %     row scale (symmetric positive definite)
            % Omega : matrix of shape (m,m)
            %     column scale (symmetric positive definite)
            %
            % returns:
            % X : matrix of shape (n,m)
            %     pseudo-random number from the matrix normal distribution

            X = M + la.cholesky_nspd(Sigma) * randn(size(M)) * la.cholesky_nspd(Omega)';
        end
        
        
        function [x] = gamma(a, b)
            
            % [x] = gamma(a, b)
            % random number generator for the Gamma distribution
            % based on algorithms d.25 and d.26
            % 
            % parameters:
            % a : float
            %     shape (positive)
            % b : float
            %     scale (positive)
            % 
            % returns:
            % x : float
            %     pseudo-random number from the Gamma distribution

            if a >= 1
                z = randg(a);
            else
                z = randg(a + 1);
                u = rand;
                z = z * u^(1 / a);
            end
            x = z * b;    
        end
        
        
        function [x] = inverse_gamma(a, b)
            
            % [x] = inverse_gamma(a, b)
            % random number generator for the inverse Gamma distribution
            % based on algorithm d.29
            % 
            % parameters:
            % a : float
            %     shape (positive)
            % b : float
            %     scale (positive)
            % 
            % returns:
            % x : float
            %     pseudo-random number from the inverse Gamma distribution

            x = 1 / rng.gamma(a, 1 / b);
        end
        
        
        function [x] = chi2(nu)
            
            % [x] = chi2(nu)
            % random number generator for the chi2 distribution
            % based on property in Table d.15
            % 
            % parameters:
            % nu : float
            %     degrees of freedom (positive)
            % 
            % returns:
            % x : float 
            %     pseudo-random number from the chi2 distribution

            x = 2 * randg(0.5 * nu);
        end
        
        
        function [X] = wishart(nu, S)
            
            % [X] = wishart(nu, S)
            % random number generator for the Wishart distribution
            % based on algorithms d.27 and d.28
            % 
            % parameters:
            % nu : float
            %     degrees of freedom (positive, nu >= n)
            % S : matrix of shape (n,n)
            %     scale (symmetric positive definite)
            % 
            % returns:
            % X : matrix of shape (n,n)
            %     pseudo-random number from the Wishart distribution

            dimension = size(S, 1);
            if nu == floor(nu) && nu <= 80 + dimension
                A = randn(dimension, nu);
            else
                degree_freedom = nu - (0:dimension-1);
                A = diag(sqrt(rng.chi2(degree_freedom)));
                index = la.lower_triangular_indices(dimension);
                index_size = dimension * (dimension - 1) / 2;
                A(index) = randn(index_size, 1);
            end
            G = la.cholesky_nspd(S);
            Z = G * A;
            X = Z * Z';
        end
        
        
        function [X] = inverse_wishart(nu, S)
            
            % [X] = inverse_wishart(nu, S)
            % random number generator for the inverse Wishart distribution
            % based on algorithms d.30 and d.31
            % 
            % parameters:
            % nu : float
            %     degrees of freedom (positive, nu >= n)
            % S : matrix of shape (n,n)
            %     scale (symmetric positive definite)
            % 
            % returns:
            % X : matrix of shape (n,n)
            %     pseudo-random number from the inverse Wishart distribution

            dimension = size(S, 1);
            G = la.cholesky_nspd(S);
            if nu == floor(nu) && nu <= 80 + dimension
                A = randn(dimension, nu);
                X = G / (A * A') * G';
            else
                degree_freedom = nu - (0:dimension-1);
                A = diag(sqrt(rng.chi2(degree_freedom)));
                index = la.upper_triangular_indices(dimension);
                index_size = dimension * (dimension - 1) / 2;
                A(index) = randn(index_size, 1);
                Z = G / A;
                X = Z * Z';
            end
        end
        
        
        function [x] = student(mu, sigma, nu)
            
            % [x] = student(mu, sigma, nu)
            % random number generator for the student distribution
            % based on algorithm d.17
            % 
            % parameters:
            % mu : float
            %     location
            % sigma : float
            %     scale (positive)
            % nu : float
            %     degrees of freedom
            % 
            % returns:
            % x : float
            %     pseudo-random number from the student distribution

            s = rng.inverse_gamma(0.5 * nu, 0.5 * nu);
            z = s^0.5 * randn;
            x = sigma^0.5 * z + mu;
        end
        
        
        function [x] = multivariate_student(mu, Sigma, nu)
            
            % [x] = multivariate_student(mu, Sigma, nu)
            % random number generator for the multivariate student distribution
            % based on algorithm d.19
            % 
            % parameters:
            % mu : matrix of shape (n,1)
            %     location
            % Sigma : matrix of shape (n,n)
            %     scale (symmetric positive definite)
            % nu : float
            %     degrees of freedom
            % 
            % returns:
            % x : matrix of shape (n,1)
            %     pseudo-random number from the multivariate student distribution

            s = rng.inverse_gamma(0.5 * nu, 0.5 * nu);
            z = s^0.5 * randn(size(Sigma,1),1);
            x = mu + la.cholesky_nspd(Sigma) * z;
        end
        
        
        function [X] = matrix_student(M, Sigma, Omega, nu)
            
            % [X] = matrix_student(M, Sigma, Omega, nu)
            % random number generator for the matrix student distribution
            % based on algorithm d.21
            % 
            % parameters:
            % M : matrix of shape(n,m)
            %     location
            % Sigma : matrix of shape (n,n)
            %     row scale (symmetric positive definite)
            % Omega : matrix of shape (m,m)
            %     column scale (symmetric positive definite)
            % nu : float
            %     degrees of freedom (nu >= 3, nu >= min(n,m))
            % 
            % returns:
            % X : matrix of shape(n,m)
            %     pseudo-random number from the matrix normal distribution

            [rows, columns] = size(M);
            Z = randn(rows, columns);
            if columns <= rows
                Phi = rng.inverse_wishart(nu + columns - 1, nu * Omega);
                G = la.cholesky_nspd(Sigma);
                H = la.cholesky_nspd(Phi);
            else
                Phi = rng.inverse_wishart(nu + rows - 1, nu * Sigma);
                G = la.cholesky_nspd(Phi);
                H = la.cholesky_nspd(Omega);
            end
            X = M + G * Z * H';
        end
        
        
        function [x] = truncated_normal(mu, sigma, a, b)
            
            % [x] = truncated_normal(mu, sigma, a, b)
            % random number generator for the truncated normal distribution
            % based on algorithm d.23
            % 
            % parameters:
            % mu : float
            %     mean of untruncated distribution
            % sigma : float
            %     variance of untruncated distribution (positive)
            % a : float
            %     lower bound of truncation
            % b : float
            %     upper bound of truncation (b > a)
            % 
            % returns:
            % x : float
            %     pseudo-random number from the truncated normal distribution

            function [x] = standard_truncated_normal(a, b)
                while true
                    x = a + (b - a) * rand;
                    if a <= 0 && b >= 0
                        w = exp(- x * x / 2);
                    elseif b < 0
                        w = exp((b * b - x * x) / 2);
                    else
                        w = exp((a * a - x * x) / 2);
                    end
                    u = rand;
                    if u <= w
                        return
                    end
                end
            end

            standard_deviation = sigma^0.5;
            a_bar = (a - mu) / standard_deviation;
            b_bar = (b - mu) / standard_deviation;
            z = standard_truncated_normal(a_bar, b_bar);
            x = mu + standard_deviation * z;
        end
        
        
        function [x] = beta(a, b)
            
            % [x] = beta(a, b)
            % random number generator for the beta distribution
            % based on algorithm d.32
            % 
            % parameters:
            % a : float
            %     shape (positive)
            % b : float
            %     shape (positive)
            % 
            % returns:
            % x : float
            %     pseudo-random number from the beta distribution

            y = randg(a);
            x = y / (y + randg(b));
        end
        
        
        function [x] = dirichlet(a)
            
            % dirichlet(a)
            % random number generator for the Dirichlet distribution
            % based on algorithm d.33
            % 
            % parameters:
            % a : matrix of shape(n,1)
            %     concentration (all positive)
            % 
            % returns:
            % x : matrix of shape(n,1)
            %     pseudo-random number from the Dirichlet distribution

            z = randg(a);
            x = z / sum(z);
        end
    
        
        function [x] = bernoulli(p)
            
            % [x] = bernoulli(p)
            % random number generator for the Bernoulli distribution
            % based on algorithm d.2
            % 
            % parameters:
            % p : float
            %     probability of success (0 <= p <= 1)
            % 
            % returns:
            % x : int
            %     pseudo-random number from the Bernoulli distribution

            x = double(rand <= p);
        end
        
        
        function [x] = binomial(n, p)
            
            % [x] = binomial(n, p)
            % random number generator for the binomial distribution
            % based on algorithm d.4
            % 
            % parameters:
            % n : int 
            %     number of trials (positive)
            % p : float
            %     probability of success (0 <= p <= 1)
            % 
            % returns:
            % x : int
            %     pseudo-random number from the binomial distribution

            x = sum(rand(n, 1) <= p);
        end
        
        
        function [x] = categorical(p)
            
            % [x] = categorical(p)
            % random number generator for the categorical distribution
            % based on algorithm d.3
            % 
            % parameters:
            % p : matrix of shape(k,1)
            %     probability of success for each category (sum = 1)
            % 
            % returns:
            % x : int
            %     pseudo-random number from the categorical distribution

            x = find(rand <= cumsum(p), 1, 'first');
        end
        
        
        function [x] = multinomial(n, p)
            
            % [x] = multinomial(n, p)
            % random number generator for the multinomial distribution
            % based on algorithm d.5
            % 
            % parameters:
            % n : int
            %     number of trials (positive)
            % p : float
            %     probability of success for each category (sum = 1)
            % 
            % returns:
            % x : matrix of shape(k,1)
            %     pseudo-random number from the multinomial distribution

            x = histc(rand(1, n), [0 cumsum(p)]);
            x = [x(1:end-2) x(end-1)+x(end)];
        end
        
        
        function [x] = discrete_uniform(a, b)
            
            % [x] = discrete_uniform(a, b)
            % random number generator for the discrete uniform distribution
            % based on algorithm d.1
            % 
            % parameters:
            % a : int 
            %     lower bound
            % b : int
            %     upper bound
            % 
            % returns:
            % x : int
            %     pseudo-random numbers from the discrete uniform distribution

            x = floor(a + (b + 1 - a) * rand);
        end
        
        
        function [x] = poisson(lambda)
            
            % [x] = poisson(lambda)
            % random number generator for the Poisson distribution
            % based on algorithms d.6 and d.7
            % 
            % parameters:
            % lambda : float
            %     intensity (positive)
            % 
            % returns:
            % x : float
            %     pseudo-random number from the Poisson distribution

            if lambda <= 30
                p = exp(-lambda);
                F = p;
                z = 0;
                u = rand;
                while u > F
                    z = z + 1;
                    p = lambda * p / z;
                    F = F + p;
                end
                x = z;

            else
                % step 2
                lambda_bar = floor(lambda + 0.5);
                alph = lambda - lambda_bar;
                c = (2 * pi * lambda_bar)^(-0.5);
                % step 4
                p_r = pr(lambda_bar, alph);
                % step 5
                u = rand;
                % step 6
                if u <= 0.5
                    step = 9;
                elseif u >= 0.5 + 7 * c / 6
                    step = 12;
                else
                    step = 7;
                end
                % step 7
                if step == 7
                    F_r = Fr(lambda_bar, alph);   
                    step = 8;
                end
                % step 8
                if step == 8
                    if u > F_r
                        step = 12;
                    else
                        step = 9;
                    end
                end
                % step 9
                if step == 9
                    if u < p_r
                        x = lambda_bar;
                        return
                    else
                        step = 10;
                    end
                end
                % step 10
                if step == 10
                    p = p_r;
                    step = 11;
                end
                % step 11
                if step == 11
                    step = 12;
                    for i = 0 : lambda_bar - 1
                        u = u - p;
                        p = (lambda_bar - i) * p / lambda;
                        if u < p
                            x = lambda_bar - i - 1;
                            return
                        end
                    end
                end
                % step 12
                if step == 12
                    u = 1 - u;
                    p = p_r;
                    step = 13;
                end
                % step 13
                if step == 13
                    lambda_max = 2 * lambda_bar + 30;
                    for i = lambda_bar + 1 : lambda_max
                        p = p * lambda / i;
                        if u < p
                            x = i;
                            return
                        end
                        u = u - p;
                    end
                end
            end

            function [x] = qr(lambda)
                x = (2 * pi * lambda)^(-0.5) * (1 - 1 / ...
                    (12 * lambda + 0.5 + 293 / (720 * lambda)));
            end

            function [x] = Gr(lambda)
                x = 0.5 + (2/3) * (2 * pi * lambda)^(-0.5) ...
                    * (1 - (23/15) / (12 * lambda + 15/14 + (30557/4508) ...
                    / (12 * lambda + 138134432 / 105880005)));
            end

            function [x]= pr(lambdabar,alpha)
                x=qr(lambdabar) * ((lambdabar + 2 * alpha / 3 - alpha^2 / 4 ...
                    - alpha^2 / (18 * lambdabar)) / (lambdabar + 2 * alpha / 3 ...
                    + alpha^2 / 4 - alpha^2 / (18 * lambdabar)));
            end

            function [x] = Fr(lambdabar,alpha)
                x = Gr(lambdabar) - alpha * qr(lambdabar)* ((lambdabar + ...
                    alpha / 2 - alpha^2 / 60 - alpha^2 / (20 * lambdabar)) ...
                    / (lambdabar + alpha / 2 + 3 * alpha^2 / 20 - alpha^2 ...
                    / (20 * lambdabar)));
            end

        end
    
    end
  
end