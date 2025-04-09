classdef la
    

    % la stands for LinearAlgebra
    % A class containing static methods for linear algebra

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
        
        
        function [G] = cholesky_nspd(A)
    
            % [G] = cholesky_nspd(A)
            % produces the Cholesky factor of either A, or the spd matrix nearest to A
            %
            % parameters:
            % A : matrix of size (n,n)
            %     possibly symmetric positive definite
            %
            % returns:
            % G : matrix of size (n,n)
            %     lower triangular Cholesky factor of A

            % test first for Cholesky factorisation
            try
                % stop if factorisation is succesful
                G = chol(A, 'lower');
            catch
                % if factorisation fails, symmetrize
                B = (A + A') / 2;
                % compute H, the symmetric polar factor of B
                [~, S, V] = svd(B);
                H = V * S * V';
                % try again to obtain an spd matrix
                C = (B + H) / 2;
                C = (C + C') / 2;
                % if matrix is only semi spd, modify it slightly until it becomes spd
                spd = false;
                k = 0;
                identity = eye(size(A,1));
                while ~ spd
                    try
                        G = chol(C, 'lower');
                        spd = true;
                    catch
                        k = k + 1;
                        mineig = min(eig(C));
                        C = C + (-mineig * k.^2 + eps) * identity;
                    end
                end
            end
        end

        
        function [X_inverse] = invert_lower_triangular_matrix(X)
    
            % [X_inverse] = invert_lower_triangular_matrix(X)
            % computes inverse of lower triangular matrix efficiently, by forward substitution
            % 
            % parameters:
            % X : matrix of shape (n,n)
            %     lower triangular matrix
            % 
            % returns:
            % X_inverse : matrix of shape (n,n)
            %     inverse of X

            % dimension of G
            dimension = size(X,1);
            % invert it by back substitution
            X_inverse = X \ eye(dimension);
        end

        
        function [A_inverse] = invert_spd_matrix(A)
    
            % [A_inverse] = invert_spd_matrix(A)
            % produces the inverse of A, where A is symmetric and positive definite
            % based on property m.35
            %
            % parameters:
            % A : matrix of size (n,n)
            %     invertible, symmetric and positive definite
            %
            % returns:
            % A_inverse : matrix of size (n,n)
            %     inverse of A

            % size of A
            dimension = size(A,1);
            % take the Cholesky factor of A and invert it
            G = la.cholesky_nspd(A);
            G_inverse = G \ eye(dimension);
            % recover inverse of A
            A_inverse = G_inverse.' * G_inverse;
        end
        
        
        function [log_determinant_A] = determinant_spd_matrix(A)

            % determinant_spd_matrix(A)
            % produces the log determinant of A, where A is symmetric and positive definite
            % based on properties m.15, m.25 and m.28
            % 
            % parameters:
            % A : matrix of size (n,n)
            %     invertible, symmetric and positive definite
            % 
            % returns:
            % log_determinant_A : float
            %     log determinant of A

            G = la.cholesky_nspd(A);
            log_determinant_A = 2 * sum(log(diag(G)));
        end
        

        function [indices] = lower_triangular_indices(dimension)
            
            % [indices] = lower_triangular_indices(dimension)
            % provides the Matlab indices of lower triangular entries
            % 
            % parameters:
            % dimension : int
            %     dimension of the square matrix
            % 
            % returns:
            % indices : 
            %     vector of Matlab array indices

            indices = ones(dimension * (dimension - 1) / 2, 1);
            indices(cumsum([1 dimension-1:-1:2])) = 2:dimension;
            indices = cumsum(indices);
        end
        
        
        function [indices] = upper_triangular_indices(dimension)
            
            % [indices] = upper_triangular_indices(dimension)
            % provides the Matlab indices of upper triangular entries
            % 
            % parameters:
            % dimension : int
            %     dimension of the square matrix
            % 
            % returns:
            % indices : 
            %     vector of Matlab array indices

            indices = ones(dimension * (dimension - 1) / 2, 1);
            indices(1 + cumsum(0:dimension-2)) = dimension+1:-1:3;
            indices = cumsum(indices);
        end
        
        
        function [F, L] = triangular_factorization(X, varargin)
            
            % [F, L] = triangular_factorization(X, varargin)
            % triangular factorization of spd matrix
            % based on property m.34
            %
            % parameters:
            % X : matrix of size (n,n) 
            %     matrix to factor (symmetric and positive definite)
            % varargin : up to one additional input, bool
            %     if true, assumes that spd matrix X is replaced by its Cholesky factor
            %
            % returns:
            % F : matrix of size (n,n) 
            %     unit lower triangular factor
            % L : matrix of size (n,1) 
            %     diagonal factor (only reports the diagonal)

            if (isempty(varargin) || varargin{1} == false)
                G = la.cholesky_nspd(X);
            elseif varargin{1} == true
                G = X;
            end 
            diagonal = diag(G)';
            F = G ./ diagonal;
            L = diagonal.^2;
        end
        
        
        function [x] = vec(X)
            
            % [x] = vec(X)
            % vectorize a matrix, column by column
            % 
            % parameters:
            % X : matrix of size (n,m) 
            %     matrix to vectorize
            % 
            % returns:
            % x : matrix of size (n*m,1)  
            %     vector resulting from vectorization

            x = X(:);
        end
        
        
        function [det] = stable_determinant(A)

            % [det] = stable_determinant(A)
            % stable (log) determinant of a matrix of the form I + A
            % uses property m.59
            %
            % parameters:
            % A : matrix of size (n,n)
            %     matrix A in the form I + A
            %
            % returns:
            % det : float
            %     log determinant of matrix I + A

            w = eig(A);
            det = real(sum(log(1 + w)));
        end
        

        function [X_present, X_lagged] = lag_matrix(X, lags)
            
            % calculates arrays of present and lagged values of X
            % 
            % parameters:
            % X : matrix of size (n,k) or (n,1)
            %     array to be lagged
            % 
            % lags : int
            %     number of lags to apply
            % 
            % returns:
            % X_present: matrix of size (n-lags,k) or (n-lags,1)
            %     array containing present values
            % 
            % X_lagged: ndarray of shape (n-lags,k*lags) or (n-lags,lags)
            %     array containing lagged values

            % array of present values
            X_present = X(lags+1:end,:);
            % first lag
            X_lagged = X(lags:end-1,:);
            % concatenate remaining lags
            for i=2:lags
                X_lagged = [X_lagged X(lags+1-i:end-i,:)];
            end
        end
        
        
        function [Y] = lag_polynomial(X, gamma)
            
            % Applies lag polynomial, defined in (3.9.56)
            % 
            % parameters:
            % X : matrix of size (n,k) or (n,1)
            %     array on which lag polynomial applies
            % 
            % gamma : matrix of size (lags,1)
            %     coefficients of the lag polynomial
            % 
            % returns:
            % Y : ndarray of shape (n-lags,k)
            %     array obtained after applying the lag polynomial

            % number of lags
            lags = size(gamma,1);
            % create present and lagged arrays
            [Y, Z] = la.lag_matrix(X, lags);
            % if Y is 1d, apply lag polynomial directly
            if isvector(Y)
                Y = Y - Z * gamma;
            % else if Y is 2d, looping over lags is necessary
            else
                columns = size(X, 2);
                for i=1:lags
                    Y = Y - gamma(i) * Z(:,columns*(i-1)+1:columns*i);
                end
            end
        end
        
        
    end
end