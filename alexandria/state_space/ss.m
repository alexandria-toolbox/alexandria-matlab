classdef ss
    

    % ss stands for state-space
    % A class containing static methods for state-space model utilities

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
        
        
        function [Z_tt Z_tt1 Ups_tt Ups_tt1] = kalman_filter(X, A, Omega, C, B, Upsilon, T, n, k)
            
            % [Z_tt Z_tt1 Ups_tt Ups_tt1] = kalman_filter(X, A, Omega, C, B, Upsilon, T, n, k)
            % Kalman filter to estimate the state variables of a general state-space model
            %
            % parameters:
            % X : matrix of size (T,n)
            %     matrix of observed variables
            % A : matrix of size (n,k,T)
            %     matrix of coefficients on observation equation
            % Omega : matrix of size (n,n,T)
            %     variance-covariance matrix of observation errors
            % C : matrix of size (T,k)
            %     intercept on observation equation  
            % B : matrix of size (k,k,T)
            %     matrix of coefficients on state equation        
            % Upsilon : matrix of size (k,k,T)
            %     variance-covariance matrix of state errors      
            % T : int
            %     number of sample periods          
            % n : int
            %     dimension of observation vector 
            % k : int
            %     dimension of state vector 
            % 
            % returns:
            % Z_tt : matrix of size (T,k)
            %     matrix of state values z_t|t          
            % Z_tt1 : matrix of size (T,k)
            %     matrix of state values z_t|t-1 
            % Ups_tt : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t       
            % Ups_tt1 : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t-1  
        
            % initiate values
            z_t1t1 = zeros(k,1);
            Upsilon_t1t1 = zeros(k,k);
            Z_tt = zeros(T,k);
            Z_tt1 = zeros(T,k);
            Ups_tt = zeros(k,k,T);
            Ups_tt1 = zeros(k,k,T);
            % Kalman recursions
            for t=1:T
                % period-specific parameters
                x_t = X(t,:)';
                A_t = A(:,:,t);
                Omega_t = Omega(:,:,t);
                c_t = C(t,:)';
                B_t = B(:,:,t);
                Upsilon_t = Upsilon(:,:,t);
                % step 1
                z_tt1 = c_t + B_t * z_t1t1;
                % step 2
                Upsilon_tt1 = B_t * Upsilon_t1t1 * B_t' + Upsilon_t;
                % step 3
                x_tt1 = A_t * z_tt1;
                % step 4
                Omega_tt1 = A_t * Upsilon_tt1 * A_t' + Omega_t;
                % Phi_t computation
                Phi_t = Upsilon_tt1 * A_t' / Omega_tt1;
                % step 5
                z_tt = z_tt1 + Phi_t * (x_t - x_tt1);
                % step 6
                Upsilon_tt = Upsilon_tt1 - Phi_t * Omega_tt1 * Phi_t';
                % record and update for incoming period
                Z_tt(t,:) = z_tt';
                Z_tt1(t,:) = z_tt1';
                z_t1t1 = z_tt;
                Ups_tt(:,:,t) = Upsilon_tt;
                Ups_tt1(:,:,t) = Upsilon_tt1;
                Upsilon_t1t1 = Upsilon_tt;
            end
        end

        
        function [Z] = backward_pass(Z_tt, Z_tt1, Ups_tt, Ups_tt1, B, T, k)
            
            % [Z] = backward_pass(Z_tt, Z_tt1, Ups_tt, Ups_tt1, B, T, k)
            % Backward pass of Carter-Kohn algorithm (algorithm k.2)
            %
            % parameters:
            % Z_tt : matrix of size (T,k)
            %     matrix of state values z_t|t          
            % Z_tt1 : matrix of size (T,k)
            %     matrix of state values z_t|t-1 
            % Ups_tt : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t       
            % Ups_tt1 : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t-1
            % B : matrix of size (k,k,T)
            %     matrix of coefficients on state equation  
            % T : int
            %     number of sample periods           
            % k : int
            %     dimension of state vector             
            % 
            % returns:
            % Z : matrix of size (T,k)
            %     matrix of sampled values for the state variables     
        
            % initiate values
            Z = zeros(T,k);
            % final period sampling
            z_TT = Z_tt(end,:)';
            Upsilon_TT = Ups_tt(:,:,end);
            Z(end,:) = rn.multivariate_normal(z_TT, Upsilon_TT);
            % backward pass, other periods
            for t=T-1:-1:1
                % period-specific parameters
                B_t1 = B(:,:,t+1);
                z_tt = Z_tt(t,:)';
                z_t1t = Z_tt1(t+1,:)';
                Upsilon_tt = Ups_tt(:,:,t);
                Upsilon_t1t = Ups_tt1(:,:,t+1);
                z_t1 = Z(t+1,:)';
                % Xi_t computation
                Xi_t = Upsilon_tt * B_t1' / Upsilon_t1t;
                % step 1
                z_bar_tt1 = z_tt + Xi_t * (z_t1 - z_t1t);
                % step 2
                Upsilon_bar_tt1 = Upsilon_tt - Xi_t * B_t1 * Upsilon_tt;
                % step 3
                Z(t,:) = rn.multivariate_normal(z_bar_tt1, Upsilon_bar_tt1);   
            end       
        end        

        
        
        function [Z_tt Z_tt1 Ups_tt Ups_tt1] = conditional_forecast_kalman_filter(X, A, Omega, C, B, Upsilon, z_00, Upsilon_00, T, n, k)
            
            % [Z_tt Z_tt1 Ups_tt Ups_tt1] = conditional_forecast_kalman_filter(X, A, Omega, C, B, Upsilon, z_00, Upsilon_00, T, n, k)
            % Kalman filter to estimate the state variables of a conditional forecast state-space model
            %
            % parameters:
            % X : matrix of size (T,n)
            %     matrix of observed variables
            % A : matrix of size (n,k)
            %     matrix of coefficients on observation equation
            % Omega : matrix of size (T,n)
            %     variance-covariance matrix of observation errors
            % C : matrix of size (T,k)
            %     intercept on observation equation  
            % B : matrix of size (k,k)
            %     matrix of coefficients on state equation        
            % Upsilon : matrix of size (k,k,T)
            %     variance-covariance matrix of state errors
            % z_00 : matrix of size (k,1)
            %     initial conditions for state variables (mean)    
            % Upsilon_00 : matrix of size (k,k)
            %     initial conditions for state variables (variance-covariance)              
            % T : int
            %     number of sample periods          
            % n : int
            %     dimension of observation vector 
            % k : int
            %     dimension of state vector 
            % 
            % returns:
            % Z_tt : matrix of size (T,k)
            %     matrix of state values z_t|t          
            % Z_tt1 : matrix of size (T,k)
            %     matrix of state values z_t|t-1 
            % Ups_tt : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t       
            % Ups_tt1 : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t-1  
        
            % initiate values
            z_t1t1 = z_00;
            Upsilon_t1t1 = Upsilon_00;
            Z_tt = zeros(T,k);
            Z_tt1 = zeros(T,k);
            Ups_tt = zeros(k,k,T);
            Ups_tt1 = zeros(k,k,T);
            % Kalman recursions
            for t=1:T
                % period-specific parameters
                x_t = X(t,:)';
                Omega_t = diag(Omega(t,:));
                c_t = C(t,:)';
                Upsilon_t = Upsilon(:,:,t);
                % step 1
                z_tt1 = c_t + B * z_t1t1;
                % step 2
                Upsilon_tt1 = B * Upsilon_t1t1 * B' + Upsilon_t;
                % step 3
                x_tt1 = A * z_tt1;
                % step 4
                Omega_tt1 = A * Upsilon_tt1 * A' + Omega_t;
                % Phi_t computation
                Phi_t = Upsilon_tt1 * A' / Omega_tt1;
                % step 5
                z_tt = z_tt1 + Phi_t * (x_t - x_tt1);
                % step 6
                Upsilon_tt = Upsilon_tt1 - Phi_t * Omega_tt1 * Phi_t';
                % record and update for incoming period
                Z_tt(t,:) = z_tt';
                Z_tt1(t,:) = z_tt1';
                z_t1t1 = z_tt;
                Ups_tt(:,:,t) = Upsilon_tt;
                Ups_tt1(:,:,t) = Upsilon_tt1;
                Upsilon_t1t1 = Upsilon_tt;
            end
        end        
        
        
        function [Z] = static_backward_pass(Z_tt, Z_tt1, Ups_tt, Ups_tt1, B, T, k)
            
            % [Z] = static_backward_pass(Z_tt, Z_tt1, Ups_tt, Ups_tt1, B, T, k)
            % Backward pass of Carter-Kohn algorithm (algorithm k.2) with static B
            %
            % parameters:
            % Z_tt : matrix of size (T,k)
            %     matrix of state values z_t|t          
            % Z_tt1 : matrix of size (T,k)
            %     matrix of state values z_t|t-1 
            % Ups_tt : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t       
            % Ups_tt1 : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t-1
            % B : matrix of size (k,k)
            %     matrix of coefficients on state equation  
            % T : int
            %     number of sample periods           
            % k : int
            %     dimension of state vector             
            % 
            % returns:
            % Z : matrix of size (T,k)
            %     matrix of sampled values for the state variables     
        
            % initiate values
            Z = zeros(T,k);
            % final period sampling
            z_TT = Z_tt(end,:)';
            Upsilon_TT = Ups_tt(:,:,end);
            Z(end,:) = rn.multivariate_normal(z_TT, Upsilon_TT);
            % backward pass, other periods
            for t=T-1:-1:1
                % period-specific parameters
                z_tt = Z_tt(t,:)';
                z_t1t = Z_tt1(t+1,:)';
                Upsilon_tt = Ups_tt(:,:,t);
                Upsilon_t1t = Ups_tt1(:,:,t+1) + 1e-10 * eye(k);
                z_t1 = Z(t+1,:)';
                % Xi_t computation
                Xi_t = Upsilon_tt * B' / Upsilon_t1t;
                % step 1
                z_bar_tt1 = z_tt + Xi_t * (z_t1 - z_t1t);
                % step 2
                Upsilon_bar_tt1 = Upsilon_tt - Xi_t * B * Upsilon_tt;
                % step 3
                Z(t,:) = rn.multivariate_normal(z_bar_tt1, Upsilon_bar_tt1);   
            end
        end         
        
        
        function [Z_tt Z_tt1 Ups_tt Ups_tt1] = varma_forward_pass(X, A, B, Upsilon, z_00, Upsilon_00, T, n, k)
            
            % varma_forward_pass(X, A, B, Upsilon, z_00, Upsilon_00, T, n, k)
            % forward pass for the state variables of a varma model
            %
            % parameters:
            % X : matrix of size (T,n)
            %     matrix of observed variables
            % A : matrix of size (n,k,T)
            %     matrix of coefficients on observation equation
            % B : matrix of size (k,k,T)
            %     matrix of coefficients on state equation        
            % Upsilon : matrix of size (k,k,T)
            %     variance-covariance matrix of state errors   
            % z_00 : matrix of size (k,1)
            %     initial conditions for state variables (mean)    
            % Upsilon_00 : matrix of size (k,k)
            %     initial conditions for state variables (variance-covariance)             
            % T : int
            %     number of sample periods          
            % n : int
            %     dimension of observation vector 
            % k : int
            %     dimension of state vector 
            % 
            % returns:
            % Z_tt : matrix of size (T,k)
            %     matrix of state values z_t|t          
            % Z_tt1 : matrix of size (T,k)
            %     matrix of state values z_t|t-1 
            % Ups_tt : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t       
            % Ups_tt1 : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t-1  
        
            % initiate values
            z_t1t1 = z_00;
            Upsilon_t1t1 = Upsilon_00;
            Z_tt = zeros(T,k);
            Z_tt1 = zeros(T,k);
            Ups_tt = zeros(k,k,T);
            Ups_tt1 = zeros(k,k,T);
            % Kalman recursions
            for t=1:T
                % period-specific parameters
                x_t = X(t,:)';
                % step 1
                z_tt1 = B * z_t1t1;
                % step 2
                Upsilon_tt1 = B * Upsilon_t1t1 * B' + Upsilon;
                % step 3
                x_tt1 = A * z_tt1;
                % step 4
                Omega_tt1 = A * Upsilon_tt1 * A';
                % Phi_t computation
                Phi_t = Upsilon_tt1 * A' / Omega_tt1;
                % step 5
                z_tt = z_tt1 + Phi_t * (x_t - x_tt1);
                % step 6
                Upsilon_tt = Upsilon_tt1 - Phi_t * Omega_tt1 * Phi_t';
                % record and update for incoming period
                Z_tt(t,:) = z_tt';
                Z_tt1(t,:) = z_tt1';
                z_t1t1 = z_tt;
                Ups_tt(:,:,t) = Upsilon_tt;
                Ups_tt1(:,:,t) = Upsilon_tt1;
                Upsilon_t1t1 = Upsilon_tt;
            end
        end        
     

        function [Z] = varma_backward_pass(Z_tt, Z_tt1, Ups_tt, Ups_tt1, B, T, k)
            
            % varma_backward_pass(Z_tt, Z_tt1, Ups_tt, Ups_tt1, B, T, k)
            % backward pass for the state variables of a varma model
            %
            % parameters:
            % Z_tt : matrix of size (T,k)
            %     matrix of state values z_t|t          
            % Z_tt1 : matrix of size (T,k)
            %     matrix of state values z_t|t-1 
            % Ups_tt : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t       
            % Ups_tt1 : matrix of size (k,k,T)
            %     matrix of state variance Upsilon_t|t-1
            % B : matrix of size (k,k)
            %     matrix of coefficients on state equation  
            % T : int
            %     number of sample periods           
            % k : int
            %     dimension of state vector             
            % 
            % returns:
            % Z : matrix of size (T,k)
            %     matrix of sampled values for the state variables     
        
            % initiate values
            Z = zeros(T,k);
            % final period sampling
            z_TT = Z_tt(end,:)';
            Upsilon_TT = Ups_tt(:,:,end);
            Z(end,:) = rn.multivariate_normal(z_TT, Upsilon_TT);
            % backward pass, other periods
            for t=T-1:-1:1
                % period-specific parameters
                z_tt = Z_tt(t,:)';
                z_t1t = Z_tt1(t+1,:)';
                Upsilon_tt = Ups_tt(:,:,t);
                Upsilon_t1t = Ups_tt1(:,:,t+1) + 1e-10 * eye(k);
                z_t1 = Z(t+1,:)';
                % Xi_t computation
                Xi_t = Upsilon_tt * B' / Upsilon_t1t;
                % step 1
                z_bar_tt1 = z_tt + Xi_t * (z_t1 - z_t1t);
                % step 2
                Upsilon_bar_tt1 = Upsilon_tt - Xi_t * B * Upsilon_tt;
                % step 3
                Z(t,:) = rn.multivariate_normal(z_bar_tt1, Upsilon_bar_tt1);   
            end
        end 


    end
    
    
end
        