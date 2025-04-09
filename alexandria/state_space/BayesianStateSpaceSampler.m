classdef BayesianStateSpaceSampler < handle
    
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
        X
        A
        Omega
        C
        B
        Upsilon
        z_00
        Upsilon_00
        kalman_type
        T
        n
        k
        Z_tt
        Z_tt1
        Ups_tt
        Ups_tt1
        Z
    end    
    
    properties (GetAccess = protected, SetAccess = protected)

    end 

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public) 
    
    
        function self = BayesianStateSpaceSampler(X, A, Omega, C, B, Upsilon, varargin)
            
            % constructor for the BayesianStateSpaceSampler class
            
            % allow for optional arguments
            parser = inputParser;
            default_z_00 = [];
            default_Upsilon_00 = [];
            default_kalman_type = 'standard';
            addRequired(parser, 'X');
            addRequired(parser, 'A');
            addRequired(parser, 'Omega');
            addRequired(parser, 'C');
            addRequired(parser, 'B');
            addRequired(parser, 'Upsilon');
            addParameter(parser, 'z_00', default_z_00);
            addParameter(parser, 'Upsilon_00', default_Upsilon_00);
            addParameter(parser, 'kalman_type', default_kalman_type);
            parse(parser, X, A, Omega, C, B, Upsilon, varargin{:});
            self.X = X;
            self.A = A;
            self.Omega = Omega;
            self.C = C;
            self.B = B;
            self.Upsilon = Upsilon;
            self.z_00 = parser.Results.z_00;
            self.Upsilon_00 = parser.Results.Upsilon_00;    
            self.kalman_type = parser.Results.kalman_type;
            self.T = size(X,1);
            self.n = size(X,2);
            self.k = size(A,2);
        end
        
        
        function carter_kohn_algorithm(self)

            % carter_kohn_algorithm()
            % Bayesian sampler for a given state-space system, using algorithm k.2
            % 
            % parameters:
            % none
            % 
            % returns:
            % none     

            self.forward_pass();
            self.backward_pass();
        end
        
        
    end
    
    
    methods (Access = protected, Hidden = true) 
        

        function forward_pass(self)
        
            % forward pass of the algorithm, based on Kalman filter
        
            if isequal(self.kalman_type, 'standard')
                [Z_tt Z_tt1 Ups_tt Ups_tt1] = ss.kalman_filter(self.X, self.A, self.Omega, ...
                                              self.C, self.B, self.Upsilon, self.T, self.n, self.k);
            elseif isequal(self.kalman_type, 'conditional_forecast')
                [Z_tt Z_tt1 Ups_tt Ups_tt1] = ss.conditional_forecast_kalman_filter(self.X, self.A, ...
                                              self.Omega, self.C, self.B, self.Upsilon, self.z_00, ...
                                              self.Upsilon_00, self.T, self.n, self.k);
            end
            self.Z_tt = Z_tt;
            self.Z_tt1 = Z_tt1;
            self.Ups_tt = Ups_tt;
            self.Ups_tt1 = Ups_tt1;    
        end
        
        
        function backward_pass(self)

            % backward pass of the algorithm, based on Kalman filter

            if isequal(self.kalman_type, 'standard')
                Z = ss.backward_pass(self.Z_tt, self.Z_tt1, self.Ups_tt, self.Ups_tt1, self.B, self.T, self.k);
            elseif isequal(self.kalman_type, 'conditional_forecast')
                Z = ss.static_backward_pass(self.Z_tt, self.Z_tt1, self.Ups_tt, self.Ups_tt1, self.B, self.T, self.k);
            end
            self.Z = Z;
        end
        
    
    end

    
end
    
    
    
    
    
    
    
    