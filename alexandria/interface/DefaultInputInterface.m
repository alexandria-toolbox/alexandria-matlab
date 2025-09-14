classdef DefaultInputInterface < handle


    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)
        

        function self = DefaultInputInterface()
        end
        
        
        function create_default_inputs(self)
        
            % overall structure
            user_inputs = struct;

            % sub structures for each of the three tabs taking user inputs
            user_inputs.tab_1 = struct;
            user_inputs.tab_2_lr = struct;
            user_inputs.tab_2_var = struct;
            user_inputs.tab_2_ext = struct;
            user_inputs.tab_3 = struct;

            % default values for tab 1
            user_inputs.tab_1.model = 1;
            user_inputs.tab_1.endogenous_variables = '';
            user_inputs.tab_1.exogenous_variables = '';
            user_inputs.tab_1.frequency = 1;
            user_inputs.tab_1.sample = '';
            user_inputs.tab_1.project_path = pwd;
            user_inputs.tab_1.data_file = '';
            user_inputs.tab_1.progress_bar = true;
            user_inputs.tab_1.create_graphics = true;
            user_inputs.tab_1.save_results = true;

            % default values for tab 2, linear regression
            user_inputs.tab_2_lr.regression_type = 1;
            user_inputs.tab_2_lr.iterations = '2000';
            user_inputs.tab_2_lr.burnin = '1000';
            user_inputs.tab_2_lr.model_credibility = '0.95';
            user_inputs.tab_2_lr.b = '0';
            user_inputs.tab_2_lr.V = '1';
            user_inputs.tab_2_lr.alpha = '0.0001';
            user_inputs.tab_2_lr.delta = '0.0001';
            user_inputs.tab_2_lr.g = '0';
            user_inputs.tab_2_lr.Q = '100';
            user_inputs.tab_2_lr.tau = '0.001';
            user_inputs.tab_2_lr.thinning = false;
            user_inputs.tab_2_lr.thinning_frequency = '10';
            user_inputs.tab_2_lr.Z_variables = '';
            user_inputs.tab_2_lr.q = '1';
            user_inputs.tab_2_lr.p = '0';
            user_inputs.tab_2_lr.H = '100';
            user_inputs.tab_2_lr.constant = true;
            user_inputs.tab_2_lr.b_constant = '0';
            user_inputs.tab_2_lr.V_constant = '1';
            user_inputs.tab_2_lr.trend = false;
            user_inputs.tab_2_lr.b_trend = '0';
            user_inputs.tab_2_lr.V_trend = '1';
            user_inputs.tab_2_lr.quadratic_trend = false;
            user_inputs.tab_2_lr.b_quadratic_trend = '0';
            user_inputs.tab_2_lr.V_quadratic_trend = '1';
            user_inputs.tab_2_lr.insample_fit = false;
            user_inputs.tab_2_lr.marginal_likelihood = false;
            user_inputs.tab_2_lr.hyperparameter_optimization = false;
            user_inputs.tab_2_lr.optimization_type = 1;

            % default values for tab 2, vector autoregression
            user_inputs.tab_2_var.var_type = 1; 
            user_inputs.tab_2_var.iterations = '2000';
            user_inputs.tab_2_var.burnin = '1000';
            user_inputs.tab_2_var.model_credibility = '0.95';
            user_inputs.tab_2_var.constant = true;
            user_inputs.tab_2_var.trend = false;
            user_inputs.tab_2_var.quadratic_trend = false;
            user_inputs.tab_2_var.lags = '4';
            user_inputs.tab_2_var.ar_coefficients = '0.9';
            user_inputs.tab_2_var.pi1 = '0.1';
            user_inputs.tab_2_var.pi2 = '0.5';
            user_inputs.tab_2_var.pi3 = '1';
            user_inputs.tab_2_var.pi4 = '100';
            user_inputs.tab_2_var.pi5 = '1';
            user_inputs.tab_2_var.pi6 = '0.1';
            user_inputs.tab_2_var.pi7 = '0.1';
            user_inputs.tab_2_var.proxy_variables = '';
            user_inputs.tab_2_var.lamda = '0.2';
            user_inputs.tab_2_var.proxy_prior = 1; 
            user_inputs.tab_2_var.insample_fit = false; 
            user_inputs.tab_2_var.constrained_coefficients = false; 
            user_inputs.tab_2_var.sums_of_coefficients = false; 
            user_inputs.tab_2_var.initial_observation = false; 
            user_inputs.tab_2_var.long_run = false; 
            user_inputs.tab_2_var.stationary = false; 
            user_inputs.tab_2_var.marginal_likelihood = false; 
            user_inputs.tab_2_var.hyperparameter_optimization = false; 
            user_inputs.tab_2_var.coefficients_file = ''; 
            user_inputs.tab_2_var.long_run_file = ''; 

            % default values for tab 2, VEC/VARMA (VAR extensions)
            user_inputs.tab_2_ext.model = 1;
            user_inputs.tab_2_ext.iterations = '3000';
            user_inputs.tab_2_ext.burnin = '1000';
            user_inputs.tab_2_ext.model_credibility = '0.95';
            user_inputs.tab_2_ext.constant = true;
            user_inputs.tab_2_ext.trend = false;
            user_inputs.tab_2_ext.quadratic_trend = false;
            user_inputs.tab_2_ext.vec_lags = '4';
            user_inputs.tab_2_ext.vec_pi1 = '0.1';
            user_inputs.tab_2_ext.vec_pi2 = '0.5';
            user_inputs.tab_2_ext.vec_pi3 = '1';
            user_inputs.tab_2_ext.vec_pi4 = '100';
            user_inputs.tab_2_ext.prior_type = 1;
            user_inputs.tab_2_ext.error_correction_type = 1;
            user_inputs.tab_2_ext.max_cointegration_rank = '1';
            user_inputs.tab_2_ext.varma_lags = '4';
            user_inputs.tab_2_ext.ar_coefficients = '0.95';
            user_inputs.tab_2_ext.varma_pi1 = '0.1';
            user_inputs.tab_2_ext.varma_pi2 = '0.5';
            user_inputs.tab_2_ext.varma_pi3 = '1';
            user_inputs.tab_2_ext.varma_pi4 = '100';
            user_inputs.tab_2_ext.residual_lags = '1';
            user_inputs.tab_2_ext.lambda1 = '0.1';
            user_inputs.tab_2_ext.lambda2 = '0.5';
            user_inputs.tab_2_ext.lambda3 = '1'; 

            % default values for tab 3
            user_inputs.tab_3.forecast = false;
            user_inputs.tab_3.conditional_forecast = false;
            user_inputs.tab_3.irf = false;
            user_inputs.tab_3.fevd = false;
            user_inputs.tab_3.hd = false;
            user_inputs.tab_3.forecast_credibility = '0.95';
            user_inputs.tab_3.conditional_forecast_credibility = '0.95';
            user_inputs.tab_3.irf_credibility = '0.95';
            user_inputs.tab_3.fevd_credibility = '0.95';
            user_inputs.tab_3.hd_credibility = '0.95';
            user_inputs.tab_3.forecast_periods = '';
            user_inputs.tab_3.conditional_forecast_type = 1;
            user_inputs.tab_3.forecast_file = '';
            user_inputs.tab_3.conditional_forecast_file = '';
            user_inputs.tab_3.forecast_evaluation = false;
            user_inputs.tab_3.irf_periods = '';
            user_inputs.tab_3.structural_identification = 1;
            user_inputs.tab_3.structural_identification_file = '';

            % save as attribute
            self.user_inputs = user_inputs;
        end
        
        
        function reset_default_inputs(self)

            % tab 1
            set(self.t1_mnu1, 'Value', self.user_inputs.tab_1.model);
            set(self.t1_edt1, 'String', self.user_inputs.tab_1.endogenous_variables);
            set(self.t1_edt2, 'String', self.user_inputs.tab_1.exogenous_variables);
            set(self.t1_mnu2, 'Value', self.user_inputs.tab_1.frequency);
            set(self.t1_edt3, 'String', self.user_inputs.tab_1.sample);
            set(self.t1_edt4, 'String', self.user_inputs.tab_1.project_path);
            set(self.t1_edt5, 'String', self.user_inputs.tab_1.data_file);
            set(self.t1_bgr1, 'SelectedObject', self.t1_rdb1);
            set(self.t1_bgr2, 'SelectedObject', self.t1_rdb3);
            set(self.t1_bgr3, 'SelectedObject', self.t1_rdb5);

            % tab 2, linear regression
            set(self.t2_lr_bgr1, 'SelectedObject', self.t2_lr_rdb1);
            set(self.t2_lr_edt1, 'String', self.user_inputs.tab_2_lr.iterations);
            set(self.t2_lr_edt2, 'String', self.user_inputs.tab_2_lr.burnin);
            set(self.t2_lr_edt3, 'String', self.user_inputs.tab_2_lr.model_credibility);
            set(self.t2_lr_edt4, 'String', self.user_inputs.tab_2_lr.b);
            set(self.t2_lr_edt5, 'String', self.user_inputs.tab_2_lr.V);
            set(self.t2_lr_edt6, 'String', self.user_inputs.tab_2_lr.alpha);
            set(self.t2_lr_edt7, 'String', self.user_inputs.tab_2_lr.delta);
            set(self.t2_lr_edt8, 'String', self.user_inputs.tab_2_lr.g);
            set(self.t2_lr_edt9, 'String', self.user_inputs.tab_2_lr.Q);
            set(self.t2_lr_edt10, 'String', self.user_inputs.tab_2_lr.tau);
            set(self.t2_lr_cbx1, 'Value', self.user_inputs.tab_2_lr.thinning);
            set(self.t2_lr_edt11, 'String', self.user_inputs.tab_2_lr.thinning_frequency);
            set(self.t2_lr_edt12, 'String', self.user_inputs.tab_2_lr.Z_variables);
            set(self.t2_lr_edt13, 'String', self.user_inputs.tab_2_lr.q);
            set(self.t2_lr_edt14, 'String', self.user_inputs.tab_2_lr.p);
            set(self.t2_lr_edt15, 'String', self.user_inputs.tab_2_lr.H);
            set(self.t2_lr_cbx2, 'Value', self.user_inputs.tab_2_lr.constant);
            set(self.t2_lr_edt16, 'String', self.user_inputs.tab_2_lr.b_constant);
            set(self.t2_lr_edt17, 'String', self.user_inputs.tab_2_lr.V_constant);
            set(self.t2_lr_cbx3, 'Value', self.user_inputs.tab_2_lr.trend);
            set(self.t2_lr_edt18, 'String', self.user_inputs.tab_2_lr.b_trend);
            set(self.t2_lr_edt19, 'String', self.user_inputs.tab_2_lr.V_trend);
            set(self.t2_lr_cbx4, 'Value', self.user_inputs.tab_2_lr.quadratic_trend);
            set(self.t2_lr_edt20, 'String', self.user_inputs.tab_2_lr.b_quadratic_trend);
            set(self.t2_lr_edt21, 'String', self.user_inputs.tab_2_lr.V_quadratic_trend);
            set(self.t2_lr_cbx5, 'Value', self.user_inputs.tab_2_lr.insample_fit);
            set(self.t2_lr_cbx6, 'Value', self.user_inputs.tab_2_lr.marginal_likelihood);
            set(self.t2_lr_cbx7, 'Value', self.user_inputs.tab_2_lr.hyperparameter_optimization);
            set(self.t2_lr_bgr2, 'SelectedObject', self.t2_lr_rdb7);
            
            % tab 3
            set(self.t3_bgr1, 'SelectedObject', self.t3_rdb2);
            set(self.t3_edt1, 'String', self.user_inputs.tab_3.forecast_credibility);
            set(self.t3_bgr2, 'SelectedObject', self.t3_rdb4);
            set(self.t3_edt2, 'String', self.user_inputs.tab_3.conditional_forecast_credibility);
            set(self.t3_bgr3, 'SelectedObject', self.t3_rdb6);
            set(self.t3_edt3, 'String', self.user_inputs.tab_3.irf_credibility);
            set(self.t3_bgr4, 'SelectedObject', self.t3_rdb8);
            set(self.t3_edt4, 'String', self.user_inputs.tab_3.fevd_credibility);
            set(self.t3_bgr5, 'SelectedObject', self.t3_rdb10);
            set(self.t3_edt5, 'String', self.user_inputs.tab_3.hd_credibility);
            set(self.t3_edt6, 'String', self.user_inputs.tab_3.forecast_periods);
            set(self.t3_mnu1, 'Value', self.user_inputs.tab_3.conditional_forecast_type);
            set(self.t3_edt7, 'String', self.user_inputs.tab_3.forecast_file);
            set(self.t3_edt8, 'String', self.user_inputs.tab_3.conditional_forecast_file);
            set(self.t3_cbx1, 'Value', self.user_inputs.tab_3.forecast_evaluation);
            set(self.t3_edt9, 'String', self.user_inputs.tab_3.irf_periods);
            set(self.t3_mnu2, 'Value', self.user_inputs.tab_3.structural_identification);
            set(self.t3_edt10, 'String', self.user_inputs.tab_3.structural_identification_file);
            
        end

    end

end






