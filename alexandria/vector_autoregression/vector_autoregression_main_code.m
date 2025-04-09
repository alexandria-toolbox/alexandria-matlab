function [var] = vector_autoregression_main_code(user_inputs)
    

    %---------------------------------------------------
    % Alexandria header and estimation start
    %--------------------------------------------------- 


    % display Alexandria header
    cu.print_alexandria_header();
    cu.print_start_message();
    
    
    %---------------------------------------------------
    % User input processing
    %---------------------------------------------------   
    
    
    % recover user inputs from input processor    
    ip = InputProcessor(user_inputs);
    ip.process_input();
    
    % recover specification parameters (interface 1)
    project_path = ip.project_path;
    progress_bar = ip.progress_bar;
    create_graphics = ip.create_graphics;
    save_results = ip.save_results; 
    
    % recover parameters specific to linear regression (interface 2)
    var_type = ip.var_type;
    iterations = ip.var_iterations;
    burnin = ip.var_burnin;
    model_credibility = ip.var_model_credibility;
    constant = ip.var_constant;
    trend = ip.var_trend;
    quadratic_trend = ip.var_quadratic_trend;
    lags = ip.lags;
    ar_coefficients = ip.ar_coefficients;
    pi1 = ip.pi1;
    pi2 = ip.pi2;
    pi3 = ip.pi3;
    pi4 = ip.pi4;
    pi5 = ip.pi5;
    pi6 = ip.pi6;
    pi7 = ip.pi7;
    proxy_variables = ip.proxy_variables;
    lamda = ip.lamda;
    proxy_prior = ip.proxy_prior;
    insample_fit = ip.var_insample_fit;
    constrained_coefficients = ip.constrained_coefficients;
    sums_of_coefficients = ip.sums_of_coefficients;
    initial_observation = ip.initial_observation;
    long_run = ip.long_run;
    stationary = ip.stationary;
    marginal_likelihood = ip.var_marginal_likelihood;
    hyperparameter_optimization = ip.var_hyperparameter_optimization;
    coefficients_file = ip.coefficients_file;
    long_run_file  = ip.long_run_file;
    
    % recover application parameters (interface 3)
    forecast = ip.forecast;
    forecast_credibility = ip.forecast_credibility;
    conditional_forecast = ip.conditional_forecast;
    conditional_forecast_credibility = ip.conditional_forecast_credibility;
    irf = ip.irf;
    irf_credibility = ip.irf_credibility;
    fevd = ip.fevd;
    fevd_credibility = ip.fevd_credibility;
    hd = ip.hd;
    hd_credibility = ip.hd_credibility;
    forecast_periods = ip.forecast_periods;
    conditional_forecast_type = ip.conditional_forecast_type;
    forecast_file = ip.forecast_file;
    conditional_forecast_file = ip.conditional_forecast_file;
    forecast_evaluation = ip.forecast_evaluation;
    irf_periods = ip.irf_periods;
    structural_identification = ip.structural_identification;
    structural_identification_file = ip.structural_identification_file;
    
    % recover remaining parameters
    endogenous = ip.var_endogenous;
    exogenous = ip.var_exogenous;
    proxys = ip.proxys;
    dates = ip.var_dates;
    forecast_dates = ip.var_forecast_dates;
    Y_p = ip.var_Y_p;
    Z_p = ip.var_Z_p;
    constrained_coefficients_table = ip.constrained_coefficients_table;
    long_run_table = ip.long_run_table;
    condition_table = ip.condition_table;
    shock_table = ip.shock_table;
    restriction_table = ip.restriction_table;

    % initialize timer
    ip.input_timer('start');

        
    %---------------------------------------------------
    % Model creation
    %---------------------------------------------------  
    
    
    % maximum likelihood VAR
    if var_type == 1
        var = MaximumLikelihoodVar(endogenous, 'exogenous', exogenous, 'structural_identification', structural_identification,...
            'lags', lags, 'constant', constant, 'trend', trend, 'quadratic_trend', quadratic_trend, ...
            'credibility_level', model_credibility, 'verbose' , progress_bar);  
        
    % Bayesian VAR with Minnesota prior
    elseif var_type == 2
        var = MinnesotaBayesianVar(endogenous, 'exogenous', exogenous, 'structural_identification', structural_identification,...
            'restriction_table', restriction_table, 'lags', lags, 'constant', constant, 'trend', trend, 'quadratic_trend', quadratic_trend, ...
            'ar_coefficients', ar_coefficients, 'pi1', pi1, 'pi2', pi2, 'pi3', pi3, 'pi4', pi4, ...
            'pi5', pi5, 'pi6', pi6, 'pi7', pi7, 'constrained_coefficients', constrained_coefficients, ...
            'constrained_coefficients_table', constrained_coefficients_table, 'sums_of_coefficients', sums_of_coefficients, ...
            'dummy_initial_observation', initial_observation, 'long_run_prior', long_run, 'long_run_table', long_run_table, ...
            'hyperparameter_optimization', hyperparameter_optimization, 'stationary_prior', stationary, ...
            'credibility_level', model_credibility, 'iterations', iterations, 'verbose' , progress_bar);  
        
    % Bayesian VAR with Normal-Wishart prior
    elseif var_type == 3
        var = NormalWishartBayesianVar(endogenous, 'exogenous', exogenous, 'structural_identification', structural_identification,...
            'restriction_table', restriction_table, 'lags', lags, 'constant', constant, 'trend', trend, 'quadratic_trend', quadratic_trend, ...
            'ar_coefficients', ar_coefficients, 'pi1', pi1, 'pi3', pi3, 'pi4', pi4, 'pi5', pi5, 'pi6', pi6, 'pi7', pi7, ...
             'sums_of_coefficients', sums_of_coefficients, 'dummy_initial_observation', initial_observation, ...
            'long_run_prior', long_run, 'long_run_table', long_run_table, ...
            'hyperparameter_optimization', hyperparameter_optimization, 'stationary_prior', stationary, ...
            'credibility_level', model_credibility, 'iterations', iterations, 'verbose' , progress_bar);          
       
    % Bayesian VAR with independent prior
    elseif var_type == 4
        var = IndependentBayesianVar(endogenous, 'exogenous', exogenous, 'structural_identification', structural_identification,...
            'restriction_table', restriction_table, 'lags', lags, 'constant', constant, 'trend', trend, 'quadratic_trend', quadratic_trend, ...
            'ar_coefficients', ar_coefficients, 'pi1', pi1, 'pi2', pi2, 'pi3', pi3, 'pi4', pi4, ...
            'pi5', pi5, 'pi6', pi6, 'pi7', pi7, 'constrained_coefficients', constrained_coefficients, ...
            'constrained_coefficients_table', constrained_coefficients_table, 'sums_of_coefficients', sums_of_coefficients, ...
            'dummy_initial_observation', initial_observation, 'long_run_prior', long_run, 'long_run_table', long_run_table, ...
            'stationary_prior', stationary, 'credibility_level', model_credibility, ...
            'iterations', iterations, 'burnin', burnin, 'verbose' , progress_bar);   
        
    % Bayesian VAR with dummy observation prior
    elseif var_type == 5
        var = DummyObservationBayesianVar(endogenous, 'exogenous', exogenous, 'structural_identification', structural_identification,...
            'restriction_table', restriction_table, 'lags', lags, 'constant', constant, 'trend', trend, 'quadratic_trend', quadratic_trend, ...
            'ar_coefficients', ar_coefficients, 'pi1', pi1, 'pi3', pi3, 'pi4', pi4, 'pi5', pi5, 'pi6', pi6, 'pi7', pi7, ...
             'sums_of_coefficients', sums_of_coefficients, 'dummy_initial_observation', initial_observation, ...
            'long_run_prior', long_run, 'long_run_table', long_run_table, 'stationary_prior', stationary, ...
            'credibility_level', model_credibility, 'iterations', iterations, 'verbose' , progress_bar);          
        
    % Bayesian VAR with large BVAR prior
    elseif var_type == 6
        var = LargeBayesianVar(endogenous, 'exogenous', exogenous, 'structural_identification', structural_identification, ...
            'restriction_table', restriction_table, 'lags', lags, 'constant', constant, 'trend', trend, 'quadratic_trend', quadratic_trend, ...
            'ar_coefficients', ar_coefficients, 'pi1', pi1, 'pi2', pi2, 'pi3', pi3, 'pi4', pi4, ...
            'pi5', pi5, 'pi6', pi6, 'pi7', pi7, 'constrained_coefficients', constrained_coefficients, ...
            'constrained_coefficients_table', constrained_coefficients_table, 'sums_of_coefficients', sums_of_coefficients, ...
            'dummy_initial_observation', initial_observation, 'long_run_prior', long_run, 'long_run_table', long_run_table, ...
            'stationary_prior', stationary, 'credibility_level', model_credibility, ...
            'iterations', iterations, 'burnin', burnin, 'verbose' , progress_bar);         

    % Bayesian VAR with proxy SVAR prior
    elseif var_type == 7
        var = BayesianProxySvar(endogenous, proxys, 'exogenous', exogenous, 'structural_identification', structural_identification, ...
              'restriction_table', restriction_table, 'lamda', lamda, 'proxy_prior', proxy_prior, 'lags', lags, 'constant', constant, ...
              'trend', trend, 'quadratic_trend', quadratic_trend, 'ar_coefficients', ar_coefficients, 'pi1', pi1, 'pi3', pi3, 'pi4', pi4, ...
              'credibility_level', model_credibility, 'iterations', iterations, 'burnin', burnin, 'verbose', progress_bar);         
    end
    

    %---------------------------------------------------
    % Model estimation
    %---------------------------------------------------
    
    
    % model estimation
    var.estimate();            
         
    
    %---------------------------------------------------
    % Model application: in-sample fit and residuals
    %--------------------------------------------------- 
    
    
    % apply if in-sample fit and residuals is selected
    if insample_fit
        var.insample_fit();
    end

    
    %---------------------------------------------------
    % Model application: marginal likelihood
    %--------------------------------------------------- 
    
    
    % apply if marginal likelihood is selected, and model is any compatible Bayesian VAR
    if marginal_likelihood && ismember(var_type, [2 3 4])
        var.marginal_likelihood();
    end
    
    
    %---------------------------------------------------
    % Model application: forecasts
    %--------------------------------------------------- 
    
    
    % estimate forecasts, if selected
    if forecast
        var.forecast(forecast_periods, forecast_credibility, Z_p);
        % estimate forecast evaluation, if selected
        if forecast_evaluation
            var.forecast_evaluation(Y_p);
        end
    end  


    %---------------------------------------------------
    % Model application: impulse response function
    %--------------------------------------------------- 


    % estimate impulse response function, if selected
    if irf
        var.impulse_response_function(irf_periods, irf_credibility);
    end
        

    %------------------------------------------------------------
    % Model application: forecast error variance decomposition
    %------------------------------------------------------------


    % estimate forecast error variance decompositionn, if selected
    if fevd
        var.forecast_error_variance_decomposition(irf_periods, fevd_credibility);
    end


    %------------------------------------------------------------
    % Model application: historical decomposition
    %------------------------------------------------------------


    % estimate historical decomposition, if selected
    if hd
        var.historical_decomposition(hd_credibility);
    end
    
    
    %------------------------------------------------------------
    % Model application: conditional forecasts
    %------------------------------------------------------------


    % estimate conditional forecast, if selected
    if var_type ~= 1 && conditional_forecast && conditional_forecast_type == 1
        var.conditional_forecast(forecast_periods, conditional_forecast_credibility, condition_table, Z_p);
    elseif var_type ~= 1 && conditional_forecast && conditional_forecast_type == 2
        var.structural_conditional_forecast(forecast_periods, ...
        conditional_forecast_credibility, condition_table, shock_table, Z_p);
    end


    %---------------------------------------------------
    % Model processor: prepare elements for results
    %---------------------------------------------------               
          
    
    % print estimation completion
    cu.print_completion_message(progress_bar);

    % end estimation timer
    ip.input_timer('end');

    % make information dictionary for result class
    ip.make_results_information();
    results_information = ip.results_information;

    % make information dictionary for graphics class
    ip.make_graphics_information();
    graphics_information = ip.graphics_information;


    %---------------------------------------------------
    % Model results: create, display and save
    %---------------------------------------------------               
            

    % recover path to result folder
    results_path = fullfile(project_path, 'results');

    % initialize results class
    res = Results(var, 'complementary_information', results_information);

    % create and save input summary if relevant
    if save_results
        res.make_input_summary();
        res.save_input_summary(results_path);
    end

    % create, show and save estimation summary if relevant
    res.make_estimation_summary();
    res.show_estimation_summary();
    if save_results
        res.save_estimation_summary(results_path);
    end

    % create and save application summary if relevant
    if save_results
        res.make_application_summary();
        res.save_application_summary(results_path);
    end


    %---------------------------------------------------
    % Model graphics: generate and save
    %---------------------------------------------------


    % if graphic selection is selected
    if create_graphics
        
        % recover path to result folder
        graphics_path = fullfile(project_path, 'graphics');

        % initialize graphics class
        grp = Graphics(var, 'complementary_information', graphics_information, ...
              'path', graphics_path, 'clear_folder', true);

        % run graphics for all applications inturn
        grp.insample_fit_graphics(false, true);
        grp.forecast_graphics(false, true);
        grp.conditional_forecast_graphics(false, true);
        grp.irf_graphics(false, true);
        grp.fevd_graphics(false, true);
        grp.hd_graphics(false, true);

        % display graphics
        gui = GraphicalUserInterface(view_graphics = true);  

    end

end