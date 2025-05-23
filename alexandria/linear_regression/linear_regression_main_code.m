function [lr] = linear_regression_main_code(user_inputs)
    

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
    regression_type = ip.regression_type;
    iterations = ip.iterations;
    burnin = ip.burnin;
    model_credibility = ip.model_credibility;
    b = ip.b;
    V = ip.V;
    alpha = ip.alpha;
    delta = ip.delta;
    g = ip.g;
    Q = ip.Q;
    tau = ip.tau;
    thinning = ip.thinning;
    thinning_frequency = ip.thinning_frequency;
    q = ip.q;
    p = ip.p;
    H = ip.H;
    constant = ip.constant;
    b_constant = ip.b_constant;
    V_constant = ip.V_constant;
    trend = ip.trend;
    b_trend = ip.b_trend;
    V_trend = ip.V_trend;
    quadratic_trend = ip.quadratic_trend;
    b_quadratic_trend = ip.b_quadratic_trend;
    V_quadratic_trend = ip.V_quadratic_trend;
    insample_fit = ip.insample_fit;
    marginal_likelihood = ip.marginal_likelihood;
    hyperparameter_optimization = ip.hyperparameter_optimization;
    optimization_type = ip.optimization_type;

    % recover application parameters (interface 3)
    forecast = ip.forecast;
    forecast_credibility = ip.forecast_credibility;
    forecast_evaluation = ip.forecast_evaluation;

    % recover remaining parameters
    endogenous = ip.endogenous;
    exogenous = ip.exogenous;
    Z = ip.Z;
    y_p = ip.y_p;
    X_p = ip.X_p;
    Z_p = ip.Z_p;

    % initialize timer
    ip.input_timer('start');


    %---------------------------------------------------
    % Model creation
    %---------------------------------------------------     
    
    
    % maximum likelihood regression
    if regression_type == 1
        lr = MaximumLikelihoodRegression(endogenous, exogenous, 'constant', constant, 'trend' , trend, ...
             'quadratic_trend' , quadratic_trend, 'credibility_level' , model_credibility, 'verbose' , progress_bar);
         
    % simple Bayesian regression
    elseif regression_type == 2
        lr = SimpleBayesianRegression(endogenous, exogenous, 'constant' , constant, 'trend' , trend, ...
             'quadratic_trend' , quadratic_trend, 'b_exogenous' , b, 'V_exogenous' , V, 'b_constant' , b_constant, ...
             'V_constant' , V_constant, 'b_trend' , b_trend, 'V_trend' , V_trend, 'b_quadratic_trend' , b_quadratic_trend, ...
             'V_quadratic_trend' , V_quadratic_trend, 'hyperparameter_optimization', hyperparameter_optimization,...
             'optimization_type', optimization_type, 'credibility_level' , model_credibility, 'verbose' , progress_bar);
         
    % hierarchical Bayesian regression
    elseif regression_type == 3
        lr = HierarchicalBayesianRegression(endogenous, exogenous, 'constant' , constant, 'trend' , trend, ...
             'quadratic_trend' , quadratic_trend, 'b_exogenous' , b, 'V_exogenous' , V, 'b_constant' , b_constant, ...
             'V_constant' , V_constant, 'b_trend' , b_trend, 'V_trend' , V_trend, 'b_quadratic_trend' , b_quadratic_trend, ...
             'V_quadratic_trend' , V_quadratic_trend, 'alpha' , alpha, 'delta' , delta, ...
             'hyperparameter_optimization', hyperparameter_optimization, 'optimization_type', optimization_type, ...
             'credibility_level', model_credibility, 'verbose' , progress_bar);
         
    % independent Bayesian regression
    elseif regression_type == 4
        lr = IndependentBayesianRegression(endogenous, exogenous, 'constant' , constant, 'trend' , trend, ...
             'quadratic_trend' , quadratic_trend, 'b_exogenous' , b, 'V_exogenous' , V, 'b_constant' , b_constant, ...
             'V_constant' , V_constant, 'b_trend' , b_trend, 'V_trend' , V_trend, 'b_quadratic_trend' , b_quadratic_trend, ...
             'V_quadratic_trend' , V_quadratic_trend, 'alpha' , alpha, 'delta' , delta, 'iterations' , iterations, ...
             'burn' , burnin, 'credibility_level' , model_credibility, 'verbose' , progress_bar);
         
    % heteroscedastic Bayesian regression
    elseif regression_type == 5
        lr = HeteroscedasticBayesianRegression(endogenous, exogenous, 'heteroscedastic' , Z, ...
            'constant' , constant, 'trend' , trend, 'quadratic_trend' , quadratic_trend, 'b_exogenous' , b, ...
            'V_exogenous' , V, 'b_constant' , b_constant, 'V_constant' , V_constant, 'b_trend' , b_trend, ...
            'V_trend' , V_trend, 'b_quadratic_trend' , b_quadratic_trend, 'V_quadratic_trend' , V_quadratic_trend, ...
            'alpha' , alpha, 'delta' , delta, 'g' , g, 'Q' , Q, 'tau' , tau, 'iterations' , iterations, ...
            'burn' , burnin, 'thinning' , thinning, 'thinning_frequency' , thinning_frequency, ...
            'credibility_level' , model_credibility, 'verbose' , progress_bar);
        
    % autocorrelated Bayesian regression
    elseif regression_type == 6
        lr = AutocorrelatedBayesianRegression(endogenous, exogenous, 'q' , q, 'constant' , constant, ...
            'trend' , trend, 'quadratic_trend' , quadratic_trend, 'b_exogenous' , b, 'V_exogenous' , V, ...
            'b_constant' , b_constant, 'V_constant' , V_constant, 'b_trend' , b_trend, ...
            'V_trend' , V_trend, 'b_quadratic_trend' , b_quadratic_trend, 'V_quadratic_trend' , V_quadratic_trend, ...
            'alpha' , alpha, 'delta' , delta, 'p' , p, 'H' , H, 'iterations' , iterations, ...
            'burn' , burnin, 'credibility_level' , model_credibility, 'verbose' , progress_bar);
    end
            

    %---------------------------------------------------
    % Model estimation
    %---------------------------------------------------
    
    
    % model estimation
    lr.estimate();
    

    %---------------------------------------------------
    % Model application: in-sample fit and residuals
    %--------------------------------------------------- 
    
    
    % apply if in-sample fit and residuals is selected
    if insample_fit
        lr.insample_fit();
    end
    
    
    %---------------------------------------------------
    % Model application: marginal likelihood
    %---------------------------------------------------         
    
    
    % apply if marginal likelihood is selected, and model is any Bayesian regression
    if marginal_likelihood && regression_type ~= 1
        lr.marginal_likelihood();
    end
    
    
    %---------------------------------------------------
    % Model application: forecasts
    %---------------------------------------------------          
        
        
    % estimate forecasts, if selected
    if forecast
        if regression_type == 5
            lr.forecast(X_p, forecast_credibility, 'Z_hat', Z_p);
        else
            lr.forecast(X_p, forecast_credibility);
        end
        % estimate forecast evaluation, if selected
        if forecast_evaluation
            lr.forecast_evaluation(y_p);
        end
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
    % res = Results(lr);
    res = Results(lr, 'complementary_information', results_information);

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
        grp = Graphics(lr, 'complementary_information', graphics_information, ...
                       'path', graphics_path, 'clear_folder', true);

        % run graphics for all applications inturn
        grp.insample_fit_graphics(false, true);
        grp.forecast_graphics(false, true);

        % display graphics
        gui = GraphicalUserInterface('view_graphics', true);
        
    end
    

end


