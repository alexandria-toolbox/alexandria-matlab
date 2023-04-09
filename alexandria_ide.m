% clear workspace and console (not to be modified)
clear all;
clc;
warning off;
close all;
addpath(genpath(pwd));

% initiate user inputs (not to be modified)
user_inputs = struct;
user_inputs.tab_1 = struct;
user_inputs.tab_2_lr = struct;
user_inputs.tab_3 = struct;


%---------------------------------------------------
% Editable part: tab 1
%--------------------------------------------------- 


% model choice (1 = linear regression)
user_inputs.tab_1.model = 1;

% endogenous variables, as string array (e.g. ["var1" "var2"])
user_inputs.tab_1.endogenous_variables = [""];

% # exogenous variables, as string array (e.g. ["var1" "var2"]; leave as empty array [] if no exogenous)
user_inputs.tab_1.exogenous_variables = [""];

% data frequency (1: cross-sectional/undated, 2: yearly, 3: quarterly, 4: monthly, 5: weekly, 6: daily)
user_inputs.tab_1.frequency = 1;

% data sample: start and end periods, as string array of timestamps (e.g. ["1990-03-31" "2020-12-31"]) or periods (e.g. ["1990Q1" "2020Q4"])
user_inputs.tab_1.sample = ["" ""];

% path to project folder, as char (e.g. 'D:\my_project')
user_inputs.tab_1.project_path = pwd();

% name of data file, as char (e.g. 'data.csv')
user_inputs.tab_1.data_file = '';

% display progress bar during estimation (true: yes, false: no)
user_inputs.tab_1.progress_bar = true;

% generate estimation graphics (true: yes, false: no)
user_inputs.tab_1.create_graphics = true;

% save estimation results (true: yes, false: no)
user_inputs.tab_1.save_results = true;


%---------------------------------------------------
% Editable part: tab 2, linear regression
%--------------------------------------------------- 


% this applies only if the selected model is linear regression (model = 1)
if user_inputs.tab_1.model == 1
    
    % choice of linear regression model (1: maximum likelihood; 2: simple Bayesian;
    % 3: hierarchical; 4: independent; 5: heteroscedastic; 6: autocorrelated)
    user_inputs.tab_2_lr.regression_type = 1;

    % post-burn iterations for MCMC algorithm (integer)
    user_inputs.tab_2_lr.iterations = 2000;

    % burnin iterations for MCMC algorithm (integer)
    user_inputs.tab_2_lr.burnin = 1000;
    
    % credibility level for model estimates (float between 0 and 1)
    user_inputs.tab_2_lr.model_credibility = 0.95;

    % prior mean for regression coefficients beta: either scalar for common mean (e.g. 0),
    % or column vector, one value for each coefficient (e.g. [0 0 0]')
    user_inputs.tab_2_lr.b = 0;

    % prior variance for regression coefficients beta: either scalar for common variance (e.g. 1),
    % or column vector, one value for each coefficient (e.g. [1 1 1]')
    user_inputs.tab_2_lr.V = 1;
    
    % prior shape for regression variance sigma (positive float)
    user_inputs.tab_2_lr.alpha = 0.0001;

    % prior scale for regression variance sigma (positive float)
    user_inputs.tab_2_lr.delta = 0.0001;

    % prior mean for heteroscedastic coefficients gamma: either scalar for common mean (e.g. 0),
    % or column vector, one value for each coefficient (e.g. [0 0 0]')
    user_inputs.tab_2_lr.g = 0;

    % prior variance for heteroscedastic coefficients gamma: either scalar for common variance (e.g. 1),
    % or column vector, one value for each coefficient (e.g. [1 1 1]')
    user_inputs.tab_2_lr.Q = 100;
    
    % variance of transition kernel for Metropolis-Hastings step in heteroscedastic model (positive float)
    user_inputs.tab_2_lr.tau = 0.001;

    % apply posterior thinning to MCMC draws in heteroscedastic model (true: yes, false: no)
    user_inputs.tab_2_lr.thinning = false;

    % frequency of posterior thinning (positive integer)
    user_inputs.tab_2_lr.thinning_frequency = 10;

    % Z variables, as string array (e.g. ["var1" "var2"]); can be empty if model is not heteroscedastic regression
    user_inputs.tab_2_lr.Z_variables = [""];
    
    % order of autoregressive process for residuals in autocorrelated models (positive integer)
    user_inputs.tab_2_lr.q = 1;
    
    % prior mean for autocorrelation coefficients phi: either scalar for common mean (e.g. 0),
    % or column vector, one value for each coefficient (e.g. [0 0 0]')
    user_inputs.tab_2_lr.p = 0;

    % prior variance for autocorrelation coefficients phi: either scalar for common variance (e.g. 1),
    % or column vector, one value for each coefficient (e.g. [1 1 1]')
    user_inputs.tab_2_lr.H = 100;
    
    % include constant in regression (true: yes, false: no)
    user_inputs.tab_2_lr.constant = true;   
    
    % prior mean for regression constant (float)
    user_inputs.tab_2_lr.b_constant = 0;

    % prior variance for constant (positive float)
    user_inputs.tab_2_lr.V_constant = 1;  
    
    % include trend in regression (true: yes, false: no)
    user_inputs.tab_2_lr.trend = false;    
    
    % prior mean for regression trend (float)
    user_inputs.tab_2_lr.b_trend = 0;

    % prior variance for trend (positive float)
    user_inputs.tab_2_lr.V_trend = 1;      
    
    % include quadratic trend in regression (true: yes, false: no)
    user_inputs.tab_2_lr.quadratic_trend = false;   
    
    % prior mean for regression quadratic trend (float)
    user_inputs.tab_2_lr.b_quadratic_trend = 0;

    % prior variance for quadratic trend (positive float)
    user_inputs.tab_2_lr.V_quadratic_trend = 1;    
    
    % estimate in-sample fit (true: yes, false: no)
    user_inputs.tab_2_lr.insample_fit = false;
    
    % estimate marginal likelihood (true: yes, false: no)
    user_inputs.tab_2_lr.marginal_likelihood = false;   
    
    % apply hyperparameter optimization (true: yes, false: no)
    user_inputs.tab_2_lr.hyperparameter_optimization = false;   
    
    % type of hyperparameter optimization (1: common variance, 2: coefficient-specific variances plus residual variance)
    user_inputs.tab_2_lr.optimization_type = 1;
end

    
%---------------------------------------------------
% Editable part: tab 3
%--------------------------------------------------- 

    
% estimate forecasts for the model (true: yes, false: no)
user_inputs.tab_3.forecast = false;

% credibility level for forecast estimates (float between 0 and 1)
user_inputs.tab_3.forecast_credibility = 0.95;

% estimate conditional forecasts for the model (true: yes, false: no)
user_inputs.tab_3.conditional_forecast = false;    

% credibility level for conditional forecast estimates (float between 0 and 1)
user_inputs.tab_3.conditional_forecast_credibility = 0.95;

% estimate impulse response functions for the model (true: yes, false: no)
user_inputs.tab_3.irf = false;

% credibility level for impulse response functions estimates (float between 0 and 1)
user_inputs.tab_3.irf_credibility = 0.95;

% estimate forecast error variance decomposition for the model (true: yes, false: no)
user_inputs.tab_3.fevd = false;

% credibility level for forecast error variance decomposition estimates (float between 0 and 1)
user_inputs.tab_3.fevd_credibility = 0.95;

% estimate historical decomposition for the model (true: yes, false: no)
user_inputs.tab_3.hd = false;

% credibility level for historical decomposition estimates (float between 0 and 1)
user_inputs.tab_3.hd_credibility = 0.95;

% number of forecast periods (positive integer)
user_inputs.tab_3.forecast_periods = 1;

% number of impulse response functions periods (positive integer)
user_inputs.tab_3.irf_periods = 1;

% type of conditional forecasts (1: all shocks, 2: shock-specific)
user_inputs.tab_3.conditional_forecast_type = 1;

% structural identification scheme (1: none, 2: Cholesky)
user_inputs.tab_3.structural_identification = 1;

% estimate forecast evaluation criteria (true: yes, false: no)
user_inputs.tab_3.forecast_evaluation = false;

% name of forecast data file, as char (e.g. 'data_forecast.csv')
user_inputs.tab_3.forecast_file = '';

% name of structural identification file, as char (e.g. 'structural_identification.csv')
user_inputs.tab_3.structural_identification_file = '';


%---------------------------------------------------
% Main code (not to be modified)
%--------------------------------------------------- 


model = user_inputs.tab_1.model;

% if model is linear regression, run main code for linear regression, and return model
if model == 1
    lr = linear_regression_main_code(user_inputs);
end

