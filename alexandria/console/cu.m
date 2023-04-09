classdef cu
    

    % cu stands for console utilities
    % A class containing static methods for console display

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
        
        
        function progress_bar(iteration, total_iterations, application_tag)
            
            % progress_bar(iteration, total_iterations, application_tag)
            % progress bar for MCMC algorithms
            %
            % parameters:
            % iteration: int
            %   current iteration of MCMC algorithm
            % total_iterations: int
            %   total number of iterations of MCMC algorithm 
            % application_tag: str
            %   string flagging the application being estimated by MCMC
            %
            % returns:
            % none

            if iteration == 1
                % print application tag (which application is being run)
                fprintf('%s\n', application_tag);
                % print iteration 1 progress bar
                string = [repmat(' ', 1, numel(num2str(total_iterations))) ...
                     '1/' num2str(total_iterations) '  [' repmat('.', 1, 33) ...
                     ']   0%'];
                 fprintf(string);
            elseif iteration == total_iterations
                % length of string to overwrite with new string
                string_length = 43 + 2 * numel(num2str(total_iterations));
                % print final iteration progress bar
                string = [repmat('\b', 1, string_length) ...
                    num2str(total_iterations) '/' num2str(total_iterations) ...
                    '  [' repmat('=',1,33) ']  —  done'];
                fprintf(string);
                fprintf('\n');
            else
                % integer defining the arrow position in the progress bar
                arrow_position = floor(iteration * 33.3333 / total_iterations);
                % integer percentage of progress
                percentage = floor(iteration * 100 / total_iterations);
                string_percentage = [repmat(' ', 1, double(percentage < 10)) ...
                    num2str(percentage) '%'];
                % length of string to overwrite with new string
                string_length = 43 + 2 * numel(num2str(total_iterations));
                string = [repmat('\b', 1, string_length) ...
                    repmat(' ', 1, numel(num2str(total_iterations)) ...
                    - numel(num2str(iteration))) num2str(iteration) '/' ...
                    num2str(total_iterations) '  [' ...
                    repmat('=', 1, arrow_position - 1) ...
                    repmat('>', 1, double(arrow_position ~= 0)) ...
                    repmat('.', 1, 33 - arrow_position) ']  ' ...
                    string_percentage '%'];
                fprintf(string);
                pause(1e-08);
            end
        end
        
        
        
        function progress_bar_complete(application_tag)
    
            % progress_bar_complete(application_tag)
            % pseudo progress bar used to show completion when no actual MCMC is run
            %
            % parameters:
            % application_tag: str
            %   string flagging the application being estimated by MCMC
            %
            % returns:
            % none

            % print application tag (which application is being run)
            fprintf('%s\n',application_tag);   
            % print pseudo progress bar
            string = ['  — / —    [' repmat('=',1,33) ']  —  done'];
            fprintf('%s\n', string);
        end
        
    
        function optimization_completion(success)
        
            % optimization_completion(success)
            % binary message that indicates success or failure of optimization
            % 
            % parameters:      
            % success: bool
            %     boolean indicating whether optimization has succeeded
            %     
            % returns:
            % none

            if success
                fprintf('Optimization conducted succesfully.\n');
            else
                fprintf('Warning! Optimization failed. Prior values may be unreliable.\n')
            end
        end
        
        
        function print_message(message)
            
            % print_message(message)
            % a simple print function that prints any message
            %
            % parameters:       
            % message: str
            %     the message to print
            %
            % returns:
            % none   

            disp(message);
        end
        
        
        function print_string_array(string_array)
        
            % print_string_list(string_list)
            % a function that prints an array of strings, line by line
            % 
            % parameters:        
            % string_array: str array
            %     array containing the strings to print
            % 
            % returns
            % none
            
            for i=1:numel(string_array)
                disp(string_array(i));
            end
        end
        

        function write_string_array(string_array, filepath)

            % print_string_list(string_list)
            % a function that writes a list of strings, line by line, in a specified path
            % 
            % parameters      
            % string_list: str list
            %     list containing the strings to print
            % filepath: str
            %     full path to file (must include file name)
            % 
            % returns
            % none
            
            file = fopen(filepath, 'a');
            for i=1:numel(string_array)
                fprintf(file, '%s\n', string_array(i));
            end
            fclose(file);
        end
            
            
        function [header] = alexandria_header()

            % [header] = alexandria_header()
            % a function that creates a list of strings constituting Alexandria's header
            % 
            % parameters      
            % none
            % 
            % returns
            % header: str array
            %     the list of strings constituting Alexandria's header
     
            header = [...
                      "                                                                               "; ...
                      "      ======================================================================== "; ...
                      "         _____   __                                   _                        "; ...
                      "        / __  / / / ___  __  __   ____    ____    ___/ / ____  (_) ____        "; ...
                      "       / /_/ / / / / _ \ \ \/ /  / __ `/ / __ \ / __  / / __/ / / / __ `/      "; ...
                      "      / __  / / / /  __/  |  |  / /_/ / / / / // /_/ / / /   / / / /_/ /       "; ...
                      "     /_/ /_/ /_/  \___/  /_/\_\ \__,_/ /_/ /_/ \____/ /_/   /_/  \__,_/        "; ...
                      "                                                                               "; ...
                      "     The library of Bayesian time-series models                                "; ...
                      "     V 0.1 - Copyright Ⓒ  Romain Legrand                                      "; ...
                      "   ========================================================================    "; ...
                      "                                                                               "; ...
                      "                                                                               "; ...
                      ];
        end
        
        
        function [formatted_number] = format_number(number)
            
            % [formatted_number] = format_number(number)
            % formats any number into a 10-character string
            %
            % parameters       
            % number: float
            %     number to be formatted
            %
            % returns
            % formatted_number: str
            %     string containing the formatted number
            
            if 0.0001 < abs(number) && abs(number) < 100000
                formatted_number = sprintf('%10.3f', number);
            else
                formatted_number = sprintf('%10.3e', number);
            end
        end
        
        
        function [formatted_name] = format_name(name)
            
            % [formatted_name] = format_name(name)
            % formats any variable name as a string containing a maximum of 20 characters
            %
            % parameters    
            % name: str
            %     variable name to be formatted
            %
            % returns
            % formatted_name: str
            %     string containing the formatted name

            if numel(name) <= 20
                formatted_name =  sprintf('%-20s', name);
            else
                formatted_name =  [name(1:16) '...'];
            end
        end
        
        
        function [string] = shorten_string(string, n)
            
            % [string] = shorten_string(string, n)
            % allows for a maximum of n characters in string; if longer, string is shortened with ' ...'
            %
            % parameters       
            % string: char
            %     string to be formatted
            %
            % returns
            % string: char
            %     string containing at most n characters

            if numel(string) > n
                string = [string(1:n-4) ' ...'];
            end
        end
        
        
        function [line] = equal_dashed_line()

            % [line] = equal_dashed_line()
            % return a line of 80 equal signs ('=')
            %
            % parameters     
            % none
            %
            % returns
            % line: char
            %     string containing the line
            
            line = convertCharsToStrings(repmat('=', 1, 80));
        end
        
        
        function [line] = hyphen_dashed_line()

            % [line] = hyphen_dashed_line()
            % return a line of 80 hyphen signs ('-')
            %
            % parameters     
            % none
            %
            % returns
            % line: char
            %     string containing the line
            
            line = convertCharsToStrings(repmat('-', 1, 80));
        end        
        
        
        function [header] = model_header(model)

            % [header] = model_header(model)
            % return a results header with model name and two wrapping equal dashed lines
            %
            % parameters     
            % model: char
            %     char containing the model name
            %
            % returns
            % header: str array
            %     array of string containing the model header

            % initiate header, add equal dashed line
            header = [cu.equal_dashed_line()];
            % center-justify model name
            header = [header; cu.center_align(model, 80)];
            % add final equal dashed line
            header = [header; cu.equal_dashed_line()];  
        end
        
        
        function [centered_string] = center_align(string, n)
            
            % [centered_string] = center_align(string)
            % returns a n-character char with the string center-aligned
            %
            % parameters:
            % string: char
            %     the char to be center-aligned
            %
            % returns
            % centered_string: char
            %     n-character char with the string center-aligned
            
            string_length = numel(string);
            left_pad = repmat(' ', 1, floor((n - string_length) / 2));
            right_pad = repmat(' ', 1, ceil((n - string_length) / 2));
            centered_string = [left_pad string right_pad];
        end
           
        
        function [header] = estimation_header(start, complete, n, endogenous, frequency, sample)

            % [header] = estimation_header(start, complete, n, endogenous, frequency, sample)
            % return an estimation header with basic model information
            %
            % parameters
            % start: datetime object
            %     datetime corresponding to date where estimation starts
            % complete: datetime object
            %     datetime corresponding to date where estimation is complete
            % n: int
            %     number of sample observations
            % endogenous: str
            %     explained variable
            % frequency: char
            %     data frequency
            % sample: str
            %     estimation sample, start and end dates
            %
            % returns
            % header: str array
            %     array of string containing the estimation header
            
            % first row: dependent variable and data frequency
            left_element = sprintf(['Dep. variable:' sprintf('%+24s', cu.shorten_string(endogenous, 20))]);
            right_element = sprintf(['Sample:' sprintf('%+31s',sample)]);
            header = [left_element  '    ' right_element];
            % second row: estimation start and sample
            left_element = sprintf(['Est. start:' sprintf('%+27s', datestr(start, 'yyyy-mm-dd hh:MM:ss'))]);
            right_element = sprintf(['Frequency:' sprintf('%+28s', frequency)]);
            header = [header; [left_element '    ' right_element]];     
            % third row: estimation complete and observations
            left_element = sprintf(['Est. complete:' sprintf('%+24s', datestr(complete, 'yyyy-mm-dd hh:MM:ss'))]);
            right_element = sprintf(['No. observations:' sprintf('%+21s', num2str(n))]);
            header = [header; [left_element  '    ' right_element]];
        end
        

        function [header] = coefficient_header(credibility_level)

            % [header] = coefficient_header(credibility_level)
            % return an estimation line with headers for coefficient, standard devations and credibility bounds
            %
            % parameters
            % credibility_level: float
            %     credibility level for model
            %
            % returns
            % line: char
            %     string containing coefficient header

            % calculate lower and upper bounds
            lower_bound = sprintf('%5.3f', 0.5 - credibility_level / 2);
            upper_bound = sprintf('%5.3f', 0.5 + credibility_level / 2);
            header = [repmat(' ', 1, 29) 'median' repmat(' ', 1, 8) ...
                      'std dev' repmat(' ', 1, 9) '[' lower_bound ...
                      repmat(' ', 1, 9) upper_bound ']'];
            % generate second line of header
            header = [header; cu.hyphen_dashed_line()];
        end
            
        
        function [line] = parameter_estimate_line(name, coefficient, standard_deviation, lower_bound, upper_bound)

            % [line] = parameter_estimate_line(name, coefficient, standard_deviation, lower_bound, upper_bound)
            % return an estimation line with formatted values for coefficient, standard devations and credibility bounds
            %
            % parameters
            % name: char
            %     name of variable to which parameter is related
            % coefficient: float
            %     coefficient value
            % standard_deviation: float (positive)
            %     standard deviation of parameter
            % lower_bound: float
            %     lower bound of parameter
            % upper_bound: float
            %     upper bound of parameter
            %
            % returns
            % line: char
            %     string containing coefficient summary

            line = [cu.format_name(name) repmat(' ', 1, 5) cu.format_number(coefficient) ...
                    repmat(' ', 1, 5) cu.format_number(standard_deviation) repmat(' ', 1, 5) ...
                    cu.format_number(lower_bound) repmat(' ', 1, 5) cu.format_number(upper_bound)];
        end
        
        
        function [line] = string_line(string)
            
            % string_line(name)
            % returns a line where name is left-aligned, and right padding is added to reach 80 characters
            %
            % parameters
            % string: char
            %     string to left-align
            %
            % returns
            % line: str
            %     80-character line with left-aligned string

            % left-justify the name, with pad to reach 80 characters
            line = sprintf('%-80s', string);
        end
     
            
        function [lines] = insample_evaluation_lines(ssr, r2, adj_r2, m_y, aic, bic)

            % [lines] = insample_evaluation_lines(ssr, r2, adj_r2, m_y, aic, bic)
            % returns the set of lines with the results for in-sample evaluation criteria
            %
            % parameters
            % ssr: float
            %     sum of squared residuals
            % r2: float
            %     coefficient of determination
            % adj_r2: float
            %     adjusted coefficient of determination
            % m_y: float
            %     log10 marginal likelihood
            % aic: float
            %     Akaike information criterion
            % bic: float
            %     Bayesian information criterion
            %
            % returns
            % lines: str array
            %     set of lines reporting in-sample evaluation criteria

            % first row: ssr and marginal likelihood
            left_element = sprintf(['ssr:' sprintf('%+34s', cu.format_number(ssr))]);
            if m_y
                right_element = sprintf(['log10 marg. lik.:' sprintf('%+21s', cu.format_number(m_y))]);
            else
                right_element = repmat(' ', 1, 38);
            end
            lines = [left_element  '    ' right_element]; 
            % second row: r2 and adjusted r2
            left_element = sprintf(['R2:' sprintf('%+35s', cu.format_number(r2))]);
            right_element = sprintf(['adj. R2:' sprintf('%+30s', cu.format_number(adj_r2))]);
            lines = [lines; [left_element  '    ' right_element]]; 
            % third row AIC and BIC
            if aic
                left_element = sprintf(['AIC:' sprintf('%+34s', cu.format_number(aic))]);
                right_element = sprintf(['BIC:' sprintf('%+34s', cu.format_number(bic))]);
            end
        end
            

        function [lines] = forecast_evaluation_lines(rmse, mae, mape, theil_u, bias, log_score, crps)

            % [lines] = forecast_evaluation_lines(rmse, mae, mape, theil_u, bias, log_score, crps)
            % returns the set of lines with the results for out-of-sample evaluation criteria
            %
            % parameters
            % rmse: float
            %     root mean squared error for predictions
            % mae: float
            %     mean absolute error for predictions
            % mape: float
            %     mean absolute percentage error for predictions
            % theil_u: float
            %     Theil's U coefficient for predictions
            % bias: float
            %     bias for predictions
            % log_score: float
            %     log score for predictions
            % crps: float
            %     continuous rank probability score for predictions
            %
            % returns
            % lines: str array
            %     set of lines reporting out-of-sample evaluation criteria

            % first row: rmse
            left_element = sprintf(['rmse:' sprintf('%+33s', cu.format_number(rmse))]);
            right_element = repmat(' ', 1, 38);
            lines = [left_element '    ' right_element];
            % second row: mae and mape
            left_element = sprintf(['mae:' sprintf('%+34s', cu.format_number(mae))]);
            right_element = sprintf(['mape:' sprintf('%+33s', cu.format_number(mape))]);
            lines = [lines; [left_element '    ' right_element]];
            % third row: Theil's U and bias
            left_element = sprintf(['Theil''s U:' sprintf('%+28s', cu.format_number(theil_u))]);
            right_element = sprintf(['bias:' sprintf('%+33s', cu.format_number(bias))]);
            lines = [lines; [left_element  '    ' right_element]]; 
            % fourth row: log score and crps
            if log_score
                left_element = sprintf(['log score:' sprintf('%+28s', cu.format_number(log_score))]);
                right_element = sprintf(['crps:' sprintf('%+33s', cu.format_number(crps))]);
                lines = [lines; [left_element  '    ' right_element]];
            end
        end
        
        
        function [lines] = tab_1_settings(model, endogenous, exogenous, frequency, ...
                           sample, path, file, progress_bar, create_graphics, save_results)

            % [lines] = tab_1_settings(model, endogenous, exogenous, frequency, ...
            %           sample, path, file, progress_bar, create_graphics, save_results)
            % returns the set of lines with the results for tab 1 settings
            %
            % parameters  
            % model: char
            %     selected model (e.g. 'linear regression')
            % endogenous: char
            %     set of endogenous variables, separated by a comma
            % exogenous: char
            %     set of exogenous variables, separated by a comma
            % frequency: char
            %     sample data frequency
            % sample: char
            %     sample dates, separated by a space
            % path: char
            %     path to project folder
            % file: char
            %     name of data file
            % progress_bar: bool
            %     user's choice for progress bar
            % create graphics: bool
            %     user's choice for graphics creation
            % save_results: bool
            %     user's choice for saving results
            %
            % returns
            % lines: string array
            %     set of lines reporting settings for tab 1

            % initiate lines
            lines = [];
            % header for tab 1
            lines = [lines;"Models"];
            lines = [lines;"---------"];
            lines = [lines;" "];
            % other elements for tab 1
            lines = [lines;['model selection: ' model]];
            lines = [lines;['endogenous variables: ' endogenous]];
            lines = [lines;['exogenous variables: ' exogenous]];            
            lines = [lines;['data frequency: ' frequency]];              
            lines = [lines;['estimation sample: ' sample]];              
            lines = [lines;['path to project folder: ' path]];              
            lines = [lines;['data file: ' file]];  
            lines = [lines;['progress bar: ' cu.bool_to_string(progress_bar)]];              
            lines = [lines;['create graphics: ' cu.bool_to_string(create_graphics)]];                 
            lines = [lines;['save_results: ' cu.bool_to_string(save_results)]]; 
            lines = [lines;" "];
            lines = [lines;" "];
        end     
        
        
        function [string] = bool_to_string(bool)

            % [string] = bool_to_string(bool)
            % converts boolean to 'yes' or 'no'
            %
            % parameters  
            % bool: bool
            %     boolean to convert
            %
            % returns
            % string: char
            %     'yes' or 'no', depending on convered boolean  

            if bool
                string = 'yes';
            else
                string = 'no';
            end
        end
        
        
    end
    
    
end
        
        
