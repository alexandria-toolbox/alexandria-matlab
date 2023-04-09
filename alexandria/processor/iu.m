classdef iu
    

    % iu stands for input utilities
    % A class containing static methods to handle user inputs such as string, dates, lists, and so on

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------
    

    methods(Static)
        
        
        function [fixed_char] = fix_char(chr)
    
            % [fixed_char] = fix_char(chr)
            % converts tabs and multiple spaces to single spaces, removes leading and final spaces
            %
            % parameters:
            % char : char
            %     char to fix
            %
            % returns:
            % fixed_char : char
            %     char with spaces made regular

            % convert tabs into spaces
            chr = regexprep(chr, '\t', ' ');
            % convert multiple spaces into single spaces
            chr = regexprep(chr,' +',' ');
            % trim initial and final spaces
            fixed_char = strtrim(chr); 
        end
        
        
        function [string_array] = char_to_array(chr)
            
            % [string_array] = char_to_array(chr)
            % converts char into array of string, with elements split at spaces
            %
            % parameters: 
            % strng: str or array of str
            %     string to split and distribute in array
            %
            % returns:
            % string_array: array of str
            %     array obtained after splitting original char         
            
            % if item is already a string array, ignore and return            
            if isstring(chr)
                string_array = chr;
            % else, item is a char: fix it in case it is incorrectly formated
            else
                chr = iu.fix_char(chr);
                % split it at spaces to convert to string array
                string_array = split(convertCharsToStrings(chr))';
            end
        end
        
        
        function [chr] = array_to_char(string_array)
            
            % [chr] = array_to_char(string_array)
            % converts list of string into char, with elements separated with commas
            %
            % parameters: 
            % string_array: array of str
            %     string array to burst into single char
            %
            % returns:
            % chr: char
            %     char obtained by bursting string array       
            
            % if array is actually a scalar, just convert to string   
            if isa(string_array, 'double') && numel(string_array) == 1
                chr = num2str(string_array);
            elseif isa(string_array, 'double') && numel(string_array) > 1
                chr = regexprep(num2str(string_array'),'\s+',', ');
            % else, item is a string array: convert to single char
            else
                chr = convertStringsToChars(strjoin(string_array,', '));
            end
        end        
        
        
        function [bool] = is_empty(x)
            
            % [bool] = is_empty(x)
            % checks whether input is empty array or empty char/string
            %
            % parameters: 
            % x : array, char or string
            %     elements to be tested for being empty
            %
            % returns:
            % bool: bool
            %     1 if x is empty, 0 otherwise
            
            bool = isempty(x) || (ischar(x) && isempty(x)) || (isstring(x) ...
                   && length(x) == 1 && strlength(x) == 0); 
        end
        
        
        function [bool] = is_digit(x)
            
            % [bool] = is_digit(x)
            % checks whether input is char representing integer value
            % 
            % parameters: 
            % x : char
            %     elements to be tested for being integer value char
            %
            % returns:
            % bool: bool
            %     1 if x is char representing integer, 0 otherwise
            
            if ischar(x) && numel(regexp(x, '\d')) == numel(x)
                bool = true;
            else
                bool = false;
            end
        end
        
        
        function [bool] = is_integer(x)
            
            % [bool] = is_integer(x)
            % checks whether input is double representing integer value
            %
            % parameters:
            % x : double
            %     element to be tested for integer value
            %
            % returns:
            % bool: bool
            %     1 if x is integer, 0 otherwise
            
            if isnumeric(x) && numel(x) == 1 && ~isnan(x) && isfinite(x) && x == floor(x)
                bool = true;
            else
                bool = false;
            end
        end
        
        
        function check_file_path(path, file)

            % check_file_path(path, file)
            % checks whether file exists at given path, and is of correct format (csv, xls, xlsx)
            % 
            % parameters: 
            % path : char
            %     path to folder containing data file
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % 
            % returns:
            % none

            % check if path is valid
            if ~isfolder(path)
                error(['Path error. Specified path ' path ' does not exist.']);
            end
            % check if file exists
            file_path = fullfile(path, file);
            if ~isfile(file_path)
                error(['File error. File ' file_path ' could not be found.']);
            end
            % check if data file is of correct format (excel or csv)
            file_type = strsplit(file, '.');
            file_type = file_type{end};
            if ~ismember(file_type, ["csv" "xls" "xlsx"])
                error(['Type error for file ' file '. File must be of type csv, xls or xlsx.']);
            end
        end
        
        
        function [data] = load_data(path, file)
            
            % [data] = load_data(path, file)
            % loads data file of type csv, xls or xlsx into Matlab table
            % 
            % parameters: 
            % path : char
            %     path to folder containing data file
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % 
            % returns:
            % data : table
            %     table containing loaded data
            
            % make sure path and file strings are properly formatted
            path = iu.fix_char(path);
            file = iu.fix_char(file);
            % prepare options for file loading
            file_path = fullfile(path, file);
            opts = detectImportOptions(file_path);
            opts.VariableTypes{1} = 'char';
            opts.VariableTypes([2:end]) = {'double'};
            % load file
            data = readtable(file_path, opts, 'ReadRowNames',true);
            data = removevars(data,{'Var1'});
            % string(data.Properties.RowNames)
        end
        
        
        function check_variables(data, file, variables, tag)

            % check_variables(data, file, variables, tag)
            % checks whether specified variables exist in data table
            %
            % parameters: 
            % data : table
            %     table to examine for data
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % variables : string array
            %     array of endogenous variables to check in table
            % tag : char
            %     tag to apply in error message, if any (e.g. 'Endogenous variables')
            %
            % returns:
            %     none
            
            % obtain the array of variables in table
            table_variables = string(data.Properties.VariableNames);
            % check that variables exist in this array
            missing_variables = setdiff(variables, table_variables);
            if ~iu.is_empty(missing_variables)
                error(['Data error for file ' file '. ' tag ' ' char(strjoin(missing_variables, ', ')) ' cannot be found.']);
            end
        end
        
        
        function check_dates(data, file, start_date, end_date)
            
            % check_dates(data, file, start_date, end_date)
            % check whether start and end dates can be found in index of dataframe
            %
            % parameters:
            % data : table
            %     table to examine for start and end dates
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % start_date : char
            %     sample start date to search in dataframe index
            % end_date : char
            %     sample end date to search in dataframe index
            %
            % returns:
            % none
            
            % first obtain the list of dates (as strings)
            dates = string(data.Properties.RowNames);
            % check for start date
            if ~ismember(start_date, dates)
                error(['Date error for file ' file '. Start date ' convertStringsToChars(start_date) ' cannot be found.']);
            end
            if ~ismember(end_date, dates)
                error(['Date error for file ' file '. End date ' convertStringsToChars(end_date) ' cannot be found.']);
            end
        end

        
        function [sample] = fetch_data(data, file, start_date, end_date, variables, tag)
            
            % [sample] = fetch_data(data, file, start_date, end_date, variables, tag)
            % fetches variables from data at given sample dates
            %
            % parameters:
            % data : table
            %     dataframe containing loaded data
            % file : char
            %     name of data file (with extension csv, xls or xlsx)
            % start_date : char
            %     sample start date to search in dataframe index
            % end_date : char
            %     sample end date to search in dataframe index    
            % variables : string array
            %     array of variables to search in table
            % tag : char
            %     tag to apply in error message, if any (e.g. 'Endogenous variables')
            %
            % returns:
            % sample : matrix
            %     matrix containing fetched data
            
            % recover sample dates
            dates = string(data.Properties.RowNames);
            start_date_index = find(all(ismember(dates, start_date), 2));
            end_date_index = find(all(ismember(dates, end_date), 2));
            dates = dates(start_date_index:end_date_index);
            % if variable array is not empty, recover sample for given variables and dates
            if ~iu.is_empty(variables)
                sample = data(dates, variables);
                % test for NaNs, and if any, raise error
                if any(any(isnan(table2array(sample))))
                    error(['Data error for file ' file '. ' tag ' contains NaN entries, which are unhandled.']);
                % else, data is valid: convert to matrix
                else
                    sample = table2array(sample);
                end
            % if no variable, return empty array
            else
                sample = [];
            end
        end        
        
        
        function [date_format] = infer_date_format(frequency, file, start_date, end_date)
            
            % [date_format] = infer_date_format(frequency, file, start_date, end_date)
            % infer date format for given data file, which can be either periods or timestamps
            %
            % parameters:
            % frequency : int
            %     data frequency, as int between 1 and 6
            % file : char
            %     name of data file (with extension csv, xls or xlsx)            
            % start_date : char
            %     sample start date to search in dataframe index
            % end_date : char
            %     sample end date to search in dataframe index
            %
            % returns:
            % date_format : char
            %     date format, either 'periods' or 'timestamps'
            
            if frequency == 1
                date_format = iu.infer_undated(file, start_date, end_date);
            elseif frequency == 2
                date_format = iu.infer_annual(file, start_date, end_date);
            elseif frequency == 3
                date_format = iu.infer_quarterly(file, start_date, end_date);                
            elseif frequency == 4
                date_format = iu.infer_monthly(file, start_date, end_date);                  
            elseif frequency == 5
                date_format = iu.infer_weekly(file, start_date, end_date);  
            elseif frequency == 6
                date_format = iu.infer_daily(file, start_date, end_date);                 
            end
        end
            
            
        function [date_format] = infer_undated(file, start_date, end_date)
            if ~iu.is_digit(start_date) || ~iu.is_digit(end_date)
                error(['Date error for file ' file '. Unrecognized format for cross-sectional/undated sample start and end dates. Should be integers.']);
            else
                date_format = 'periods';
            end
        end
            
 
        function [date_format] = infer_annual(file, start_date, end_date)
            if iu.is_digit(start_date) && iu.is_digit(end_date)
                date_format = 'periods';
            elseif ~iu.is_empty(regexp(start_date, '\d{4}[-]\d{1,2}[-]\d{1,2}')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}[-]\d{1,2}[-]\d{1,2}'))
                date_format = 'timestamps';
            else
                error(['Date error for file ' file '. Unrecognized format for annual sample start and end dates. Should be integers (e.g. 1990) or timestamp (e.g. 1990-12-31).']);
            end
        end

  
        function [date_format] = infer_quarterly(file, start_date, end_date)
            if ~iu.is_empty(regexp(start_date, '\d{4}Q[1-4]')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}Q[1-4]'))
                date_format = 'periods';
            elseif ~iu.is_empty(regexp(start_date, '\d{4}[-]\d{1,2}[-]\d{1,2}')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}[-]\d{1,2}[-]\d{1,2}'))
                date_format = 'timestamps';
            else
                error(['Date error for file ' file '. Unrecognized format for quarterly sample start and end dates. Should be period (e.g. 1990Q1) or timestamp (e.g. 1990-03-31).']);
            end
        end       
        
 
        function [date_format] = infer_monthly(file, start_date, end_date)
            if ~iu.is_empty(regexp(start_date, '\d{4}M(0?[1-9]|[1][0-2])$')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}M(0?[1-9]|[1][0-2])$'))
                date_format = 'periods';
            elseif ~iu.is_empty(regexp(start_date, '\d{4}[-]\d{1,2}[-]\d{1,2}')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}[-]\d{1,2}[-]\d{1,2}'))
                date_format = 'timestamps';
            else
                error(['Date error for file ' file '. Unrecognized format for monthly sample start and end dates. Should be period (e.g. 1990M1) or timestamp (e.g. 1990-01-31).']);
            end
        end  

        
        function [date_format] = infer_weekly(file, start_date, end_date)
            if ~iu.is_empty(regexp(start_date, '\d{4}W(0?[1-9]|[1-4][0-9]|[5][0-3])$')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}W(0?[1-9]|[1-4][0-9]|[5][0-3])$'))
                date_format = 'periods';
            elseif ~iu.is_empty(regexp(start_date, '\d{4}[-]\d{1,2}[-]\d{1,2}')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}[-]\d{1,2}[-]\d{1,2}'))
                date_format = 'timestamps';
            else
                error(['Date error for file ' file '. Unrecognized format for weekly sample start and end dates. Should be period (e.g. 1990W1) or timestamp (e.g. 1990-01-05).']);
            end
        end  
        
        
        function [date_format] = infer_daily(file, start_date, end_date)
            if ~iu.is_empty(regexp(start_date, '\d{4}D(0?0?[1-9]|0?[1-9]\d|[1-2]\d{2}|3[0-5][0-9]|36[0-6])$')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}D(0?0?[1-9]|0?[1-9]\d|[1-2]\d{2}|3[0-5][0-9]|36[0-6])$'))
                date_format = 'periods';
            elseif ~iu.is_empty(regexp(start_date, '\d{4}[-]\d{1,2}[-]\d{1,2}')) ...
                    && ~iu.is_empty(regexp(end_date, '\d{4}[-]\d{1,2}[-]\d{1,2}'))
                date_format = 'timestamps';
            else
                error(['Date error for file ' file '. Unrecognized format for daily sample start and end dates. Should be period (e.g. 1990D10) or timestamp (e.g. 1990-01-10).']);
            end
        end
        
        
        function [dates] = generate_dates(data, date_format, frequency, file, start_date, end_date)
            
            % [dates] = generate_dates(data, date_format, frequency, file, start_date, end_date)
            % generates date series, under the form of a date array
            %
            % parameters:
            % data : table
            %     table to examine for start and end dates
            % date_format : char
            %     date format, either 'periods' or 'timestamps'        
            % frequency : int
            %     data frequency, as int between 1 and 6
            % file : char
            %     name of data file (with extension csv, xls or xlsx)            
            % start_date : char
            %     sample start date to search in dataframe index
            % end_date : char
            %     sample end date to search in dataframe index
            %
            % returns:
            % dates : date array
            %     array of datetime entries
            
            % recover sample dates
            dates = string(data.Properties.RowNames);
            start_date_index = find(all(ismember(dates, start_date), 2));
            end_date_index = find(all(ismember(dates, end_date), 2));
            dates = dates(start_date_index:end_date_index);
            % convert to dates, depending on frequency
            if frequency == 1
                dates = iu.generate_undated_dates(dates, date_format, file);
            elseif frequency == 2
                dates = iu.generate_annual_dates(dates, date_format, file);
            elseif frequency == 3
                dates = iu.generate_quarterly_dates(dates, date_format, file);  
            elseif frequency == 4
                dates = iu.generate_monthly_dates(dates, date_format, file);
            elseif frequency == 5
                dates = iu.generate_weekly_dates(dates, date_format, file);                 
            elseif frequency == 6
                dates = iu.generate_daily_dates(dates, date_format, file);                 
            end
        end
        
        
        function [dates] = generate_undated_dates(dates, date_format, file)
            if strcmp(date_format, 'periods')
                try
                    dates = str2num(char(dates));
                catch
                    error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
                end
            else
                error(['Date error for file ' file '. Date format should be periods, not timestamps.']);
            end
        end
                    

        function [dates] = generate_annual_dates(dates, date_format, file)
            try
                if strcmp(date_format, 'periods')
                    dates = datetime(str2num(char(dates))+ 1, 1, 0);
                elseif strcmp(date_format, 'timestamps')
                    dates = dateshift(datetime(dates),'end','year');
                end  
            catch
                error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
            end
        end
        
        
        function [dates] = generate_quarterly_dates(dates, date_format, file)
            try
                if strcmp(date_format, 'periods')
                    dates = char(dates);
                    years = str2num(dates(:,1:4));
                    months = str2num(dates(:,6:end)) * 3 + 1;
                    dates = datetime(years, months, 0);
                elseif strcmp(date_format, 'timestamps')
                    dates = dateshift(datetime(dates),'end','quarter');
                end  
            catch
                error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
            end
        end        
        
        
        function [dates] = generate_monthly_dates(dates, date_format, file)
            try
                if strcmp(date_format, 'periods')
                    dates = char(dates);
                    years = str2num(dates(:,1:4));
                    months = str2num(dates(:,6:end)) + 1;
                    dates = datetime(years, months, 0);
                elseif strcmp(date_format, 'timestamps')
                    dates = dateshift(datetime(dates),'end','month');
                end  
            catch
                error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
            end
        end         
        

        function [dates] = generate_weekly_dates(dates, date_format, file)
            try
                if strcmp(date_format, 'periods')
                    dates = char(dates);
                    years = str2num(dates(:,1:4));
                    days = str2num(dates(:,6:end)) * 7;
                    dates = dateshift(datetime(years, 1, days),'dayofweek','Friday');
                elseif strcmp(date_format, 'timestamps')
                    dates = dateshift(datetime(dates),'dayofweek','Friday');
                end  
            catch
                error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
            end
        end 

        
        function [dates] = generate_daily_dates(dates, date_format, file)
            try
                if strcmp(date_format, 'periods')
                    dates = char(dates);
                    years = str2num(dates(:,1:4));
                    days = str2num(dates(:,6:end));
                    dates = datetime(years, 1, days);
                elseif strcmp(date_format, 'timestamps')
                    dates = datetime(dates);
                end  
            catch
                error(['Date error for file ' file '. Some sample periods seem to be improperly formatted. Please verify the dates in the data file.'])
            end
        end
        
        
        function [sample_p] = fetch_forecast_data(data, insample_data, ...
                              variables, file, required, periods, tag)
                          
            % [sample_p] = fetch_forecast_data(data, insample_data, ...
            %              variables, file, required, periods, tag)
            % fetches predictor variables from forecast data file
            %
            % parameters:
            % data : table
            %     table containing loaded data
            % insample_data : matrix
            %     matrix containing in-sample counterparts of forecast variables        
            % variables : string array
            %     array of variables to search in table
            % exogenous_variables : string array
            %     array of variables to extract from table  
            % file : char
            %     name of data file (with extension csv, xls or xlsx) 
            % required : bool
            %     if true, checks that variables are provided in file
            % periods : int
            %     number of forecast periods
            %
            % returns:
            % sample_p : matrix
            %     matrix containing forecast data
            
            
            % find any missing variable
            data_variables = string(data.Properties.VariableNames);
            missing_endogenous = setdiff(variables, data_variables);            
            % if some variables are missing
            if ~iu.is_empty(missing_endogenous)            
                % if variables are required and no in-sample data is provided, raise error
                if required && iu.is_empty(insample_data)     
                    error(['Data error for file ' file '. ' tag ' ' char(strjoin(missing_variables, ', ')) ' required to conduct forecast (or forecast evaluation), but cannot be found.'])
                % if variables are required but some in-sample data is provided, replicate final in-sample value
                elseif required && ~iu.is_empty(insample_data)
                    sample_p = repmat(insample_data(end,:), periods, 1);
                % if not required, return empty array
                else
                    sample_p = [];
                end
            % if no variable is missing, recover prediction data 
            else
                sample_p = table2array(data(:,variables));
                % if too few periods of data are provided, raise error
                if size(data, 1) < periods
                    error(['Data error for file ' file '. Forecasts must be conducted for ' str2num(periods) ' periods, but ' tag ' is provided for fewer periods.'])
                % else, reduce data to number of periods
                else
                    sample_p = sample_p(1:periods,:);
                end
                % test for NaNs, and if any, raise error
                if any(isnan(sample_p(:)))
                    error(['Data error for file ' file '. ' tag ' contains NaN entries, which are unhandled.'])
                end
            end
        end
                
        
        function [forecast_dates] = generate_forecast_dates(end_date, periods, frequency)
            
            % [forecast_dates] = generate_forecast_dates(end_date, periods, frequency)
            % generates date series for forecasts, under the form of a pandas index, for a given number of periods
            %
            % parameters:
            % end_date : char
            %     sample end date to search in dataframe index            
            % periods : int
            %     number of forecast periods            
            % frequency : int
            %     data frequency, as int between 1 and 6
            %
            % returns:
            % forecast_dates : date array
            %     array of datetime entries
            
            % if frequency is undated, create a range from 1 to periods as forecast dates
            if frequency == 1
                forecast_dates = (1:periods)';
            % if frequency is annual, expand sample by years equal to periods
            elseif frequency == 2
                forecast_dates = (dateshift(end_date + calyears(1:periods),'end','year'))';
            % if frequency is quarterly, expand sample by quarters equal to periods            
            elseif frequency == 3
                forecast_dates = (dateshift(end_date + calquarters(1:periods),'end','quarter'))';
            % if frequency is monthly, expand sample by months equal to periods
            elseif frequency == 4
                forecast_dates = (dateshift(end_date + calmonths(1:periods),'end','month'))';
            % if frequency is weekly, expand sample by weeks equal to periods
            elseif frequency == 5
                forecast_dates = (dateshift(end_date + calweeks(1:periods),'dayofweek','Friday'))';            
            % if frequency is daily, expand sample by days equal to periods
            elseif frequency == 6
                forecast_dates = (end_date + caldays(1:periods))';
            end
        end
        

    end
        
        
end
    
    
    