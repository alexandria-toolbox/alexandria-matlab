classdef Results < handle
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------    
    
    
    properties (GetAccess = public, SetAccess = protected)
        settings
    end    

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------

    
    methods (Access = public) 

        
        function self = Results()
        end 
        
                
        function create_result_file(self, project_path, save_results)
            % if save_results is activated
            if save_results
                % create path to results folder
                results_folder_path = fullfile(project_path, 'results');
                results_file_path = fullfile(results_folder_path, 'results.txt');
                % if results folder already exists, delete it 
                if exist(results_folder_path) == 7
                    rmdir(results_folder_path, 's');
                end
                % create results folder
                mkdir(project_path, 'results');
                % get Alexandria header
                header = cu.alexandria_header();
                % write header in results file
                cu.write_string_array(header, results_file_path);
            end
        end       
        
        
    end
    

    methods (Access = protected, Hidden = true)
        
   
        function print_alexandria_header(self)
            % get Alexandria header
            header = cu.alexandria_header();
            % print header
            cu.print_string_array(header);     
        end    

        
        function print_start_message(self)
            cu.print_message('Starting estimation of your model...');
            cu.print_message(' ');
        end     

        
        function print_completion_message(self)
            if self.progress_bar
                cu.print_message(' ');
            end
            cu.print_message('Estimation completed successfully.');
            cu.print_message(' ');
            cu.print_message(' ');       
        end
        
        
        function print_and_save_summary(self)
            % display result summary on console
            cu.print_string_array(self.summary);
            % if save_results is activated, save summary in file
            if self.save_results
                results_file_path = fullfile(self.project_path, 'results', 'results.txt');
                cu.write_string_array(self.summary, results_file_path);
            end
        end

        
        function add_settings_header(self)
            % Alexandria header
            lines = cu.alexandria_header();
            lines = [lines; 'Estimation date:  ' datestr(self.estimation_start, 'yyyy-mm-dd hh:MM:ss')];
            lines = [lines; ' '];
            lines = [lines; ' '];
            self.settings = [self.settings;lines];
        end  
                    
        
        function add_tab_1_settings(self)
            % recover tab 1 elements
            endogenous = convertStringsToChars(strjoin(self.endogenous,', '));
            exogenous = convertStringsToChars(strjoin(self.exogenous,', '));
            if self.frequency == 1
                frequency = 'cross-sectional/undated';
            elseif self.frequency == 2
                frequency = 'annual';           
            elseif self.frequency == 3
                frequency = 'quarterly';  
            elseif self.frequency == 4
                frequency = 'monthly';   
            elseif self.frequency == 5
                frequency = 'weekly';   
            elseif self.frequency == 6
                frequency = 'daily';
            end
            sample = [convertStringsToChars(self.start_date) ' ' ...
                      convertStringsToChars(self.end_date)];
            project_path = self.project_path;
            data_file = self.data_file;
            progress_bar = self.progress_bar;
            create_graphics = self.create_graphics;
            save_results = self.save_results;
            lines = cu.tab_1_settings('linear regression', endogenous, exogenous, ...
                                      frequency, sample, project_path, data_file, ...
                                      progress_bar, create_graphics, save_results);
            self.settings = [self.settings;lines];
        end
        
        function save_settings(self)
            % if save_results is activated, save settings in file
            if self.save_results
                settings_file_path = fullfile(self.project_path, 'results', 'settings.txt');
                cu.write_string_array(self.settings, settings_file_path);
            end
        end
                

        

    end
    
    
end