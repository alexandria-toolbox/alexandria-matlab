classdef GraphicalUserInterface < handle & DefaultInputInterface & Tab1Interface ...
                                  & Tab2RegressionInterface ...
                                  & Tab3Interface & Tab4Interface & Tab5Interface
    

    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess = public)
        % property indicating that the 'Run' button has not been pressed yet
        press_run = false
        % property indicating that user has closed interface manually
        user_interrupt = false
        % main structure, storing user inputs
        user_inputs
        % properties for lazy creation of tabs 2
        created_tab_2_lr = false
        % interface and tabs properties
        background_color
        backtabs_color
        interface
        current_tab
        tab_pbt1
        tab_pbt2
        tab_pbt3
        tab_pbt4
        tab_pbt5
        interface_path
        % property to decide if interface is used to navigate graphics
        view_graphics
    end
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------     
    
    
    methods (Access = public)
        

        function self = GraphicalUserInterface(varargin)
            % save graphics view as attribute
            parser = inputParser;
            default_view_graphics = false;
            addParameter(parser, 'view_graphics', default_view_graphics);
            parse(parser, varargin{:});
            self.view_graphics = parser.Results.view_graphics;            
            % initiate structure of user inputs
            self.initiate_inputs();
            % generate overall interface and tabs
            self.create_interface_and_tabs();
            % create elements of tab 1
            self.create_tab_1();
            % create elements of tab 2
            self.create_tab_2();
            % create elements of tab 3
            self.create_tab_3();
            % create elements of tab 4
            self.create_tab_4();
            % create elements of tab 5
            self.create_tab_5();
            % check if GUI is switched back to tab 4
            self.switch_tab_4();
        end


        function initiate_inputs(self)
            % get path to interface folder
            self.interface_path = fileparts(which(mfilename));
            % if previous user inputs have been saved, load them
            if exist(fullfile(self.interface_path, 'user_inputs.mat'), 'file') == 2
                load(fullfile(self.interface_path, 'user_inputs.mat'));
                self.user_inputs = user_inputs;
            % otherwise, implement default inputs
            else
                self.create_default_inputs();      
            end
        end
        
        
        function create_interface_and_tabs(self)
            % define colors, depending on operating system
            if ispc
                self.background_color = [1 1 0.77];
                self.backtabs_color = [0.95 0.87 0.55];                
            else    
                self.background_color = [1 0.97 0.88];
                self.backtabs_color = [0.85 0.77 0.60];
            end
            % get screen resolution
            set(0, 'units', 'pixels');
            screen_dimensions = get(0, 'ScreenSize');
            screen_width = screen_dimensions(3);
            screen_heigth = screen_dimensions(4);
            % calculate interface position and create main window
            interface_width = 1000;
            interface_heigth = 650;
            left_shift = (screen_width - interface_width) / 2;
            up_shift = (screen_heigth - interface_heigth) / 2;        
            self.interface = figure('CloseRequestFcn', @self.cb_user_interrupt);
            set(self.interface, 'units', 'pixels', 'position', ...
                   [left_shift up_shift interface_width interface_heigth]);
            set(self.interface, 'name', 'Alexandria', 'NumberTitle', 'off');
            set(self.interface, 'ToolBar', 'none', 'MenuBar', 'none');
            set(self.interface, 'Color', self.background_color);
            % first tab button            
            self.tab_pbt1 = uicontrol('style', 'pushbutton');
            set(self.tab_pbt1, 'unit', 'pixels', 'position', [1 610 200 40]);
            set(self.tab_pbt1, 'String', 'Models');
            set(self.tab_pbt1, 'FontSize', 16);
            set(self.tab_pbt1, 'FontName', 'ZapfDingbats');
            set(self.tab_pbt1, 'BackgroundColor', self.background_color);
            set(self.tab_pbt1, 'CallBack', @self.cb_tab_pbt1);
            % second tab button            
            self.tab_pbt2 = uicontrol('style', 'pushbutton');
            set(self.tab_pbt2, 'unit', 'pixels', 'position', [202 610 200 40]);
            set(self.tab_pbt2, 'String', 'Specifications');
            set(self.tab_pbt2, 'FontSize', 16);
            set(self.tab_pbt2, 'FontName', 'ZapfDingbats');
            set(self.tab_pbt2, 'BackgroundColor', self.backtabs_color);
            set(self.tab_pbt2, 'CallBack', @self.cb_tab_pbt2);
            % third tab button            
            self.tab_pbt3 = uicontrol('style', 'pushbutton');
            set(self.tab_pbt3, 'unit', 'pixels', 'position', [403 610 200 40]);
            set(self.tab_pbt3, 'String', 'Applications');
            set(self.tab_pbt3, 'FontSize', 16);
            set(self.tab_pbt3, 'FontName', 'ZapfDingbats');
            set(self.tab_pbt3, 'BackgroundColor', self.backtabs_color);
            set(self.tab_pbt3, 'CallBack', @self.cb_tab_pbt3);
            % fourth tab button            
            self.tab_pbt4 = uicontrol('style', 'pushbutton');
            set(self.tab_pbt4, 'unit', 'pixels', 'position', [604 610 200 40]);
            set(self.tab_pbt4, 'String', 'Graphics');
            set(self.tab_pbt4, 'FontSize', 16);
            set(self.tab_pbt4, 'FontName', 'ZapfDingbats');
            set(self.tab_pbt4, 'BackgroundColor', self.backtabs_color);
            set(self.tab_pbt4, 'CallBack', @self.cb_tab_pbt4);
            % fifth tab button            
            self.tab_pbt5 = uicontrol('style', 'pushbutton');
            set(self.tab_pbt5, 'unit', 'pixels', 'position', [805 610 200 40]);
            set(self.tab_pbt5, 'String', 'Credits');
            set(self.tab_pbt5, 'FontSize', 16);
            set(self.tab_pbt5, 'FontName', 'ZapfDingbats');
            set(self.tab_pbt5, 'BackgroundColor', self.backtabs_color);
            set(self.tab_pbt5, 'CallBack', @self.cb_tab_pbt5);
        end


        function cb_tab_pbt1(self, hObject, callbackdata)
            % hide current tab, show tab 1
            self.hide_current_tab();
            self.show_tab_1();
            % set current tab as tab 1, update tab button color
            self.current_tab = 'tab_1';
            set(self.tab_pbt1, 'BackgroundColor', self.background_color);
        end  
        
        
        function cb_tab_pbt2(self, hObject, callbackdata)
            % hide current tab
            self.hide_current_tab();
            % tab 2 is created in a lazy fashion: create the tab only when it is called
            % if tab2 is called for linear regression:
            if self.user_inputs.tab_1.model == 1
                % if tab 2 for linear regression does not exist, create it
                if self.created_tab_2_lr == false
                    self.create_tab_2_lr();
                end
                % show tab 2 for linear regression
                self.show_tab_2_lr();
            end
            % set current tab as tab 2, linear regression
            self.current_tab = 'tab_2_lr';
            % update tab button color
            set(self.tab_pbt2, 'BackgroundColor', self.background_color);
        end         
        
        
        function cb_tab_pbt3(self, hObject, callbackdata)
            % hide current tab, show tab 3
            self.hide_current_tab();
            self.show_tab_3();
            % set current tab as tab 3, update tab button color
            self.current_tab = 'tab_3';
            set(self.tab_pbt3, 'BackgroundColor', self.background_color);
        end          

        
        function cb_tab_pbt4(self, hObject, callbackdata)
            % hide current tab, show tab 4
            self.hide_current_tab();
            self.show_tab_4();
            % set current tab as tab 3, update tab button color
            self.current_tab = 'tab_4';
            set(self.tab_pbt4, 'BackgroundColor', self.background_color);
        end         
        
        
        function cb_tab_pbt5(self, hObject, callbackdata)
            % hide current tab, show tab 5
            self.hide_current_tab();
            self.show_tab_5();
            % set current tab as tab 5, update tab button color
            self.current_tab = 'tab_5';
            set(self.tab_pbt5, 'BackgroundColor', self.background_color);
        end   

        
        function create_tab_2(self)
            if self.user_inputs.tab_1.model == 1
                self.create_tab_2_lr();
            end
        end
        
        
        function hide_current_tab(self)
            % if current tab is tab 1, hide it
            if strcmp(self.current_tab, 'tab_1')
                self.hide_tab_1();
            % if current tab is tab 2 for regression, hide it
            elseif strcmp(self.current_tab, 'tab_2_lr')
                self.hide_tab_2_lr();                
            % if current tab is tab 3, hide it
            elseif strcmp(self.current_tab, 'tab_3')
                self.hide_tab_3();
            % if current tab is tab 4, hide it
            elseif strcmp(self.current_tab, 'tab_4')
                self.hide_tab_4(); 
            % if current tab is tab 5, hide it
            elseif strcmp(self.current_tab, 'tab_5')
                self.hide_tab_5();
            end
        end
        
        
        function switch_tab_4(self)
            % if view graphics is True, move back to tab 4
            if self.view_graphics
                self.cb_tab_pbt4();
            end
        end
        

        function validate_interface(self)
            % save user inputs to drive
            user_inputs = self.user_inputs;
            save(fullfile(self.interface_path, 'user_inputs.mat'), 'user_inputs');
            % close interface and set press_run to true; this terminates waitfor
            delete(self.interface);
            self.press_run = true;
        end


    end        
    

end
    
    
    

