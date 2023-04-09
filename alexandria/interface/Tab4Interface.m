classdef Tab4Interface < handle


    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess = public)
        % tab 4 properties
        t4_txt1
        t4_frm1
        t4_txt2
        t4_sld1
        t4_mnu1
        t4_txt3
        t4_frm2
        t4_bgr1
        t4_rdb1
        t4_rdb2
        t4_txt4
        t4_frm3
        t4_txt5
        t4_sld2
        t4_mnu2
        t4_txt6
        t4_frm4
        t4_txt7
        t4_sld3
        t4_mnu3
        t4_img1
        graphics_folder_path
        split_image_string
        rollmenus
        applications        
        current_application
        joint_plot
        current_variables
        current_variable
        current_responses
        current_response
        image_name
    end
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)
        

        function self = Tab4Interface()
        end
        
        
        function create_tab_4(self)

            % application label
            self.t4_txt1 = uicontrol('style', 'text');
            set(self.t4_txt1, 'unit', 'pixels', 'position', [30 550 250 30]);
            set(self.t4_txt1, 'String', ' Application');
            set(self.t4_txt1, 'HorizontalAlignment', 'left');
            set(self.t4_txt1, 'FontName', 'Serif');
            set(self.t4_txt1, 'FontSize', 14);
            set(self.t4_txt1, 'FontWeight', 'bold');
            set(self.t4_txt1, 'FontAngle', 'italic');
            set(self.t4_txt1, 'BackgroundColor', self.background_color);
            set(self.t4_txt1, 'Visible', 'off');

            % frame around application
            self.t4_frm1 = uicontrol('style','frame');
            set(self.t4_frm1, 'unit', 'pixels', 'position', [20 445 280 110]);
            set(self.t4_frm1, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t4_frm1, 'BackgroundColor', self.background_color);
            set(self.t4_frm1, 'Visible', 'off');

            % application selection label
            self.t4_txt2 = uicontrol('style', 'text');
            set(self.t4_txt2, 'unit', 'pixels', 'position', [30 515 200 30]);
            set(self.t4_txt2, 'String', ' none');
            set(self.t4_txt2, 'HorizontalAlignment', 'left');
            set(self.t4_txt2, 'FontName', 'Serif');
            set(self.t4_txt2, 'FontSize', 12);
            set(self.t4_txt2, 'FontWeight', 'bold');
            set(self.t4_txt2, 'BackgroundColor', self.background_color); 
            set(self.t4_txt2, 'Visible', 'off');

            % application slider
            self.t4_sld1 = uicontrol('style','slider');
            set(self.t4_sld1, 'unit', 'pixels', 'position', [40 500 240 20]);
            set(self.t4_sld1, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t4_sld1, 'BackgroundColor', self.background_color);
            set(self.t4_sld1, 'Min', 1, 'Max', 10, 'Value', 1);
            set(self.t4_sld1, 'Visible', 'off');

            % application menu
            self.t4_mnu1 = uicontrol('style', 'popupmenu');
            set(self.t4_mnu1, 'position',[40 455 240 30]);            
            set(self.t4_mnu1, 'String', "none"); 
            set(self.t4_mnu1, 'Visible', 'off');

            % display label
            self.t4_txt3 = uicontrol('style', 'text');
            set(self.t4_txt3, 'unit', 'pixels', 'position', [30 400 250 30]);
            set(self.t4_txt3, 'String', ' Display type');
            set(self.t4_txt3, 'HorizontalAlignment', 'left');
            set(self.t4_txt3, 'FontName', 'Serif');
            set(self.t4_txt3, 'FontSize', 14);
            set(self.t4_txt3, 'FontWeight', 'bold');
            set(self.t4_txt3, 'FontAngle', 'italic');
            set(self.t4_txt3, 'BackgroundColor', self.background_color);  
            set(self.t4_txt3, 'Visible', 'off');

            % frame around display type
            self.t4_frm2 = uicontrol('style','frame');
            set(self.t4_frm2, 'unit', 'pixels', 'position', [20 330 280 70]);
            set(self.t4_frm2, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t4_frm2, 'BackgroundColor', self.background_color);
            set(self.t4_frm2, 'Visible', 'off');

            % display radiobuttons
            self.t4_bgr1 = uibuttongroup('unit','pixels', 'Position',[30 340 200 50]);
            set(self.t4_bgr1, 'BorderType', 'none');
            set(self.t4_bgr1, 'BackgroundColor', self.background_color);            
            self.t4_rdb1 = uicontrol(self.t4_bgr1,'Style','radiobutton');
            set(self.t4_rdb1, 'Position',[0 30 150 25]);
            set(self.t4_rdb1, 'String',' joint plot');
            set(self.t4_rdb1, 'FontName', 'Serif');
            set(self.t4_rdb1, 'FontSize', 12);
            set(self.t4_rdb1, 'FontWeight', 'bold');
            set(self.t4_rdb1, 'BackgroundColor', self.background_color);
            self.t4_rdb2 = uicontrol(self.t4_bgr1,'Style','radiobutton');
            set(self.t4_rdb2, 'Position',[0 0 1500 25]);
            set(self.t4_rdb2, 'String',' individual graphs');
            set(self.t4_rdb2, 'FontName', 'Serif');
            set(self.t4_rdb2, 'FontSize', 12);
            set(self.t4_rdb2, 'FontWeight', 'bold');
            set(self.t4_rdb2, 'BackgroundColor', self.background_color);
            set(self.t4_bgr1, 'Visible', 'off');

            % variable label
            self.t4_txt4 = uicontrol('style', 'text');
            set(self.t4_txt4, 'unit', 'pixels', 'position', [30 280 250 30]);
            set(self.t4_txt4, 'String', ' Variable');
            set(self.t4_txt4, 'HorizontalAlignment', 'left');
            set(self.t4_txt4, 'FontName', 'Serif');
            set(self.t4_txt4, 'FontSize', 14);
            set(self.t4_txt4, 'FontWeight', 'bold');
            set(self.t4_txt4, 'FontAngle', 'italic');
            set(self.t4_txt4, 'BackgroundColor', self.background_color); 
            set(self.t4_txt4, 'Visible', 'off');

            % Frame around variable
            self.t4_frm3 = uicontrol('style','frame');
            set(self.t4_frm3, 'unit', 'pixels', 'position', [20 175 280 110]);
            set(self.t4_frm3, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t4_frm3, 'BackgroundColor', self.background_color);
            set(self.t4_frm3, 'Visible', 'off');

            % variable label
            self.t4_txt5 = uicontrol('style', 'text');
            set(self.t4_txt5, 'unit', 'pixels', 'position', [30 245 200 30]);
            set(self.t4_txt5, 'String', ' none');
            set(self.t4_txt5, 'HorizontalAlignment', 'left');
            set(self.t4_txt5, 'FontName', 'Serif');
            set(self.t4_txt5, 'FontSize', 12);
            set(self.t4_txt5, 'FontWeight', 'bold');
            set(self.t4_txt5, 'BackgroundColor', self.background_color); 
            set(self.t4_txt5, 'Visible', 'off');

            % variable slider
            self.t4_sld2 = uicontrol('style','slider');
            set(self.t4_sld2, 'unit', 'pixels', 'position', [40 230 240 20]);
            set(self.t4_sld2, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t4_sld2, 'BackgroundColor', self.background_color);
            set(self.t4_sld2, 'Min', 1, 'Max', 10, 'Value', 1);
            set(self.t4_sld2, 'Visible', 'off');

            % variable menu
            self.t4_mnu2 = uicontrol('style', 'popupmenu');
            set(self.t4_mnu2, 'position',[40 185 240 30]);            
            set(self.t4_mnu2, 'String', "none");
            set(self.t4_mnu2, 'Visible', 'off');

            % response label
            self.t4_txt6 = uicontrol('style', 'text');
            set(self.t4_txt6, 'unit', 'pixels', 'position', [30 125 250 30]);
            set(self.t4_txt6, 'String', ' Responding to');
            set(self.t4_txt6, 'HorizontalAlignment', 'left');
            set(self.t4_txt6, 'FontName', 'Serif');
            set(self.t4_txt6, 'FontSize', 14);
            set(self.t4_txt6, 'FontWeight', 'bold');
            set(self.t4_txt6, 'FontAngle', 'italic');
            set(self.t4_txt6, 'BackgroundColor', self.background_color);  
            set(self.t4_txt6, 'Visible', 'off');

            % frame around response
            self.t4_frm4 = uicontrol('style','frame');
            set(self.t4_frm4, 'unit', 'pixels', 'position', [20 20 280 110]);
            set(self.t4_frm4, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t4_frm4, 'BackgroundColor', self.background_color);
            set(self.t4_frm4, 'Visible', 'off');

            % response selection label
            self.t4_txt7 = uicontrol('style', 'text');
            set(self.t4_txt7, 'unit', 'pixels', 'position', [30 90 200 30]);
            set(self.t4_txt7, 'String', ' none');
            set(self.t4_txt7, 'HorizontalAlignment', 'left');
            set(self.t4_txt7, 'FontName', 'Serif');
            set(self.t4_txt7, 'FontSize', 12);
            set(self.t4_txt7, 'FontWeight', 'bold');
            set(self.t4_txt7, 'BackgroundColor', self.background_color); 
            set(self.t4_txt7, 'Visible', 'off');

            % response slider
            self.t4_sld3 = uicontrol('style','slider');
            set(self.t4_sld3, 'unit', 'pixels', 'position', [40 75 240 20]);
            set(self.t4_sld3, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t4_sld3, 'BackgroundColor', self.background_color);
            set(self.t4_sld3, 'Min', 1, 'Max', 10, 'Value', 1); 
            set(self.t4_sld3, 'Visible', 'off');

            % response menu
            self.t4_mnu3 = uicontrol('style', 'popupmenu');
            set(self.t4_mnu3, 'position',[40 30 240 30]);            
            set(self.t4_mnu3, 'String', "none");
            set(self.t4_mnu3, 'Visible', 'off');

            % default no-graphic image
            self.t4_img1 = axes('units', 'pixels', 'position', [320 20 660 570]);
            [img, map, alphachannel] = imread(fullfile(self.interface_path, 'no_graphic.png'));
            image(img);
            set(self.t4_img1, 'xtick', [], 'ytick', []);
            set(self.t4_img1, 'Visible', 'off');
            set(get(self.t4_img1,'children'),'visible','off');
            
            % if graphic view is activated, update tab 4 with graphics to turn it into a graphics navigator
            self.update_tab_4();
        end
        
        
        function hide_tab_4(self)

            % hide all controls            
            set(self.t4_txt1, 'Visible', 'off');
            set(self.t4_frm1, 'Visible', 'off');
            set(self.t4_txt2, 'Visible', 'off');
            set(self.t4_sld1, 'Visible', 'off');
            set(self.t4_mnu1, 'Visible', 'off');
            set(self.t4_txt3, 'Visible', 'off');
            set(self.t4_frm2, 'Visible', 'off');
            set(self.t4_bgr1, 'Visible', 'off');
            set(self.t4_rdb1, 'Visible', 'off');
            set(self.t4_rdb2, 'Visible', 'off');
            set(self.t4_txt4, 'Visible', 'off');
            set(self.t4_frm3, 'Visible', 'off');
            set(self.t4_txt5, 'Visible', 'off');
            set(self.t4_sld2, 'Visible', 'off');
            set(self.t4_mnu2, 'Visible', 'off');
            set(self.t4_txt6, 'Visible', 'off');
            set(self.t4_frm4, 'Visible', 'off');
            set(self.t4_txt7, 'Visible', 'off');
            set(self.t4_sld3, 'Visible', 'off');
            set(self.t4_mnu3, 'Visible', 'off');
            set(self.t4_img1, 'Visible', 'off');
            set(get(self.t4_img1,'children'),'visible','off');

            % update tab color
            set(self.tab_pbt4, 'BackgroundColor', self.backtabs_color);  
        end
        
        
        function show_tab_4(self)

            % show all controls
            set(self.t4_txt1, 'Visible', 'on');
            set(self.t4_frm1, 'Visible', 'on');
            set(self.t4_txt2, 'Visible', 'on');
            set(self.t4_sld1, 'Visible', 'on');
            set(self.t4_mnu1, 'Visible', 'on');
            set(self.t4_txt3, 'Visible', 'on');
            set(self.t4_frm2, 'Visible', 'on');
            set(self.t4_bgr1, 'Visible', 'on');
            set(self.t4_rdb1, 'Visible', 'on');
            set(self.t4_rdb2, 'Visible', 'on');
            set(self.t4_txt4, 'Visible', 'on');
            set(self.t4_frm3, 'Visible', 'on');
            set(self.t4_txt5, 'Visible', 'on');
            set(self.t4_sld2, 'Visible', 'on');
            set(self.t4_mnu2, 'Visible', 'on');
            set(self.t4_txt6, 'Visible', 'on');
            set(self.t4_frm4, 'Visible', 'on');
            set(self.t4_txt7, 'Visible', 'on');
            set(self.t4_sld3, 'Visible', 'on');
            set(self.t4_mnu3, 'Visible', 'on');
            set(self.t4_img1, 'Visible', 'on');
            set(get(self.t4_img1, 'children'), 'visible', 'on'); 
        end
        
        
        function update_tab_4(self)

            if self.view_graphics
                % recover path to graphics folder
                self.get_graphics_folder_path();               
                % recover the names of all images in graphics folder
                self.get_image_names();
                % generate rollmenus from the list of splitted images
                self.get_rollmenus();
                % initiate application, variable, and response values
                self.initiate_values();
                % update all controls
                self.update_application_slider();
                self.update_application_rollmenu();
                self.update_variable_slider();
                self.update_variable_rollmenu(); 
                self.update_response_slider();
                self.update_response_rollmenu();              
                % set all controls to initial positions
                self.set_application_label();
                self.set_application_slider();         
                self.set_application_rollmenu(); 
                self.set_variable_label();
                self.set_variable_slider();        
                self.set_variable_rollmenu();
                self.set_response_label();
                self.set_response_slider();
                self.set_response_rollmenu();              
                % initiate callbacks
                self.initiate_callbacks();
                % initiate image
                self.update_image();
            end
        end
        
        
        function cb_t4_sld1(self, hObject, callbackdata)
            % update current application
            items = get(self.t4_mnu1,'string');
            self.current_application = string(items(get(self.t4_sld1, 'value')+1));
            % set application controls
            self.set_application_label();
            self.set_application_rollmenu();
            % update variables 
            self.update_current_variables();
            self.update_current_variable();
            % update variable controls
            self.update_variable_slider();
            self.update_variable_rollmenu();
            % set variable controls
            self.set_variable_label();
            self.set_variable_slider();       
            self.set_variable_rollmenu();
            % update responses
            self.update_current_responses();
            self.update_current_response();
            % update response controls
            self.update_response_slider();
            self.update_response_rollmenu();
            % set response controls
            self.set_response_label();
            self.set_response_slider();
            self.set_response_rollmenu();       
            % update image
            self.update_image();
        end
        
        
        function cb_t4_mnu1(self, hObject, callbackdata)
            % update current application
            index = get(self.t4_mnu1, 'value');
            items = get(self.t4_mnu1,'string');
            if index ~= 1
                self.current_application = string(items(index));
            end 
            % set application controls
            self.set_application_label();
            self.set_application_slider();      
            self.set_application_rollmenu();
            % update variables 
            self.update_current_variables();
            self.update_current_variable();
            % update variable controls
            self.update_variable_slider();
            self.update_variable_rollmenu();
            % set variable controls
            self.set_variable_label();
            self.set_variable_slider();       
            self.set_variable_rollmenu();
            % update responses
            self.update_current_responses();
            self.update_current_response();
            % update response controls
            self.update_response_slider();
            self.update_response_rollmenu();
            % set response controls
            self.set_response_label();
            self.set_response_slider();
            self.set_response_rollmenu();  
            % update image
            self.update_image();   
        end
     
                
        function cb_t4_bgr1(self, hObject, callbackdata)
            % if joint plot is True, update and set current variable to all
            if get(self.t4_bgr1, 'SelectedObject') == self.t4_rdb1
                self.joint_plot = true;
            elseif get(self.t4_bgr1, 'SelectedObject') == self.t4_rdb2
                self.joint_plot = false;
            end
            % update variables
            self.update_current_variable();
            % update variable controls
            self.update_variable_slider();
            self.update_variable_rollmenu();   
            % set variable controls
            self.set_variable_label();
            self.set_variable_slider();
            self.set_variable_rollmenu();   
            self.switch_variable_controls();       
            % update responses
            self.update_current_responses();
            self.update_current_response();
            % update response controls
            self.update_response_slider();
            self.update_response_rollmenu();
            % set response controls
            self.switch_response_controls();
            self.set_response_label();
            self.set_response_slider();
            self.set_response_rollmenu();       
            % update image
            self.update_image();
        end
   
        
        function cb_t4_sld2(self, hObject, callbackdata)
            % update current variable
            items = get(self.t4_mnu2,'string');
            self.current_variable = string(items(get(self.t4_sld2, 'value')+1)); 
            % set variable controls
            self.set_variable_label();      
            self.set_variable_rollmenu();
            % update responses
            self.update_current_responses();
            self.update_current_response();
            % update response controls
            self.update_response_slider();
            self.update_response_rollmenu();
            % set response controls
            self.set_response_label();
            self.set_response_slider();
            self.set_response_rollmenu();       
            % update image
            self.update_image();
        end        
       
        
        function cb_t4_mnu2(self, hObject, callbackdata)
            % update current variable
            index = get(self.t4_mnu2, 'value');
            items = get(self.t4_mnu2,'string');
            if index ~= 1
                self.current_variable = string(items(index));
            end 
            % set variable controls
            self.set_variable_label();
            self.set_variable_slider();       
            self.set_variable_rollmenu();
            % update responses
            self.update_current_responses();
            self.update_current_response();
            % update response controls
            self.update_response_slider();
            self.update_response_rollmenu();
            % set response controls
            self.set_response_label();
            self.set_response_slider();
            self.set_response_rollmenu();  
            % update image
            self.update_image();   
        end        
        
        
        function cb_t4_sld3(self, hObject, callbackdata)
            % update current response
            items = get(self.t4_mnu3,'string');
            self.current_response = string(items(get(self.t4_sld3, 'value')+1)); 
            % set response controls
            self.set_response_label();
            self.set_response_rollmenu();       
            % update image
            self.update_image();
        end         

        
        function cb_t4_mnu3(self, hObject, callbackdata)
            % update current response
            index = get(self.t4_mnu3, 'value');
            items = get(self.t4_mnu3,'string');
            if index ~= 1
                self.current_response = string(items(index));
            end 
            % set response controls
            self.set_response_label();
            self.set_response_slider();
            self.set_response_rollmenu();  
            % update image
            self.update_image();   
        end         
        

    end
    
    
    methods (Access = protected, Hidden = true)
        
        
        function get_graphics_folder_path(self)
            project_path = self.user_inputs.tab_1.project_path;
            graphics_folder_path = fullfile(project_path, 'graphics');
            self.graphics_folder_path = graphics_folder_path;
        end
        
        
        function get_image_names(self)
            % get list of all image files in graphics folder
            graphics_folder_files = dir(self.graphics_folder_path);
            file_list = string({graphics_folder_files.name});
            image_list = file_list(3:end);
            n = numel(image_list);
            % initiate list of split names
            split_image_string = strings([n,3]);         
            for i=1:n
                image = image_list(i);
                % application is first part, before underscore
                temp = split(image,'_');
                application = temp(1);
                % variables and response are parts after underscore, removing png extentions
                temp = extractAfter(image, "_");
                variable_and_response = extractBefore(temp,'.png');
                % variable is part before @, response part after @ (if any)
                split_variable_and_response = split(variable_and_response, "@");              
                variable = split_variable_and_response(1);
                if numel(split_variable_and_response) == 2
                    response = split_variable_and_response(2);
                else
                    response = "none";
                end
                split_image_string(i,1) = application;
                split_image_string(i,2) = variable;
                split_image_string(i,3) = response;
            end
            self.split_image_string = split_image_string;
            % also, get list of all applications (obtained in arbitrary order for now)
            applications = unique(split_image_string(:,1));
            % reorganize applications so that they are in the right order
            sorted_applications = [];
            possible_applications = ["fit", "residuals", "forecasts", ...
                                     "conditional_forecasts", "irf", "fevd", "hd"];
            for i=1:numel(possible_applications)
                if any(contains(applications, possible_applications(i)))
                    sorted_applications = [sorted_applications possible_applications(i)];
                end
            end
            self.applications = sorted_applications;
            % check whether there are no images to display
            if iu.is_empty(self.applications)
                error(['Image error: graphics and figures are selected, but there are no images to display. Select applications producing figures.']);
            end
        end
        
        
        function get_rollmenus(self)
            % initiate structure of rollmenus
            rollmenus = struct;
            split_image_string = self.split_image_string;
            applications = self.applications;
            rollmenus.applications = applications;
            % loop over applications
            for i=1:numel(applications)
                application = applications(i);
                % add applications to structure of rollmenus
                rollmenus.(application) = struct;                
                % list of all files matching the application
                application_files = split_image_string(split_image_string(:,1) == application,:);
                % list of all variables
                variables = sort(unique(application_files(:,2)'));
                if any(contains(variables, "all"))
                    variables = variables(variables~="all");
                    variables = ["all" variables];
                end
                % add variables to application structure
                rollmenus.(application).variables = variables;
                % loop over variables
                for j=1:numel(variables)
                    variable = variables(j);
                    % list of all files matching the variable
                    variable_files = application_files(application_files(:,2) == variable,:);
                    % list of all corresponding responses
                    responses = sort(variable_files(:,3)');
                    if any(contains(responses, "all"))
                        responses = responses(responses~="all");
                        responses = ["all" responses];
                    end
                    % add responses to variable structure
                    rollmenus.(application).(variable) = responses;
                end
            end
            self.rollmenus = rollmenus;
        end
        
        
        function initiate_values(self)
            self.current_application = self.applications(1);
            self.joint_plot = true;
            self.current_variables = self.rollmenus.(self.current_application).variables;
            self.current_variable = self.current_variables(1);
            self.current_responses = self.rollmenus.(self.current_application).(self.current_variable);
            self.current_response = self.current_responses(1);
        end

            
        function initiate_callbacks(self)
            set(self.t4_sld1, 'CallBack', @self.cb_t4_sld1);
            set(self.t4_mnu1, 'CallBack', @self.cb_t4_mnu1);
            set(self.t4_bgr1, 'SelectionChangeFcn', @self.cb_t4_bgr1);
            set(self.t4_sld2, 'CallBack', @self.cb_t4_sld2);
            set(self.t4_mnu2, 'CallBack', @self.cb_t4_mnu2);
            set(self.t4_sld3, 'CallBack', @self.cb_t4_sld3);
            set(self.t4_mnu3, 'CallBack', @self.cb_t4_mnu3);
        end
        
 
        function switch_variable_controls(self)
            if self.joint_plot
                set(self.t4_sld2, 'Enable', 'off');
                set(self.t4_mnu2, 'Enable', 'off');
            else
                set(self.t4_sld2, 'Enable', 'on');
                set(self.t4_mnu2, 'Enable', 'on');
            end
        end
                           
            
        function switch_response_controls(self)
            if self.joint_plot
                set(self.t4_sld3, 'Enable', 'off');
                set(self.t4_mnu3, 'Enable', 'off');
            else
                set(self.t4_sld3, 'Enable', 'on');
                set(self.t4_mnu3, 'Enable', 'on');
            end
        end
        
        
        function update_current_application(self)
            self.current_application = self.rollmenus.applications(1);
        end
        
        
        function update_current_variables(self)
            self.current_variables = self.rollmenus.(self.current_application).variables;
        end
        
        
        function update_current_variable(self)
            % recover list of all variables of current application
            current_variables = self.current_variables;
            % if joint plot, current variable must be all by default
            if self.joint_plot
                self.current_variable = "all";
            % if not joint plot, make sure all is not in the list of current variables
            else
                if any(contains(current_variables, "all"))
                    current_variables = current_variables(current_variables~="all");
                end
            end
            % now, if current variable is not in current variables, use first variable by default
            if ~any(contains(current_variables, self.current_variable))
                self.current_variable = current_variables(1);
            end
        end
                

        function update_current_responses(self)
            self.current_responses = self.rollmenus.(self.current_application).(self.current_variable);
        end


        function update_current_response(self)
            if ~any(contains(self.current_responses, self.current_response))
                self.current_response = self.current_responses(1);
            end
        end
        
        
        function update_application_slider(self)
            total_steps = numel(self.applications);
            set(self.t4_sld1, 'Max', total_steps);
            if total_steps == 1
                set(self.t4_sld1, 'Max', 1.000001);
                set(self.t4_sld1, 'SliderStep', [1 , 10000]);
            else
                set(self.t4_sld1, 'Max', total_steps);
                set(self.t4_sld1, 'SliderStep', [1/(total_steps-1) , 3/(total_steps-1)]);
            end  
        end
            
            
        function update_application_rollmenu(self)
            items = ["select" self.applications];
            set(self.t4_mnu1, 'String', items);
            set(self.t4_mnu1, 'Value', 1);
        end
        
        
        function update_variable_slider(self)
            if self.joint_plot
                set(self.t4_sld2, 'Enable', 'off');
                total_steps = 1;
            else
                set(self.t4_sld2, 'Enable', 'on');
                current_variables = self.current_variables;
                if any(contains(current_variables, "all"))
                    current_variables = current_variables(current_variables~="all");
                end
                total_steps = numel(current_variables);
            end
            if total_steps == 1
                set(self.t4_sld2, 'Max', 1.000001);
                set(self.t4_sld2, 'SliderStep', [1 , 10000]);
            else
                set(self.t4_sld2, 'Max', total_steps);
                set(self.t4_sld2, 'SliderStep', [1/(total_steps-1) , 3/(total_steps-1)]);
            end            
        end
                     
            
        function update_variable_rollmenu(self)
            if self.joint_plot
                set(self.t4_mnu2, 'Enable', 'off');
                set(self.t1_mnu2, 'String', "none");
            else
                set(self.t4_mnu2, 'Enable', 'on');
                current_variables = self.current_variables;
                if any(contains(current_variables, "all"))
                    current_variables = current_variables(current_variables~="all");
                end
                items = ["select" current_variables];
                set(self.t4_mnu2, 'String', items);
            end
        end
        
        
        function update_response_slider(self)
            if self.joint_plot
                set(self.t4_sld3, 'Enable', 'off');
                total_steps = 1;
            else
                set(self.t4_sld3, 'Enable', 'on');
                total_steps = numel(self.current_responses);
            end
            if total_steps == 1
                set(self.t4_sld3, 'Max', 1.000001);
                set(self.t4_sld3, 'SliderStep', [1 , 10000]);
            else
                set(self.t4_sld3, 'Max', total_steps);
                set(self.t4_sld3, 'SliderStep', [1/(total_steps-1) , 3/(total_steps-1)]);
            end
        end
          
        
        function update_response_rollmenu(self)
            if self.joint_plot
                set(self.t4_mnu3, 'Enable', 'off');
                set(self.t4_mnu3, 'String', "none");
            else
                set(self.t4_mnu3, 'Enable', 'on');
                current_variables = self.current_variables;
                items = ["select" self.current_responses];
                set(self.t4_mnu3, 'String', items);
            end
        end            
            
        
        function set_application_label(self)
            set(self.t4_txt2, 'String', [' ' convertStringsToChars(self.current_application)]);
        end
            
            
        function set_application_slider(self)
            index = find(get(self.t4_mnu1,'String') == self.current_application) - 1;
            set(self.t4_sld1, 'Value', index);
        end
        
        
        function set_application_rollmenu(self)
            set(self.t4_mnu1, 'Value', 1);
        end
        
        
        function set_variable_label(self)
            if self.current_variable == 'all'
                set(self.t4_txt5, 'String', ' none');
            else
                set(self.t4_txt5, 'String', [' ' convertStringsToChars(self.current_variable)]);
            end
        end
                
            
        function set_variable_slider(self)
            if self.current_variable == 'all' || self.joint_plot
                index = 1;
            else
                index = find(get(self.t4_mnu2,'String') == self.current_variable) - 1;
            end
            set(self.t4_sld2, 'Value', index);
        end
        
        
        function set_variable_rollmenu(self)
            set(self.t4_mnu2, 'Value', 1);
        end
        
        
        function set_response_label(self)
            if self.current_variable == 'none'
                set(self.t4_txt7, 'String', ' none');
            else
                set(self.t4_txt7, 'String', [' ' convertStringsToChars(self.current_response)]);
            end
        end        
        
        
        function set_response_slider(self)
            if self.current_response == 'none'
                index = 1;
            else
                 index = find(get(self.t4_mnu3,'String') == self.current_response) - 1;
            end
            set(self.t4_sld3, 'Value', index);
        end   
        
        
        function set_response_rollmenu(self)
            set(self.t4_mnu3, 'Value', 1);
        end    
        
        
        function update_image(self)
            image_name = [convertStringsToChars(self.current_application)...
                          '_' convertStringsToChars(self.current_variable)];
            if self.current_response ~= 'none'
                image_name = [image_name '@' convertStringsToChars(self.current_response)];
            end
            image_name = [image_name '.png'];     
            self.image_name = image_name;
            [img, map, alphachannel] = imread(fullfile(self.graphics_folder_path, image_name));
            image = imresize(img, [570 660]);
            image_handle = imhandles(self.t4_img1);
            set(image_handle, 'CData', image);
        end
        
           
    end   
    
    
end



