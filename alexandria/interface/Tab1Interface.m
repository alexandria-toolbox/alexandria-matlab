classdef Tab1Interface < handle
    
    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess = public)
        % tab 1 properties
        t1_txt1
        t1_txt2
        t1_txt3
        t1_img1
        t1_txt4
        t1_frm1
        t1_txt5
        t1_mnu1
        t1_txt6
        t1_edt1
        t1_txt7
        t1_edt2
        t1_txt8
        t1_mnu2
        t1_txt9
        t1_edt3
        t1_txt10
        t1_frm2
        t1_txt11
        t1_edt4
        t1_txt12
        t1_edt5
        t1_txt13
        t1_bgr1
        t1_rdb1
        t1_rdb2
        t1_txt14
        t1_bgr2
        t1_rdb3
        t1_rdb4
        t1_txt15
        t1_bgr3
        t1_rdb5
        t1_rdb6
        t1_pbt1
        t1_pbt2        
    end

    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)
        

        function self = Tab1Interface()
        end
        
        
        function create_tab_1(self)

            % main title
            self.t1_txt1 = uicontrol('style', 'text');
            set(self.t1_txt1, 'unit', 'pixels', 'position', [10 535 600 70]);
            set(self.t1_txt1, 'String', ' Alexandria');
            set(self.t1_txt1, 'HorizontalAlignment', 'left');
            set(self.t1_txt1, 'FontName', 'Serif');
            set(self.t1_txt1, 'FontSize', 34);
            set(self.t1_txt1, 'FontWeight', 'bold');
            set(self.t1_txt1, 'FontAngle', 'italic');
            set(self.t1_txt1, 'BackgroundColor', self.background_color);

            % subtitle
            self.t1_txt2 = uicontrol('style', 'text');
            set(self.t1_txt2, 'unit', 'pixels', 'position', [15 520 600 30]);
            set(self.t1_txt2, 'String', ' The library of Bayesian time-series models');
            set(self.t1_txt2, 'HorizontalAlignment', 'left');
            set(self.t1_txt2, 'FontName', 'Serif');
            set(self.t1_txt2, 'FontSize', 16);
            set(self.t1_txt2, 'FontWeight', 'bold');
            set(self.t1_txt2, 'FontAngle', 'italic');
            set(self.t1_txt2, 'BackgroundColor', self.background_color); 

            % version
            self.t1_txt3 = uicontrol('style', 'text');
            set(self.t1_txt3, 'unit', 'pixels', 'position', [15 490 600 30]);
            set(self.t1_txt3, 'String', ' V 1.0 - Matlab edition');
            set(self.t1_txt3, 'HorizontalAlignment', 'left');
            set(self.t1_txt3, 'FontName', 'Serif');
            set(self.t1_txt3, 'FontSize', 16);
            set(self.t1_txt3, 'FontWeight', 'bold');
            set(self.t1_txt3, 'FontAngle', 'italic');
            set(self.t1_txt3, 'BackgroundColor', self.background_color);   

            % Matlab logo
            self.t1_img1 = axes('units', 'pixels', 'position', [850 490 125 100]);
            [img, map, alphachannel] = imread(fullfile(self.interface_path, 'matlab.png'));
            image(img, 'AlphaData', alphachannel);
            set(self.t1_img1, 'color', self.background_color);
            set(self.t1_img1, 'XColor', self.background_color, 'YColor', self.background_color);
            set(self.t1_img1, 'xtick', [], 'ytick', []);

            % model label
            self.t1_txt4 = uicontrol('style', 'text');
            set(self.t1_txt4, 'unit', 'pixels', 'position', [30 430 600 30]);
            set(self.t1_txt4, 'String', ' Model');
            set(self.t1_txt4, 'HorizontalAlignment', 'left');
            set(self.t1_txt4, 'FontName', 'Serif');
            set(self.t1_txt4, 'FontSize', 16);
            set(self.t1_txt4, 'FontWeight', 'bold');
            set(self.t1_txt4, 'FontAngle', 'italic');
            set(self.t1_txt4, 'BackgroundColor', self.background_color);

            % frame around model
            self.t1_frm1 = uicontrol('style','frame');
            set(self.t1_frm1, 'unit', 'pixels', 'position', [20 20 380 410]);
            set(self.t1_frm1, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t1_frm1, 'BackgroundColor', self.background_color);

            % model selection label
            self.t1_txt5 = uicontrol('style', 'text');
            set(self.t1_txt5, 'unit', 'pixels', 'position', [30 390 300 30]);
            set(self.t1_txt5, 'String', ' model selection');
            set(self.t1_txt5, 'HorizontalAlignment', 'left');
            set(self.t1_txt5, 'FontName', 'Serif');
            set(self.t1_txt5, 'FontSize', 12);
            set(self.t1_txt5, 'FontWeight', 'bold');
            set(self.t1_txt5, 'BackgroundColor', self.background_color); 

            % model selection menu
            self.t1_mnu1 = uicontrol('style', 'popupmenu');
            set(self.t1_mnu1, 'position',[38 345 250 50]);            
            set(self.t1_mnu1, 'String', {'1. linear regression', '2. vector autoregression'});    
            set(self.t1_mnu1, 'Value', self.user_inputs.tab_1.model);
            set(self.t1_mnu1, 'CallBack', @self.cb_t1_mnu1);

            % endogenous label
            self.t1_txt6 = uicontrol('style', 'text');
            set(self.t1_txt6, 'unit', 'pixels', 'position', [30 320 300 30]);
            set(self.t1_txt6, 'String', ' endogenous variables');
            set(self.t1_txt6, 'HorizontalAlignment', 'left');
            set(self.t1_txt6, 'FontName', 'Serif');
            set(self.t1_txt6, 'FontSize', 12);
            set(self.t1_txt6, 'FontWeight', 'bold');
            set(self.t1_txt6, 'BackgroundColor', self.background_color);

            % endogenous edit
            self.t1_edt1 = uicontrol('style','edit');
            set(self.t1_edt1, 'unit', 'pixels', 'position', [35 270 340 60]);
            set(self.t1_edt1, 'HorizontalAlignment', 'left', 'Max', 2);     
            set(self.t1_edt1, 'String', self.user_inputs.tab_1.endogenous_variables);
            set(self.t1_edt1, 'CallBack', @self.cb_t1_edt1);

            % exogenous label
            self.t1_txt7 = uicontrol('style', 'text');
            set(self.t1_txt7, 'unit', 'pixels', 'position', [30 225 300 30]);
            set(self.t1_txt7, 'String', ' exogenous variables');
            set(self.t1_txt7, 'HorizontalAlignment', 'left');
            set(self.t1_txt7, 'FontName', 'Serif');
            set(self.t1_txt7, 'FontSize', 12);
            set(self.t1_txt7, 'FontWeight', 'bold');
            set(self.t1_txt7, 'BackgroundColor', self.background_color);

            % exogenous edit
            self.t1_edt2 = uicontrol('style','edit');
            set(self.t1_edt2, 'unit', 'pixels', 'position', [35 175 340 60]);
            set(self.t1_edt2, 'HorizontalAlignment', 'left', 'Max', 2);           
            set(self.t1_edt2, 'String', self.user_inputs.tab_1.exogenous_variables);
            set(self.t1_edt2, 'CallBack', @self.cb_t1_edt2);

            % frequency label
            self.t1_txt8 = uicontrol('style', 'text');
            set(self.t1_txt8, 'unit', 'pixels', 'position', [30 125 300 30]);
            set(self.t1_txt8, 'String', ' data frequency');
            set(self.t1_txt8, 'HorizontalAlignment', 'left');
            set(self.t1_txt8, 'FontName', 'Serif');
            set(self.t1_txt8, 'FontSize', 12);
            set(self.t1_txt8, 'FontWeight', 'bold');
            set(self.t1_txt8, 'BackgroundColor', self.background_color);  

            % data frequency menu
            self.t1_mnu2 = uicontrol('style', 'popupmenu');
            set(self.t1_mnu2, 'position',[38 80 250 50]);            
            set(self.t1_mnu2, 'String', ...
                {'1. cross-sectional/undated', '2. annual', '3. quarterly', ...
                '4. monthly', '5. weekly', '6. daily'});
            set(self.t1_mnu2, 'Value', self.user_inputs.tab_1.frequency);
            set(self.t1_mnu2, 'CallBack', @self.cb_t1_mnu2);

            % sample label
            self.t1_txt9 = uicontrol('style', 'text');
            set(self.t1_txt9, 'unit', 'pixels', 'position', [30 55 300 30]);
            set(self.t1_txt9, 'String', ' estimation sample');
            set(self.t1_txt9, 'HorizontalAlignment', 'left');
            set(self.t1_txt9, 'FontName', 'Serif');
            set(self.t1_txt9, 'FontSize', 12);
            set(self.t1_txt9, 'FontWeight', 'bold');
            set(self.t1_txt9, 'BackgroundColor', self.background_color);   

            % sample edit
            self.t1_edt3 = uicontrol('style','edit');
            set(self.t1_edt3, 'unit', 'pixels', 'position', [35 35 340 25]);
            set(self.t1_edt3, 'HorizontalAlignment', 'left');              
            set(self.t1_edt3, 'String', self.user_inputs.tab_1.sample);
            set(self.t1_edt3, 'CallBack', @self.cb_t1_edt3);

            % settings label
            self.t1_txt10 = uicontrol('style', 'text');
            set(self.t1_txt10, 'unit', 'pixels', 'position', [430 430 600 30]);
            set(self.t1_txt10, 'String', 'Settings');
            set(self.t1_txt10, 'HorizontalAlignment', 'left');
            set(self.t1_txt10, 'FontName', 'Serif');
            set(self.t1_txt10, 'FontSize', 16);
            set(self.t1_txt10, 'FontWeight', 'bold');
            set(self.t1_txt10, 'FontAngle', 'italic');
            set(self.t1_txt10, 'BackgroundColor', self.background_color); 

            % frame around settings
            self.t1_frm2 = uicontrol('style','frame');
            set(self.t1_frm2, 'unit', 'pixels', 'position', [420 20 380 410]);
            set(self.t1_frm2, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t1_frm2, 'BackgroundColor', self.background_color);

            % project path label
            self.t1_txt11 = uicontrol('style', 'text');
            set(self.t1_txt11, 'unit', 'pixels', 'position', [430 390 300 30]);
            set(self.t1_txt11, 'String', ' path to project folder');
            set(self.t1_txt11, 'HorizontalAlignment', 'left');
            set(self.t1_txt11, 'FontName', 'Serif');
            set(self.t1_txt11, 'FontSize', 12);
            set(self.t1_txt11, 'FontWeight', 'bold');
            set(self.t1_txt11, 'BackgroundColor', self.background_color); 

            % project path edit
            self.t1_edt4 = uicontrol('style','edit');
            set(self.t1_edt4, 'position',[438 370 340 25]);            
            set(self.t1_edt4, 'HorizontalAlignment', 'left'); 
            set(self.t1_edt4, 'String', self.user_inputs.tab_1.project_path);
            set(self.t1_edt4, 'CallBack', @self.cb_t1_edt4);

            % data file label
            self.t1_txt12 = uicontrol('style', 'text');
            set(self.t1_txt12, 'unit', 'pixels', 'position', [430 320 300 30]);
            set(self.t1_txt12, 'String', ' data file');
            set(self.t1_txt12, 'HorizontalAlignment', 'left');
            set(self.t1_txt12, 'FontName', 'Serif');
            set(self.t1_txt12, 'FontSize', 12);
            set(self.t1_txt12, 'FontWeight', 'bold');
            set(self.t1_txt12, 'BackgroundColor', self.background_color);

            % data file edit
            self.t1_edt5 = uicontrol('style','edit');
            set(self.t1_edt5, 'unit', 'pixels', 'position', [438 305 340 25]);
            set(self.t1_edt5, 'HorizontalAlignment', 'left');             
            set(self.t1_edt5, 'String', self.user_inputs.tab_1.data_file);
            set(self.t1_edt5, 'CallBack', @self.cb_t1_edt5);

            % progress bar label
            self.t1_txt13 = uicontrol('style', 'text');
            set(self.t1_txt13, 'unit', 'pixels', 'position', [430 255 300 30]);
            set(self.t1_txt13, 'String', ' progress bar');
            set(self.t1_txt13, 'HorizontalAlignment', 'left');
            set(self.t1_txt13, 'FontName', 'Serif');
            set(self.t1_txt13, 'FontSize', 12);
            set(self.t1_txt13, 'FontWeight', 'bold');
            set(self.t1_txt13, 'BackgroundColor', self.background_color);

            % progress bar radiobuttons
            self.t1_bgr1 = uibuttongroup('unit','pixels', 'Position',[430 230 200 30]);
            set(self.t1_bgr1, 'BorderType', 'none');
            set(self.t1_bgr1, 'BackgroundColor', self.background_color);            
            self.t1_rdb1 = uicontrol(self.t1_bgr1,'Style','radiobutton');
            set(self.t1_rdb1, 'Position',[8 10 80 25]);
            set(self.t1_rdb1, 'String',' yes');
            set(self.t1_rdb1, 'FontName', 'Serif');
            set(self.t1_rdb1, 'FontSize', 12);
            set(self.t1_rdb1, 'FontWeight', 'bold');
            set(self.t1_rdb1, 'BackgroundColor', self.background_color);
            self.t1_rdb2 = uicontrol(self.t1_bgr1,'Style','radiobutton');
            set(self.t1_rdb2, 'Position',[88 10 80 25]);
            set(self.t1_rdb2, 'String',' no');
            set(self.t1_rdb2, 'FontName', 'Serif');
            set(self.t1_rdb2, 'FontSize', 12);
            set(self.t1_rdb2, 'FontWeight', 'bold');
            set(self.t1_rdb2, 'BackgroundColor', self.background_color);
            if self.user_inputs.tab_1.progress_bar
                set(self.t1_bgr1, 'SelectedObject', self.t1_rdb1);
            else
                set(self.t1_bgr1, 'SelectedObject', self.t1_rdb2);
            end
            set(self.t1_bgr1, 'SelectionChangeFcn', @self.cb_t1_bgr1);

            % graphics label
            self.t1_txt14 = uicontrol('style', 'text');
            set(self.t1_txt14, 'unit', 'pixels', 'position', [430 190 300 30]);
            set(self.t1_txt14, 'String', ' graphics and figures');
            set(self.t1_txt14, 'HorizontalAlignment', 'left');
            set(self.t1_txt14, 'FontName', 'Serif');
            set(self.t1_txt14, 'FontSize', 12);
            set(self.t1_txt14, 'FontWeight', 'bold');
            set(self.t1_txt14, 'BackgroundColor', self.background_color);

            % graphics radiobuttons
            self.t1_bgr2 = uibuttongroup('unit','pixels', 'Position',[430 165 200 30]);
            set(self.t1_bgr2, 'BorderType', 'none');
            set(self.t1_bgr2, 'BackgroundColor', self.background_color);            
            self.t1_rdb3 = uicontrol(self.t1_bgr2,'Style','radiobutton');
            set(self.t1_rdb3, 'Position',[8 10 80 25]);
            set(self.t1_rdb3, 'String',' yes');
            set(self.t1_rdb3, 'FontName', 'Serif');
            set(self.t1_rdb3, 'FontSize', 12);
            set(self.t1_rdb3, 'FontWeight', 'bold');
            set(self.t1_rdb3, 'BackgroundColor', self.background_color);
            self.t1_rdb4 = uicontrol(self.t1_bgr2,'Style','radiobutton');
            set(self.t1_rdb4, 'Position',[88 10 80 25]);
            set(self.t1_rdb4, 'String',' no');
            set(self.t1_rdb4, 'FontName', 'Serif');
            set(self.t1_rdb4, 'FontSize', 12);
            set(self.t1_rdb4, 'FontWeight', 'bold');
            set(self.t1_rdb4, 'BackgroundColor', self.background_color); 
            if self.user_inputs.tab_1.create_graphics
                set(self.t1_bgr2, 'SelectedObject', self.t1_rdb3);
            else
                set(self.t1_bgr2, 'SelectedObject', self.t1_rdb4);
            end
            set(self.t1_bgr2, 'SelectionChangeFcn', @self.cb_t1_bgr2);

            % save label
            self.t1_txt15 = uicontrol('style', 'text');
            set(self.t1_txt15, 'unit', 'pixels', 'position', [430 125 300 30]);
            set(self.t1_txt15, 'String', ' save results in project folder');
            set(self.t1_txt15, 'HorizontalAlignment', 'left');
            set(self.t1_txt15, 'FontName', 'Serif');
            set(self.t1_txt15, 'FontSize', 12);
            set(self.t1_txt15, 'FontWeight', 'bold');
            set(self.t1_txt15, 'BackgroundColor', self.background_color);

            % save radiobuttons
            self.t1_bgr3 = uibuttongroup('unit','pixels', 'Position',[430 100 200 30]);
            set(self.t1_bgr3, 'BorderType', 'none');
            set(self.t1_bgr3, 'BackgroundColor', self.background_color);            
            self.t1_rdb5 = uicontrol(self.t1_bgr3,'Style','radiobutton');
            set(self.t1_rdb5, 'Position',[8 10 80 25]);
            set(self.t1_rdb5, 'String',' yes');
            set(self.t1_rdb5, 'FontName', 'Serif');
            set(self.t1_rdb5, 'FontSize', 12);
            set(self.t1_rdb5, 'FontWeight', 'bold');
            set(self.t1_rdb5, 'BackgroundColor', self.background_color);
            self.t1_rdb6 = uicontrol(self.t1_bgr3,'Style','radiobutton');
            set(self.t1_rdb6, 'Position',[88 10 80 25]);
            set(self.t1_rdb6, 'String',' no');
            set(self.t1_rdb6, 'FontName', 'Serif');
            set(self.t1_rdb6, 'FontSize', 12);
            set(self.t1_rdb6, 'FontWeight', 'bold');
            set(self.t1_rdb6, 'BackgroundColor', self.background_color);
            if self.user_inputs.tab_1.save_results
                set(self.t1_bgr3, 'SelectedObject', self.t1_rdb5);
            else
                set(self.t1_bgr3, 'SelectedObject', self.t1_rdb6);
            end
            set(self.t1_bgr3, 'SelectionChangeFcn', @self.cb_t1_bgr3);

            % reset pushbutton
            self.t1_pbt1 = uicontrol('style', 'pushbutton');
            set(self.t1_pbt1, 'unit', 'pixels', 'position', [510 50 200 30]);
            set(self.t1_pbt1, 'String', 'Reset all');
            set(self.t1_pbt1, 'FontSize', 14);
            set(self.t1_pbt1, 'FontName', 'Serif');      
            set(self.t1_pbt1, 'CallBack', @self.cb_t1_pbt1);

            % run pushbutton
            self.t1_pbt2 = uicontrol('style', 'pushbutton');
            set(self.t1_pbt2, 'unit', 'pixels', 'position', [820 50 160 260]);
            run_image = imread(fullfile(self.interface_path, 'run_button.png'));
            set(self.t1_pbt2,'cdata', run_image);
            set(self.t1_pbt2, 'CallBack', @self.cb_t1_pbt2);

            % initiate current tab
            self.current_tab = 'tab_1';
        end
        
        
        function hide_tab_1(self)
        
            % hide all controls
            set(self.t1_txt1, 'Visible', 'off');
            set(self.t1_txt2, 'Visible', 'off');
            set(self.t1_txt3, 'Visible', 'off');
            set(self.t1_img1, 'Visible', 'off');
            set(get(self.t1_img1,'children'),'visible','off');
            set(self.t1_txt4, 'Visible', 'off');
            set(self.t1_frm1, 'Visible', 'off');
            set(self.t1_txt5, 'Visible', 'off');
            set(self.t1_mnu1, 'Visible', 'off');
            set(self.t1_txt6, 'Visible', 'off');
            set(self.t1_edt1, 'Visible', 'off');
            set(self.t1_txt7, 'Visible', 'off');
            set(self.t1_edt2, 'Visible', 'off');
            set(self.t1_txt8, 'Visible', 'off');
            set(self.t1_mnu2, 'Visible', 'off');
            set(self.t1_txt9, 'Visible', 'off');
            set(self.t1_edt3, 'Visible', 'off');
            set(self.t1_txt10, 'Visible', 'off');
            set(self.t1_frm2, 'Visible', 'off');
            set(self.t1_txt11, 'Visible', 'off');
            set(self.t1_edt4, 'Visible', 'off');
            set(self.t1_txt12, 'Visible', 'off');
            set(self.t1_edt5, 'Visible', 'off');
            set(self.t1_txt13, 'Visible', 'off');
            set(self.t1_bgr1, 'Visible', 'off');
            set(self.t1_rdb1, 'Visible', 'off');
            set(self.t1_rdb2, 'Visible', 'off');
            set(self.t1_txt14, 'Visible', 'off');
            set(self.t1_bgr2, 'Visible', 'off');
            set(self.t1_rdb3, 'Visible', 'off');
            set(self.t1_rdb4, 'Visible', 'off');
            set(self.t1_txt15, 'Visible', 'off');
            set(self.t1_bgr3, 'Visible', 'off');
            set(self.t1_rdb5, 'Visible', 'off');
            set(self.t1_rdb6, 'Visible', 'off');
            set(self.t1_pbt1, 'Visible', 'off');
            set(self.t1_pbt2, 'Visible', 'off');

            % update tab color
            set(self.tab_pbt1, 'BackgroundColor', self.backtabs_color);
        end
        
        
        function show_tab_1(self)
            
            % show all controls
            set(self.t1_txt1, 'Visible', 'on');
            set(self.t1_txt2, 'Visible', 'on');
            set(self.t1_txt3, 'Visible', 'on');
            set(self.t1_img1, 'Visible', 'on');
            set(get(self.t1_img1,'children'), 'Visible', 'on');
            set(self.t1_txt4, 'Visible', 'on');
            set(self.t1_frm1, 'Visible', 'on');
            set(self.t1_txt5, 'Visible', 'on');
            set(self.t1_mnu1, 'Visible', 'on');
            set(self.t1_txt6, 'Visible', 'on');
            set(self.t1_edt1, 'Visible', 'on');
            set(self.t1_txt7, 'Visible', 'on');
            set(self.t1_edt2, 'Visible', 'on');
            set(self.t1_txt8, 'Visible', 'on');
            set(self.t1_mnu2, 'Visible', 'on');
            set(self.t1_txt9, 'Visible', 'on');
            set(self.t1_edt3, 'Visible', 'on');
            set(self.t1_txt10, 'Visible', 'on');
            set(self.t1_frm2, 'Visible', 'on');
            set(self.t1_txt11, 'Visible', 'on');
            set(self.t1_edt4, 'Visible', 'on');
            set(self.t1_txt12, 'Visible', 'on');
            set(self.t1_edt5, 'Visible', 'on');
            set(self.t1_txt13, 'Visible', 'on');
            set(self.t1_bgr1, 'Visible', 'on');
            set(self.t1_rdb1, 'Visible', 'on');
            set(self.t1_rdb2, 'Visible', 'on');
            set(self.t1_txt14, 'Visible', 'on');
            set(self.t1_bgr2, 'Visible', 'on');
            set(self.t1_rdb3, 'Visible', 'on');
            set(self.t1_rdb4, 'Visible', 'on');
            set(self.t1_txt15, 'Visible', 'on');
            set(self.t1_bgr3, 'Visible', 'on');
            set(self.t1_rdb5, 'Visible', 'on');
            set(self.t1_rdb6, 'Visible', 'on');
            set(self.t1_pbt1, 'Visible', 'on');
            set(self.t1_pbt2, 'Visible', 'on');
        end


        function cb_user_interrupt(self, src, callbackdata)
            % indicate that interface has been manually closed by user
            self.user_interrupt = true;
            % close interface and set press_run to true; this terminates waitfor
            delete(self.interface);
            self.press_run = true;
        end  
               
        
        function cb_t1_mnu1(self, hObject, callbackdata)
            self.user_inputs.tab_1.model = get(self.t1_mnu1, 'Value');
        end
        
        
        function cb_t1_edt1(self, hObject, callbackdata)
            self.user_inputs.tab_1.endogenous_variables = get(self.t1_edt1, 'String');
        end        
        
        
        function cb_t1_edt2(self, hObject, callbackdata)
            self.user_inputs.tab_1.exogenous_variables = get(self.t1_edt2, 'String');
        end
        
        
        function cb_t1_mnu2(self, hObject, callbackdata)
            self.user_inputs.tab_1.frequency = get(self.t1_mnu2, 'Value');
        end   
        
        
        function cb_t1_edt3(self, hObject, callbackdata)
            self.user_inputs.tab_1.sample = get(self.t1_edt3, 'String');
        end        
        
        
        function cb_t1_edt4(self, hObject, callbackdata)
            self.user_inputs.tab_1.project_path = get(self.t1_edt4, 'String');
        end         
        
        
        function cb_t1_edt5(self, hObject, callbackdata)
            self.user_inputs.tab_1.data_file = get(self.t1_edt5, 'String');
        end        
        
        
        function cb_t1_bgr1(self, hObject, callbackdata)
           if get(self.t1_bgr1, 'SelectedObject') == self.t1_rdb1
               self.user_inputs.tab_1.progress_bar = true;
           elseif get(self.t1_bgr1, 'SelectedObject') == self.t1_rdb2
               self.user_inputs.tab_1.progress_bar = false;
           end
        end        
        
        
        function cb_t1_bgr2(self, hObject, callbackdata)
           if get(self.t1_bgr2, 'SelectedObject') == self.t1_rdb3
               self.user_inputs.tab_1.create_graphics = true;
           elseif get(self.t1_bgr2, 'SelectedObject') == self.t1_rdb4
               self.user_inputs.tab_1.create_graphics = false;
           end
        end         
        
        
        function cb_t1_bgr3(self, hObject, callbackdata)
           if get(self.t1_bgr3, 'SelectedObject') == self.t1_rdb5
               self.user_inputs.tab_1.save_results = true;
           elseif get(self.t1_bgr3, 'SelectedObject') == self.t1_rdb6
               self.user_inputs.tab_1.save_results = false;
           end
        end         
        
        
        function cb_t1_pbt1(self, hObject, callbackdata)
            if exist(fullfile(self.interface_path, 'user_inputs.mat'), 'file') == 2
                delete(fullfile(self.interface_path, 'user_inputs.mat'));
            end
            self.create_default_inputs();
            self.create_tab_2_lr;
            self.reset_default_inputs();
            self.created_tab_2_var = false;
        end          
        

        function cb_t1_pbt2(self, hObject, callbackdata)
            self.validate_interface();
        end         
        

    end
  
    
end