classdef Tab3Interface < handle


    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess = public)
        % tab 3 properties
        t3_txt1
        t3_frm1
        t3_txt2
        t3_bgr1
        t3_rdb1
        t3_rdb2
        t3_edt1
        t3_txt3
        t3_bgr2
        t3_rdb3
        t3_rdb4
        t3_edt2
        t3_txt4
        t3_bgr3
        t3_rdb5
        t3_rdb6
        t3_edt3
        t3_txt5
        t3_bgr4
        t3_rdb7
        t3_rdb8
        t3_edt4
        t3_txt6
        t3_bgr5
        t3_rdb9
        t3_rdb10
        t3_edt5
        t3_txt7
        t3_frm2
        t3_txt8
        t3_txt9
        t3_edt6
        t3_txt10
        t3_mnu1
        t3_txt11
        t3_edt7
        t3_txt12
        t3_cbx1
        t3_txt13
        t3_txt14
        t3_edt8
        t3_txt15
        t3_mnu2
        t3_txt16
        t3_edt9
        t3_pbt1
        t3_txt17
        t3_txt18
    end
    
    
    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)
        

        function self = Tab3Interface()
        end
        
        
        function create_tab_3(self)
        
            % applications label
            self.t3_txt1 = uicontrol('style', 'text');
            set(self.t3_txt1, 'unit', 'pixels', 'position', [30 560 450 30]);
            set(self.t3_txt1, 'String', ' Applications and credibility levels');
            set(self.t3_txt1, 'HorizontalAlignment', 'left');
            set(self.t3_txt1, 'FontName', 'Serif');
            set(self.t3_txt1, 'FontSize', 16);
            set(self.t3_txt1, 'FontWeight', 'bold');
            set(self.t3_txt1, 'FontAngle', 'italic');
            set(self.t3_txt1, 'BackgroundColor', self.background_color);
            set(self.t3_txt1, 'Visible', 'off');

            % frame around applications
            self.t3_frm1 = uicontrol('style','frame');
            set(self.t3_frm1, 'unit', 'pixels', 'position', [20 370 780 190]);
            set(self.t3_frm1, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t3_frm1, 'BackgroundColor', self.background_color);
            set(self.t3_frm1, 'Visible', 'off');

            % forecast activation label
            self.t3_txt2 = uicontrol('style', 'text');
            set(self.t3_txt2, 'unit', 'pixels', 'position', [30 520 380 30]);
            set(self.t3_txt2, 'String', ' forecasts');
            set(self.t3_txt2, 'HorizontalAlignment', 'left');
            set(self.t3_txt2, 'FontName', 'Serif');
            set(self.t3_txt2, 'FontSize', 12);
            set(self.t3_txt2, 'FontWeight', 'bold');
            set(self.t3_txt2, 'BackgroundColor', self.background_color);
            set(self.t3_txt2, 'Visible', 'off');

            % forecast radiobuttons
            self.t3_bgr1 = uibuttongroup('unit','pixels', 'Position',[500 520 160 30]);
            set(self.t3_bgr1, 'BorderType', 'none');
            set(self.t3_bgr1, 'BackgroundColor', self.background_color);            
            self.t3_rdb1 = uicontrol(self.t3_bgr1,'Style','radiobutton');
            set(self.t3_rdb1, 'Position',[8 10 80 25]);
            set(self.t3_rdb1, 'String',' yes');
            set(self.t3_rdb1, 'FontName', 'Serif');
            set(self.t3_rdb1, 'FontSize', 12);
            set(self.t3_rdb1, 'FontWeight', 'bold');
            set(self.t3_rdb1, 'BackgroundColor', self.background_color);
            self.t3_rdb2 = uicontrol(self.t3_bgr1,'Style','radiobutton');
            set(self.t3_rdb2, 'Position',[88 10 80 25]);
            set(self.t3_rdb2, 'String',' no');
            set(self.t3_rdb2, 'FontName', 'Serif');
            set(self.t3_rdb2, 'FontSize', 12);
            set(self.t3_rdb2, 'FontWeight', 'bold');
            set(self.t3_rdb2, 'BackgroundColor', self.background_color);
            set(self.t3_bgr1, 'Visible', 'off');
            if self.user_inputs.tab_3.forecast
                set(self.t3_bgr1, 'SelectedObject', self.t3_rdb1);
            else
                set(self.t3_bgr1, 'SelectedObject', self.t3_rdb2);
            end
            set(self.t3_bgr1, 'SelectionChangeFcn', @self.cb_t3_bgr1);

            % forecast credibility edit
            self.t3_edt1 = uicontrol('style','edit');
            set(self.t3_edt1, 'unit', 'pixels', 'position', [690 525 70 25]);
            set(self.t3_edt1, 'HorizontalAlignment', 'center');
            set(self.t3_edt1, 'Visible', 'off');
            set(self.t3_edt1, 'String', self.user_inputs.tab_3.forecast_credibility);
            set(self.t3_edt1, 'CallBack', @self.cb_t3_edt1);

            % conditional forecast activation label
            self.t3_txt3 = uicontrol('style', 'text');
            set(self.t3_txt3, 'unit', 'pixels', 'position', [30 485 380 30]);
            set(self.t3_txt3, 'String', ' conditional forecasts');
            set(self.t3_txt3, 'HorizontalAlignment', 'left');
            set(self.t3_txt3, 'FontName', 'Serif');
            set(self.t3_txt3, 'FontSize', 12);
            set(self.t3_txt3, 'FontWeight', 'bold');
            set(self.t3_txt3, 'BackgroundColor', self.background_color);
            set(self.t3_txt3, 'Visible', 'off');

            % conditional forecast radiobuttons
            self.t3_bgr2 = uibuttongroup('unit','pixels', 'Position',[500 485 160 30]);
            set(self.t3_bgr2, 'BorderType', 'none');
            set(self.t3_bgr2, 'BackgroundColor', self.background_color);            
            self.t3_rdb3 = uicontrol(self.t3_bgr2,'Style','radiobutton');
            set(self.t3_rdb3, 'Position',[8 10 80 25]);
            set(self.t3_rdb3, 'String',' yes');
            set(self.t3_rdb3, 'FontName', 'Serif');
            set(self.t3_rdb3, 'FontSize', 12);
            set(self.t3_rdb3, 'FontWeight', 'bold');
            set(self.t3_rdb3, 'BackgroundColor', self.background_color);
            set(self.t3_rdb3, 'Enable', 'off');
            self.t3_rdb4 = uicontrol(self.t3_bgr2,'Style','radiobutton');
            set(self.t3_rdb4, 'Position',[88 10 80 25]);
            set(self.t3_rdb4, 'String',' no');
            set(self.t3_rdb4, 'FontName', 'Serif');
            set(self.t3_rdb4, 'FontSize', 12);
            set(self.t3_rdb4, 'FontWeight', 'bold');
            set(self.t3_rdb4, 'BackgroundColor', self.background_color);
            set(self.t3_rdb4, 'Enable', 'off');
            set(self.t3_bgr2, 'Visible', 'off');
            if self.user_inputs.tab_3.conditional_forecast
                set(self.t3_bgr2, 'SelectedObject', self.t3_rdb3);
            else
                set(self.t3_bgr2, 'SelectedObject', self.t3_rdb4);
            end
            set(self.t3_bgr2, 'SelectionChangeFcn', @self.cb_t3_bgr2);

            % conditional forecast credibility edit
            self.t3_edt2 = uicontrol('style','edit');
            set(self.t3_edt2, 'unit', 'pixels', 'position', [690 490 70 25]);
            set(self.t3_edt2, 'HorizontalAlignment', 'center');
            set(self.t3_edt2, 'Visible', 'off');
            set(self.t3_edt2, 'String', self.user_inputs.tab_3.conditional_forecast_credibility);
            set(self.t3_edt2, 'CallBack', @self.cb_t3_edt2);
            set(self.t3_edt2, 'Enable', 'off');

            % irf activation label
            self.t3_txt4 = uicontrol('style', 'text');
            set(self.t3_txt4, 'unit', 'pixels', 'position', [30 450 380 30]);
            set(self.t3_txt4, 'String', ' impulse response functions');
            set(self.t3_txt4, 'HorizontalAlignment', 'left');
            set(self.t3_txt4, 'FontName', 'Serif');
            set(self.t3_txt4, 'FontSize', 12);
            set(self.t3_txt4, 'FontWeight', 'bold');
            set(self.t3_txt4, 'BackgroundColor', self.background_color);
            set(self.t3_txt4, 'Visible', 'off');

            % irf radiobuttons
            self.t3_bgr3 = uibuttongroup('unit','pixels', 'Position',[500 450 160 30]);
            set(self.t3_bgr3, 'BorderType', 'none');
            set(self.t3_bgr3, 'BackgroundColor', self.background_color);            
            self.t3_rdb5 = uicontrol(self.t3_bgr3,'Style','radiobutton');
            set(self.t3_rdb5, 'Position',[8 10 80 25]);
            set(self.t3_rdb5, 'String',' yes');
            set(self.t3_rdb5, 'FontName', 'Serif');
            set(self.t3_rdb5, 'FontSize', 12);
            set(self.t3_rdb5, 'FontWeight', 'bold');
            set(self.t3_rdb5, 'BackgroundColor', self.background_color);
            set(self.t3_rdb5, 'Enable', 'off');
            self.t3_rdb6 = uicontrol(self.t3_bgr3,'Style','radiobutton');
            set(self.t3_rdb6, 'Position',[88 10 80 25]);
            set(self.t3_rdb6, 'String',' no');
            set(self.t3_rdb6, 'FontName', 'Serif');
            set(self.t3_rdb6, 'FontSize', 12);
            set(self.t3_rdb6, 'FontWeight', 'bold');
            set(self.t3_rdb6, 'BackgroundColor', self.background_color);
            set(self.t3_rdb6, 'Enable', 'off');
            set(self.t3_bgr3, 'Visible', 'off');
            if self.user_inputs.tab_3.irf
                set(self.t3_bgr3, 'SelectedObject', self.t3_rdb5);
            else
                set(self.t3_bgr3, 'SelectedObject', self.t3_rdb6);
            end
            set(self.t3_bgr3, 'SelectionChangeFcn', @self.cb_t3_bgr3);

            % irf credibility edit
            self.t3_edt3 = uicontrol('style','edit');
            set(self.t3_edt3, 'unit', 'pixels', 'position', [690 455 70 25]);
            set(self.t3_edt3, 'HorizontalAlignment', 'center');
            set(self.t3_edt3, 'Visible', 'off');
            set(self.t3_edt3, 'String', self.user_inputs.tab_3.irf_credibility);
            set(self.t3_edt3, 'CallBack', @self.cb_t3_edt3);
            set(self.t3_edt3, 'Enable', 'off');

            % fevd activation label
            self.t3_txt5 = uicontrol('style', 'text');
            set(self.t3_txt5, 'unit', 'pixels', 'position', [30 415 380 30]);
            set(self.t3_txt5, 'String', ' forecast error variance decomposition');
            set(self.t3_txt5, 'HorizontalAlignment', 'left');
            set(self.t3_txt5, 'FontName', 'Serif');
            set(self.t3_txt5, 'FontSize', 12);
            set(self.t3_txt5, 'FontWeight', 'bold');
            set(self.t3_txt5, 'BackgroundColor', self.background_color);
            set(self.t3_txt5, 'Visible', 'off');

            % fevd radiobuttons
            self.t3_bgr4 = uibuttongroup('unit','pixels', 'Position',[500 415 160 30]);
            set(self.t3_bgr4, 'BorderType', 'none');
            set(self.t3_bgr4, 'BackgroundColor', self.background_color);            
            self.t3_rdb7 = uicontrol(self.t3_bgr4,'Style','radiobutton');
            set(self.t3_rdb7, 'Position',[8 10 80 25]);
            set(self.t3_rdb7, 'String',' yes');
            set(self.t3_rdb7, 'FontName', 'Serif');
            set(self.t3_rdb7, 'FontSize', 12);
            set(self.t3_rdb7, 'FontWeight', 'bold');
            set(self.t3_rdb7, 'BackgroundColor', self.background_color);
            set(self.t3_rdb7, 'Enable', 'off');
            self.t3_rdb8 = uicontrol(self.t3_bgr4,'Style','radiobutton');
            set(self.t3_rdb8, 'Position',[88 10 80 25]);
            set(self.t3_rdb8, 'String',' no');
            set(self.t3_rdb8, 'FontName', 'Serif');
            set(self.t3_rdb8, 'FontSize', 12);
            set(self.t3_rdb8, 'FontWeight', 'bold');
            set(self.t3_rdb8, 'BackgroundColor', self.background_color);
            set(self.t3_rdb8, 'Enable', 'off');
            set(self.t3_bgr4, 'Visible', 'off');
            if self.user_inputs.tab_3.fevd
                set(self.t3_bgr4, 'SelectedObject', self.t3_rdb7);
            else
                set(self.t3_bgr4, 'SelectedObject', self.t3_rdb8);
            end
            set(self.t3_bgr4, 'SelectionChangeFcn', @self.cb_t3_bgr4);

            % fevd credibility edit
            self.t3_edt4 = uicontrol('style','edit');
            set(self.t3_edt4, 'unit', 'pixels', 'position', [690 420 70 25]);
            set(self.t3_edt4, 'HorizontalAlignment', 'center');
            set(self.t3_edt4, 'Visible', 'off');
            set(self.t3_edt4, 'String', self.user_inputs.tab_3.fevd_credibility);
            set(self.t3_edt4, 'CallBack', @self.cb_t3_edt4);
            set(self.t3_edt4, 'Enable', 'off');

            % historical decomposition activation label
            self.t3_txt6 = uicontrol('style', 'text');
            set(self.t3_txt6, 'unit', 'pixels', 'position', [30 380 380 30]);
            set(self.t3_txt6, 'String', ' historical decomposition');
            set(self.t3_txt6, 'HorizontalAlignment', 'left');
            set(self.t3_txt6, 'FontName', 'Serif');
            set(self.t3_txt6, 'FontSize', 12);
            set(self.t3_txt6, 'FontWeight', 'bold');
            set(self.t3_txt6, 'BackgroundColor', self.background_color);
            set(self.t3_txt6, 'Visible', 'off');

            % historical decomposition radiobuttons
            self.t3_bgr5 = uibuttongroup('unit','pixels', 'Position',[500 380 160 30]);
            set(self.t3_bgr5, 'BorderType', 'none');
            set(self.t3_bgr5, 'BackgroundColor', self.background_color);            
            self.t3_rdb9 = uicontrol(self.t3_bgr5,'Style','radiobutton');
            set(self.t3_rdb9, 'Position',[8 10 80 25]);
            set(self.t3_rdb9, 'String',' yes');
            set(self.t3_rdb9, 'FontName', 'Serif');
            set(self.t3_rdb9, 'FontSize', 12);
            set(self.t3_rdb9, 'FontWeight', 'bold');
            set(self.t3_rdb9, 'BackgroundColor', self.background_color);
            set(self.t3_rdb9, 'Enable', 'off');
            self.t3_rdb10 = uicontrol(self.t3_bgr5,'Style','radiobutton');
            set(self.t3_rdb10, 'Position',[88 10 80 25]);
            set(self.t3_rdb10, 'String',' no');
            set(self.t3_rdb10, 'FontName', 'Serif');
            set(self.t3_rdb10, 'FontSize', 12);
            set(self.t3_rdb10, 'FontWeight', 'bold');
            set(self.t3_rdb10, 'BackgroundColor', self.background_color);
            set(self.t3_rdb10, 'Enable', 'off');
            set(self.t3_bgr5, 'Visible', 'off');
            if self.user_inputs.tab_3.hd
                set(self.t3_bgr5, 'SelectedObject', self.t3_rdb9);
            else
                set(self.t3_bgr5, 'SelectedObject', self.t3_rdb10);
            end
            set(self.t3_bgr5, 'SelectionChangeFcn', @self.cb_t3_bgr5);

            % historical decomposition edit
            self.t3_edt5 = uicontrol('style','edit');
            set(self.t3_edt5, 'unit', 'pixels', 'position', [690 385 70 25]);
            set(self.t3_edt5, 'HorizontalAlignment', 'center');
            set(self.t3_edt5, 'Visible', 'off');
            set(self.t3_edt5, 'String', self.user_inputs.tab_3.hd_credibility);
            set(self.t3_edt5, 'CallBack', @self.cb_t3_edt5);
            set(self.t3_edt5, 'Enable', 'off');

            % application settings label
            self.t3_txt7 = uicontrol('style', 'text');
            set(self.t3_txt7, 'unit', 'pixels', 'position', [30 320 300 30]);
            set(self.t3_txt7, 'String', ' Application settings');
            set(self.t3_txt7, 'HorizontalAlignment', 'left');
            set(self.t3_txt7, 'FontName', 'Serif');
            set(self.t3_txt7, 'FontSize', 16);
            set(self.t3_txt7, 'FontWeight', 'bold');
            set(self.t3_txt7, 'FontAngle', 'italic');
            set(self.t3_txt7, 'BackgroundColor', self.background_color); 
            set(self.t3_txt7, 'Visible', 'off');

            % frame around options
            self.t3_frm2 = uicontrol('style','frame');
            set(self.t3_frm2, 'unit', 'pixels', 'position', [20 20 780 300]);
            set(self.t3_frm2, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t3_frm2, 'BackgroundColor', self.background_color);
            set(self.t3_frm2, 'Visible', 'off');

            % periods option label
            self.t3_txt8 = uicontrol('style', 'text');
            set(self.t3_txt8, 'unit', 'pixels', 'position', [30 275 300 30]);
            set(self.t3_txt8, 'String', ' Periods');
            set(self.t3_txt8, 'HorizontalAlignment', 'left');
            set(self.t3_txt8, 'FontName', 'Serif');
            set(self.t3_txt8, 'FontSize', 14);
            set(self.t3_txt8, 'FontAngle', 'italic');
            set(self.t3_txt8, 'BackgroundColor', self.background_color); 
            set(self.t3_txt8, 'Visible', 'off');

            % forecast periods label
            self.t3_txt9 = uicontrol('style', 'text');
            set(self.t3_txt9, 'unit', 'pixels', 'position', [30 235 200 30]);
            set(self.t3_txt9, 'String', ' forecast periods');
            set(self.t3_txt9, 'HorizontalAlignment', 'left');
            set(self.t3_txt9, 'FontName', 'Serif');
            set(self.t3_txt9, 'FontSize', 12);
            set(self.t3_txt9, 'FontWeight', 'bold');
            set(self.t3_txt9, 'BackgroundColor', self.background_color);
            set(self.t3_txt9, 'Visible', 'off');

            % forecast periods edit
            self.t3_edt6 = uicontrol('style','edit');
            set(self.t3_edt6, 'unit', 'pixels', 'position', [300 245 70 25]);
            set(self.t3_edt6, 'HorizontalAlignment', 'center');
            set(self.t3_edt6, 'Visible', 'off');
            set(self.t3_edt6, 'String', self.user_inputs.tab_3.forecast_periods);
            set(self.t3_edt6, 'CallBack', @self.cb_t3_edt6); 
            set(self.t3_edt6, 'Enable', 'off');
            
            % irf periods label
            self.t3_txt10 = uicontrol('style', 'text');
            set(self.t3_txt10, 'unit', 'pixels', 'position', [30 200 200 30]);
            set(self.t3_txt10, 'String', ' IRF periods');
            set(self.t3_txt10, 'HorizontalAlignment', 'left');
            set(self.t3_txt10, 'FontName', 'Serif');
            set(self.t3_txt10, 'FontSize', 12);
            set(self.t3_txt10, 'FontWeight', 'bold');
            set(self.t3_txt10, 'BackgroundColor', self.background_color); 
            set(self.t3_txt10, 'Visible', 'off');

            % irf periods edit
            self.t3_edt7 = uicontrol('style','edit');
            set(self.t3_edt7, 'unit', 'pixels', 'position', [300 210 70 25]);
            set(self.t3_edt7, 'HorizontalAlignment', 'center');
            set(self.t3_edt7, 'Visible', 'off');
            set(self.t3_edt7, 'String', self.user_inputs.tab_3.irf_periods);
            set(self.t3_edt7, 'CallBack', @self.cb_t3_edt7); 
            set(self.t3_edt7, 'Enable', 'off');
            
            % identification scheme label
            self.t3_txt11 = uicontrol('style', 'text');
            set(self.t3_txt11, 'unit', 'pixels', 'position', [30 155 300 30]);
            set(self.t3_txt11, 'String', ' Identification scheme');
            set(self.t3_txt11, 'HorizontalAlignment', 'left');
            set(self.t3_txt11, 'FontName', 'Serif');
            set(self.t3_txt11, 'FontSize', 14);
            set(self.t3_txt11, 'FontAngle', 'italic');
            set(self.t3_txt11, 'BackgroundColor', self.background_color); 
            set(self.t3_txt11, 'Visible', 'off');      
            
            % conditional forecast type label
            self.t3_txt12 = uicontrol('style', 'text');
            set(self.t3_txt12, 'unit', 'pixels', 'position', [30 115 250 30]);
            set(self.t3_txt12, 'String', ' conditional forecast type');
            set(self.t3_txt12, 'HorizontalAlignment', 'left');
            set(self.t3_txt12, 'FontName', 'Serif');
            set(self.t3_txt12, 'FontSize', 12);
            set(self.t3_txt12, 'FontWeight', 'bold');
            set(self.t3_txt12, 'BackgroundColor', self.background_color); 
            set(self.t3_txt12, 'Visible', 'off');

            % conditional forecast type menu
            self.t3_mnu1 = uicontrol('style', 'popupmenu');
            set(self.t3_mnu1, 'position',[35 90 200 30]);            
            set(self.t3_mnu1, 'String', {'1. general', '2. all shocks', '3. shock-specific'}); 
            set(self.t3_mnu1, 'Visible', 'off');
            set(self.t3_mnu1, 'Value', self.user_inputs.tab_3.conditional_forecast_type);
            set(self.t3_mnu1, 'CallBack', @self.cb_t3_mnu1);
            set(self.t3_mnu1, 'Enable', 'off');

            % structural identification label
            self.t3_txt13 = uicontrol('style', 'text');
            set(self.t3_txt13, 'unit', 'pixels', 'position', [30 55 250 30]);
            set(self.t3_txt13, 'String', ' structural identification');
            set(self.t3_txt13, 'HorizontalAlignment', 'left');
            set(self.t3_txt13, 'FontName', 'Serif');
            set(self.t3_txt13, 'FontSize', 12);
            set(self.t3_txt13, 'FontWeight', 'bold');
            set(self.t3_txt13, 'BackgroundColor', self.background_color); 
            set(self.t3_txt13, 'Visible', 'off');

            % structural identification menu
            self.t3_mnu2 = uicontrol('style', 'popupmenu');
            set(self.t3_mnu2, 'position',[35 30 200 30]);
            set(self.t3_mnu2, 'String', {'1. none', '2. Cholesky', '3. triangular'}); 
            set(self.t3_mnu2, 'Visible', 'off');
            set(self.t3_mnu2, 'Value', self.user_inputs.tab_3.structural_identification);
            set(self.t3_mnu2, 'CallBack', @self.cb_t3_mnu2);
            set(self.t3_mnu2, 'Enable', 'off');
            
            % options label
            self.t3_txt14 = uicontrol('style', 'text');
            set(self.t3_txt14, 'unit', 'pixels', 'position', [430 275 300 30]);
            set(self.t3_txt14, 'String', ' Options');
            set(self.t3_txt14, 'HorizontalAlignment', 'left');
            set(self.t3_txt14, 'FontName', 'Serif');
            set(self.t3_txt14, 'FontSize', 14);
            set(self.t3_txt14, 'FontAngle', 'italic');
            set(self.t3_txt14, 'BackgroundColor', self.background_color); 
            set(self.t3_txt14, 'Visible', 'off');            
            
            % forecast evaluation label
            self.t3_txt15 = uicontrol('style', 'text');
            set(self.t3_txt15, 'unit', 'pixels', 'position', [430 235 200 30]);
            set(self.t3_txt15, 'String', ' forecast evaluation');
            set(self.t3_txt15, 'HorizontalAlignment', 'left');
            set(self.t3_txt15, 'FontName', 'Serif');
            set(self.t3_txt15, 'FontSize', 12);
            set(self.t3_txt15, 'FontWeight', 'bold');
            set(self.t3_txt15, 'BackgroundColor', self.background_color); 
            set(self.t3_txt15, 'Visible', 'off');

            % forecast evaluation checkbox
            self.t3_cbx1 = uicontrol('style', 'checkbox');
            set(self.t3_cbx1, 'unit', 'pixels', 'position', [745 250 20 20]);
            set(self.t3_cbx1, 'BackgroundColor', self.background_color);
            set(self.t3_cbx1, 'Visible', 'off');
            set(self.t3_cbx1, 'Value', self.user_inputs.tab_3.forecast_evaluation);
            set(self.t3_cbx1, 'CallBack', @self.cb_t3_cbx1);

            % input files label
            self.t3_txt16 = uicontrol('style', 'text');
            set(self.t3_txt16, 'unit', 'pixels', 'position', [430 155 300 30]);
            set(self.t3_txt16, 'String', ' Input files');
            set(self.t3_txt16, 'HorizontalAlignment', 'left');
            set(self.t3_txt16, 'FontName', 'Serif');
            set(self.t3_txt16, 'FontSize', 14);
            set(self.t3_txt16, 'FontAngle', 'italic');
            set(self.t3_txt16, 'BackgroundColor', self.background_color); 
            set(self.t3_txt16, 'Visible', 'off');                

            % forecast input file label
            self.t3_txt17 = uicontrol('style', 'text');
            set(self.t3_txt17, 'unit', 'pixels', 'position', [430 115 250 30]);
            set(self.t3_txt17, 'String', ' forecast input file');
            set(self.t3_txt17, 'HorizontalAlignment', 'left');
            set(self.t3_txt17, 'FontName', 'Serif');
            set(self.t3_txt17, 'FontSize', 12);
            set(self.t3_txt17, 'FontWeight', 'bold');
            set(self.t3_txt17, 'BackgroundColor', self.background_color);
            set(self.t3_txt17, 'Visible', 'off');

            % forecast file edit
            self.t3_edt8 = uicontrol('style','edit');
            set(self.t3_edt8, 'unit', 'pixels', 'position', [435 95 325 25]);
            set(self.t3_edt8, 'HorizontalAlignment', 'left');
            set(self.t3_edt8, 'Visible', 'off');
            set(self.t3_edt8, 'String', self.user_inputs.tab_3.forecast_file);
            set(self.t3_edt8, 'CallBack', @self.cb_t3_edt8);            

            % structural identification file label
            self.t3_txt18 = uicontrol('style', 'text');
            set(self.t3_txt18, 'unit', 'pixels', 'position', [430 55 300 30]);
            set(self.t3_txt18, 'String', ' structural identification file');
            set(self.t3_txt18, 'HorizontalAlignment', 'left');
            set(self.t3_txt18, 'FontName', 'Serif');
            set(self.t3_txt18, 'FontSize', 12);
            set(self.t3_txt18, 'FontWeight', 'bold');
            set(self.t3_txt18, 'BackgroundColor', self.background_color);
            set(self.t3_txt18, 'Visible', 'off');

            % structural identification file edit
            self.t3_edt9 = uicontrol('style','edit');
            set(self.t3_edt9, 'unit', 'pixels', 'position', [435 35 325 25]);
            set(self.t3_edt9, 'HorizontalAlignment', 'left');
            set(self.t3_edt9, 'Visible', 'off');
            set(self.t3_edt9, 'String', self.user_inputs.tab_3.structural_identification_file);
            set(self.t3_edt9, 'CallBack', @self.cb_t3_edt9);
            set(self.t3_edt9, 'Enable', 'off');

            % run pushbutton
            self.t3_pbt1 = uicontrol('style', 'pushbutton');
            set(self.t3_pbt1, 'unit', 'pixels', 'position', [820 50 160 260]);
            run_image = imread(fullfile(self.interface_path, 'run_button.png'));
            set(self.t3_pbt1,'cdata', run_image);
            set(self.t3_pbt1, 'Visible', 'off');
            set(self.t3_pbt1, 'CallBack', @self.cb_t3_pbt1);
        end
        
        
        function hide_tab_3(self)

            % hide all controls   
            set(self.t3_txt1, 'Visible', 'off');
            set(self.t3_frm1, 'Visible', 'off');
            set(self.t3_txt2, 'Visible', 'off');
            set(self.t3_bgr1, 'Visible', 'off');
            set(self.t3_rdb1, 'Visible', 'off');
            set(self.t3_rdb2, 'Visible', 'off');
            set(self.t3_edt1, 'Visible', 'off');
            set(self.t3_txt3, 'Visible', 'off');
            set(self.t3_bgr2, 'Visible', 'off');
            set(self.t3_rdb3, 'Visible', 'off');
            set(self.t3_rdb4, 'Visible', 'off');
            set(self.t3_edt2, 'Visible', 'off');
            set(self.t3_txt4, 'Visible', 'off');
            set(self.t3_bgr3, 'Visible', 'off');
            set(self.t3_rdb5, 'Visible', 'off');
            set(self.t3_rdb6, 'Visible', 'off');
            set(self.t3_edt3, 'Visible', 'off');
            set(self.t3_txt5, 'Visible', 'off');
            set(self.t3_bgr4, 'Visible', 'off');
            set(self.t3_rdb7, 'Visible', 'off');
            set(self.t3_rdb8, 'Visible', 'off');
            set(self.t3_edt4, 'Visible', 'off');
            set(self.t3_txt6, 'Visible', 'off');
            set(self.t3_bgr5, 'Visible', 'off');
            set(self.t3_rdb9, 'Visible', 'off');
            set(self.t3_rdb10, 'Visible', 'off');
            set(self.t3_edt5, 'Visible', 'off');
            set(self.t3_txt7, 'Visible', 'off');
            set(self.t3_frm2, 'Visible', 'off');
            set(self.t3_txt8, 'Visible', 'off');
            set(self.t3_txt9, 'Visible', 'off');
            set(self.t3_edt6, 'Visible', 'off');
            set(self.t3_txt10, 'Visible', 'off');
            set(self.t3_mnu1, 'Visible', 'off');
            set(self.t3_txt11, 'Visible', 'off');
            set(self.t3_edt7, 'Visible', 'off');
            set(self.t3_txt12, 'Visible', 'off');
            set(self.t3_cbx1, 'Visible', 'off');
            set(self.t3_txt13, 'Visible', 'off');
            set(self.t3_txt14, 'Visible', 'off');
            set(self.t3_edt8, 'Visible', 'off');
            set(self.t3_txt15, 'Visible', 'off');
            set(self.t3_mnu2, 'Visible', 'off');
            set(self.t3_txt16, 'Visible', 'off');
            set(self.t3_edt9, 'Visible', 'off');
            set(self.t3_pbt1, 'Visible', 'off');
            set(self.t3_txt17, 'Visible', 'off');
            set(self.t3_txt18, 'Visible', 'off');

            % update tab color
            set(self.tab_pbt3, 'BackgroundColor', self.backtabs_color); 
        end
        
        
        function show_tab_3(self)

            % show all controls
            set(self.t3_txt1, 'Visible', 'on');
            set(self.t3_frm1, 'Visible', 'on');
            set(self.t3_txt2, 'Visible', 'on');
            set(self.t3_bgr1, 'Visible', 'on');
            set(self.t3_rdb1, 'Visible', 'on');
            set(self.t3_rdb2, 'Visible', 'on');
            set(self.t3_edt1, 'Visible', 'on');
            set(self.t3_txt3, 'Visible', 'on');
            set(self.t3_bgr2, 'Visible', 'on');
            set(self.t3_rdb3, 'Visible', 'on');
            set(self.t3_rdb4, 'Visible', 'on');
            set(self.t3_edt2, 'Visible', 'on');
            set(self.t3_txt4, 'Visible', 'on');
            set(self.t3_bgr3, 'Visible', 'on');
            set(self.t3_rdb5, 'Visible', 'on');
            set(self.t3_rdb6, 'Visible', 'on');
            set(self.t3_edt3, 'Visible', 'on');
            set(self.t3_txt5, 'Visible', 'on');
            set(self.t3_bgr4, 'Visible', 'on');
            set(self.t3_rdb7, 'Visible', 'on');
            set(self.t3_rdb8, 'Visible', 'on');
            set(self.t3_edt4, 'Visible', 'on');
            set(self.t3_txt6, 'Visible', 'on');
            set(self.t3_bgr5, 'Visible', 'on');
            set(self.t3_rdb9, 'Visible', 'on');
            set(self.t3_rdb10, 'Visible', 'on');
            set(self.t3_edt5, 'Visible', 'on');
            set(self.t3_txt7, 'Visible', 'on');
            set(self.t3_frm2, 'Visible', 'on');
            set(self.t3_txt8, 'Visible', 'on');
            set(self.t3_txt9, 'Visible', 'on');
            set(self.t3_edt6, 'Visible', 'on');
            set(self.t3_txt10, 'Visible', 'on');
            set(self.t3_mnu1, 'Visible', 'on');
            set(self.t3_txt11, 'Visible', 'on');
            set(self.t3_edt7, 'Visible', 'on');
            set(self.t3_txt12, 'Visible', 'on');
            set(self.t3_cbx1, 'Visible', 'on');
            set(self.t3_txt13, 'Visible', 'on');
            set(self.t3_txt14, 'Visible', 'on');
            set(self.t3_edt8, 'Visible', 'on');
            set(self.t3_txt15, 'Visible', 'on');
            set(self.t3_mnu2, 'Visible', 'on');
            set(self.t3_txt16, 'Visible', 'on');
            set(self.t3_edt9, 'Visible', 'on');
            set(self.t3_pbt1, 'Visible', 'on');
            set(self.t3_txt17, 'Visible', 'on');
            set(self.t3_txt18, 'Visible', 'on');
        end


        function cb_t3_bgr1(self, hObject, callbackdata)
           if get(self.t3_bgr1, 'SelectedObject') == self.t3_rdb1
               self.user_inputs.tab_3.forecast = true;
           elseif get(self.t3_bgr1, 'SelectedObject') == self.t3_rdb2
               self.user_inputs.tab_3.forecast = false;
           end
        end 

        
        function cb_t3_edt1(self, hObject, callbackdata)
            self.user_inputs.tab_3.forecast_credibility = get(self.t3_edt1, 'String');
        end 
        
        
        function cb_t3_bgr2(self, hObject, callbackdata)
           if get(self.t3_bgr2, 'SelectedObject') == self.t3_rdb3
               self.user_inputs.tab_3.conditional_forecast = true;
           elseif get(self.t3_bgr2, 'SelectedObject') == self.t3_rdb4
               self.user_inputs.tab_3.conditional_forecast = false;
           end
        end
        
        
        function cb_t3_edt2(self, hObject, callbackdata)
            self.user_inputs.tab_3.conditional_forecast_credibility = get(self.t3_edt2, 'String');
        end
        
        
        function cb_t3_bgr3(self, hObject, callbackdata)
           if get(self.t3_bgr3, 'SelectedObject') == self.t3_rdb5
               self.user_inputs.tab_3.irf = true;
           elseif get(self.t3_bgr3, 'SelectedObject') == self.t3_rdb6
               self.user_inputs.tab_3.irf = false;
           end
        end
        
        
        function cb_t3_edt3(self, hObject, callbackdata)
            self.user_inputs.tab_3.irf_credibility = get(self.t3_edt3, 'String');
        end
        
        
        function cb_t3_bgr4(self, hObject, callbackdata)
           if get(self.t3_bgr4, 'SelectedObject') == self.t3_rdb7
               self.user_inputs.tab_3.fevd = true;
           elseif get(self.t3_bgr4, 'SelectedObject') == self.t3_rdb8
               self.user_inputs.tab_3.fevd = false;
           end
        end
        
        
        function cb_t3_edt4(self, hObject, callbackdata)
            self.user_inputs.tab_3.fevd_credibility = get(self.t3_edt4, 'String');
        end
        
        
        function cb_t3_bgr5(self, hObject, callbackdata)
           if get(self.t3_bgr5, 'SelectedObject') == self.t3_rdb9
               self.user_inputs.tab_3.hd = true;
           elseif get(self.t3_bgr5, 'SelectedObject') == self.t3_rdb10
               self.user_inputs.tab_3.hd = false;
           end
        end
        
        
        function cb_t3_edt5(self, hObject, callbackdata)
            self.user_inputs.tab_3.hd_credibility = get(self.t3_edt5, 'String');
        end
        
        
        function cb_t3_edt6(self, hObject, callbackdata)
            self.user_inputs.tab_3.forecast_periods = get(self.t3_edt6, 'String');
        end
        
        
        function cb_t3_mnu1(self, hObject, callbackdata)
            self.user_inputs.tab_3.conditional_forecast_type = get(self.t3_mnu1, 'Value');
        end
        
        
        function cb_t3_edt7(self, hObject, callbackdata)
            self.user_inputs.tab_3.irf_periods = get(self.t3_edt7, 'String');
        end
        
        
        function cb_t3_cbx1(self, hObject, callbackdata)
            self.user_inputs.tab_3.forecast_evaluation = logical(get(self.t3_cbx1, 'Value'));
        end
        
        
        function cb_t3_edt8(self, hObject, callbackdata)
            self.user_inputs.tab_3.forecast_file = get(self.t3_edt8, 'String');
        end
        
        
        function cb_t3_mnu2(self, hObject, callbackdata)
            self.user_inputs.tab_3.structural_identification = get(self.t3_mnu2, 'Value');
        end
        
        
        function cb_t3_edt9(self, hObject, callbackdata)
            self.user_inputs.tab_3.structural_identification_file = get(self.t3_edt9, 'String');
        end
        
        
        function cb_t3_pbt1(self, hObject, callbackdata)
            self.validate_interface();
        end 
        
        
    end
    
    
end
        