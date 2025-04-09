classdef Tab2VectorAutoregressionInterface < handle

    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess = public)
        % tab 2 properties (vector autoregression)
        t2_var_txt1
        t2_var_txt2        
        t2_var_txt3
        t2_var_txt4        
        t2_var_txt5        
        t2_var_txt6
        t2_var_txt7
        t2_var_txt8        
        t2_var_txt9   
        t2_var_txt10
        t2_var_txt11
        t2_var_txt12
        t2_var_txt13
        t2_var_txt14
        t2_var_txt15
        t2_var_txt16
        t2_var_txt17
        t2_var_txt18
        t2_var_txt19
        t2_var_txt20
        t2_var_txt21
        t2_var_txt22
        t2_var_txt23        
        t2_var_txt24
        t2_var_txt25        
        t2_var_txt26        
        t2_var_txt27
        t2_var_txt28
        t2_var_txt29
        t2_var_txt30
        t2_var_txt31
        t2_var_txt32
        t2_var_txt33
        t2_var_txt34        
        t2_var_txt35
        t2_var_txt36
        t2_var_txt37     
        t2_var_txt38
        t2_var_frm1
        t2_var_frm2        
        t2_var_frm3
        t2_var_frm4        
        t2_var_bgr1
        t2_var_bgr2        
        t2_var_bgr3        
        t2_var_bgr4
        t2_var_bgr5
        t2_var_bgr6
        t2_var_bgr7
        t2_var_bgr8
        t2_var_bgr9
        t2_var_bgr10
        t2_var_rdb1
        t2_var_rdb2
        t2_var_rdb3
        t2_var_rdb4
        t2_var_rdb5
        t2_var_rdb6
        t2_var_rdb7
        t2_var_rdb8
        t2_var_rdb9
        t2_var_rdb10
        t2_var_rdb11
        t2_var_rdb12
        t2_var_rdb13
        t2_var_rdb14
        t2_var_rdb15
        t2_var_rdb16
        t2_var_rdb17
        t2_var_rdb18
        t2_var_rdb19
        t2_var_rdb20
        t2_var_rdb21
        t2_var_rdb22
        t2_var_rdb23
        t2_var_rdb24
        t2_var_rdb25
        t2_var_edt1
        t2_var_edt2
        t2_var_edt3
        t2_var_edt4
        t2_var_edt5
        t2_var_edt6
        t2_var_edt7
        t2_var_edt8
        t2_var_edt9
        t2_var_edt10
        t2_var_edt11
        t2_var_edt12
        t2_var_edt13
        t2_var_edt14
        t2_var_edt15
        t2_var_edt16         
        t2_var_cbx1
        t2_var_cbx2
        t2_var_cbx3
    end
    

    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)    
    
    
        function self = Tab2VectorAutoregressionInterface()
        end

        
        function create_tab_2_var(self)
    
            % vector autoregression label
            self.t2_var_txt1 = uicontrol('style', 'text');
            set(self.t2_var_txt1, 'unit', 'pixels', 'position', [30 560 300 30]);
            set(self.t2_var_txt1, 'String', ' VAR type');
            set(self.t2_var_txt1, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt1, 'FontName', 'Serif');
            set(self.t2_var_txt1, 'FontSize', 16);
            set(self.t2_var_txt1, 'FontWeight', 'bold');
            set(self.t2_var_txt1, 'FontAngle', 'italic');
            set(self.t2_var_txt1, 'BackgroundColor', self.background_color);  
            set(self.t2_var_txt1, 'Visible', 'off');

            % frame around VAR type
            self.t2_var_frm1 = uicontrol('style','frame');
            set(self.t2_var_frm1, 'unit', 'pixels', 'position', [20 450 470 110]);
            set(self.t2_var_frm1, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_var_frm1, 'BackgroundColor', self.background_color);
            set(self.t2_var_frm1, 'Visible', 'off');
            
            % VAR type radiobuttons
            self.t2_var_bgr1 = uibuttongroup('unit','pixels', 'Position',[25 455 440 100]);
            set(self.t2_var_bgr1, 'BorderType', 'none');
            set(self.t2_var_bgr1, 'BackgroundColor', self.background_color); 
            self.t2_var_rdb1 = uicontrol(self.t2_var_bgr1,'Style','radiobutton');
            set(self.t2_var_rdb1, 'Position',[5 75 220 25]);
            set(self.t2_var_rdb1, 'String',' maximum likelihood');
            set(self.t2_var_rdb1, 'FontName', 'Serif');
            set(self.t2_var_rdb1, 'FontSize', 12);
            set(self.t2_var_rdb1, 'FontWeight', 'bold');
            set(self.t2_var_rdb1, 'BackgroundColor', self.background_color);
            self.t2_var_rdb2 = uicontrol(self.t2_var_bgr1,'Style','radiobutton');
            set(self.t2_var_rdb2, 'Position',[235 75 220 25]);
            set(self.t2_var_rdb2, 'String',' Minnesota');
            set(self.t2_var_rdb2, 'FontName', 'Serif');
            set(self.t2_var_rdb2, 'FontSize', 12);
            set(self.t2_var_rdb2, 'FontWeight', 'bold');
            set(self.t2_var_rdb2, 'BackgroundColor', self.background_color);
            self.t2_var_rdb3 = uicontrol(self.t2_var_bgr1,'Style','radiobutton');
            set(self.t2_var_rdb3, 'Position',[5 50 220 25]);
            set(self.t2_var_rdb3, 'String',' normal-Wishart');
            set(self.t2_var_rdb3, 'FontName', 'Serif');
            set(self.t2_var_rdb3, 'FontSize', 12);
            set(self.t2_var_rdb3, 'FontWeight', 'bold');
            set(self.t2_var_rdb3, 'BackgroundColor', self.background_color);
            self.t2_var_rdb4 = uicontrol(self.t2_var_bgr1,'Style','radiobutton');
            set(self.t2_var_rdb4, 'Position',[235 50 220 25]);
            set(self.t2_var_rdb4, 'String',' independent');
            set(self.t2_var_rdb4, 'FontName', 'Serif');
            set(self.t2_var_rdb4, 'FontSize', 12);
            set(self.t2_var_rdb4, 'FontWeight', 'bold');
            set(self.t2_var_rdb4, 'BackgroundColor', self.background_color);
            self.t2_var_rdb5 = uicontrol(self.t2_var_bgr1,'Style','radiobutton');
            set(self.t2_var_rdb5, 'Position',[5 25 220 25]);
            set(self.t2_var_rdb5, 'String',' dummy observations');
            set(self.t2_var_rdb5, 'FontName', 'Serif');
            set(self.t2_var_rdb5, 'FontSize', 12);
            set(self.t2_var_rdb5, 'FontWeight', 'bold');
            set(self.t2_var_rdb5, 'BackgroundColor', self.background_color);
            self.t2_var_rdb6 = uicontrol(self.t2_var_bgr1,'Style','radiobutton');
            set(self.t2_var_rdb6, 'Position',[235 25 250 25]);
            set(self.t2_var_rdb6, 'String',' large Bayesian VAR');
            set(self.t2_var_rdb6, 'FontName', 'Serif');
            set(self.t2_var_rdb6, 'FontSize', 12);
            set(self.t2_var_rdb6, 'FontWeight', 'bold');
            set(self.t2_var_rdb6, 'BackgroundColor', self.background_color);
            self.t2_var_rdb7 = uicontrol(self.t2_var_bgr1,'Style','radiobutton');
            set(self.t2_var_rdb7, 'Position',[5 0 220 25]);
            set(self.t2_var_rdb7, 'String',' proxy-SVAR');
            set(self.t2_var_rdb7, 'FontName', 'Serif');
            set(self.t2_var_rdb7, 'FontSize', 12);
            set(self.t2_var_rdb7, 'FontWeight', 'bold');
            set(self.t2_var_rdb7, 'BackgroundColor', self.background_color);
            set(self.t2_var_bgr1, 'Visible', 'off');
            if self.user_inputs.tab_2_var.var_type == 1
                set(self.t2_var_bgr1, 'SelectedObject', self.t2_var_rdb1);
            elseif self.user_inputs.tab_2_var.var_type == 2
                set(self.t2_var_bgr1, 'SelectedObject', self.t2_var_rdb2);
            elseif self.user_inputs.tab_2_var.var_type == 3
                set(self.t2_var_bgr1, 'SelectedObject', self.t2_var_rdb3);
            elseif self.user_inputs.tab_2_var.var_type == 4
                set(self.t2_var_bgr1, 'SelectedObject', self.t2_var_rdb4);
            elseif self.user_inputs.tab_2_var.var_type == 5
                set(self.t2_var_bgr1, 'SelectedObject', self.t2_var_rdb5);
            elseif self.user_inputs.tab_2_var.var_type == 6
                set(self.t2_var_bgr1, 'SelectedObject', self.t2_var_rdb6);
            elseif self.user_inputs.tab_2_var.var_type == 7
                set(self.t2_var_bgr1, 'SelectedObject', self.t2_var_rdb7); 
            end
            set(self.t2_var_bgr1, 'SelectionChangeFcn', @self.cb_t2_var_bgr1);
            
            % estimation label
            self.t2_var_txt2 = uicontrol('style', 'text');
            set(self.t2_var_txt2, 'unit', 'pixels', 'position', [520 560 300 30]);
            set(self.t2_var_txt2, 'String', ' Estimation');
            set(self.t2_var_txt2, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt2, 'FontName', 'Serif');
            set(self.t2_var_txt2, 'FontSize', 16);
            set(self.t2_var_txt2, 'FontWeight', 'bold');
            set(self.t2_var_txt2, 'FontAngle', 'italic');
            set(self.t2_var_txt2, 'BackgroundColor', self.background_color);  
            set(self.t2_var_txt2, 'Visible', 'off');
            
            % frame around estimation
            self.t2_var_frm2 = uicontrol('style','frame');
            set(self.t2_var_frm2, 'unit', 'pixels', 'position', [510 450 470 110]);
            set(self.t2_var_frm2, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_var_frm2, 'BackgroundColor', self.background_color);
            set(self.t2_var_frm2, 'Visible', 'off');
            
            % Gibbs sampling label
            self.t2_var_txt3 = uicontrol('style', 'text');
            set(self.t2_var_txt3, 'unit', 'pixels', 'position', [520 525 200 30]);
            set(self.t2_var_txt3, 'String', ' Gibbs sampling');
            set(self.t2_var_txt3, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt3, 'FontName', 'Serif');
            set(self.t2_var_txt3, 'FontSize', 14);
            set(self.t2_var_txt3, 'FontAngle', 'italic');
            set(self.t2_var_txt3, 'BackgroundColor', self.background_color);
            set(self.t2_var_txt3, 'Visible', 'off');            
            
            % iteration label
            self.t2_var_txt4 = uicontrol('style', 'text');
            set(self.t2_var_txt4, 'unit', 'pixels', 'position', [520 500 200 25]);
            set(self.t2_var_txt4, 'String', ' iterations');
            set(self.t2_var_txt4, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt4, 'FontName', 'Serif');
            set(self.t2_var_txt4, 'FontSize', 11);
            set(self.t2_var_txt4, 'FontWeight', 'bold');
            set(self.t2_var_txt4, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt4, 'Visible', 'off');  
            
            % iteration edit
            self.t2_var_edt1 = uicontrol('style','edit');
            set(self.t2_var_edt1, 'unit', 'pixels', 'position', [670 505 70 23]);
            set(self.t2_var_edt1, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt1, 'Visible', 'off');
            set(self.t2_var_edt1, 'String', self.user_inputs.tab_2_var.iterations);
            set(self.t2_var_edt1, 'CallBack', @self.cb_t2_var_edt1);
            
            % burn-in label
            self.t2_var_txt5 = uicontrol('style', 'text');
            set(self.t2_var_txt5, 'unit', 'pixels', 'position', [520 475 200 25]);
            set(self.t2_var_txt5, 'String', ' burn-in');
            set(self.t2_var_txt5, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt5, 'FontName', 'Serif');
            set(self.t2_var_txt5, 'FontSize', 11);
            set(self.t2_var_txt5, 'FontWeight', 'bold');
            set(self.t2_var_txt5, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt5, 'Visible', 'off');   
            
            % burn-in edit
            self.t2_var_edt2 = uicontrol('style','edit');
            set(self.t2_var_edt2, 'unit', 'pixels', 'position', [670 480 70 23]);
            set(self.t2_var_edt2, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt2, 'Visible', 'off');
            set(self.t2_var_edt2, 'String', self.user_inputs.tab_2_var.burnin);
            set(self.t2_var_edt2, 'CallBack', @self.cb_t2_var_edt2);
            
            % credibility label
            self.t2_var_txt6 = uicontrol('style', 'text');
            set(self.t2_var_txt6, 'unit', 'pixels', 'position', [520 455 200 20]);
            set(self.t2_var_txt6, 'String', ' credibility level');
            set(self.t2_var_txt6, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt6, 'FontName', 'Serif');
            set(self.t2_var_txt6, 'FontSize', 11);
            set(self.t2_var_txt6, 'FontWeight', 'bold');
            set(self.t2_var_txt6, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt6, 'Visible', 'off');         
            
            % credibility edit
            self.t2_var_edt3 = uicontrol('style','edit');
            set(self.t2_var_edt3, 'unit', 'pixels', 'position', [670 455 70 23]);
            set(self.t2_var_edt3, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt3, 'Visible', 'off');
            set(self.t2_var_edt3, 'String', self.user_inputs.tab_2_var.model_credibility);
            set(self.t2_var_edt3, 'CallBack', @self.cb_t2_var_edt3);

            % Exogenous label
            self.t2_var_txt7 = uicontrol('style', 'text');
            set(self.t2_var_txt7, 'unit', 'pixels', 'position', [770 525 200 30]);
            set(self.t2_var_txt7, 'String', ' Exogenous');
            set(self.t2_var_txt7, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt7, 'FontName', 'Serif');
            set(self.t2_var_txt7, 'FontSize', 14);
            set(self.t2_var_txt7, 'FontAngle', 'italic');
            set(self.t2_var_txt7, 'BackgroundColor', self.background_color);
            set(self.t2_var_txt7, 'Visible', 'off');              

            % constant label
            self.t2_var_txt8 = uicontrol('style', 'text');
            set(self.t2_var_txt8, 'unit', 'pixels', 'position', [770 500 200 25]);
            set(self.t2_var_txt8, 'String', ' constant');
            set(self.t2_var_txt8, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt8, 'FontName', 'Serif');
            set(self.t2_var_txt8, 'FontSize', 11);
            set(self.t2_var_txt8, 'FontWeight', 'bold');
            set(self.t2_var_txt8, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt8, 'Visible', 'off');  
            
            % constant checkbox
            self.t2_var_cbx1 = uicontrol('style', 'checkbox');
            set(self.t2_var_cbx1, 'unit', 'pixels', 'position', [950 505 20 20]);
            set(self.t2_var_cbx1, 'BackgroundColor', self.background_color);
            set(self.t2_var_cbx1, 'Visible', 'off');
            set(self.t2_var_cbx1, 'Value', self.user_inputs.tab_2_var.constant);
            set(self.t2_var_cbx1, 'CallBack', @self.cb_t2_var_cbx1);
            
            % linear trend label
            self.t2_var_txt9 = uicontrol('style', 'text');
            set(self.t2_var_txt9, 'unit', 'pixels', 'position', [770 475 200 25]);
            set(self.t2_var_txt9, 'String', ' linear trend');
            set(self.t2_var_txt9, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt9, 'FontName', 'Serif');
            set(self.t2_var_txt9, 'FontSize', 11);
            set(self.t2_var_txt9, 'FontWeight', 'bold');
            set(self.t2_var_txt9, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt9, 'Visible', 'off');             
            
            % linear trend checkbox
            self.t2_var_cbx2 = uicontrol('style', 'checkbox');
            set(self.t2_var_cbx2, 'unit', 'pixels', 'position', [950 480 20 20]);
            set(self.t2_var_cbx2, 'BackgroundColor', self.background_color);
            set(self.t2_var_cbx2, 'Visible', 'off');
            set(self.t2_var_cbx2, 'Value', self.user_inputs.tab_2_var.trend);
            set(self.t2_var_cbx2, 'CallBack', @self.cb_t2_var_cbx2);
            
            % quadratic trend label
            self.t2_var_txt10 = uicontrol('style', 'text');
            set(self.t2_var_txt10, 'unit', 'pixels', 'position', [770 455 200 20]);
            set(self.t2_var_txt10, 'String', ' quadratic trend');
            set(self.t2_var_txt10, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt10, 'FontName', 'Serif');
            set(self.t2_var_txt10, 'FontSize', 11);
            set(self.t2_var_txt10, 'FontWeight', 'bold');
            set(self.t2_var_txt10, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt10, 'Visible', 'off');             
            
            % quadratic trend checkbox
            self.t2_var_cbx3 = uicontrol('style', 'checkbox');
            set(self.t2_var_cbx3, 'unit', 'pixels', 'position', [950 455 20 20]);
            set(self.t2_var_cbx3, 'BackgroundColor', self.background_color);
            set(self.t2_var_cbx3, 'Visible', 'off');
            set(self.t2_var_cbx3, 'Value', self.user_inputs.tab_2_var.quadratic_trend);
            set(self.t2_var_cbx3, 'CallBack', @self.cb_t2_var_cbx3);
            
            % hyperparameter label
            self.t2_var_txt11 = uicontrol('style', 'text');
            set(self.t2_var_txt11, 'unit', 'pixels', 'position', [30 400 300 30]);
            set(self.t2_var_txt11, 'String', ' Hyperparameters');
            set(self.t2_var_txt11, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt11, 'FontName', 'Serif');
            set(self.t2_var_txt11, 'FontSize', 16);
            set(self.t2_var_txt11, 'FontWeight', 'bold');
            set(self.t2_var_txt11, 'FontAngle', 'italic');
            set(self.t2_var_txt11, 'BackgroundColor', self.background_color);  
            set(self.t2_var_txt11, 'Visible', 'off');
            
            % frame around hyperparameters
            self.t2_var_frm3 = uicontrol('style','frame');
            set(self.t2_var_frm3, 'unit', 'pixels', 'position', [20 20 470 380]);
            set(self.t2_var_frm3, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_var_frm3, 'BackgroundColor', self.background_color);
            set(self.t2_var_frm3, 'Visible', 'off');    
            
            % specification label
            self.t2_var_txt12 = uicontrol('style', 'text');
            set(self.t2_var_txt12, 'unit', 'pixels', 'position', [30 370 300 25]);
            set(self.t2_var_txt12, 'String', ' VAR specification');
            set(self.t2_var_txt12, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt12, 'FontName', 'Serif');
            set(self.t2_var_txt12, 'FontSize', 14);
            set(self.t2_var_txt12, 'FontAngle', 'italic');
            set(self.t2_var_txt12, 'BackgroundColor', self.background_color);
            set(self.t2_var_txt12, 'Visible', 'off'); 
            
            % lag label
            self.t2_var_txt13 = uicontrol('style', 'text');
            set(self.t2_var_txt13, 'unit', 'pixels', 'position', [30 342 300 25]);
            set(self.t2_var_txt13, 'String', ' p:    lags');
            set(self.t2_var_txt13, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt13, 'FontName', 'Serif');
            set(self.t2_var_txt13, 'FontSize', 11);
            set(self.t2_var_txt13, 'FontWeight', 'bold');
            set(self.t2_var_txt13, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt13, 'Visible', 'off');    
            
            % lag edit
            self.t2_var_edt4 = uicontrol('style','edit');
            set(self.t2_var_edt4, 'unit', 'pixels', 'position', [330 347 140 22]);
            set(self.t2_var_edt4, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt4, 'Visible', 'off');
            set(self.t2_var_edt4, 'String', self.user_inputs.tab_2_var.lags);
            set(self.t2_var_edt4, 'CallBack', @self.cb_t2_var_edt4);            
            
            % AR coefficients label
            self.t2_var_txt14 = uicontrol('style', 'text');
            set(self.t2_var_txt14, 'unit', 'pixels', 'position', [30 318 300 25]);
            set(self.t2_var_txt14, 'String', ' δ:    AR coefficients');
            set(self.t2_var_txt14, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt14, 'FontName', 'Serif');
            set(self.t2_var_txt14, 'FontSize', 11);
            set(self.t2_var_txt14, 'FontWeight', 'bold');
            set(self.t2_var_txt14, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt14, 'Visible', 'off'); 
            
            % AR coefficients edit
            self.t2_var_edt5 = uicontrol('style','edit');
            set(self.t2_var_edt5, 'unit', 'pixels', 'position', [330 323 140 22]);
            set(self.t2_var_edt5, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt5, 'Visible', 'off');
            set(self.t2_var_edt5, 'String', self.user_inputs.tab_2_var.ar_coefficients);
            set(self.t2_var_edt5, 'CallBack', @self.cb_t2_var_edt5); 
            
            % pi1 label
            self.t2_var_txt15 = uicontrol('style', 'text');
            set(self.t2_var_txt15, 'unit', 'pixels', 'position', [30 294 300 25]);
            set(self.t2_var_txt15, 'String', ' π₁:  overall tightness');
            set(self.t2_var_txt15, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt15, 'FontName', 'Serif');
            set(self.t2_var_txt15, 'FontSize', 11);
            set(self.t2_var_txt15, 'FontWeight', 'bold');
            set(self.t2_var_txt15, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt15, 'Visible', 'off'); 
            
            % pi1 edit
            self.t2_var_edt6 = uicontrol('style','edit');
            set(self.t2_var_edt6, 'unit', 'pixels', 'position', [330 299 140 22]);
            set(self.t2_var_edt6, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt6, 'Visible', 'off');
            set(self.t2_var_edt6, 'String', self.user_inputs.tab_2_var.pi1);
            set(self.t2_var_edt6, 'CallBack', @self.cb_t2_var_edt6); 
            
            % pi2 label
            self.t2_var_txt16 = uicontrol('style', 'text');
            set(self.t2_var_txt16, 'unit', 'pixels', 'position', [30 270 300 25]);
            set(self.t2_var_txt16, 'String', ' π₂:  cross-variable shrinkage');
            set(self.t2_var_txt16, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt16, 'FontName', 'Serif');
            set(self.t2_var_txt16, 'FontSize', 11);
            set(self.t2_var_txt16, 'FontWeight', 'bold');
            set(self.t2_var_txt16, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt16, 'Visible', 'off'); 
            
            % pi2 edit
            self.t2_var_edt7 = uicontrol('style','edit');
            set(self.t2_var_edt7, 'unit', 'pixels', 'position', [330 275 140 22]);
            set(self.t2_var_edt7, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt7, 'Visible', 'off');
            set(self.t2_var_edt7, 'String', self.user_inputs.tab_2_var.pi2);
            set(self.t2_var_edt7, 'CallBack', @self.cb_t2_var_edt7); 
            
            % pi3 label
            self.t2_var_txt17 = uicontrol('style', 'text');
            set(self.t2_var_txt17, 'unit', 'pixels', 'position', [30 246 300 25]);
            set(self.t2_var_txt17, 'String', ' π₃:  lag decay');
            set(self.t2_var_txt17, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt17, 'FontName', 'Serif');
            set(self.t2_var_txt17, 'FontSize', 11);
            set(self.t2_var_txt17, 'FontWeight', 'bold');
            set(self.t2_var_txt17, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt17, 'Visible', 'off'); 
            
            % pi3 edit
            self.t2_var_edt8 = uicontrol('style','edit');
            set(self.t2_var_edt8, 'unit', 'pixels', 'position', [330 251 140 22]);
            set(self.t2_var_edt8, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt8, 'Visible', 'off');
            set(self.t2_var_edt8, 'String', self.user_inputs.tab_2_var.pi3);
            set(self.t2_var_edt8, 'CallBack', @self.cb_t2_var_edt8); 
            
            % pi4 label
            self.t2_var_txt18 = uicontrol('style', 'text');
            set(self.t2_var_txt18, 'unit', 'pixels', 'position', [30 222 300 25]);
            set(self.t2_var_txt18, 'String', ' π₄:  exogenous slackness');
            set(self.t2_var_txt18, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt18, 'FontName', 'Serif');
            set(self.t2_var_txt18, 'FontSize', 11);
            set(self.t2_var_txt18, 'FontWeight', 'bold');
            set(self.t2_var_txt18, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt18, 'Visible', 'off'); 
            
            % pi4 edit
            self.t2_var_edt9 = uicontrol('style','edit');
            set(self.t2_var_edt9, 'unit', 'pixels', 'position', [330 227 140 22]);
            set(self.t2_var_edt9, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt9, 'Visible', 'off');
            set(self.t2_var_edt9, 'String', self.user_inputs.tab_2_var.pi4);
            set(self.t2_var_edt9, 'CallBack', @self.cb_t2_var_edt9); 
            
            % pi5 label
            self.t2_var_txt19 = uicontrol('style', 'text');
            set(self.t2_var_txt19, 'unit', 'pixels', 'position', [30 198 300 25]);
            set(self.t2_var_txt19, 'String', ' π₅:  sums-of-coefficients tightness');
            set(self.t2_var_txt19, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt19, 'FontName', 'Serif');
            set(self.t2_var_txt19, 'FontSize', 11);
            set(self.t2_var_txt19, 'FontWeight', 'bold');
            set(self.t2_var_txt19, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt19, 'Visible', 'off'); 
            
            % pi5 edit
            self.t2_var_edt10 = uicontrol('style','edit');
            set(self.t2_var_edt10, 'unit', 'pixels', 'position', [330 203 140 22]);
            set(self.t2_var_edt10, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt10, 'Visible', 'off');
            set(self.t2_var_edt10, 'String', self.user_inputs.tab_2_var.pi5);
            set(self.t2_var_edt10, 'CallBack', @self.cb_t2_var_edt10); 
            
            % pi6 label
            self.t2_var_txt20 = uicontrol('style', 'text');
            set(self.t2_var_txt20, 'unit', 'pixels', 'position', [30 174 300 25]);
            set(self.t2_var_txt20, 'String', ' π₆:  initial observation tightness');
            set(self.t2_var_txt20, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt20, 'FontName', 'Serif');
            set(self.t2_var_txt20, 'FontSize', 11);
            set(self.t2_var_txt20, 'FontWeight', 'bold');
            set(self.t2_var_txt20, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt20, 'Visible', 'off'); 
            
            % pi6 edit
            self.t2_var_edt11 = uicontrol('style','edit');
            set(self.t2_var_edt11, 'unit', 'pixels', 'position', [330 179 140 22]);
            set(self.t2_var_edt11, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt11, 'Visible', 'off');
            set(self.t2_var_edt11, 'String', self.user_inputs.tab_2_var.pi6);
            set(self.t2_var_edt11, 'CallBack', @self.cb_t2_var_edt11); 
            
            % pi7 label
            self.t2_var_txt21 = uicontrol('style', 'text');
            set(self.t2_var_txt21, 'unit', 'pixels', 'position', [30 150 300 25]);
            set(self.t2_var_txt21, 'String', ' π₇:  long-run tightness');
            set(self.t2_var_txt21, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt21, 'FontName', 'Serif');
            set(self.t2_var_txt21, 'FontSize', 11);
            set(self.t2_var_txt21, 'FontWeight', 'bold');
            set(self.t2_var_txt21, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt21, 'Visible', 'off'); 
            
            % pi7 edit
            self.t2_var_edt12 = uicontrol('style','edit');
            set(self.t2_var_edt12, 'unit', 'pixels', 'position', [330 155 140 22]);
            set(self.t2_var_edt12, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt12, 'Visible', 'off');
            set(self.t2_var_edt12, 'String', self.user_inputs.tab_2_var.pi7);
            set(self.t2_var_edt12, 'CallBack', @self.cb_t2_var_edt12);
            
            % proxy SVAR label
            self.t2_var_txt22 = uicontrol('style', 'text');
            set(self.t2_var_txt22, 'unit', 'pixels', 'position', [30 120 300 25]);
            set(self.t2_var_txt22, 'String', ' Proxy-SVAR');
            set(self.t2_var_txt22, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt22, 'FontName', 'Serif');
            set(self.t2_var_txt22, 'FontSize', 14);
            set(self.t2_var_txt22, 'FontAngle', 'italic');
            set(self.t2_var_txt22, 'BackgroundColor', self.background_color);
            set(self.t2_var_txt22, 'Visible', 'off'); 
            
            % proxy label
            self.t2_var_txt23 = uicontrol('style', 'text');
            set(self.t2_var_txt23, 'unit', 'pixels', 'position', [30 92 300 25]);
            set(self.t2_var_txt23, 'String', ' proxys');
            set(self.t2_var_txt23, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt23, 'FontName', 'Serif');
            set(self.t2_var_txt23, 'FontSize', 11);
            set(self.t2_var_txt23, 'FontWeight', 'bold');
            set(self.t2_var_txt23, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt23, 'Visible', 'off'); 
            
            % proxy edit
            self.t2_var_edt13 = uicontrol('style','edit');
            set(self.t2_var_edt13, 'unit', 'pixels', 'position', [330 97 140 22]);
            set(self.t2_var_edt13, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt13, 'Visible', 'off');
            set(self.t2_var_edt13, 'String', self.user_inputs.tab_2_var.proxy_variables);
            set(self.t2_var_edt13, 'CallBack', @self.cb_t2_var_edt13);
            
            % relevance label
            self.t2_var_txt24 = uicontrol('style', 'text');
            set(self.t2_var_txt24, 'unit', 'pixels', 'position', [30 68 300 25]);
            set(self.t2_var_txt24, 'String', ' λ:  relevance');
            set(self.t2_var_txt24, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt24, 'FontName', 'Serif');
            set(self.t2_var_txt24, 'FontSize', 11);
            set(self.t2_var_txt24, 'FontWeight', 'bold');
            set(self.t2_var_txt24, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt24, 'Visible', 'off'); 
            
            % relevance edit
            self.t2_var_edt14 = uicontrol('style','edit');
            set(self.t2_var_edt14, 'unit', 'pixels', 'position', [330 73 140 22]);
            set(self.t2_var_edt14, 'HorizontalAlignment', 'center');  
            set(self.t2_var_edt14, 'Visible', 'off');
            set(self.t2_var_edt14, 'String', self.user_inputs.tab_2_var.lamda);
            set(self.t2_var_edt14, 'CallBack', @self.cb_t2_var_edt14);
            
            % prior scheme label
            self.t2_var_txt25 = uicontrol('style', 'text');
            set(self.t2_var_txt25, 'unit', 'pixels', 'position', [30 43 300 25]);
            set(self.t2_var_txt25, 'String', ' prior scheme');
            set(self.t2_var_txt25, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt25, 'FontName', 'Serif');
            set(self.t2_var_txt25, 'FontSize', 11);
            set(self.t2_var_txt25, 'FontWeight', 'bold');
            set(self.t2_var_txt25, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt25, 'Visible', 'off'); 

            % prior radiobuttons
            self.t2_var_bgr2 = uibuttongroup('unit','pixels', 'Position',[325 23 160 50]);
            set(self.t2_var_bgr2, 'BorderType', 'none');
            set(self.t2_var_bgr2, 'BackgroundColor', self.background_color); 
            self.t2_var_rdb8 = uicontrol(self.t2_var_bgr2,'Style','radiobutton');
            set(self.t2_var_rdb8, 'Position',[5 25 220 25]);
            set(self.t2_var_rdb8, 'String',' uninformative');
            set(self.t2_var_rdb8, 'FontName', 'Serif');
            set(self.t2_var_rdb8, 'FontSize', 11);
            set(self.t2_var_rdb8, 'FontWeight', 'bold');
            set(self.t2_var_rdb8, 'BackgroundColor', self.background_color);
            self.t2_var_rdb9 = uicontrol(self.t2_var_bgr2,'Style','radiobutton');
            set(self.t2_var_rdb9, 'Position',[5 0 220 25]);
            set(self.t2_var_rdb9, 'String',' Minnesota');
            set(self.t2_var_rdb9, 'FontName', 'Serif');
            set(self.t2_var_rdb9, 'FontSize', 11);
            set(self.t2_var_rdb9, 'FontWeight', 'bold');
            set(self.t2_var_rdb9, 'BackgroundColor', self.background_color);
            set(self.t2_var_bgr2, 'Visible', 'off');
            if self.user_inputs.tab_2_var.proxy_prior == 1
                set(self.t2_var_bgr2, 'SelectedObject', self.t2_var_rdb8);
            elseif self.user_inputs.tab_2_var.proxy_prior == 2
                set(self.t2_var_bgr2, 'SelectedObject', self.t2_var_rdb9);
            end
            set(self.t2_var_bgr2, 'SelectionChangeFcn', @self.cb_t2_var_bgr2);
            
            % options label
            self.t2_var_txt26 = uicontrol('style', 'text');
            set(self.t2_var_txt26, 'unit', 'pixels', 'position', [520 400 300 30]);
            set(self.t2_var_txt26, 'String', ' Options');
            set(self.t2_var_txt26, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt26, 'FontName', 'Serif');
            set(self.t2_var_txt26, 'FontSize', 16);
            set(self.t2_var_txt26, 'FontWeight', 'bold');
            set(self.t2_var_txt26, 'FontAngle', 'italic');
            set(self.t2_var_txt26, 'BackgroundColor', self.background_color);  
            set(self.t2_var_txt26, 'Visible', 'off');
            
            % frame around options
            self.t2_var_frm4 = uicontrol('style','frame');
            set(self.t2_var_frm4, 'unit', 'pixels', 'position', [510 20 470 380]);
            set(self.t2_var_frm4, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_var_frm4, 'BackgroundColor', self.background_color);
            set(self.t2_var_frm4, 'Visible', 'off');   
            
            % applications label
            self.t2_var_txt27 = uicontrol('style', 'text');
            set(self.t2_var_txt27, 'unit', 'pixels', 'position', [520 370 300 25]);
            set(self.t2_var_txt27, 'String', ' Applications');
            set(self.t2_var_txt27, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt27, 'FontName', 'Serif');
            set(self.t2_var_txt27, 'FontSize', 14);
            set(self.t2_var_txt27, 'FontAngle', 'italic');
            set(self.t2_var_txt27, 'BackgroundColor', self.background_color);
            set(self.t2_var_txt27, 'Visible', 'off'); 

            % in-sample label
            self.t2_var_txt38 = uicontrol('style', 'text');
            set(self.t2_var_txt38, 'unit', 'pixels', 'position', [520 342 300 25]);
            set(self.t2_var_txt38, 'String', ' in-sample fit');
            set(self.t2_var_txt38, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt38, 'FontName', 'Serif');
            set(self.t2_var_txt38, 'FontSize', 11);
            set(self.t2_var_txt38, 'FontWeight', 'bold');
            set(self.t2_var_txt38, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt38, 'Visible', 'off');   
            
            % in-sample radiobuttons
            self.t2_var_bgr10 = uibuttongroup('unit','pixels', 'Position',[825 340 140 30]);
            set(self.t2_var_bgr10, 'BorderType', 'none');
            set(self.t2_var_bgr10, 'BackgroundColor', self.background_color); 
            self.t2_var_rdb24 = uicontrol(self.t2_var_bgr10,'Style','radiobutton');
            set(self.t2_var_rdb24, 'Position',[5 5 60 25]);
            set(self.t2_var_rdb24, 'String',' yes');
            set(self.t2_var_rdb24, 'FontName', 'Serif');
            set(self.t2_var_rdb24, 'FontSize', 12);
            set(self.t2_var_rdb24, 'FontWeight', 'bold');
            set(self.t2_var_rdb24, 'BackgroundColor', self.background_color);
            self.t2_var_rdb25 = uicontrol(self.t2_var_bgr10,'Style','radiobutton');
            set(self.t2_var_rdb25, 'Position',[85 5 60 25]);
            set(self.t2_var_rdb25, 'String',' no');
            set(self.t2_var_rdb25, 'FontName', 'Serif');
            set(self.t2_var_rdb25, 'FontSize', 12);
            set(self.t2_var_rdb25, 'FontWeight', 'bold');
            set(self.t2_var_rdb25, 'BackgroundColor', self.background_color);
            set(self.t2_var_bgr10, 'Visible', 'off');
            if self.user_inputs.tab_2_var.insample_fit
                set(self.t2_var_bgr10, 'SelectedObject', self.t2_var_rdb24);
            else
                set(self.t2_var_bgr10, 'SelectedObject', self.t2_var_rdb25);
            end
            set(self.t2_var_bgr10, 'SelectionChangeFcn', @self.cb_t2_var_bgr10);             

            % constrained coefficients label
            self.t2_var_txt28 = uicontrol('style', 'text');
            set(self.t2_var_txt28, 'unit', 'pixels', 'position', [520 318 300 25]);
            set(self.t2_var_txt28, 'String', ' constrained coefficients');
            set(self.t2_var_txt28, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt28, 'FontName', 'Serif');
            set(self.t2_var_txt28, 'FontSize', 11);
            set(self.t2_var_txt28, 'FontWeight', 'bold');
            set(self.t2_var_txt28, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt28, 'Visible', 'off');   
            
            % constrained coefficients radiobuttons
            self.t2_var_bgr3 = uibuttongroup('unit','pixels', 'Position',[825 316 140 30]);
            set(self.t2_var_bgr3, 'BorderType', 'none');
            set(self.t2_var_bgr3, 'BackgroundColor', self.background_color); 
            self.t2_var_rdb10 = uicontrol(self.t2_var_bgr3,'Style','radiobutton');
            set(self.t2_var_rdb10, 'Position',[5 5 60 25]);
            set(self.t2_var_rdb10, 'String',' yes');
            set(self.t2_var_rdb10, 'FontName', 'Serif');
            set(self.t2_var_rdb10, 'FontSize', 12);
            set(self.t2_var_rdb10, 'FontWeight', 'bold');
            set(self.t2_var_rdb10, 'BackgroundColor', self.background_color);
            self.t2_var_rdb11 = uicontrol(self.t2_var_bgr3,'Style','radiobutton');
            set(self.t2_var_rdb11, 'Position',[85 5 60 25]);
            set(self.t2_var_rdb11, 'String',' no');
            set(self.t2_var_rdb11, 'FontName', 'Serif');
            set(self.t2_var_rdb11, 'FontSize', 12);
            set(self.t2_var_rdb11, 'FontWeight', 'bold');
            set(self.t2_var_rdb11, 'BackgroundColor', self.background_color);
            set(self.t2_var_bgr3, 'Visible', 'off');
            if self.user_inputs.tab_2_var.constrained_coefficients
                set(self.t2_var_bgr3, 'SelectedObject', self.t2_var_rdb10);
            else
                set(self.t2_var_bgr3, 'SelectedObject', self.t2_var_rdb11);
            end
            set(self.t2_var_bgr3, 'SelectionChangeFcn', @self.cb_t2_var_bgr3);         

            % sums-of-coefficients label
            self.t2_var_txt29 = uicontrol('style', 'text');
            set(self.t2_var_txt29, 'unit', 'pixels', 'position', [520 294 300 25]);
            set(self.t2_var_txt29, 'String', ' sums-of-coefficients');
            set(self.t2_var_txt29, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt29, 'FontName', 'Serif');
            set(self.t2_var_txt29, 'FontSize', 11);
            set(self.t2_var_txt29, 'FontWeight', 'bold');
            set(self.t2_var_txt29, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt29, 'Visible', 'off');  
            
            % sums-of-coefficients radiobuttons
            self.t2_var_bgr4 = uibuttongroup('unit','pixels', 'Position',[825 292 140 30]);
            set(self.t2_var_bgr4, 'BorderType', 'none');
            set(self.t2_var_bgr4, 'BackgroundColor', self.background_color); 
            self.t2_var_rdb12 = uicontrol(self.t2_var_bgr4,'Style','radiobutton');
            set(self.t2_var_rdb12, 'Position',[5 5 60 25]);
            set(self.t2_var_rdb12, 'String',' yes');
            set(self.t2_var_rdb12, 'FontName', 'Serif');
            set(self.t2_var_rdb12, 'FontSize', 12);
            set(self.t2_var_rdb12, 'FontWeight', 'bold');
            set(self.t2_var_rdb12, 'BackgroundColor', self.background_color);
            self.t2_var_rdb13 = uicontrol(self.t2_var_bgr4,'Style','radiobutton');
            set(self.t2_var_rdb13, 'Position',[85 5 60 25]);
            set(self.t2_var_rdb13, 'String',' no');
            set(self.t2_var_rdb13, 'FontName', 'Serif');
            set(self.t2_var_rdb13, 'FontSize', 12);
            set(self.t2_var_rdb13, 'FontWeight', 'bold');
            set(self.t2_var_rdb13, 'BackgroundColor', self.background_color);
            set(self.t2_var_bgr4, 'Visible', 'off');
            if self.user_inputs.tab_2_var.sums_of_coefficients
                set(self.t2_var_bgr4, 'SelectedObject', self.t2_var_rdb12);
            else
                set(self.t2_var_bgr4, 'SelectedObject', self.t2_var_rdb13);
            end
            set(self.t2_var_bgr4, 'SelectionChangeFcn', @self.cb_t2_var_bgr4);  
            
            % dummy initial observation label
            self.t2_var_txt30 = uicontrol('style', 'text');
            set(self.t2_var_txt30, 'unit', 'pixels', 'position', [520 270 300 25]);
            set(self.t2_var_txt30, 'String', ' dummy initial observation');
            set(self.t2_var_txt30, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt30, 'FontName', 'Serif');
            set(self.t2_var_txt30, 'FontSize', 11);
            set(self.t2_var_txt30, 'FontWeight', 'bold');
            set(self.t2_var_txt30, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt30, 'Visible', 'off'); 
            
            % dummy initial observation radiobuttons
            self.t2_var_bgr5 = uibuttongroup('unit','pixels', 'Position',[825 268 140 30]);
            set(self.t2_var_bgr5, 'BorderType', 'none');
            set(self.t2_var_bgr5, 'BackgroundColor', self.background_color); 
            self.t2_var_rdb14 = uicontrol(self.t2_var_bgr5,'Style','radiobutton');
            set(self.t2_var_rdb14, 'Position',[5 5 60 25]);
            set(self.t2_var_rdb14, 'String',' yes');
            set(self.t2_var_rdb14, 'FontName', 'Serif');
            set(self.t2_var_rdb14, 'FontSize', 12);
            set(self.t2_var_rdb14, 'FontWeight', 'bold');
            set(self.t2_var_rdb14, 'BackgroundColor', self.background_color);
            self.t2_var_rdb15 = uicontrol(self.t2_var_bgr5,'Style','radiobutton');
            set(self.t2_var_rdb15, 'Position',[85 5 60 25]);
            set(self.t2_var_rdb15, 'String',' no');
            set(self.t2_var_rdb15, 'FontName', 'Serif');
            set(self.t2_var_rdb15, 'FontSize', 12);
            set(self.t2_var_rdb15, 'FontWeight', 'bold');
            set(self.t2_var_rdb15, 'BackgroundColor', self.background_color);
            set(self.t2_var_bgr5, 'Visible', 'off');
            if self.user_inputs.tab_2_var.initial_observation
                set(self.t2_var_bgr5, 'SelectedObject', self.t2_var_rdb14);
            else
                set(self.t2_var_bgr5, 'SelectedObject', self.t2_var_rdb15);
            end
            set(self.t2_var_bgr5, 'SelectionChangeFcn', @self.cb_t2_var_bgr5);             
   
            % long-run label
            self.t2_var_txt31 = uicontrol('style', 'text');
            set(self.t2_var_txt31, 'unit', 'pixels', 'position', [520 246 300 25]);
            set(self.t2_var_txt31, 'String', ' long-run prior');
            set(self.t2_var_txt31, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt31, 'FontName', 'Serif');
            set(self.t2_var_txt31, 'FontSize', 11);
            set(self.t2_var_txt31, 'FontWeight', 'bold');
            set(self.t2_var_txt31, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt31, 'Visible', 'off');  
            
            % long-run radiobuttons
            self.t2_var_bgr6 = uibuttongroup('unit','pixels', 'Position',[825 244 140 30]);
            set(self.t2_var_bgr6, 'BorderType', 'none');
            set(self.t2_var_bgr6, 'BackgroundColor', self.background_color); 
            self.t2_var_rdb16 = uicontrol(self.t2_var_bgr6,'Style','radiobutton');
            set(self.t2_var_rdb16, 'Position',[5 5 60 25]);
            set(self.t2_var_rdb16, 'String',' yes');
            set(self.t2_var_rdb16, 'FontName', 'Serif');
            set(self.t2_var_rdb16, 'FontSize', 12);
            set(self.t2_var_rdb16, 'FontWeight', 'bold');
            set(self.t2_var_rdb16, 'BackgroundColor', self.background_color);
            self.t2_var_rdb17 = uicontrol(self.t2_var_bgr6,'Style','radiobutton');
            set(self.t2_var_rdb17, 'Position',[85 5 60 25]);
            set(self.t2_var_rdb17, 'String',' no');
            set(self.t2_var_rdb17, 'FontName', 'Serif');
            set(self.t2_var_rdb17, 'FontSize', 12);
            set(self.t2_var_rdb17, 'FontWeight', 'bold');
            set(self.t2_var_rdb17, 'BackgroundColor', self.background_color);
            set(self.t2_var_bgr6, 'Visible', 'off');
            if self.user_inputs.tab_2_var.long_run
                set(self.t2_var_bgr6, 'SelectedObject', self.t2_var_rdb16);
            else
                set(self.t2_var_bgr6, 'SelectedObject', self.t2_var_rdb17);
            end
            set(self.t2_var_bgr6, 'SelectionChangeFcn', @self.cb_t2_var_bgr6);             

            % stationary prior label
            self.t2_var_txt32 = uicontrol('style', 'text');
            set(self.t2_var_txt32, 'unit', 'pixels', 'position', [520 222 300 25]);
            set(self.t2_var_txt32, 'String', ' stationary prior');
            set(self.t2_var_txt32, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt32, 'FontName', 'Serif');
            set(self.t2_var_txt32, 'FontSize', 11);
            set(self.t2_var_txt32, 'FontWeight', 'bold');
            set(self.t2_var_txt32, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt32, 'Visible', 'off'); 
            
            % stationary prior radiobuttons
            self.t2_var_bgr7 = uibuttongroup('unit','pixels', 'Position',[825 220 140 30]);
            set(self.t2_var_bgr7, 'BorderType', 'none');
            set(self.t2_var_bgr7, 'BackgroundColor', self.background_color); 
            self.t2_var_rdb18 = uicontrol(self.t2_var_bgr7,'Style','radiobutton');
            set(self.t2_var_rdb18, 'Position',[5 5 60 25]);
            set(self.t2_var_rdb18, 'String',' yes');
            set(self.t2_var_rdb18, 'FontName', 'Serif');
            set(self.t2_var_rdb18, 'FontSize', 12);
            set(self.t2_var_rdb18, 'FontWeight', 'bold');
            set(self.t2_var_rdb18, 'BackgroundColor', self.background_color);
            self.t2_var_rdb19 = uicontrol(self.t2_var_bgr7,'Style','radiobutton');
            set(self.t2_var_rdb19, 'Position',[85 5 60 25]);
            set(self.t2_var_rdb19, 'String',' no');
            set(self.t2_var_rdb19, 'FontName', 'Serif');
            set(self.t2_var_rdb19, 'FontSize', 12);
            set(self.t2_var_rdb19, 'FontWeight', 'bold');
            set(self.t2_var_rdb19, 'BackgroundColor', self.background_color);
            set(self.t2_var_bgr7, 'Visible', 'off');
            if self.user_inputs.tab_2_var.stationary
                set(self.t2_var_bgr7, 'SelectedObject', self.t2_var_rdb18);
            else
                set(self.t2_var_bgr7, 'SelectedObject', self.t2_var_rdb19);
            end
            set(self.t2_var_bgr7, 'SelectionChangeFcn', @self.cb_t2_var_bgr7);             

            % marginal likelihood label
            self.t2_var_txt33 = uicontrol('style', 'text');
            set(self.t2_var_txt33, 'unit', 'pixels', 'position', [520 198 300 25]);
            set(self.t2_var_txt33, 'String', ' marginal likelihood');
            set(self.t2_var_txt33, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt33, 'FontName', 'Serif');
            set(self.t2_var_txt33, 'FontSize', 11);
            set(self.t2_var_txt33, 'FontWeight', 'bold');
            set(self.t2_var_txt33, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt33, 'Visible', 'off'); 
            
            % marginal likelihood radiobuttons
            self.t2_var_bgr8 = uibuttongroup('unit','pixels', 'Position',[825 196 140 30]);
            set(self.t2_var_bgr8, 'BorderType', 'none');
            set(self.t2_var_bgr8, 'BackgroundColor', self.background_color); 
            self.t2_var_rdb20 = uicontrol(self.t2_var_bgr8,'Style','radiobutton');
            set(self.t2_var_rdb20, 'Position',[5 5 60 25]);
            set(self.t2_var_rdb20, 'String',' yes');
            set(self.t2_var_rdb20, 'FontName', 'Serif');
            set(self.t2_var_rdb20, 'FontSize', 12);
            set(self.t2_var_rdb20, 'FontWeight', 'bold');
            set(self.t2_var_rdb20, 'BackgroundColor', self.background_color);
            self.t2_var_rdb21 = uicontrol(self.t2_var_bgr8,'Style','radiobutton');
            set(self.t2_var_rdb21, 'Position',[85 5 60 25]);
            set(self.t2_var_rdb21, 'String',' no');
            set(self.t2_var_rdb21, 'FontName', 'Serif');
            set(self.t2_var_rdb21, 'FontSize', 12);
            set(self.t2_var_rdb21, 'FontWeight', 'bold');
            set(self.t2_var_rdb21, 'BackgroundColor', self.background_color);
            set(self.t2_var_bgr8, 'Visible', 'off');
            if self.user_inputs.tab_2_var.marginal_likelihood
                set(self.t2_var_bgr8, 'SelectedObject', self.t2_var_rdb20);
            else
                set(self.t2_var_bgr8, 'SelectedObject', self.t2_var_rdb21);
            end
            set(self.t2_var_bgr8, 'SelectionChangeFcn', @self.cb_t2_var_bgr8);             

            % optimization label
            self.t2_var_txt34 = uicontrol('style', 'text');
            set(self.t2_var_txt34, 'unit', 'pixels', 'position', [520 174 300 25]);
            set(self.t2_var_txt34, 'String', ' hyperparameter optimization');
            set(self.t2_var_txt34, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt34, 'FontName', 'Serif');
            set(self.t2_var_txt34, 'FontSize', 11);
            set(self.t2_var_txt34, 'FontWeight', 'bold');
            set(self.t2_var_txt34, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt34, 'Visible', 'off'); 
            
            % optimization radiobuttons
            self.t2_var_bgr9 = uibuttongroup('unit','pixels', 'Position',[825 172 140 30]);
            set(self.t2_var_bgr9, 'BorderType', 'none');
            set(self.t2_var_bgr9, 'BackgroundColor', self.background_color); 
            self.t2_var_rdb22 = uicontrol(self.t2_var_bgr9,'Style','radiobutton');
            set(self.t2_var_rdb22, 'Position',[5 5 60 25]);
            set(self.t2_var_rdb22, 'String',' yes');
            set(self.t2_var_rdb22, 'FontName', 'Serif');
            set(self.t2_var_rdb22, 'FontSize', 12);
            set(self.t2_var_rdb22, 'FontWeight', 'bold');
            set(self.t2_var_rdb22, 'BackgroundColor', self.background_color);
            self.t2_var_rdb23 = uicontrol(self.t2_var_bgr9,'Style','radiobutton');
            set(self.t2_var_rdb23, 'Position',[85 5 60 25]);
            set(self.t2_var_rdb23, 'String',' no');
            set(self.t2_var_rdb23, 'FontName', 'Serif');
            set(self.t2_var_rdb23, 'FontSize', 12);
            set(self.t2_var_rdb23, 'FontWeight', 'bold');
            set(self.t2_var_rdb23, 'BackgroundColor', self.background_color);
            set(self.t2_var_bgr9, 'Visible', 'off');
            if self.user_inputs.tab_2_var.hyperparameter_optimization
                set(self.t2_var_bgr9, 'SelectedObject', self.t2_var_rdb22);
            else
                set(self.t2_var_bgr9, 'SelectedObject', self.t2_var_rdb23);
            end
            set(self.t2_var_bgr9, 'SelectionChangeFcn', @self.cb_t2_var_bgr9); 
            
            % files label
            self.t2_var_txt35 = uicontrol('style', 'text');
            set(self.t2_var_txt35, 'unit', 'pixels', 'position', [520 120 300 25]);
            set(self.t2_var_txt35, 'String', ' Files');
            set(self.t2_var_txt35, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt35, 'FontName', 'Serif');
            set(self.t2_var_txt35, 'FontSize', 14);
            set(self.t2_var_txt35, 'FontAngle', 'italic');
            set(self.t2_var_txt35, 'BackgroundColor', self.background_color);
            set(self.t2_var_txt35, 'Visible', 'off'); 
            
            % constrained coefficients label
            self.t2_var_txt36 = uicontrol('style', 'text');
            set(self.t2_var_txt36, 'unit', 'pixels', 'position', [520 92 300 25]);
            set(self.t2_var_txt36, 'String', ' constrained coefficients');
            set(self.t2_var_txt36, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt36, 'FontName', 'Serif');
            set(self.t2_var_txt36, 'FontSize', 11);
            set(self.t2_var_txt36, 'FontWeight', 'bold');
            set(self.t2_var_txt36, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt36, 'Visible', 'off');     
            
            % constrained coefficients edit
            self.t2_var_edt15 = uicontrol('style','edit');
            set(self.t2_var_edt15, 'unit', 'pixels', 'position', [525 73 430 22]);
            set(self.t2_var_edt15, 'HorizontalAlignment', 'left');  
            set(self.t2_var_edt15, 'Visible', 'off');
            set(self.t2_var_edt15, 'String', self.user_inputs.tab_2_var.coefficients_file);
            set(self.t2_var_edt15, 'CallBack', @self.cb_t2_var_edt15);
            
            % long run prior label
            self.t2_var_txt37 = uicontrol('style', 'text');
            set(self.t2_var_txt37, 'unit', 'pixels', 'position', [520 43 300 25]);
            set(self.t2_var_txt37, 'String', ' long-run prior');
            set(self.t2_var_txt37, 'HorizontalAlignment', 'left');
            set(self.t2_var_txt37, 'FontName', 'Serif');
            set(self.t2_var_txt37, 'FontSize', 11);
            set(self.t2_var_txt37, 'FontWeight', 'bold');
            set(self.t2_var_txt37, 'BackgroundColor', self.background_color); 
            set(self.t2_var_txt37, 'Visible', 'off');             
            
            % long run prior edit
            self.t2_var_edt16 = uicontrol('style','edit');
            set(self.t2_var_edt16, 'unit', 'pixels', 'position', [525 24 430 22]);
            set(self.t2_var_edt16, 'HorizontalAlignment', 'left');  
            set(self.t2_var_edt16, 'Visible', 'off');
            set(self.t2_var_edt16, 'String', self.user_inputs.tab_2_var.long_run_file);
            set(self.t2_var_edt16, 'CallBack', @self.cb_t2_var_edt16);            
            
            % indicate that tab 2 for vector autoregression is now created
            self.created_tab_2_var = true;            
        end

        
        function hide_tab_2_var(self)

            % hide all controls            
            set(self.t2_var_txt1, 'Visible', 'off'); 
            set(self.t2_var_txt2, 'Visible', 'off'); 
            set(self.t2_var_txt3, 'Visible', 'off'); 
            set(self.t2_var_txt4, 'Visible', 'off'); 
            set(self.t2_var_txt5, 'Visible', 'off'); 
            set(self.t2_var_txt6, 'Visible', 'off'); 
            set(self.t2_var_txt7, 'Visible', 'off'); 
            set(self.t2_var_txt8, 'Visible', 'off'); 
            set(self.t2_var_txt9, 'Visible', 'off'); 
            set(self.t2_var_txt10, 'Visible', 'off'); 
            set(self.t2_var_txt11, 'Visible', 'off'); 
            set(self.t2_var_txt12, 'Visible', 'off'); 
            set(self.t2_var_txt13, 'Visible', 'off'); 
            set(self.t2_var_txt14, 'Visible', 'off'); 
            set(self.t2_var_txt15, 'Visible', 'off'); 
            set(self.t2_var_txt16, 'Visible', 'off'); 
            set(self.t2_var_txt17, 'Visible', 'off'); 
            set(self.t2_var_txt18, 'Visible', 'off'); 
            set(self.t2_var_txt19, 'Visible', 'off'); 
            set(self.t2_var_txt20, 'Visible', 'off'); 
            set(self.t2_var_txt21, 'Visible', 'off'); 
            set(self.t2_var_txt22, 'Visible', 'off'); 
            set(self.t2_var_txt23, 'Visible', 'off'); 
            set(self.t2_var_txt24, 'Visible', 'off'); 
            set(self.t2_var_txt25, 'Visible', 'off'); 
            set(self.t2_var_txt26, 'Visible', 'off'); 
            set(self.t2_var_txt27, 'Visible', 'off'); 
            set(self.t2_var_txt28, 'Visible', 'off'); 
            set(self.t2_var_txt29, 'Visible', 'off'); 
            set(self.t2_var_txt30, 'Visible', 'off'); 
            set(self.t2_var_txt31, 'Visible', 'off'); 
            set(self.t2_var_txt32, 'Visible', 'off'); 
            set(self.t2_var_txt33, 'Visible', 'off'); 
            set(self.t2_var_txt34, 'Visible', 'off'); 
            set(self.t2_var_txt35, 'Visible', 'off'); 
            set(self.t2_var_txt36, 'Visible', 'off'); 
            set(self.t2_var_txt37, 'Visible', 'off'); 
            set(self.t2_var_txt38, 'Visible', 'off'); 
            set(self.t2_var_frm1, 'Visible', 'off');
            set(self.t2_var_frm2, 'Visible', 'off');
            set(self.t2_var_frm3, 'Visible', 'off');
            set(self.t2_var_frm4, 'Visible', 'off');
            set(self.t2_var_bgr1, 'Visible', 'off');
            set(self.t2_var_bgr2, 'Visible', 'off');            
            set(self.t2_var_bgr3, 'Visible', 'off');            
            set(self.t2_var_bgr4, 'Visible', 'off');            
            set(self.t2_var_bgr5, 'Visible', 'off');
            set(self.t2_var_bgr6, 'Visible', 'off');            
            set(self.t2_var_bgr7, 'Visible', 'off');            
            set(self.t2_var_bgr8, 'Visible', 'off');
            set(self.t2_var_bgr9, 'Visible', 'off');    
            set(self.t2_var_bgr10, 'Visible', 'off');    
            set(self.t2_var_rdb1, 'Visible', 'off');
            set(self.t2_var_rdb2, 'Visible', 'off');
            set(self.t2_var_rdb3, 'Visible', 'off');
            set(self.t2_var_rdb4, 'Visible', 'off');
            set(self.t2_var_rdb5, 'Visible', 'off');
            set(self.t2_var_rdb6, 'Visible', 'off');
            set(self.t2_var_rdb7, 'Visible', 'off');
            set(self.t2_var_rdb8, 'Visible', 'off');
            set(self.t2_var_rdb9, 'Visible', 'off');
            set(self.t2_var_rdb8, 'Visible', 'off');
            set(self.t2_var_rdb9, 'Visible', 'off');
            set(self.t2_var_rdb10, 'Visible', 'off');
            set(self.t2_var_rdb11, 'Visible', 'off');
            set(self.t2_var_rdb12, 'Visible', 'off');
            set(self.t2_var_rdb13, 'Visible', 'off');
            set(self.t2_var_rdb14, 'Visible', 'off');
            set(self.t2_var_rdb15, 'Visible', 'off');
            set(self.t2_var_rdb16, 'Visible', 'off');
            set(self.t2_var_rdb17, 'Visible', 'off');
            set(self.t2_var_rdb18, 'Visible', 'off');
            set(self.t2_var_rdb19, 'Visible', 'off');
            set(self.t2_var_rdb20, 'Visible', 'off');
            set(self.t2_var_rdb21, 'Visible', 'off');
            set(self.t2_var_rdb22, 'Visible', 'off');
            set(self.t2_var_rdb23, 'Visible', 'off');
            set(self.t2_var_rdb24, 'Visible', 'off');       
            set(self.t2_var_rdb25, 'Visible', 'off');       
            set(self.t2_var_edt1, 'Visible', 'off'); 
            set(self.t2_var_edt2, 'Visible', 'off'); 
            set(self.t2_var_edt3, 'Visible', 'off'); 
            set(self.t2_var_edt4, 'Visible', 'off'); 
            set(self.t2_var_edt5, 'Visible', 'off'); 
            set(self.t2_var_edt6, 'Visible', 'off'); 
            set(self.t2_var_edt7, 'Visible', 'off'); 
            set(self.t2_var_edt8, 'Visible', 'off'); 
            set(self.t2_var_edt9, 'Visible', 'off'); 
            set(self.t2_var_edt10, 'Visible', 'off'); 
            set(self.t2_var_edt11, 'Visible', 'off'); 
            set(self.t2_var_edt12, 'Visible', 'off'); 
            set(self.t2_var_edt13, 'Visible', 'off'); 
            set(self.t2_var_edt14, 'Visible', 'off');
            set(self.t2_var_edt15, 'Visible', 'off'); 
            set(self.t2_var_edt16, 'Visible', 'off'); 
            set(self.t2_var_cbx1, 'Visible', 'off');
            set(self.t2_var_cbx2, 'Visible', 'off');
            set(self.t2_var_cbx3, 'Visible', 'off'); 
            
            % update tab color
            set(self.tab_pbt2, 'BackgroundColor', self.backtabs_color);
        end
        
        
        function show_tab_2_var(self)

            % show all controls
            set(self.t2_var_txt1, 'Visible', 'on');     
            set(self.t2_var_txt2, 'Visible', 'on');   
            set(self.t2_var_txt3, 'Visible', 'on');    
            set(self.t2_var_txt4, 'Visible', 'on');    
            set(self.t2_var_txt5, 'Visible', 'on');    
            set(self.t2_var_txt6, 'Visible', 'on');  
            set(self.t2_var_txt7, 'Visible', 'on');  
            set(self.t2_var_txt8, 'Visible', 'on');  
            set(self.t2_var_txt9, 'Visible', 'on');  
            set(self.t2_var_txt10, 'Visible', 'on');  
            set(self.t2_var_txt11, 'Visible', 'on'); 
            set(self.t2_var_txt12, 'Visible', 'on'); 
            set(self.t2_var_txt13, 'Visible', 'on'); 
            set(self.t2_var_txt14, 'Visible', 'on'); 
            set(self.t2_var_txt15, 'Visible', 'on'); 
            set(self.t2_var_txt16, 'Visible', 'on'); 
            set(self.t2_var_txt17, 'Visible', 'on'); 
            set(self.t2_var_txt18, 'Visible', 'on'); 
            set(self.t2_var_txt19, 'Visible', 'on'); 
            set(self.t2_var_txt20, 'Visible', 'on'); 
            set(self.t2_var_txt21, 'Visible', 'on'); 
            set(self.t2_var_txt22, 'Visible', 'on'); 
            set(self.t2_var_txt23, 'Visible', 'on'); 
            set(self.t2_var_txt24, 'Visible', 'on'); 
            set(self.t2_var_txt25, 'Visible', 'on'); 
            set(self.t2_var_txt26, 'Visible', 'on'); 
            set(self.t2_var_txt27, 'Visible', 'on'); 
            set(self.t2_var_txt28, 'Visible', 'on'); 
            set(self.t2_var_txt29, 'Visible', 'on'); 
            set(self.t2_var_txt30, 'Visible', 'on'); 
            set(self.t2_var_txt31, 'Visible', 'on'); 
            set(self.t2_var_txt32, 'Visible', 'on'); 
            set(self.t2_var_txt33, 'Visible', 'on'); 
            set(self.t2_var_txt34, 'Visible', 'on'); 
            set(self.t2_var_txt35, 'Visible', 'on'); 
            set(self.t2_var_txt36, 'Visible', 'on'); 
            set(self.t2_var_txt37, 'Visible', 'on'); 
            set(self.t2_var_txt38, 'Visible', 'on'); 
            set(self.t2_var_frm1, 'Visible', 'on');
            set(self.t2_var_frm2, 'Visible', 'on');
            set(self.t2_var_frm3, 'Visible', 'on');
            set(self.t2_var_frm4, 'Visible', 'on');
            set(self.t2_var_bgr1, 'Visible', 'on');
            set(self.t2_var_bgr2, 'Visible', 'on');            
            set(self.t2_var_bgr3, 'Visible', 'on');            
            set(self.t2_var_bgr4, 'Visible', 'on');            
            set(self.t2_var_bgr5, 'Visible', 'on');
            set(self.t2_var_bgr6, 'Visible', 'on');            
            set(self.t2_var_bgr7, 'Visible', 'on');            
            set(self.t2_var_bgr8, 'Visible', 'on');
            set(self.t2_var_bgr9, 'Visible', 'on');  
            set(self.t2_var_bgr10, 'Visible', 'on');   
            set(self.t2_var_rdb1, 'Visible', 'on');
            set(self.t2_var_rdb2, 'Visible', 'on');
            set(self.t2_var_rdb3, 'Visible', 'on');
            set(self.t2_var_rdb4, 'Visible', 'on');
            set(self.t2_var_rdb5, 'Visible', 'on');
            set(self.t2_var_rdb6, 'Visible', 'on');
            set(self.t2_var_rdb7, 'Visible', 'on');
            set(self.t2_var_rdb8, 'Visible', 'on');
            set(self.t2_var_rdb9, 'Visible', 'on');
            set(self.t2_var_rdb8, 'Visible', 'on');
            set(self.t2_var_rdb9, 'Visible', 'on');
            set(self.t2_var_rdb10, 'Visible', 'on');
            set(self.t2_var_rdb11, 'Visible', 'on');
            set(self.t2_var_rdb12, 'Visible', 'on');
            set(self.t2_var_rdb13, 'Visible', 'on');
            set(self.t2_var_rdb14, 'Visible', 'on');
            set(self.t2_var_rdb15, 'Visible', 'on');
            set(self.t2_var_rdb16, 'Visible', 'on');
            set(self.t2_var_rdb17, 'Visible', 'on');
            set(self.t2_var_rdb18, 'Visible', 'on');
            set(self.t2_var_rdb19, 'Visible', 'on');
            set(self.t2_var_rdb20, 'Visible', 'on');
            set(self.t2_var_rdb21, 'Visible', 'on');
            set(self.t2_var_rdb22, 'Visible', 'on');
            set(self.t2_var_rdb23, 'Visible', 'on'); 
            set(self.t2_var_rdb24, 'Visible', 'on'); 
            set(self.t2_var_rdb25, 'Visible', 'on'); 
            set(self.t2_var_edt1, 'Visible', 'on'); 
            set(self.t2_var_edt2, 'Visible', 'on'); 
            set(self.t2_var_edt3, 'Visible', 'on');
            set(self.t2_var_edt4, 'Visible', 'on');
            set(self.t2_var_edt5, 'Visible', 'on');
            set(self.t2_var_edt6, 'Visible', 'on');
            set(self.t2_var_edt7, 'Visible', 'on');
            set(self.t2_var_edt8, 'Visible', 'on');
            set(self.t2_var_edt9, 'Visible', 'on');
            set(self.t2_var_edt10, 'Visible', 'on');
            set(self.t2_var_edt11, 'Visible', 'on');
            set(self.t2_var_edt12, 'Visible', 'on');
            set(self.t2_var_edt13, 'Visible', 'on');
            set(self.t2_var_edt14, 'Visible', 'on');
            set(self.t2_var_edt15, 'Visible', 'on');
            set(self.t2_var_edt16, 'Visible', 'on');
            set(self.t2_var_cbx1, 'Visible', 'on');
            set(self.t2_var_cbx2, 'Visible', 'on');
            set(self.t2_var_cbx3, 'Visible', 'on');
        end
        
        
        function cb_t2_var_bgr1(self, hObject, callbackdata)
           if get(self.t2_var_bgr1, 'SelectedObject') == self.t2_var_rdb1
               self.user_inputs.tab_2_var.var_type = 1;
           elseif get(self.t2_var_bgr1, 'SelectedObject') == self.t2_var_rdb2
               self.user_inputs.tab_2_var.var_type = 2;
           elseif get(self.t2_var_bgr1, 'SelectedObject') == self.t2_var_rdb3
               self.user_inputs.tab_2_var.var_type = 3;               
           elseif get(self.t2_var_bgr1, 'SelectedObject') == self.t2_var_rdb4
               self.user_inputs.tab_2_var.var_type = 4;               
           elseif get(self.t2_var_bgr1, 'SelectedObject') == self.t2_var_rdb5
               self.user_inputs.tab_2_var.var_type = 5;               
           elseif get(self.t2_var_bgr1, 'SelectedObject') == self.t2_var_rdb6
               self.user_inputs.tab_2_var.var_type = 6;   
           elseif get(self.t2_var_bgr1, 'SelectedObject') == self.t2_var_rdb7
               self.user_inputs.tab_2_var.var_type = 7; 
           end
        end
        
        function cb_t2_var_edt1(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.iterations = get(self.t2_var_edt1, 'String');
        end 
        
        function cb_t2_var_edt2(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.burnin = get(self.t2_var_edt2, 'String');
        end 
        
        function cb_t2_var_edt3(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.model_credibility = get(self.t2_var_edt3, 'String');
        end 
        
        function cb_t2_var_cbx1(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.constant = logical(get(self.t2_var_cbx1, 'Value'));
        end 
        
        function cb_t2_var_cbx2(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.trend = logical(get(self.t2_var_cbx2, 'Value'));
        end 
        
        function cb_t2_var_cbx3(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.quadratic_trend = logical(get(self.t2_var_cbx3, 'Value'));
        end 
        
        function cb_t2_var_edt4(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.lags = get(self.t2_var_edt4, 'String');
        end 
        
        function cb_t2_var_edt5(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.ar_coefficients = get(self.t2_var_edt5, 'String');
        end
        
        function cb_t2_var_edt6(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.pi1 = get(self.t2_var_edt6, 'String');
        end 
        
        function cb_t2_var_edt7(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.pi2 = get(self.t2_var_edt7, 'String');
        end 
        
        function cb_t2_var_edt8(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.pi3 = get(self.t2_var_edt8, 'String');
        end 
        
        function cb_t2_var_edt9(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.pi4 = get(self.t2_var_edt9, 'String');
        end 
        
        function cb_t2_var_edt10(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.pi5 = get(self.t2_var_edt10, 'String');
        end 
        
        function cb_t2_var_edt11(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.pi6 = get(self.t2_var_edt11, 'String');
        end 
        
        function cb_t2_var_edt12(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.pi7 = get(self.t2_var_edt12, 'String');
        end 
        
        function cb_t2_var_edt13(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.proxy_variables = get(self.t2_var_edt13, 'String');
        end 
        
        function cb_t2_var_edt14(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.lamda = get(self.t2_var_edt14, 'String');
        end 
        
        function cb_t2_var_bgr2(self, hObject, callbackdata)
           if get(self.t2_var_bgr2, 'SelectedObject') == self.t2_var_rdb8
               self.user_inputs.tab_2_var.proxy_prior = 1;
           elseif get(self.t2_var_bgr2, 'SelectedObject') == self.t2_var_rdb9
               self.user_inputs.tab_2_var.proxy_prior = 2;
           end
        end
        
        function cb_t2_var_bgr3(self, hObject, callbackdata)
           if get(self.t2_var_bgr3, 'SelectedObject') == self.t2_var_rdb10
               self.user_inputs.tab_2_var.constrained_coefficients = true;
           elseif get(self.t2_var_bgr3, 'SelectedObject') == self.t2_var_rdb11
               self.user_inputs.tab_2_var.constrained_coefficients = false;
           end
        end
        
        function cb_t2_var_bgr4(self, hObject, callbackdata)
           if get(self.t2_var_bgr4, 'SelectedObject') == self.t2_var_rdb12
               self.user_inputs.tab_2_var.sums_of_coefficients = true;
           elseif get(self.t2_var_bgr4, 'SelectedObject') == self.t2_var_rdb13
               self.user_inputs.tab_2_var.sums_of_coefficients = false;
           end
        end
        
        function cb_t2_var_bgr5(self, hObject, callbackdata)
           if get(self.t2_var_bgr5, 'SelectedObject') == self.t2_var_rdb14
               self.user_inputs.tab_2_var.initial_observation = true;
           elseif get(self.t2_var_bgr5, 'SelectedObject') == self.t2_var_rdb15
               self.user_inputs.tab_2_var.initial_observation = false;
           end
        end
        
        function cb_t2_var_bgr6(self, hObject, callbackdata)
           if get(self.t2_var_bgr6, 'SelectedObject') == self.t2_var_rdb16
               self.user_inputs.tab_2_var.long_run = true;
           elseif get(self.t2_var_bgr6, 'SelectedObject') == self.t2_var_rdb17
               self.user_inputs.tab_2_var.long_run = false;
           end
        end
        
        function cb_t2_var_bgr7(self, hObject, callbackdata)
           if get(self.t2_var_bgr7, 'SelectedObject') == self.t2_var_rdb18
               self.user_inputs.tab_2_var.stationary = true;
           elseif get(self.t2_var_bgr7, 'SelectedObject') == self.t2_var_rdb19
               self.user_inputs.tab_2_var.stationary = false;
           end
        end
        
        function cb_t2_var_bgr8(self, hObject, callbackdata)
           if get(self.t2_var_bgr8, 'SelectedObject') == self.t2_var_rdb20
               self.user_inputs.tab_2_var.marginal_likelihood = true;
           elseif get(self.t2_var_bgr8, 'SelectedObject') == self.t2_var_rdb21
               self.user_inputs.tab_2_var.marginal_likelihood = false;
           end
        end
        
        function cb_t2_var_bgr9(self, hObject, callbackdata)
           if get(self.t2_var_bgr9, 'SelectedObject') == self.t2_var_rdb22
               self.user_inputs.tab_2_var.hyperparameter_optimization = true;
           elseif get(self.t2_var_bgr9, 'SelectedObject') == self.t2_var_rdb23
               self.user_inputs.tab_2_var.hyperparameter_optimization = false;
           end
        end
        
        function cb_t2_var_bgr10(self, hObject, callbackdata)
           if get(self.t2_var_bgr10, 'SelectedObject') == self.t2_var_rdb24
               self.user_inputs.tab_2_var.insample_fit = true;
           elseif get(self.t2_var_bgr10, 'SelectedObject') == self.t2_var_rdb25
               self.user_inputs.tab_2_var.insample_fit = false;
           end
        end

        function cb_t2_var_edt15(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.coefficients_file = get(self.t2_var_edt15, 'String');
        end
        
        function cb_t2_var_edt16(self, hObject, callbackdata)
            self.user_inputs.tab_2_var.long_run_file = get(self.t2_var_edt16, 'String');
        end 
        
    end  
     
end