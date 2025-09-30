classdef Tab2VecVarmaInterface < handle

    
    
    %---------------------------------------------------
    % Properties
    %---------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess = public)
        % tab 2 properties (var extension)
        t2_ext_txt1
        t2_ext_txt2        
        t2_ext_txt3
        t2_ext_txt4        
        t2_ext_txt5        
        t2_ext_txt6
        t2_ext_txt7
        t2_ext_txt8        
        t2_ext_txt9   
        t2_ext_txt10
        t2_ext_txt11
        t2_ext_txt12
        t2_ext_txt13
        t2_ext_txt14
        t2_ext_txt15
        t2_ext_txt16
        t2_ext_txt17
        t2_ext_txt18
        t2_ext_txt19
        t2_ext_txt20
        t2_ext_txt21
        t2_ext_txt22
        t2_ext_txt23        
        t2_ext_txt24
        t2_ext_txt25        
        t2_ext_txt26        
        t2_ext_txt27
        t2_ext_txt28
        t2_ext_txt29
        t2_ext_txt30
        t2_ext_txt31
        t2_ext_txt32
        t2_ext_txt33
        t2_ext_txt34        
        t2_ext_txt35
        t2_ext_txt36
        t2_ext_frm1
        t2_ext_frm2        
        t2_ext_frm3
        t2_ext_frm4        
        t2_ext_bgr1
        t2_ext_bgr2        
        t2_ext_bgr3        
        t2_ext_rdb1
        t2_ext_rdb2
        t2_ext_rdb3
        t2_ext_rdb4
        t2_ext_rdb5
        t2_ext_rdb6
        t2_ext_rdb7
        t2_ext_edt1
        t2_ext_edt2
        t2_ext_edt3
        t2_ext_edt4
        t2_ext_edt5
        t2_ext_edt6
        t2_ext_edt7
        t2_ext_edt8
        t2_ext_edt9
        t2_ext_edt10
        t2_ext_edt11
        t2_ext_edt12
        t2_ext_edt13
        t2_ext_edt14
        t2_ext_edt15
        t2_ext_edt16   
        t2_ext_edt17   
        t2_ext_edt18   
        t2_ext_edt19   
        t2_ext_cbx1
        t2_ext_cbx2
        t2_ext_cbx3
    end
    

    %---------------------------------------------------
    % Methods
    %---------------------------------------------------


    methods (Access = public)    
    
    
        function self = Tab2VecVarmaInterface()
        end

        
        function create_tab_2_ext(self)
    
            % vector autoregression label
            self.t2_ext_txt1 = uicontrol('style', 'text');
            set(self.t2_ext_txt1, 'unit', 'pixels', 'position', [30 560 300 30]);
            set(self.t2_ext_txt1, 'String', ' Model');
            set(self.t2_ext_txt1, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt1, 'FontName', 'Serif');
            set(self.t2_ext_txt1, 'FontSize', 16);
            set(self.t2_ext_txt1, 'FontWeight', 'bold');
            set(self.t2_ext_txt1, 'FontAngle', 'italic');
            set(self.t2_ext_txt1, 'BackgroundColor', self.background_color);  
            set(self.t2_ext_txt1, 'Visible', 'off');

            % frame around VAR type
            self.t2_ext_frm1 = uicontrol('style','frame');
            set(self.t2_ext_frm1, 'unit', 'pixels', 'position', [20 450 470 110]);
            set(self.t2_ext_frm1, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_ext_frm1, 'BackgroundColor', self.background_color);
            set(self.t2_ext_frm1, 'Visible', 'off');
            
            % model radiobuttons
            self.t2_ext_bgr1 = uibuttongroup('unit','pixels', 'Position',[25 455 440 100]);
            set(self.t2_ext_bgr1, 'BorderType', 'none');
            set(self.t2_ext_bgr1, 'BackgroundColor', self.background_color); 
            self.t2_ext_rdb1 = uicontrol(self.t2_ext_bgr1,'Style','radiobutton');
            set(self.t2_ext_rdb1, 'Position',[5 75 220 25]);
            set(self.t2_ext_rdb1, 'String',' VEC');
            set(self.t2_ext_rdb1, 'FontName', 'Serif');
            set(self.t2_ext_rdb1, 'FontSize', 12);
            set(self.t2_ext_rdb1, 'FontWeight', 'bold');
            set(self.t2_ext_rdb1, 'BackgroundColor', self.background_color);
            self.t2_ext_rdb2 = uicontrol(self.t2_ext_bgr1,'Style','radiobutton');
            set(self.t2_ext_rdb2, 'Position',[5 25 220 25]);
            set(self.t2_ext_rdb2, 'String',' VARMA');
            set(self.t2_ext_rdb2, 'FontName', 'Serif');
            set(self.t2_ext_rdb2, 'FontSize', 12);
            set(self.t2_ext_rdb2, 'FontWeight', 'bold');
            set(self.t2_ext_rdb2, 'BackgroundColor', self.background_color);
            set(self.t2_ext_bgr1, 'Visible', 'off');
            if self.user_inputs.tab_2_ext.model == 1
                set(self.t2_ext_bgr1, 'SelectedObject', self.t2_ext_rdb1);
            elseif self.user_inputs.tab_2_ext.model == 2
                set(self.t2_ext_bgr1, 'SelectedObject', self.t2_ext_rdb2);
            end
            set(self.t2_ext_bgr1, 'SelectionChangeFcn', @self.cb_t2_ext_bgr1);

            % VEC label
            self.t2_ext_txt2 = uicontrol('style', 'text');
            set(self.t2_ext_txt2, 'unit', 'pixels', 'position', [50 501 400 25]);
            set(self.t2_ext_txt2, 'String', '(Bayesian Vector Error Correction)');
            set(self.t2_ext_txt2, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt2, 'FontName', 'Serif');
            set(self.t2_ext_txt2, 'FontSize', 11);
            set(self.t2_ext_txt2, 'FontWeight', 'bold');
            set(self.t2_ext_txt2, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt2, 'Visible', 'off'); 

            % VARMA label
            self.t2_ext_txt3 = uicontrol('style', 'text');
            set(self.t2_ext_txt3, 'unit', 'pixels', 'position', [50 451 430 25]);
            set(self.t2_ext_txt3, 'String', '(Bayesian Vector Autoregressive Moving Average)');
            set(self.t2_ext_txt3, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt3, 'FontName', 'Serif');
            set(self.t2_ext_txt3, 'FontSize', 11);
            set(self.t2_ext_txt3, 'FontWeight', 'bold');
            set(self.t2_ext_txt3, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt3, 'Visible', 'off'); 

            % estimation label
            self.t2_ext_txt4 = uicontrol('style', 'text');
            set(self.t2_ext_txt4, 'unit', 'pixels', 'position', [520 560 300 30]);
            set(self.t2_ext_txt4, 'String', ' Estimation');
            set(self.t2_ext_txt4, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt4, 'FontName', 'Serif');
            set(self.t2_ext_txt4, 'FontSize', 16);
            set(self.t2_ext_txt4, 'FontWeight', 'bold');
            set(self.t2_ext_txt4, 'FontAngle', 'italic');
            set(self.t2_ext_txt4, 'BackgroundColor', self.background_color);  
            set(self.t2_ext_txt4, 'Visible', 'off');

            % frame around estimation
            self.t2_ext_frm2 = uicontrol('style','frame');
            set(self.t2_ext_frm2, 'unit', 'pixels', 'position', [510 450 470 110]);
            set(self.t2_ext_frm2, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_ext_frm2, 'BackgroundColor', self.background_color);
            set(self.t2_ext_frm2, 'Visible', 'off');

            % Gibbs sampling label
            self.t2_ext_txt5 = uicontrol('style', 'text');
            set(self.t2_ext_txt5, 'unit', 'pixels', 'position', [520 525 200 30]);
            set(self.t2_ext_txt5, 'String', ' Gibbs sampling');
            set(self.t2_ext_txt5, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt5, 'FontName', 'Serif');
            set(self.t2_ext_txt5, 'FontSize', 14);
            set(self.t2_ext_txt5, 'FontAngle', 'italic');
            set(self.t2_ext_txt5, 'BackgroundColor', self.background_color);
            set(self.t2_ext_txt5, 'Visible', 'off');            

            % iteration label
            self.t2_ext_txt6 = uicontrol('style', 'text');
            set(self.t2_ext_txt6, 'unit', 'pixels', 'position', [520 500 200 25]);
            set(self.t2_ext_txt6, 'String', ' iterations');
            set(self.t2_ext_txt6, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt6, 'FontName', 'Serif');
            set(self.t2_ext_txt6, 'FontSize', 11);
            set(self.t2_ext_txt6, 'FontWeight', 'bold');
            set(self.t2_ext_txt6, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt6, 'Visible', 'off');  

            % iteration edit
            self.t2_ext_edt1 = uicontrol('style','edit');
            set(self.t2_ext_edt1, 'unit', 'pixels', 'position', [670 505 70 23]);
            set(self.t2_ext_edt1, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt1, 'Visible', 'off');
            set(self.t2_ext_edt1, 'String', self.user_inputs.tab_2_ext.iterations);
            set(self.t2_ext_edt1, 'CallBack', @self.cb_t2_ext_edt1);

            % burn-in label
            self.t2_ext_txt7 = uicontrol('style', 'text');
            set(self.t2_ext_txt7, 'unit', 'pixels', 'position', [520 475 200 25]);
            set(self.t2_ext_txt7, 'String', ' burn-in');
            set(self.t2_ext_txt7, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt7, 'FontName', 'Serif');
            set(self.t2_ext_txt7, 'FontSize', 11);
            set(self.t2_ext_txt7, 'FontWeight', 'bold');
            set(self.t2_ext_txt7, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt7, 'Visible', 'off');   

            % burn-in edit
            self.t2_ext_edt2 = uicontrol('style','edit');
            set(self.t2_ext_edt2, 'unit', 'pixels', 'position', [670 480 70 23]);
            set(self.t2_ext_edt2, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt2, 'Visible', 'off');
            set(self.t2_ext_edt2, 'String', self.user_inputs.tab_2_ext.burnin);
            set(self.t2_ext_edt2, 'CallBack', @self.cb_t2_ext_edt2);

            % credibility label
            self.t2_ext_txt8 = uicontrol('style', 'text');
            set(self.t2_ext_txt8, 'unit', 'pixels', 'position', [520 455 200 20]);
            set(self.t2_ext_txt8, 'String', ' credibility level');
            set(self.t2_ext_txt8, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt8, 'FontName', 'Serif');
            set(self.t2_ext_txt8, 'FontSize', 11);
            set(self.t2_ext_txt8, 'FontWeight', 'bold');
            set(self.t2_ext_txt8, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt8, 'Visible', 'off');         

            % credibility edit
            self.t2_ext_edt3 = uicontrol('style','edit');
            set(self.t2_ext_edt3, 'unit', 'pixels', 'position', [670 455 70 23]);
            set(self.t2_ext_edt3, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt3, 'Visible', 'off');
            set(self.t2_ext_edt3, 'String', self.user_inputs.tab_2_ext.model_credibility);
            set(self.t2_ext_edt3, 'CallBack', @self.cb_t2_ext_edt3);

            % Exogenous label
            self.t2_ext_txt9 = uicontrol('style', 'text');
            set(self.t2_ext_txt9, 'unit', 'pixels', 'position', [770 525 200 30]);
            set(self.t2_ext_txt9, 'String', ' Exogenous');
            set(self.t2_ext_txt9, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt9, 'FontName', 'Serif');
            set(self.t2_ext_txt9, 'FontSize', 14);
            set(self.t2_ext_txt9, 'FontAngle', 'italic');
            set(self.t2_ext_txt9, 'BackgroundColor', self.background_color);
            set(self.t2_ext_txt9, 'Visible', 'off');              

            % constant label
            self.t2_ext_txt10 = uicontrol('style', 'text');
            set(self.t2_ext_txt10, 'unit', 'pixels', 'position', [770 500 200 25]);
            set(self.t2_ext_txt10, 'String', ' constant');
            set(self.t2_ext_txt10, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt10, 'FontName', 'Serif');
            set(self.t2_ext_txt10, 'FontSize', 11);
            set(self.t2_ext_txt10, 'FontWeight', 'bold');
            set(self.t2_ext_txt10, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt10, 'Visible', 'off');  

            % constant checkbox
            self.t2_ext_cbx1 = uicontrol('style', 'checkbox');
            set(self.t2_ext_cbx1, 'unit', 'pixels', 'position', [950 505 20 20]);
            set(self.t2_ext_cbx1, 'BackgroundColor', self.background_color);
            set(self.t2_ext_cbx1, 'Visible', 'off');
            set(self.t2_ext_cbx1, 'Value', self.user_inputs.tab_2_ext.constant);
            set(self.t2_ext_cbx1, 'CallBack', @self.cb_t2_ext_cbx1);

            % linear trend label
            self.t2_ext_txt11 = uicontrol('style', 'text');
            set(self.t2_ext_txt11, 'unit', 'pixels', 'position', [770 475 200 25]);
            set(self.t2_ext_txt11, 'String', ' linear trend');
            set(self.t2_ext_txt11, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt11, 'FontName', 'Serif');
            set(self.t2_ext_txt11, 'FontSize', 11);
            set(self.t2_ext_txt11, 'FontWeight', 'bold');
            set(self.t2_ext_txt11, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt11, 'Visible', 'off');             

            % linear trend checkbox
            self.t2_ext_cbx2 = uicontrol('style', 'checkbox');
            set(self.t2_ext_cbx2, 'unit', 'pixels', 'position', [950 480 20 20]);
            set(self.t2_ext_cbx2, 'BackgroundColor', self.background_color);
            set(self.t2_ext_cbx2, 'Visible', 'off');
            set(self.t2_ext_cbx2, 'Value', self.user_inputs.tab_2_ext.trend);
            set(self.t2_ext_cbx2, 'CallBack', @self.cb_t2_ext_cbx2);

            % quadratic trend label
            self.t2_ext_txt12 = uicontrol('style', 'text');
            set(self.t2_ext_txt12, 'unit', 'pixels', 'position', [770 455 200 20]);
            set(self.t2_ext_txt12, 'String', ' quadratic trend');
            set(self.t2_ext_txt12, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt12, 'FontName', 'Serif');
            set(self.t2_ext_txt12, 'FontSize', 11);
            set(self.t2_ext_txt12, 'FontWeight', 'bold');
            set(self.t2_ext_txt12, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt12, 'Visible', 'off');             

            % quadratic trend checkbox
            self.t2_ext_cbx3 = uicontrol('style', 'checkbox');
            set(self.t2_ext_cbx3, 'unit', 'pixels', 'position', [950 455 20 20]);
            set(self.t2_ext_cbx3, 'BackgroundColor', self.background_color);
            set(self.t2_ext_cbx3, 'Visible', 'off');
            set(self.t2_ext_cbx3, 'Value', self.user_inputs.tab_2_ext.quadratic_trend);
            set(self.t2_ext_cbx3, 'CallBack', @self.cb_t2_ext_cbx3);

            % vec label
            self.t2_ext_txt13 = uicontrol('style', 'text');
            set(self.t2_ext_txt13, 'unit', 'pixels', 'position', [30 400 300 30]);
            set(self.t2_ext_txt13, 'String', ' Hyperparameters');
            set(self.t2_ext_txt13, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt13, 'FontName', 'Serif');
            set(self.t2_ext_txt13, 'FontSize', 16);
            set(self.t2_ext_txt13, 'FontWeight', 'bold');
            set(self.t2_ext_txt13, 'FontAngle', 'italic');
            set(self.t2_ext_txt13, 'BackgroundColor', self.background_color);  
            set(self.t2_ext_txt13, 'Visible', 'off');

            % frame around vec
            self.t2_ext_frm3 = uicontrol('style','frame');
            set(self.t2_ext_frm3, 'unit', 'pixels', 'position', [20 20 470 380]);
            set(self.t2_ext_frm3, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_ext_frm3, 'BackgroundColor', self.background_color);
            set(self.t2_ext_frm3, 'Visible', 'off');   

            % specification label
            self.t2_ext_txt14 = uicontrol('style', 'text');
            set(self.t2_ext_txt14, 'unit', 'pixels', 'position', [30 365 300 25]);
            set(self.t2_ext_txt14, 'String', ' Autoregressive specification');
            set(self.t2_ext_txt14, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt14, 'FontName', 'Serif');
            set(self.t2_ext_txt14, 'FontSize', 14);
            set(self.t2_ext_txt14, 'FontAngle', 'italic');
            set(self.t2_ext_txt14, 'BackgroundColor', self.background_color);
            set(self.t2_ext_txt14, 'Visible', 'off'); 

            % lag label
            self.t2_ext_txt15 = uicontrol('style', 'text');
            set(self.t2_ext_txt15, 'unit', 'pixels', 'position', [30 337 300 25]);
            set(self.t2_ext_txt15, 'String', ' p:    lags');
            set(self.t2_ext_txt15, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt15, 'FontName', 'Serif');
            set(self.t2_ext_txt15, 'FontSize', 12);
            set(self.t2_ext_txt15, 'FontWeight', 'bold');
            set(self.t2_ext_txt15, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt15, 'Visible', 'off');    

            % lag edit
            self.t2_ext_edt4 = uicontrol('style','edit');
            set(self.t2_ext_edt4, 'unit', 'pixels', 'position', [330 342 140 22]);
            set(self.t2_ext_edt4, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt4, 'Visible', 'off');
            set(self.t2_ext_edt4, 'String', self.user_inputs.tab_2_ext.vec_lags);
            set(self.t2_ext_edt4, 'CallBack', @self.cb_t2_ext_edt4);           

            % pi1 label
            self.t2_ext_txt16 = uicontrol('style', 'text');
            set(self.t2_ext_txt16, 'unit', 'pixels', 'position', [30 287 300 25]);
            set(self.t2_ext_txt16, 'String', ' π₁:  overall tightness');
            set(self.t2_ext_txt16, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt16, 'FontName', 'Serif');
            set(self.t2_ext_txt16, 'FontSize', 12);
            set(self.t2_ext_txt16, 'FontWeight', 'bold');
            set(self.t2_ext_txt16, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt16, 'Visible', 'off'); 

            % pi1 edit
            self.t2_ext_edt5 = uicontrol('style','edit');
            set(self.t2_ext_edt5, 'unit', 'pixels', 'position', [330 292 140 22]);
            set(self.t2_ext_edt5, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt5, 'Visible', 'off');
            set(self.t2_ext_edt5, 'String', self.user_inputs.tab_2_ext.vec_pi1);
            set(self.t2_ext_edt5, 'CallBack', @self.cb_t2_ext_edt5); 

            % pi2 label
            self.t2_ext_txt17 = uicontrol('style', 'text');
            set(self.t2_ext_txt17, 'unit', 'pixels', 'position', [30 262 300 25]);
            set(self.t2_ext_txt17, 'String', ' π₂:  cross-variable shrinkage');
            set(self.t2_ext_txt17, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt17, 'FontName', 'Serif');
            set(self.t2_ext_txt17, 'FontSize', 12);
            set(self.t2_ext_txt17, 'FontWeight', 'bold');
            set(self.t2_ext_txt17, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt17, 'Visible', 'off'); 

            % pi2 edit
            self.t2_ext_edt6 = uicontrol('style','edit');
            set(self.t2_ext_edt6, 'unit', 'pixels', 'position', [330 267 140 22]);
            set(self.t2_ext_edt6, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt6, 'Visible', 'off');
            set(self.t2_ext_edt6, 'String', self.user_inputs.tab_2_ext.vec_pi2);
            set(self.t2_ext_edt6, 'CallBack', @self.cb_t2_ext_edt6); 

            % pi3 label
            self.t2_ext_txt18 = uicontrol('style', 'text');
            set(self.t2_ext_txt18, 'unit', 'pixels', 'position', [30 237 300 25]);
            set(self.t2_ext_txt18, 'String', ' π₃:  lag decay');
            set(self.t2_ext_txt18, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt18, 'FontName', 'Serif');
            set(self.t2_ext_txt18, 'FontSize', 12);
            set(self.t2_ext_txt18, 'FontWeight', 'bold');
            set(self.t2_ext_txt18, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt18, 'Visible', 'off'); 

            % pi3 edit
            self.t2_ext_edt7 = uicontrol('style','edit');
            set(self.t2_ext_edt7, 'unit', 'pixels', 'position', [330 242 140 22]);
            set(self.t2_ext_edt7, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt7, 'Visible', 'off');
            set(self.t2_ext_edt7, 'String', self.user_inputs.tab_2_ext.vec_pi3);
            set(self.t2_ext_edt7, 'CallBack', @self.cb_t2_ext_edt7); 

            % pi4 label
            self.t2_ext_txt19 = uicontrol('style', 'text');
            set(self.t2_ext_txt19, 'unit', 'pixels', 'position', [30 212 300 25]);
            set(self.t2_ext_txt19, 'String', ' π₄:  exogenous slackness');
            set(self.t2_ext_txt19, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt19, 'FontName', 'Serif');
            set(self.t2_ext_txt19, 'FontSize', 12);
            set(self.t2_ext_txt19, 'FontWeight', 'bold');
            set(self.t2_ext_txt19, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt19, 'Visible', 'off'); 

            % pi4 edit
            self.t2_ext_edt8 = uicontrol('style','edit');
            set(self.t2_ext_edt8, 'unit', 'pixels', 'position', [330 217 140 22]);
            set(self.t2_ext_edt8, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt8, 'Visible', 'off');
            set(self.t2_ext_edt8, 'String', self.user_inputs.tab_2_ext.vec_pi4);
            set(self.t2_ext_edt8, 'CallBack', @self.cb_t2_ext_edt8);  
            
            % error correction label
            self.t2_ext_txt20 = uicontrol('style', 'text');
            set(self.t2_ext_txt20, 'unit', 'pixels', 'position', [30 182 300 25]);
            set(self.t2_ext_txt20, 'String', ' Error correction specification');
            set(self.t2_ext_txt20, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt20, 'FontName', 'Serif');
            set(self.t2_ext_txt20, 'FontSize', 14);
            set(self.t2_ext_txt20, 'FontAngle', 'italic');
            set(self.t2_ext_txt20, 'BackgroundColor', self.background_color);
            set(self.t2_ext_txt20, 'Visible', 'off'); 

            % prior label
            self.t2_ext_txt21 = uicontrol('style', 'text');
            set(self.t2_ext_txt21, 'unit', 'pixels', 'position', [30 153 100 25]);
            set(self.t2_ext_txt21, 'String', ' prior:');
            set(self.t2_ext_txt21, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt21, 'FontName', 'Serif');
            set(self.t2_ext_txt21, 'FontSize', 12);
            set(self.t2_ext_txt21, 'FontWeight', 'bold');
            set(self.t2_ext_txt21, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt21, 'Visible', 'off'); 

            % prior radiobuttons
            self.t2_ext_bgr2 = uibuttongroup('unit','pixels', 'Position',[225 80 260 100]);
            set(self.t2_ext_bgr2, 'BorderType', 'none');
            set(self.t2_ext_bgr2, 'BackgroundColor', self.background_color); 
            self.t2_ext_rdb3 = uicontrol(self.t2_ext_bgr2,'Style','radiobutton');
            set(self.t2_ext_rdb3, 'Position',[105 75 220 25]);
            set(self.t2_ext_rdb3, 'String','uninformative');
            set(self.t2_ext_rdb3, 'FontName', 'Serif');
            set(self.t2_ext_rdb3, 'FontSize', 12);
            set(self.t2_ext_rdb3, 'FontWeight', 'bold');
            set(self.t2_ext_rdb3, 'BackgroundColor', self.background_color);
            self.t2_ext_rdb4 = uicontrol(self.t2_ext_bgr2,'Style','radiobutton');
            set(self.t2_ext_rdb4, 'Position',[105 50 220 25]);
            set(self.t2_ext_rdb4, 'String','horseshoe');
            set(self.t2_ext_rdb4, 'FontName', 'Serif');
            set(self.t2_ext_rdb4, 'FontSize', 12);
            set(self.t2_ext_rdb4, 'FontWeight', 'bold');
            set(self.t2_ext_rdb4, 'BackgroundColor', self.background_color);
            self.t2_ext_rdb5 = uicontrol(self.t2_ext_bgr2,'Style','radiobutton');
            set(self.t2_ext_rdb5, 'Position',[105 25 220 25]);
            set(self.t2_ext_rdb5, 'String','selection');
            set(self.t2_ext_rdb5, 'FontName', 'Serif');
            set(self.t2_ext_rdb5, 'FontSize', 12);
            set(self.t2_ext_rdb5, 'FontWeight', 'bold');
            set(self.t2_ext_rdb5, 'BackgroundColor', self.background_color);
            set(self.t2_ext_bgr2, 'Visible', 'off');
            if self.user_inputs.tab_2_ext.prior_type == 1
                set(self.t2_ext_bgr2, 'SelectedObject', self.t2_ext_rdb3);
            elseif self.user_inputs.tab_2_ext.prior_type == 2
                set(self.t2_ext_bgr2, 'SelectedObject', self.t2_ext_rdb4);
            elseif self.user_inputs.tab_2_ext.prior_type == 3
                set(self.t2_ext_bgr2, 'SelectedObject', self.t2_ext_rdb5);                
            end
            set(self.t2_ext_bgr2, 'SelectionChangeFcn', @self.cb_t2_ext_bgr2);

            % correction type label
            self.t2_ext_txt22 = uicontrol('style', 'text');
            set(self.t2_ext_txt22, 'unit', 'pixels', 'position', [30 78 200 25]);
            set(self.t2_ext_txt22, 'String', ' correction type:');
            set(self.t2_ext_txt22, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt22, 'FontName', 'Serif');
            set(self.t2_ext_txt22, 'FontSize', 12);
            set(self.t2_ext_txt22, 'FontWeight', 'bold');
            set(self.t2_ext_txt22, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt22, 'Visible', 'off'); 

            % prior radiobuttons
            self.t2_ext_bgr3 = uibuttongroup('unit','pixels', 'Position',[225 30 260 75]);
            set(self.t2_ext_bgr3, 'BorderType', 'none');
            set(self.t2_ext_bgr3, 'BackgroundColor', self.background_color); 
            self.t2_ext_rdb6 = uicontrol(self.t2_ext_bgr3,'Style','radiobutton');
            set(self.t2_ext_rdb6, 'Position',[105 50 220 25]);
            set(self.t2_ext_rdb6, 'String','general');
            set(self.t2_ext_rdb6, 'FontName', 'Serif');
            set(self.t2_ext_rdb6, 'FontSize', 12);
            set(self.t2_ext_rdb6, 'FontWeight', 'bold');
            set(self.t2_ext_rdb6, 'BackgroundColor', self.background_color);
            self.t2_ext_rdb7 = uicontrol(self.t2_ext_bgr3,'Style','radiobutton');
            set(self.t2_ext_rdb7, 'Position',[105 25 220 25]);
            set(self.t2_ext_rdb7, 'String','reduced-rank');
            set(self.t2_ext_rdb7, 'FontName', 'Serif');
            set(self.t2_ext_rdb7, 'FontSize', 12);
            set(self.t2_ext_rdb7, 'FontWeight', 'bold');
            set(self.t2_ext_rdb7, 'BackgroundColor', self.background_color);
            set(self.t2_ext_bgr3, 'Visible', 'off');
            if self.user_inputs.tab_2_ext.error_correction_type == 1
                set(self.t2_ext_bgr3, 'SelectedObject', self.t2_ext_rdb6);
            elseif self.user_inputs.tab_2_ext.error_correction_type == 2
                set(self.t2_ext_bgr3, 'SelectedObject', self.t2_ext_rdb7);              
            end
            set(self.t2_ext_bgr3, 'SelectionChangeFcn', @self.cb_t2_ext_bgr3);

            % correction type label
            self.t2_ext_txt23 = uicontrol('style', 'text');
            set(self.t2_ext_txt23, 'unit', 'pixels', 'position', [30 25 250 25]);
            set(self.t2_ext_txt23, 'String', ' r: max cointegration rank');
            set(self.t2_ext_txt23, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt23, 'FontName', 'Serif');
            set(self.t2_ext_txt23, 'FontSize', 12);
            set(self.t2_ext_txt23, 'FontWeight', 'bold');
            set(self.t2_ext_txt23, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt23, 'Visible', 'off'); 

            % cointegration rank edit
            self.t2_ext_edt9 = uicontrol('style','edit');
            set(self.t2_ext_edt9, 'unit', 'pixels', 'position', [330 30 140 22]);
            set(self.t2_ext_edt9, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt9, 'Visible', 'off');
            set(self.t2_ext_edt9, 'String', self.user_inputs.tab_2_ext.max_cointegration_rank);
            set(self.t2_ext_edt9, 'CallBack', @self.cb_t2_ext_edt9);

            % varma label
            self.t2_ext_txt24 = uicontrol('style', 'text');
            set(self.t2_ext_txt24, 'unit', 'pixels', 'position', [520 400 300 30]);
            set(self.t2_ext_txt24, 'String', ' VARMA');
            set(self.t2_ext_txt24, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt24, 'FontName', 'Serif');
            set(self.t2_ext_txt24, 'FontSize', 16);
            set(self.t2_ext_txt24, 'FontWeight', 'bold');
            set(self.t2_ext_txt24, 'FontAngle', 'italic');
            set(self.t2_ext_txt24, 'BackgroundColor', self.background_color);  
            set(self.t2_ext_txt24, 'Visible', 'off');

            % frame around varma
            self.t2_ext_frm4 = uicontrol('style','frame');
            set(self.t2_ext_frm4, 'unit', 'pixels', 'position', [510 20 470 380]);
            set(self.t2_ext_frm4, 'ForegroundColor', [0.7 0.7 0.7]);
            set(self.t2_ext_frm4, 'BackgroundColor', self.background_color);
            set(self.t2_ext_frm4, 'Visible', 'off');  

            % specification label
            self.t2_ext_txt25 = uicontrol('style', 'text');
            set(self.t2_ext_txt25, 'unit', 'pixels', 'position', [520 365 300 25]);
            set(self.t2_ext_txt25, 'String', ' Autoregressive specification');
            set(self.t2_ext_txt25, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt25, 'FontName', 'Serif');
            set(self.t2_ext_txt25, 'FontSize', 14);
            set(self.t2_ext_txt25, 'FontAngle', 'italic');
            set(self.t2_ext_txt25, 'BackgroundColor', self.background_color);
            set(self.t2_ext_txt25, 'Visible', 'off'); 

            % lag label
            self.t2_ext_txt26 = uicontrol('style', 'text');
            set(self.t2_ext_txt26, 'unit', 'pixels', 'position', [520 337 300 25]);
            set(self.t2_ext_txt26, 'String', ' p:    lags');
            set(self.t2_ext_txt26, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt26, 'FontName', 'Serif');
            set(self.t2_ext_txt26, 'FontSize', 12);
            set(self.t2_ext_txt26, 'FontWeight', 'bold');
            set(self.t2_ext_txt26, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt26, 'Visible', 'off');    

            % lag edit
            self.t2_ext_edt10 = uicontrol('style','edit');
            set(self.t2_ext_edt10, 'unit', 'pixels', 'position', [820 342 140 22]);
            set(self.t2_ext_edt10, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt10, 'Visible', 'off');
            set(self.t2_ext_edt10, 'String', self.user_inputs.tab_2_ext.varma_lags);
            set(self.t2_ext_edt10, 'CallBack', @self.cb_t2_ext_edt10);     

            % AR coefficients label
            self.t2_ext_txt27 = uicontrol('style', 'text');
            set(self.t2_ext_txt27, 'unit', 'pixels', 'position', [520 312 300 25]);
            set(self.t2_ext_txt27, 'String', ' δ:    AR coefficients');
            set(self.t2_ext_txt27, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt27, 'FontName', 'Serif');
            set(self.t2_ext_txt27, 'FontSize', 12);
            set(self.t2_ext_txt27, 'FontWeight', 'bold');
            set(self.t2_ext_txt27, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt27, 'Visible', 'off');    

            % AR coefficients edit
            self.t2_ext_edt11 = uicontrol('style','edit');
            set(self.t2_ext_edt11, 'unit', 'pixels', 'position', [820 317 140 22]);
            set(self.t2_ext_edt11, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt11, 'Visible', 'off');
            set(self.t2_ext_edt11, 'String', self.user_inputs.tab_2_ext.ar_coefficients);
            set(self.t2_ext_edt11, 'CallBack', @self.cb_t2_ext_edt11); 

            % pi1 label
            self.t2_ext_txt28 = uicontrol('style', 'text');
            set(self.t2_ext_txt28, 'unit', 'pixels', 'position', [520 287 300 25]);
            set(self.t2_ext_txt28, 'String', ' π₁:  overall tightness');
            set(self.t2_ext_txt28, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt28, 'FontName', 'Serif');
            set(self.t2_ext_txt28, 'FontSize', 12);
            set(self.t2_ext_txt28, 'FontWeight', 'bold');
            set(self.t2_ext_txt28, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt28, 'Visible', 'off'); 

            % pi1 edit
            self.t2_ext_edt12 = uicontrol('style','edit');
            set(self.t2_ext_edt12, 'unit', 'pixels', 'position', [820 292 140 22]);
            set(self.t2_ext_edt12, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt12, 'Visible', 'off');
            set(self.t2_ext_edt12, 'String', self.user_inputs.tab_2_ext.varma_pi1);
            set(self.t2_ext_edt12, 'CallBack', @self.cb_t2_ext_edt12); 

            % pi2 label
            self.t2_ext_txt29 = uicontrol('style', 'text');
            set(self.t2_ext_txt29, 'unit', 'pixels', 'position', [520 262 300 25]);
            set(self.t2_ext_txt29, 'String', ' π₂:  cross-variable shrinkage');
            set(self.t2_ext_txt29, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt29, 'FontName', 'Serif');
            set(self.t2_ext_txt29, 'FontSize', 12);
            set(self.t2_ext_txt29, 'FontWeight', 'bold');
            set(self.t2_ext_txt29, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt29, 'Visible', 'off'); 

            % pi2 edit
            self.t2_ext_edt13 = uicontrol('style','edit');
            set(self.t2_ext_edt13, 'unit', 'pixels', 'position', [820 267 140 22]);
            set(self.t2_ext_edt13, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt13, 'Visible', 'off');
            set(self.t2_ext_edt13, 'String', self.user_inputs.tab_2_ext.varma_pi2);
            set(self.t2_ext_edt13, 'CallBack', @self.cb_t2_ext_edt13); 

            % pi3 label
            self.t2_ext_txt30 = uicontrol('style', 'text');
            set(self.t2_ext_txt30, 'unit', 'pixels', 'position', [520 237 300 25]);
            set(self.t2_ext_txt30, 'String', ' π₃:  lag decay');
            set(self.t2_ext_txt30, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt30, 'FontName', 'Serif');
            set(self.t2_ext_txt30, 'FontSize', 12);
            set(self.t2_ext_txt30, 'FontWeight', 'bold');
            set(self.t2_ext_txt30, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt30, 'Visible', 'off'); 

            % pi3 edit
            self.t2_ext_edt14 = uicontrol('style','edit');
            set(self.t2_ext_edt14, 'unit', 'pixels', 'position', [820 242 140 22]);
            set(self.t2_ext_edt14, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt14, 'Visible', 'off');
            set(self.t2_ext_edt14, 'String', self.user_inputs.tab_2_ext.varma_pi3);
            set(self.t2_ext_edt14, 'CallBack', @self.cb_t2_ext_edt14); 

            % pi4 label
            self.t2_ext_txt31 = uicontrol('style', 'text');
            set(self.t2_ext_txt31, 'unit', 'pixels', 'position', [520 212 300 25]);
            set(self.t2_ext_txt31, 'String', ' π₄:  exogenous slackness');
            set(self.t2_ext_txt31, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt31, 'FontName', 'Serif');
            set(self.t2_ext_txt31, 'FontSize', 12);
            set(self.t2_ext_txt31, 'FontWeight', 'bold');
            set(self.t2_ext_txt31, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt31, 'Visible', 'off'); 

            % pi4 edit
            self.t2_ext_edt15 = uicontrol('style','edit');
            set(self.t2_ext_edt15, 'unit', 'pixels', 'position', [820 217 140 22]);
            set(self.t2_ext_edt15, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt15, 'Visible', 'off');
            set(self.t2_ext_edt15, 'String', self.user_inputs.tab_2_ext.varma_pi4);
            set(self.t2_ext_edt15, 'CallBack', @self.cb_t2_ext_edt15);  
            
            % moving average label
            self.t2_ext_txt32 = uicontrol('style', 'text');
            set(self.t2_ext_txt32, 'unit', 'pixels', 'position', [520 182 300 25]);
            set(self.t2_ext_txt32, 'String', 'Moving average specification');
            set(self.t2_ext_txt32, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt32, 'FontName', 'Serif');
            set(self.t2_ext_txt32, 'FontSize', 14);
            set(self.t2_ext_txt32, 'FontAngle', 'italic');
            set(self.t2_ext_txt32, 'BackgroundColor', self.background_color);
            set(self.t2_ext_txt32, 'Visible', 'off'); 

            % residual lag label
            self.t2_ext_txt33 = uicontrol('style', 'text');
            set(self.t2_ext_txt33, 'unit', 'pixels', 'position', [520 153 200 25]);
            set(self.t2_ext_txt33, 'String', ' q:    residual lags');
            set(self.t2_ext_txt33, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt33, 'FontName', 'Serif');
            set(self.t2_ext_txt33, 'FontSize', 12);
            set(self.t2_ext_txt33, 'FontWeight', 'bold');
            set(self.t2_ext_txt33, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt33, 'Visible', 'off'); 

            % residual lag edit
            self.t2_ext_edt16 = uicontrol('style','edit');
            set(self.t2_ext_edt16, 'unit', 'pixels', 'position', [820 158 140 22]);
            set(self.t2_ext_edt16, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt16, 'Visible', 'off');
            set(self.t2_ext_edt16, 'String', self.user_inputs.tab_2_ext.residual_lags);
            set(self.t2_ext_edt16, 'CallBack', @self.cb_t2_ext_edt16); 

            % lambda1 label
            self.t2_ext_txt34 = uicontrol('style', 'text');
            set(self.t2_ext_txt34, 'unit', 'pixels', 'position', [520 128 300 25]);
            set(self.t2_ext_txt34, 'String', ' λ₁:  overall tightness');
            set(self.t2_ext_txt34, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt34, 'FontName', 'Serif');
            set(self.t2_ext_txt34, 'FontSize', 12);
            set(self.t2_ext_txt34, 'FontWeight', 'bold');
            set(self.t2_ext_txt34, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt34, 'Visible', 'off'); 

            % lambda1 edit
            self.t2_ext_edt17 = uicontrol('style','edit');
            set(self.t2_ext_edt17, 'unit', 'pixels', 'position', [820 133 140 22]);
            set(self.t2_ext_edt17, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt17, 'Visible', 'off');
            set(self.t2_ext_edt17, 'String', self.user_inputs.tab_2_ext.lambda1);
            set(self.t2_ext_edt17, 'CallBack', @self.cb_t2_ext_edt17); 

            % lambda2 label
            self.t2_ext_txt35 = uicontrol('style', 'text');
            set(self.t2_ext_txt35, 'unit', 'pixels', 'position', [520 103 300 25]);
            set(self.t2_ext_txt35, 'String', ' λ₂:  cross-variable shrinkage');
            set(self.t2_ext_txt35, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt35, 'FontName', 'Serif');
            set(self.t2_ext_txt35, 'FontSize', 12);
            set(self.t2_ext_txt35, 'FontWeight', 'bold');
            set(self.t2_ext_txt35, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt35, 'Visible', 'off'); 

            % lambda2 edit
            self.t2_ext_edt18 = uicontrol('style','edit');
            set(self.t2_ext_edt18, 'unit', 'pixels', 'position', [820 108 140 22]);
            set(self.t2_ext_edt18, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt18, 'Visible', 'off');
            set(self.t2_ext_edt18, 'String', self.user_inputs.tab_2_ext.lambda2);
            set(self.t2_ext_edt18, 'CallBack', @self.cb_t2_ext_edt18); 

            % lambda3 label
            self.t2_ext_txt36 = uicontrol('style', 'text');
            set(self.t2_ext_txt36, 'unit', 'pixels', 'position', [520 78 300 25]);
            set(self.t2_ext_txt36, 'String', ' λ₃:  lag decay');
            set(self.t2_ext_txt36, 'HorizontalAlignment', 'left');
            set(self.t2_ext_txt36, 'FontName', 'Serif');
            set(self.t2_ext_txt36, 'FontSize', 12);
            set(self.t2_ext_txt36, 'FontWeight', 'bold');
            set(self.t2_ext_txt36, 'BackgroundColor', self.background_color); 
            set(self.t2_ext_txt36, 'Visible', 'off'); 

            % lambda3 edit
            self.t2_ext_edt19 = uicontrol('style','edit');
            set(self.t2_ext_edt19, 'unit', 'pixels', 'position', [820 83 140 22]);
            set(self.t2_ext_edt19, 'HorizontalAlignment', 'center');  
            set(self.t2_ext_edt19, 'Visible', 'off');
            set(self.t2_ext_edt19, 'String', self.user_inputs.tab_2_ext.lambda3);
            set(self.t2_ext_edt19, 'CallBack', @self.cb_t2_ext_edt19);        
            
            % indicate that tab 2 for vec/varma is now created
            self.created_tab_2_ext = true;            
        end

        
        function hide_tab_2_ext(self)

            % hide all controls            
            set(self.t2_ext_txt1, 'Visible', 'off'); 
            set(self.t2_ext_txt2, 'Visible', 'off'); 
            set(self.t2_ext_txt3, 'Visible', 'off'); 
            set(self.t2_ext_txt4, 'Visible', 'off'); 
            set(self.t2_ext_txt5, 'Visible', 'off'); 
            set(self.t2_ext_txt6, 'Visible', 'off'); 
            set(self.t2_ext_txt7, 'Visible', 'off'); 
            set(self.t2_ext_txt8, 'Visible', 'off'); 
            set(self.t2_ext_txt9, 'Visible', 'off'); 
            set(self.t2_ext_txt10, 'Visible', 'off'); 
            set(self.t2_ext_txt11, 'Visible', 'off'); 
            set(self.t2_ext_txt12, 'Visible', 'off'); 
            set(self.t2_ext_txt13, 'Visible', 'off'); 
            set(self.t2_ext_txt14, 'Visible', 'off'); 
            set(self.t2_ext_txt15, 'Visible', 'off'); 
            set(self.t2_ext_txt16, 'Visible', 'off'); 
            set(self.t2_ext_txt17, 'Visible', 'off'); 
            set(self.t2_ext_txt18, 'Visible', 'off'); 
            set(self.t2_ext_txt19, 'Visible', 'off'); 
            set(self.t2_ext_txt20, 'Visible', 'off'); 
            set(self.t2_ext_txt21, 'Visible', 'off'); 
            set(self.t2_ext_txt22, 'Visible', 'off'); 
            set(self.t2_ext_txt23, 'Visible', 'off'); 
            set(self.t2_ext_txt24, 'Visible', 'off'); 
            set(self.t2_ext_txt25, 'Visible', 'off'); 
            set(self.t2_ext_txt26, 'Visible', 'off'); 
            set(self.t2_ext_txt27, 'Visible', 'off'); 
            set(self.t2_ext_txt28, 'Visible', 'off'); 
            set(self.t2_ext_txt29, 'Visible', 'off'); 
            set(self.t2_ext_txt30, 'Visible', 'off'); 
            set(self.t2_ext_txt31, 'Visible', 'off'); 
            set(self.t2_ext_txt32, 'Visible', 'off'); 
            set(self.t2_ext_txt33, 'Visible', 'off'); 
            set(self.t2_ext_txt34, 'Visible', 'off'); 
            set(self.t2_ext_txt35, 'Visible', 'off'); 
            set(self.t2_ext_txt36, 'Visible', 'off'); 
            set(self.t2_ext_frm1, 'Visible', 'off');
            set(self.t2_ext_frm2, 'Visible', 'off');
            set(self.t2_ext_frm3, 'Visible', 'off');
            set(self.t2_ext_frm4, 'Visible', 'off');
            set(self.t2_ext_bgr1, 'Visible', 'off');
            set(self.t2_ext_bgr2, 'Visible', 'off');            
            set(self.t2_ext_bgr3, 'Visible', 'off');            
            set(self.t2_ext_rdb1, 'Visible', 'off');
            set(self.t2_ext_rdb2, 'Visible', 'off');
            set(self.t2_ext_rdb3, 'Visible', 'off');
            set(self.t2_ext_rdb4, 'Visible', 'off');
            set(self.t2_ext_rdb5, 'Visible', 'off');
            set(self.t2_ext_rdb6, 'Visible', 'off');
            set(self.t2_ext_rdb7, 'Visible', 'off');
            set(self.t2_ext_edt1, 'Visible', 'off'); 
            set(self.t2_ext_edt2, 'Visible', 'off'); 
            set(self.t2_ext_edt3, 'Visible', 'off'); 
            set(self.t2_ext_edt4, 'Visible', 'off'); 
            set(self.t2_ext_edt5, 'Visible', 'off'); 
            set(self.t2_ext_edt6, 'Visible', 'off'); 
            set(self.t2_ext_edt7, 'Visible', 'off'); 
            set(self.t2_ext_edt8, 'Visible', 'off'); 
            set(self.t2_ext_edt9, 'Visible', 'off'); 
            set(self.t2_ext_edt10, 'Visible', 'off'); 
            set(self.t2_ext_edt11, 'Visible', 'off'); 
            set(self.t2_ext_edt12, 'Visible', 'off'); 
            set(self.t2_ext_edt13, 'Visible', 'off'); 
            set(self.t2_ext_edt14, 'Visible', 'off');
            set(self.t2_ext_edt15, 'Visible', 'off'); 
            set(self.t2_ext_edt16, 'Visible', 'off'); 
            set(self.t2_ext_edt17, 'Visible', 'off'); 
            set(self.t2_ext_edt18, 'Visible', 'off'); 
            set(self.t2_ext_edt19, 'Visible', 'off'); 
            set(self.t2_ext_cbx1, 'Visible', 'off');
            set(self.t2_ext_cbx2, 'Visible', 'off');
            set(self.t2_ext_cbx3, 'Visible', 'off'); 
            
            % update tab color
            set(self.tab_pbt2, 'BackgroundColor', self.backtabs_color);
        end
        
        
        function show_tab_2_ext(self)

            % show all controls
            set(self.t2_ext_txt1, 'Visible', 'on');     
            set(self.t2_ext_txt2, 'Visible', 'on');   
            set(self.t2_ext_txt3, 'Visible', 'on');    
            set(self.t2_ext_txt4, 'Visible', 'on');    
            set(self.t2_ext_txt5, 'Visible', 'on');    
            set(self.t2_ext_txt6, 'Visible', 'on');  
            set(self.t2_ext_txt7, 'Visible', 'on');  
            set(self.t2_ext_txt8, 'Visible', 'on');  
            set(self.t2_ext_txt9, 'Visible', 'on');  
            set(self.t2_ext_txt10, 'Visible', 'on');  
            set(self.t2_ext_txt11, 'Visible', 'on'); 
            set(self.t2_ext_txt12, 'Visible', 'on'); 
            set(self.t2_ext_txt13, 'Visible', 'on'); 
            set(self.t2_ext_txt14, 'Visible', 'on'); 
            set(self.t2_ext_txt15, 'Visible', 'on'); 
            set(self.t2_ext_txt16, 'Visible', 'on'); 
            set(self.t2_ext_txt17, 'Visible', 'on'); 
            set(self.t2_ext_txt18, 'Visible', 'on'); 
            set(self.t2_ext_txt19, 'Visible', 'on'); 
            set(self.t2_ext_txt20, 'Visible', 'on'); 
            set(self.t2_ext_txt21, 'Visible', 'on'); 
            set(self.t2_ext_txt22, 'Visible', 'on'); 
            set(self.t2_ext_txt23, 'Visible', 'on'); 
            set(self.t2_ext_txt24, 'Visible', 'on'); 
            set(self.t2_ext_txt25, 'Visible', 'on'); 
            set(self.t2_ext_txt26, 'Visible', 'on'); 
            set(self.t2_ext_txt27, 'Visible', 'on'); 
            set(self.t2_ext_txt28, 'Visible', 'on'); 
            set(self.t2_ext_txt29, 'Visible', 'on'); 
            set(self.t2_ext_txt30, 'Visible', 'on'); 
            set(self.t2_ext_txt31, 'Visible', 'on'); 
            set(self.t2_ext_txt32, 'Visible', 'on'); 
            set(self.t2_ext_txt33, 'Visible', 'on'); 
            set(self.t2_ext_txt34, 'Visible', 'on'); 
            set(self.t2_ext_txt35, 'Visible', 'on'); 
            set(self.t2_ext_txt36, 'Visible', 'on'); 
            set(self.t2_ext_frm1, 'Visible', 'on');
            set(self.t2_ext_frm2, 'Visible', 'on');
            set(self.t2_ext_frm3, 'Visible', 'on');
            set(self.t2_ext_frm4, 'Visible', 'on');
            set(self.t2_ext_bgr1, 'Visible', 'on');
            set(self.t2_ext_bgr2, 'Visible', 'on');            
            set(self.t2_ext_bgr3, 'Visible', 'on');             
            set(self.t2_ext_rdb1, 'Visible', 'on');
            set(self.t2_ext_rdb2, 'Visible', 'on');
            set(self.t2_ext_rdb3, 'Visible', 'on');
            set(self.t2_ext_rdb4, 'Visible', 'on');
            set(self.t2_ext_rdb5, 'Visible', 'on');
            set(self.t2_ext_rdb6, 'Visible', 'on');
            set(self.t2_ext_rdb7, 'Visible', 'on');
            set(self.t2_ext_edt1, 'Visible', 'on'); 
            set(self.t2_ext_edt2, 'Visible', 'on'); 
            set(self.t2_ext_edt3, 'Visible', 'on');
            set(self.t2_ext_edt4, 'Visible', 'on');
            set(self.t2_ext_edt5, 'Visible', 'on');
            set(self.t2_ext_edt6, 'Visible', 'on');
            set(self.t2_ext_edt7, 'Visible', 'on');
            set(self.t2_ext_edt8, 'Visible', 'on');
            set(self.t2_ext_edt9, 'Visible', 'on');
            set(self.t2_ext_edt10, 'Visible', 'on');
            set(self.t2_ext_edt11, 'Visible', 'on');
            set(self.t2_ext_edt12, 'Visible', 'on');
            set(self.t2_ext_edt13, 'Visible', 'on');
            set(self.t2_ext_edt14, 'Visible', 'on');
            set(self.t2_ext_edt15, 'Visible', 'on');
            set(self.t2_ext_edt16, 'Visible', 'on');
            set(self.t2_ext_edt17, 'Visible', 'on');
            set(self.t2_ext_edt18, 'Visible', 'on');
            set(self.t2_ext_edt19, 'Visible', 'on');
            set(self.t2_ext_cbx1, 'Visible', 'on');
            set(self.t2_ext_cbx2, 'Visible', 'on');
            set(self.t2_ext_cbx3, 'Visible', 'on');
        end
        
        
        function cb_t2_ext_bgr1(self, hObject, callbackdata)
           if get(self.t2_ext_bgr1, 'SelectedObject') == self.t2_ext_rdb1
               self.user_inputs.tab_2_ext.model = 1;
           elseif get(self.t2_ext_bgr1, 'SelectedObject') == self.t2_ext_rdb2
               self.user_inputs.tab_2_ext.model = 2;
           end
        end
        
        function cb_t2_ext_edt1(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.iterations = get(self.t2_ext_edt1, 'String');
        end 
        
        function cb_t2_ext_edt2(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.burnin = get(self.t2_ext_edt2, 'String');
        end 
        
        function cb_t2_ext_edt3(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.model_credibility = get(self.t2_ext_edt3, 'String');
        end 
        
        function cb_t2_ext_cbx1(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.constant = logical(get(self.t2_ext_cbx1, 'Value'));
        end 
        
        function cb_t2_ext_cbx2(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.trend = logical(get(self.t2_ext_cbx2, 'Value'));
        end 
        
        function cb_t2_ext_cbx3(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.quadratic_trend = logical(get(self.t2_ext_cbx3, 'Value'));
        end 
        
        function cb_t2_ext_edt4(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.vec_lags = get(self.t2_ext_edt4, 'String');
        end 

        function cb_t2_ext_edt5(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.vec_pi1 = get(self.t2_ext_edt5, 'String');
        end 

        function cb_t2_ext_edt6(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.vec_pi2 = get(self.t2_ext_edt6, 'String');
        end 

        function cb_t2_ext_edt7(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.vec_pi3 = get(self.t2_ext_edt7, 'String');
        end 

        function cb_t2_ext_edt8(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.vec_pi4 = get(self.t2_ext_edt8, 'String');
        end 
        
        function cb_t2_ext_bgr2(self, hObject, callbackdata)
           if get(self.t2_ext_bgr2, 'SelectedObject') == self.t2_ext_rdb3
               self.user_inputs.tab_2_ext.prior_type = 1;
           elseif get(self.t2_ext_bgr2, 'SelectedObject') == self.t2_ext_rdb4
               self.user_inputs.tab_2_ext.prior_type = 2;
           elseif get(self.t2_ext_bgr2, 'SelectedObject') == self.t2_ext_rdb5
               self.user_inputs.tab_2_ext.prior_type = 3;               
           end
        end

        function cb_t2_ext_bgr3(self, hObject, callbackdata)
           if get(self.t2_ext_bgr3, 'SelectedObject') == self.t2_ext_rdb6
               self.user_inputs.tab_2_ext.error_correction_type = 1;
           elseif get(self.t2_ext_bgr3, 'SelectedObject') == self.t2_ext_rdb7
               self.user_inputs.tab_2_ext.error_correction_type = 2;             
           end
        end

        function cb_t2_ext_edt9(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.max_cointegration_rank = get(self.t2_ext_edt9, 'String');
        end 
        
        function cb_t2_ext_edt10(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.varma_lags = get(self.t2_ext_edt10, 'String');
        end 

        function cb_t2_ext_edt11(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.ar_coefficients = get(self.t2_ext_edt11, 'String');
        end 

        function cb_t2_ext_edt12(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.varma_pi1 = get(self.t2_ext_edt12, 'String');
        end 

        function cb_t2_ext_edt13(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.varma_pi2 = get(self.t2_ext_edt13, 'String');
        end 

        function cb_t2_ext_edt14(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.varma_pi3 = get(self.t2_ext_edt14, 'String');
        end 

        function cb_t2_ext_edt15(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.varma_pi4 = get(self.t2_ext_edt15, 'String');
        end 

        function cb_t2_ext_edt16(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.residual_lags = get(self.t2_ext_edt16, 'String');
        end 

        function cb_t2_ext_edt17(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.lambda1 = get(self.t2_ext_edt17, 'String');
        end 

        function cb_t2_ext_edt18(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.lambda2 = get(self.t2_ext_edt18, 'String');
        end 

        function cb_t2_ext_edt19(self, hObject, callbackdata)
            self.user_inputs.tab_2_ext.lambda3 = get(self.t2_ext_edt19, 'String');
        end         
        
    end  
     
end