% EXAMPLE Nonlinear Limit state
clear probdata femodel analysisopt gfundata randomfield systems results output_filename

addpath(genpath(strcat(pwd,'\ferumcore')))

% Probability Data
probdata.marg(1,:) =  [ 1  5   5*0.3   5  0 0 0 0 0]; % dist of Var 1
probdata.marg(2,:) =  [ 1  3   3*0.3   3  0 0 0 0 0]; % dist of Var 2
probdata.correlation = eye(2); % Correlation matrix                     
probdata.parameter = distribution_parameter(probdata.marg);
               
% Analysis Options
analysisopt.ig_max    = 100;
analysisopt.il_max    = 5;
analysisopt.e1        = 0.001;
analysisopt.e2        = 0.001; 
analysisopt.step_code = 0;
analysisopt.grad_flag = 'FFD';
analysisopt.sim_point = 'dspt';
analysisopt.stdv_sim  = 1;
analysisopt.num_sim   = 100000;
analysisopt.target_cov = 0.0125;

% Limit State function (gfun) data
gfundata(1).thetag = d; % design variables
gfundata(1).evaluator = 'basic';
gfundata(1).type = 'expression';
gfundata(1).parameter = 'yes';
gfundata(1).expression = 'g_func(x(1),x(2),gfundata(1).thetag(1),gfundata(1).thetag(2))';
gfundata(1).dgthetag = {'1','1'};

% Other Data
femodel = 0;
randomfield.mesh = 0;

% Run modified FORM analysis from FERUM
ferum_form;