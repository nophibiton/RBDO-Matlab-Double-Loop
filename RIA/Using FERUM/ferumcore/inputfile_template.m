%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a template inputfile for FERUM. The data are given in 6 data structures with  various data fields.               %                            %
%                                                                                                                          %
% The 6 data structures are:                                                                                               %
%  - probdata    (uncertainty characterization)                                                                            %
%  - analysisopt (analysis options)                                                                                        %
%  - femodel     (finite element model data)                                                                               % 
%        - linear and non linear code                                                                                      %
%        - fedeas code                                                                                                     %                                                                                                                         %
%  - gfundata    (limit-state function information)                                                                        %
%  - randomfield (random field characterization)                                                                           %
%  - systems     (systems reliability analysis)                                                                            %
%                                                                                                                          %
%                                                                                                                          %  
%----------------------------------------------------                                                                      %
% THE REQUIRED RULES TO PERFORM ANALYSIS WITH FERUM %                                                                      %
%----------------------------------------------------                                                                      %
% The syntax used here is as follows:                                                                                      %
% (something)  indicates that a number is required.                                                                        %
% 'something'  indicates that a string is required.                                                                        %
%                                                                                                                          %
% Some of the data fields in these 6 structures may not have to be specified. It depends on the problem which is analysed. %
%                                                                                                                          %
% 1.Clear all possible data in the MATLAB workspace                                                                        %
%                                                                                                                          %
% 2.Define a file to store the output results                                                                              %
%    - If not defined the output of the analysis won't be saved.                                                           %
%    - If the file already exists, the results of the analysis will be written at the end of the results obtained          %
%      by a previous analysis.                                                                                             %
%                                                                                                                          %
% 3.Data fields in Probdata                                                                                                %
%    - Random variables can be defined by their marginal distribution and their correlation matrix                         %
%          * Available distributions are listed below.                                                                     %
%          * Distributions can be defined either giving the mean and standard deviation                                    % 
%            either giving their associated parameters.                                                                    %
%    - If the random field option is used, number of random variables should here be equal to 'nterm' from                 %
%      the random field option defined in data fields in random fields. Further, the random variables should               %
%      then have the N(0,1) distribution, that is, normal with zero mean and unit standard deviation.                      %
%                                                                                                                          %
% 4.Data fields in Analysisopt                                                                                             %
% The user can define:                                                                                                     %                                                                                              %
%    - Parameters in search algorithm (necessary for FORM, SORM and SYSTEM analysis).                                      %
%    - Parameters for simulation analysis (necessary for simulation analysis).                                             %
%    - Parameters for inverse reliability analysis (necessary for inverse reliability analysis).                           %
% If these parameters are NOT defined the associated analysis will failed.                                                 %
%                                                                                                                          %
% 5.Data fields in Femodel                                                                                                 %
%   - Set femodel = 0 if femodel option is NOT used.                                                                       %
%   - Available element types and their corresponding node and parameters are listed below.                                %
%   - Random quantities (load magnitude, element parameter...) are DUMMY when the user define "femodel.el",                % 
%     "femodel.loading" or "femodel.nodal_spring".                                                                         %
%   - There is a special structure "fem.id" to define the random quantities associating the random quantity with           % 
%     a random variable defined in Probdata.                                                                               %
%                                                                                                                          %
% 6.Data fields in Gfundata                                                                                                %
% It depends on the analysis type:                                                                                         %
%   - Reliability analysis                                                                                                 %
%   - Inverse reliability analysis: the search deterministic parameter must be named "theta"                               %
%                                                                                                                          %
% 7.Data fields in randomfield                                                                                             %
%   - Set randomfield.mesh = 0 if random field option is NOT used.                                                         %
%                                                                                                                          %
% 8.Data fields in systems                                                                                                 %
%     system.system                                                                                                        %
%         - Cutset Formulation:   (1)'-': compliment event                                                                 %
%                                 (2)'0': divider between cutsets                                                            %
%         - Series system: system.system = [-no.of events] e.g. [-12] : series system with 12 events                       % 
%         - Parallel system: system.system = [no.of events] e.g. [13] : parallel system with 13 events                     %
%         - system.system=[1],[-1],[0],[] will be considered as component reliability problem                              %
%                                                                                                                          %
%                                                                                                                          %
%                                                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR POSSIBLE OLD DATA IN THE MATLAB WORKSPACE
clear probdata femodel analysisopt gfundata randomfield systems results output_filename


% Define the name and type of the file to store the output of the analysis:
output_filename = 'filename.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FIELDS IN 'PROBDATA':                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Marginal distributions for each random variable:
probdata.marg(1,:) =  [ (type)  (mean)   (std.dev.)  (startpoint) (p1) (p2) (p3) (p4) (Input_type)]; 																			  
probdata.marg(2,:) =  [ (type)  (mean)   (std.dev.)  (startpoint) (p1) (p2) (p3) (p4) (Input_type)];
...
   
%Notes:
%    - Each field (mean, std.dev., p1, p2, p3, p4, Input_type) must be fill in. If not used input a dummy value.
%
%    - Input_type = 0 when distribution defined thanks to the mean and std.dev.
%      Input_type = 1 when distribution defined thanks to the distribution parameters pi
%
%    - For Type III Smallest value marginal distribution
%      You must give the value of epsilon as p3 when using the mean  and std.dev. input  
%
%    - For Beta marginal distribution , you have to give the value of a as p3 and b as p4 when using 
%      the mean and std. dev. input. 
%      User will not be able to use this beta distribution while using a student version of MATLAB.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              FERUM distributions library                                          %
% Type: 1 = Normal distribution                                                                                     %
%       2 = Lognormal distribution                                                                                  %
%       3 = Gamma                                                                                                   %
%       4 = Shifted Exponential marginal distribution                                                               %
%       5 = Shifted Rayleigh marginal distribution                                                                  %
%       6 = Uniform distribution                                                                                    %
%       7 = Beta                                                                                                    %
%       8 = Chi-square                                                                                              %
%                                                                                                                   %
%      11 = Type I Largest Value marginal distribution                                                              %
%      12 = Type I Smallest Value marginal distribution                                                             %
%      13 = Type II Largest Value marginal distribution                                                             %
%      14 = Type III Smallest Value marginal distribution                                                           %
%      15 = Gumbel (same as type I largest value)                                                                   %  
%      16 = Weibull marginal distribution (same as Type III Smallest Value marginal distribution with epsilon = 0 ) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the parameters,the mean and standard deviation associated with the distribution of each random variable
probdata.parameter = distribution_parameter(probdata.marg);
% Correlation matrix (square matrix with dimension equal to number of r.v.'s)
probdata.correlation=[1.0  0.2  0.3  ;
                      0.2  1.0  0.4  ;
                      0.3  0.4  1.0 ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FIELDS IN 'ANALYSISOPT':                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters in search algorithm
analysisopt.ig_max = 100;       % Maximum number of global iterations allowed in the search algorithm
analysisopt.il_max = 5;         % Maximum number of line iterations allowed in the search algorithm
analysisopt.e1 = 0.001;         % Tolerance on how close design point is to limit-state surface
analysisopt.e2 = 0.001;         % Tolerance on how accurately the gradient points towards the origin
analysisopt.step_code = 0;      % 0: step size by Armijo rule, otherwise: given value (0 < s <= 1) is the step size.
analysisopt.grad_flag = 'DDM';  % 'DDM': direct differentiation, 'FFD': forward finite difference

% Simulation analysis
analysisopt.sim_point = 'dspt'; % 'dspt': design point, 'origin': origin in standard normal space
analysisopt.stdv_sim  = 1;      % Standard deviation of sampling distribution
analysisopt.num_sim   = 1000;   % Number of simulations
analysisopt.target_cov = 0.05;  % Target coefficient of variation of failure probability estimate

% Inverse FORM analysis
analysisopt.beta_target = 2;    % Target reliability index
analysisopt.e3 = 0.001;         % Tolerance on target beta 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FIELDS IN 'FEMODEL' for any of FERUMs FE codes: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of degrees of freedom per node
femodel.ndf = 2;

% Nodal information
femodel.node(1,:) = [  (x-coord)   (y-coord)  ];
femodel.node(2,:) = [  (x-coord)   (y-coord)  ];
...

% Element information
femodel.el(1,:) = [ (type) (node1) (node2) (node3) .. (parameter1) (parameter2) .. ];
femodel.el(2,:) = [ (type) (node1) (node2) (node3) .. (parameter1) (parameter2) .. ];
...
% Available element types and their corresponding node and parameter list:
% Elastic truss:    (1) (node1) (node2) (Elastic Modulus) (Area)
% Elastic beam:     (2) (node1) (node2) (E) (I) (A)
% Elastic quad4:    (3) (node1) (node2) (node3) (node4) (E) (nu) (t) (PS)
% J2 plastic truss: (4) (node1) (node2) (E) (A)  (nu) (sy) (Hi) (Hk) 
% J2 plastic quad4: (5) (node1) (node2) (node3) (node4) (E) (nu) (t) (sy) (Hi) (Hk) (PS)
%                   only PS=1 option available for J2 plastic quad4 at this time.
%
% Note: nodei=global node number coinciding element node i,
%       E=elastic modulus, A=area, I=moment of inertia, nu=Poisson's ratio, t=thickness,
%       sy=yield stress, Hi=isotropic hardening modulus, Hk=kinematic hardening modulus,
%       PS=1 for plane-strain analysis, PS=2 for plane stress analysis

% Nodal loads for linear FE analysis 
femodel.loading(1,:) = [ (nodenumber)  (magnitude)  (direction)  (factor)  ];
femodel.loading(2,:) = [ (nodenumber)  (magnitude)  (direction)  (factor)  ];
...
% Nodal loads for nonlinear FE analysis
femodel.loading(1,:) = [ (nodenumber)  (magnitude)  (direction)  (factor)  (time) (loadfactor)  (time) (loadfactor)  ..  ];
femodel.loading(2,:) = [ (nodenumber)  (magnitude)  (direction)  (factor)  (time) (loadfactor)  (time) (loadfactor)  ..  ];
...

% Direction:  1 : x-direction
%             2 : y-direction
%             3 : clockwise moment
% A load in NEGATIVE direction is specified by a minus sign on the direction specification.
% The 'factor' variable can be used when several loads (with different magnitude) are based 
% on the same random variable.


% Concentrated nodal spring stiffnesses
% Set nodal_spring=0 if no springs are included in the model
femodel.nodal_spring(1,:) = [ (nodenumber)   (magnitude)   (direction)  ];
femodel.nodal_spring(2,:) = [ (nodenumber)   (magnitude)   (direction) ];
...

% Fixed degrees of freedoms;  0=free, 1=fixed
femodel.fixed_dof(1,:) = [ (node)  (0/1) .. (0/1) ];
femodel.fixed_dof(2,:) = [ (node)  (0/1) .. (0/1) ];
...

% ID array to identify where the r.v.'s enter into the FE model:
femodel.id(1,:) = [  (r.v.number) (phys.meaning) (load/node/el number)   ];
femodel.id(2,:) = [  (r.v.number) (phys.meaning) (load/node/el number)   ];
...
% Phys.meaning:  - 1:  Nodal load
%                - 2:  Young's modulus (E)
%                - 3:  Moment of inertia (I)
%                - 4:  Cross sectional area (A)
%                - 5:  Nodal spring stiffness
%                - 6:  Poisson's ratio (nu)
%                - 7:  Thickness of 2D element
%                - 8:  Yield limit (s0)
%                - 9:  Isotropic hardening parameter (Hi)
%                - 10: Kinematic hardening parameter (Hk)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FIELDS IN 'FEMODEL' for FEDEAS gfun evaluator:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
femodel.Model;
femodel.ElData;
femodel.Load;
femodel.SolStrat;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FIELDS IN 'GFUNDATA' (one structure per gfun):  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------
%Reliability analysis -
%----------------------

% User can define limit-state function parameters
% if there is no parameter in the limit-state function
gfundata(1).parameter = 'no';
% if there is parameters
gfundata(1).parameter = 'yes';
    % define the parameter values
    gfundata(1).thetag = [ (par1) (par2) ..];

% Type of limit-state function evaluator (Alternatives: 'basic', 'FERUMlinearfecode', 'FERUMnonlinearfecode', 'fedeas')

gfundata(1).evaluator = 'basic';
% Type of limit-state function
% (Alternatives: 'expression','matlabfile')
    %In case of 'expression'
    gfundata(1).type = 'expression';
        % if there is no parameter in the limit-state definition
	    gfundata(1).expression = '1.0 - x(2)/(1000*x(3)) - (x(1)/(200*x(3)))^2';
        % if there is parameters in the limit-state definition
        gfundata(1).expression = 'gfundata(1).thetag(1) - x(2)/(1000*x(3)) - (x(1)/(200*x(3)))^2';
        % Give explicit gradient expressions with respect to the involved quantities (in the order x(1), x(2), ...) if DDM is used:
        gfundata(1).dgdq = { '-x(1)/(20000*x(3)^2)' ;
                             '-1/(1000*x(3))' ;
                             '(20*x(2)*x(3)+x(1)^2)/(20000*x(3)^2)'};
        % Give explicit gradient expressions with respect to the limit-state function parameters 
        %(in the order thetag(1), thetag(2), ...) if DDM is used:
        gfundata(1).dgthetag = {'1'};
    
    %In case of 'matlabfile'
    %NOTE THAT THIS OPTION IS NOT AVAILABLE YET%
    gfundata(1).type = 'matlabfile';
    %user have to define is limit-state function using the user_lsf.m file and follow the instructions
 
gfundata(1).evaluator = 'FERUMlinearfecode'; % gfundata(1).evaluator = 'FERUMnonlinearfecode'; gfundata(1).evaluator = 'fedeas';
% Type of limit-state function
% (Alternatives: 'displacementlimit','matlabfile')
    %In case of 'displacementlimit'
    gfundata(1).type = 'displacementlimit';
        % Identification of response quantity
        gfundata(1).resp = [  (nodenumber)    (dofnumber)  ];
        % Design limit that should not be exceeded by the specified response quantity
        gfundata(1).lim = (displacementlimit);
        
    %In case of 'matlabfile'
    %NOTE THAT THIS OPTION IS NOT AVAILABLE YET%
    gfundata(1).type = 'matlabfile';
    %user have to define is limit-state function using the user_lsf.m file and follow the instructions


%------------------------------
%Inverse reliability analysis -
%------------------------------
% The deterministic parameter which has to be found must be defined as theta or theta(1)

% Define starting value for the deterministic parameter
gfundata(1).deterministic_parameter_start_value = 2000 ;

% Type of limit-state function evaluator
% (Alternatives: 'basic', 'FERUMlinearfecode', 'FERUMnonlinearfecode', 'fedeas')
gfundata(1).evaluator = 'basic';

% Type of limit-state function
% (Alternatives: 'expression', 'matlabfile')
    %In case of 'expression'
        gfundata(1).type = 'expression';
    
        gfundata(1).expression = '1.0 - theta(1)/(1000*x(2)) - (x(1)/(200*x(2)))^2';
        % If DDM is used, give explicit derivative expressions with respect to :  - the involved quantities
        gfundata(1).dgdq = { '-1/(1000*x(3))' ;
						     '(20*x(2)*x(3)+x(1)^2)/(20000*x(3)^2)'};
        %                                                                         - the deterministic parameter
        gfundata(1).dgdthetadeter = { '-1/(1000*x(2))'};

    %In case of 'matlabfile'   
    %NOTE THAT THIS OPTION IS NOT AVAILABLE YET%
        gfundata(1).type = 'matlabfile';
        %user have to define is limit-state function using the user_lsf.m file and follow the instructions
    
 gfundata(1).evaluator = 'FERUMlinearfecode'; % gfundata(1).evaluator = 'FERUMnonlinearfecode'; gfundata(1).evaluator = 'fedeas';
% Type of limit-state function
% (Alternatives: 'displacementlimit','matlabfile')
    %In case of 'displacementlimit'
    gfundata(1).type = 'displacementlimit';
        % Identification of response quantity
        gfundata(1).resp = [  (nodenumber)    (dofnumber)  ];
        % Design limit that should not be exceeded by the specified response quantity
        gfundata(1).lim = 'theta(1)';
        
    %In case of 'matlabfile'
    %NOTE THAT THIS OPTION IS NOT AVAILABLE YET%
    gfundata(1).type = 'matlabfile';
    %user have to define is limit-state function using the user_lsf.m file and follow the instructions

     



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FIELDS IN 'RANDOMFIELD':                             %
% Input for random field analysis (random Young's modulus). %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randomfield.mesh = [ (x1) (x2) (x3) .. ]; % Nodes in random field discretization (should be more coarse than the finite element mesh)
randomfield.meanfunc = '200000';          % Mean function
randomfield.covariancefunc = '20000^2*exp(-(abs(xn(i)-x(j)))^2/1500^2)'; % Covariance function in terms of element
                                                                         % midpoints 'x' and random field nodes 'xn'. 
randomfield.nterm      = 5;               % Number of eigenvalues to include in the realization (max: xn)
randomfield.fieldtype  = 1;               % 1=Gaussian, 2=Lognormal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FIELDS IN 'SYSTEMS':                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note:
% (1) For series system: system.system=[-no.of events]
%      e.g. [-12] : series system with 12 events
% (2) For parallel system: system.system=[no.of events]
%      e.g. [13] : parallel system with 13 events
% (3) system.system=[1],[-1],[0],[] will be considered as component reliability problem
% (4) For a cutset formulation
system.system = [1 4 0 1 5 0 1 6 0 2 7 0 0 2 -5 0 2 8 0 3 9 0]; % Cutset Formulation
% This 'system' corresonds to (e1e4)U(e1e5)U(e1e6)U(e2e7)U(e2e5c)U(e2e8)U(e3e9)
% ([0 -2 7 0 3 2 5 0...]: C1(-e2,e7), C2(e3,e2,e5)...
% ([1 3 -5 0 0 2 3 9...]: C1(e1,e3,-e5), C2(e2,e3,e9)...)
% 0's are just separating different cutsets.
% 0's can be repeated and be placed in the first or last place.

system.scis_max = 20000;  % maximum no. of simulation for each scis (Default : 20000)
system.scis_min = 1000;   % minimum no. of simulation for each scis (Default : 1000)
system.cov_max  = 0.05;   % target c.o.v. of each scis (Default : 0.05)

