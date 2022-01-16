%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite Element Reliability using Matlab, FERUM, Version 3.0       %
%                                                                   %
% This program is free software; you can redistribute it and/or     %
% modify it under the terms of the GNU General Public License       %
% as published by the Free Software Foundation; either version 2    %
% of the License, or (at your option) any later version.            %
%                                                                   %
% This program is distributed in the hope that it will be useful,   %
% but WITHOUT ANY WARRANTY; without even the implied warranty of    %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     %
% GNU General Public License for more details.                      %
%                                                                   %
% A copy of the GNU General Public License is found in the file     %
% <gpl.txt> following this collection of program files.             %
%                                                                   %
% Developed under the sponsorship of the Pacific                    %
% Earthquake Engineering (PEER) Center.                             %
%                                                                   %
% For more information, visit: http://www.ce.berkeley.edu/~haukaas  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc
% disp('   ____________________________________________');
% disp('  | Welcome to FERUM Version 3.1 beta          |');
% disp('  | (Finite Element Reliability Using Matlab)  |');
% disp('  | For more information, visit:               |');
% disp('  | http://www.ce.berkeley.edu/~haukaas.       |');
% disp('  | Note: All the analysis options below       |');
% disp('  | assumes that necessary data are            |');
% disp('  | available in the current Matlab workspace. |');
% disp('   ____________________________________________');
% disp(' ');
% disp(' ');

% test_ouputfile = exist('output_filename','var');
test_ouputfile = 0;
% if test_ouputfile == 1
%     filename = output_filename;
%     disp(['  The results of the analysis will be store in your current Matlab directory in the file: ',filename,''])
% else
%     disp('   WARNING: No defined output file, the results of the analysis will not be a store in a file')
% end

% disp(' ');
% disp('  0: Exit');
% disp('  1: Help');
% disp('  2: Input Check');
% disp('  3: FORM Analysis');
% disp('  4: SORM Analysis');
% disp('  5: Importance Sampling Simulation Analysis');
% disp('  6: Systems Analysis');
% disp('  7: Inverse FORM Analysis');
% disp(' ');
% disp(' ');
% 
% analysistype = input('  CHOOSE OPTION FROM THE LIST ABOVE: ');

analysistype = 3;
switch analysistype
     
   
case 0 % ---- EXIT ----------------------------------------------------------------------------
   disp(' ');
   disp('  Bye, bye.');
   disp(' ');
   
      
case 1 % ---- HELP ----------------------------------------------------------------------------
   clc
   disp(' ');
   disp('  FERUM HELP');
   disp('  If you are new to FERUM, you are recommended to visit the web page:');
   disp('  http://www.ce.berkeley.edu/~haukaas/matlab_toolbox');
   disp(' ');
   disp('  TO RUN E.G. A FINITE ELEMENT RELIABILITY ANALYSIS, DO THE FOLLOWING:');
   disp('  1. Specify necessary parameters in your current Matlab workspace. ');
   disp('     (The format of the input can be found in the template.m file.');
   disp('     If you want to try one of the provided example inputfiles,');
   disp('     simply read the file into your workspace by writing the file');
   disp('     name and press enter. ');
   disp('  2. Start FERUM (the shell program) by issuing the command "ferum".');
   disp('  3. Choose the alternative that fits your purpose.');
   disp('  4. If the analysis converged; view the results as indicated on the screen.');
   disp(' ');
   
   
   
case 2 % ---- INPUT CHECK ----------------------------------------------------------------------------
   disp(' ');
   disp('  Currently this option only checks the random field.');
   disp(' ');
   positionrandomfield(0,femodel,randomfield,'yes');
   
   
   
case 3 % ---- FORM ----------------------------------------------------------------------------
   
   % Clear screen and display message
%    disp('FORM analysis is running, please wait... (Ctrl+C breaks)')
   nrv = size(probdata.marg,1);
   % Run the analysis
   t = clock;
   [formresults] = form(1,probdata,analysisopt,gfundata,femodel,randomfield);
   results.form = formresults;
   time = etime(clock,t); 
   % Display results
%    clc
   if test_ouputfile == 1
      diary(filename)
   end
%    disp(' ');
%    disp(['          ################################################################################################'])
%    disp(['          #                         RESULTS FROM RUNNING FORM RELIABILITY ANALYSIS                       #'])
%    disp(['          ################################################################################################'])
%    disp(' ');
%    disp(['Number of iterations: ',num2str(formresults.iter),''])
%    disp([' '])
%    disp(['Time to complete the analysis:    ',num2str(time),''])
%    disp(' ');
%    disp(['Reliability index beta1: ',num2str(formresults.beta1),''])
%    disp([' '])
%    disp(['Failure probability pf1: ',num2str(formresults.pf1,'%0.5e'),''])
%    disp([' '])
%    diary off
   %formresults.sensi = disp_sensitivities (1,formresults,probdata,gfundata);
   %disp_sensitivities (1,formresults,probdata,gfundata)
   
   if nrv < 20
   
   if test_ouputfile == 1
      diary(filename)
   end
%    disp_sensitivities (1,formresults,probdata,gfundata)
%    diary off
%    
%    disp(['.................................................'])
%    disp('The following parameters are now available in your current workspace:')
%    disp('   formresults.iter  = Number of iterations')
%    disp('   formresults.beta1  = Reliability index beta from FORM analysis')
%    disp('   formresults.pf1   = Failure probability pf1')
%    disp('   formresults.dsptu = Design point u_star')
%    disp('   formresults.dsptx = Design point in original space')
%    disp('   formresults.alpha = Alpha vector')
%    disp('   formresults.imptg = Importance vector gamma')
%    disp('   formresults.gfcn  = Recorded  values of the limit-state function during search')
%    disp('   formresults.stpsz = Recorded step size values during search')
%    disp([' '])
%    disp('   formresults.beta_sensi_thetaf  = Beta sensitivities with respect to distribution parameters')
%    disp('   formresults.pf_sensi_thetaf    = Probability of failure sensitivities with respect to distribution parameters')
   switch lower(gfundata(1).parameter)
   case 'yes' 
%         disp('   formresults.beta_sensi_thetag  = Beta sensitivities with respect to limit-state function parameters')
%         disp('   formresults.pf_sensi_thetag    = Probability of failure sensitivities with respect to limit-state function parameters')
   case 'no'
   otherwise
   end
%    disp(['.................................................'])
%    disp([' '])
   
   else
%    disp(['.................................................'])
%    disp('The following parameters are now available in your current workspace:')
%    disp('   formresults.iter  = Number of iterations')
%    disp('   formresults.beta  = Reliability index beta from FORM analysis')
%    disp('   formresults.pf1   = Failure probability pf1')
%    disp('   formresults.dsptu = Design point u_star')
%    disp('   formresults.alpha = Alpha vector')
%    disp('   formresults.dsptx = Design point in original space')
%    disp('   formresults.imptg = Importance vector gamma')
%    disp('   formresults.gfcn  = Recorded  values of the limit-state function')
%    disp('   formresults.stpsz = Recorded step size values')
%    disp([' '])
%   	disp('  Do you want to see the sensitivities ? ');
%    disp('  1: YES');
%    disp('  2: NO');
%    sensitivities_display = input('  CHOOSE YES or NO ? ');
%    disp(['..............................................................................'])
%    disp(' ');
%    switch sensitivities_display
%    case 1  %----- YES -----------------------
%        if test_ouputfile == 1
%          diary(filename)
%        end
%        disp_sensitivities (1,formresults,probdata,gfundata)
%        diary off
%    case 2  %----- NO -----------------------
%    end
%    disp([' '])
%    disp('   formresults.beta_sensi         = Beta sensitivities with respect to the distribution parameters')
%    disp('   formresults.pf_sensi           = Probability of failure sensitivities with respect to the distribution parameters')
%    switch lower(gfundata(1).parameter)
%    case 'yes' 
%    disp('   formresults.beta_sensi_thetag  = Beta sensitivities with respect to the limit-state function parameters')
%    disp('   formresults.pf_sensi_thetag    = Probability of failure sensitivities with respect to the limit-state function parameters')
%    case 'no'
%    otherwise
%    end
%    disp(['.................................................'])
%    disp([' '])
end

   
case 4 % ---- SORM ----------------------------------------------------------------------------
   
   disp(' ');
   disp(' ');
   disp(' ');
   disp(' ');
   disp('Choose the SORM method you want to use ');
	disp('  0: Exit');
	disp('  1: Curvature-fitted SORM');
	disp('  2: Improved Point-fitted SORM');
    disp('  3: Gradient-based SORM');
    %disp('  4: All the methods');
	disp(' ');
	sorm_analysistype = input('  CHOOSE OPTION FROM THE LIST ABOVE: ');

	switch sorm_analysistype
   
   case 0 % ---- EXIT ----------------------------------------------------------------------------
   disp(' ');
   disp(' ');
   disp(' ');
   disp('  Bye Bye...  : ) ');
   disp(' ');

	case 1 % ---- CURVATURE-FITTED SORM ----------------------------------------------------------- 
   % Clear screen and display message
   disp('SORM analysis is running, please wait... (Ctrl+C breaks)')
   nrv = size(probdata.marg,1);
   % Run the analysis
   t = clock;
   [sorm_curfit_results] = sorm_curvature_fitting(1,formresults,probdata,analysisopt,gfundata,femodel,randomfield);
   results.sorm_curve_fitting = sorm_curfit_results;
   time = etime(clock,t); 
   clc
   if test_ouputfile == 1
      diary(filename)
   end
   disp(' ');
   disp(['               #################################################################################'])
   disp(['               #        RESULTS FROM RUNNING CURVATURE-FITTED SORM RELIABILITY ANALYSIS        #'])
   disp(['               #################################################################################'])
   disp(' ');
   diary off
   disp(['------------------------------------------------------------------------------'])
   disp(['                       Results from FORM analysis                             '])
   disp(['------------------------------------------------------------------------------'])
   disp(['Number of iterations:',num2str(sorm_curfit_results.iter),''])
   disp(['Reliability index beta1: ',num2str(sorm_curfit_results.beta1),'']) 
   disp(['Failure probability pf1: ',num2str(sorm_curfit_results.pf1,'%0.5e'),''])
   disp(' ');
   disp(['------------------------------------------------------------------------------'])
   disp(['                Curvature-fitted SORM and Breitung formula                    '])
   disp(['------------------------------------------------------------------------------'])
   if test_ouputfile == 1
      diary(filename)
   end
   disp(['Time to complete the analysis:    ',num2str(time),''])
   disp(' ');
   disp(['Main curvatures in (n-1)x(n-1) space:    ',num2str(sorm_curfit_results.kappa','%15.5e'),''])
   disp(' ');
   disp(['                                           Breitung formula'])
   disp(['Reliability index beta2 :                  ',num2str(sorm_curfit_results.beta2_breitung),''])
   disp(['Failure probability  pf2:                  ',num2str(sorm_curfit_results.pf2_breitung,'%0.5e'),''])
   disp(' ');
   disp(['                                           Improved Breitung (Hohenbichler / Rackwitz)'])
   disp(['Reliability index beta2 :                  ',num2str(sorm_curfit_results.beta2_breitung_mod),''])
   disp(['Failure probability  pf2:                  ',num2str(sorm_curfit_results.pf2_breitung_mod,'%0.5e'),''])
   disp(' ');
   disp(['                                           Tvedt Exact Integral'])
   disp(['Reliability index beta2 :                  ',num2str(sorm_curfit_results.beta2_tvedt_EI),''])
   disp(['Failure probability  pf2:                  ',num2str(sorm_curfit_results.pf2_tvedt_EI,'%0.5e'),''])
   disp(' ');
   diary off
   disp(['               #################################################################################'])
   disp(['               #     The following parameters are now available in your current workspace      #'])
   disp(['               #################################################################################'])
   disp(' ');
   disp('   sorm_curfit_results.beta2_breitung     = Reliability index beta SORM curvature fitting analysis and Breitung formula')
   disp('   sorm_curfit_results.pf2_breitung       = Failure probability pf2 SORM curvature fitting analysis and Breitung formula')
   disp('   sorm_curfit_results.beta2_breitung_mod = Reliability index beta SORM curvature fitting analysis and Breitung modified formula')
   disp('   sorm_curfit_results.pf2_breitung_mod   = Failure probability pf2 SORM curvature fitting analysis and Breitung formula')
   disp('   sorm_curfit_results.beta2_tvedt_EI     = Reliability index beta2  according to Tvedt exact integral formula')
   disp('   sorm_curfit_results.pf2_tvedt_EI       = Failure probability pf2  according to Tvedt exact integral formula')
   disp('   sorm_curfit_results.hess_G             = Hessian of G evaluated at the design point')
   disp('   sorm_curfit_results.kappa              = Curvatures in the (n-1)(n-1) space')
   disp('   sorm_curfit_results.eig_vectors        = Matrix of eigenvectors')
   disp('   sorm_curfit_results.R1                 = Rotation matrix')
   disp(' ');
   
   
      
   case 2 % ---- IMPROVED POINT-FITTED SORM ----------------------------------------------------------- 
   nrv = size(probdata.marg,1);
   disp(['..............................................................................'])
   disp('  Which method do you want to use to find the fitting points? ');
   disp('  1: Secant iteration');
   disp('  2: Newton iteration');
   improved_point_fitting = input('  CHOOSE OPTION # FROM THE LIST ABOVE ? ');
   disp(['..............................................................................'])
   disp(' ');
  
   switch improved_point_fitting
    case 1  %---- Secant method ----------   
      % Clear screen and display message
   	  disp('SORM analysis is running, please wait... (Ctrl+C breaks)')
      % Run the analysis
      t = clock;
   	  [sorm_point_fitted_secant_results] = sorm_point_fitted_mod_secant(1,formresults,probdata,analysisopt,gfundata,femodel,randomfield);
   	  time = etime(clock,t);
      results.sorm_point_fitted_secant = sorm_point_fitted_secant_results;
      
   % Display results
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   clc
   if test_ouputfile == 1
      diary(filename)
   end
   disp(' ');
   disp(['             ###################################################################################################'])
   disp(['             #      RESULTS FROM RUNNING IMPROVED POINT-FITTED SORM RELIABILITY ANALYSIS USING SECANT METHOD   #' ])
   disp(['             ###################################################################################################'])
   disp(' ');
   diary off
   disp(['------------------------------------------------------------------------------'])
   disp(['                       Results from FORM analysis                             '])
   disp(['------------------------------------------------------------------------------'])
   disp(['Number of iterations:     ',num2str(sorm_point_fitted_secant_results.iter),''])
   disp(['Reliability index beta1 : ',num2str(sorm_point_fitted_secant_results.beta_form),''])
   disp(['Failure probability pf1:  ',num2str(sorm_point_fitted_secant_results.pf1_form,'%0.5e'),''])
   disp(' ');
   disp(['------------------------------------------------------------------------------'])
   disp(['              Improved Point-fitted SORM and Breitung formula                 '])
   disp(['------------------------------------------------------------------------------'])
   if test_ouputfile == 1
      diary(filename)
   end
   disp(['Time to complete the analysis:    ',num2str(time),''])
   disp([' '])
   disp(' ');
   disp(['                                            Breitung formula'])
   disp(['Reliability index beta2 :                   ',num2str(sorm_point_fitted_secant_results.beta2_breitung),''])
   disp(['Failure probability pf2:                    ',num2str(sorm_point_fitted_secant_results.pf2_breitung,'%0.5e'),''])
   disp(' ');
   disp(['                                            Improved Breitung (Hohenbichler / Rackwitz)'])
   disp(['Reliability index beta2 :                   ',num2str(sorm_point_fitted_secant_results.beta2_breitung_mod),''])
   disp(['Failure probability pf2:                    ',num2str(sorm_point_fitted_secant_results.pf2_breitung_mod,'%0.5e'),''])
   disp(' ');
   disp(['                                            Tvedt Exact Integral'])
   disp(['Reliability index beta2 :                   ',num2str(sorm_point_fitted_secant_results.beta2_tvedt_EI),''])
   disp(['Failure probability  pf2:                   ',num2str(sorm_point_fitted_secant_results.pf2_tvedt_EI,'%0.5e'),''])
   disp(' ');
   
   
   if nrv < 20
       
   disp(sprintf('FITTING POINTS ALONG POSITIVE SEMI-AXIS IN THE ROTATED NORMAL SPACE'))
   disp(sprintf('-----------------------------------------------------------------------'))
   disp(sprintf('Axis      abscissa       ordinate    Value of G      Curvatures'))
   for k = 1:nrv-1
    disp(sprintf('%3.d%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n',k,sorm_point_fitted_secant_results.U_prime_plus(k,2),sorm_point_fitted_secant_results.U_prime_plus(k,3),...
                                                                      sorm_point_fitted_secant_results.U_prime_plus(k,4),sorm_point_fitted_secant_results.U_prime_plus(k,5)))
   end
   disp(sprintf('-----------------------------------------------------------------------\n'))
   disp(sprintf('FITTING POINTS ALONG NEGATIVE SEMI-AXIS IN THE ROTATED NORMAL SPACE'))
   disp(sprintf('-----------------------------------------------------------------------'))
   disp(sprintf('Axis      abscissa       ordinate    Value of G      Curvatures'))
   for k = 1:nrv-1
    disp(sprintf('%3.d%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n',k,sorm_point_fitted_secant_results.U_prime_minus(k,2),sorm_point_fitted_secant_results.U_prime_minus(k,3),...
                                                                      sorm_point_fitted_secant_results.U_prime_minus(k,4),sorm_point_fitted_secant_results.U_prime_minus(k,5)))
   end
   disp(sprintf('-----------------------------------------------------------------------\n'))
   end
   diary off
   
   disp(['..............................................................................'])
   disp('  Do you want to see the parameters and the results now available ? ');
   disp('  1: YES');
   disp('  2: NO');
   sorm_all_results = input('  CHOOSE YES or NO ? ');
   disp(['..............................................................................'])
   disp(' ');
   
   switch sorm_all_results
      case 1  %----- YES -----------------------
   disp(' ');
   disp(['               #################################################################################'])
   disp(['               #     The following parameters are now available in your current workspace      #'])
   disp(['               #################################################################################'])
   disp(' ');
   disp('   sorm_point_fitted_secant_results.beta2_breitung     = Reliability index beta Breitung from SORM point-fitted analysis')
   disp('   sorm_point_fitted_secant_results.pf2_breitung       = Failure probability pf2_Breitung from SORM point-fitted analysis')
   disp('   sorm_point_fitted_secant_results.beta2_breitung_mod = Reliability index beta modified Breitung from SORM point-fitted analysis')
   disp('   sorm_point_fitted_secant_results.pf2_breitung_mod   = Failure probability pf2_modified_ Breitung from SORM point-fitted analysis')
   disp('   sorm_point_fitted_secant_results.beta2_tvedt_EI     = Reliability index beta2  according to Tvedt exact integral formula')
   disp('   sorm_point_fitted_secant_results.pf2_tvedt_EI       = Failure probability pf2  according to Tvedt exact integral formula')
   disp('   sorm_point_fitted_secant_results.U_prime_plus       = Coordinates, g values, and curvatures for the u_prime + points')
   disp('   sorm_point_fitted_secant_results.U_prime_minus      = Coordinates, g values, and curvatures for the u_prime - points')
   disp('   sorm_point_fitted_secant_results.kapa_plus          = Curvatures along positive semi-axis')
   disp('   sorm_point_fitted_secant_results.kapa_minus         = Curvatures along negative semi-axis')
   disp(' ');   
		case 2  %----------- NO ----------------------------
   disp(' ');      
   disp(' The analysis is over ')
   disp(' ');
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
          
   	case 2  %---- Newton iteration method ----------
      % Clear screen and display message
   	  disp('SORM analysis is running, please wait... (Ctrl+C breaks)')
      % Run the analysis
      t = clock;
   	  [sorm_point_fitted_newton_results] = sorm_point_fitted_mod(1,formresults,probdata,analysisopt,gfundata,femodel,randomfield);     
      time = etime(clock,t);
      results.sorm_point_fitted_newton = sorm_point_fitted_newton_results;
      
   % Display results
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   clc
   if test_ouputfile == 1
      diary(filename)
   end
   disp(' ');
   disp(['        #######################################################################################################'])
   disp(['        #      RESULTS FROM RUNNING IMPROVED POINT-FITTED SORM RELIABILITY ANALYSIS USING NEWTON ITERATION    #' ])
   disp(['        #######################################################################################################'])
   disp(' ');
   diary off
   disp(['------------------------------------------------------------------------------'])
   disp(['                       Results from FORM analysis                             '])
   disp(['------------------------------------------------------------------------------'])
   disp(['Number of iterations:     ',num2str(sorm_point_fitted_newton_results.iter),''])
   disp(['Reliability index beta1 : ',num2str(sorm_point_fitted_newton_results.beta1),''])
   disp(['Failure probability pf1:  ',num2str(sorm_point_fitted_newton_results.pf1,'%0.5e'),''])
   disp(' ');
   disp(['------------------------------------------------------------------------------'])
   disp(['              Improved Point-fitted SORM and Breitung formula                 '])
   disp(['------------------------------------------------------------------------------'])
   if test_ouputfile == 1
      diary(filename)
   end
   disp(['Time to complete the analysis:    ',num2str(time),''])
   disp([' '])
   disp(' ');
   disp(['                                            Breitung formula'])
   disp(['Reliability index beta2 :                   ',num2str(sorm_point_fitted_newton_results.beta2_breitung),''])
   disp(['Failure probability pf2:                    ',num2str(sorm_point_fitted_newton_results.pf2_breitung,'%0.5e'),''])
   disp(' ');
   disp(['                                            Improved Breitung (Hohenbichler / Rackwitz)'])
   disp(['Reliability index beta2 :                   ',num2str(sorm_point_fitted_newton_results.beta2_breitung_mod),''])
   disp(['Failure probability pf2:                    ',num2str(sorm_point_fitted_newton_results.pf2_breitung_mod,'%0.5e'),''])
   disp(' ');
   disp(['                                            Tvedt Exact Integral'])
   disp(['Reliability index beta2 :                   ',num2str(sorm_point_fitted_newton_results.beta2_tvedt_EI),''])
   disp(['Failure probability  pf2:                   ',num2str(sorm_point_fitted_newton_results.pf2_tvedt_EI,'%0.5e'),''])
   disp(' ');
   
   if nrv < 20
   disp(sprintf('FITTING POINTS ALONG POSITIVE SEMI-AXIS IN THE ROTATED NORMAL SPACE'))
   disp(sprintf('-----------------------------------------------------------------------'))
   disp(sprintf('Axis      abscissa       ordinate    Value of G      Curvatures'))
   for k = 1:nrv-1
    disp(sprintf('%3.d%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n',k,sorm_point_fitted_newton_results.U_prime_plus(k,2),sorm_point_fitted_newton_results.U_prime_plus(k,3),...
                                                                      sorm_point_fitted_newton_results.U_prime_plus(k,4),sorm_point_fitted_newton_results.U_prime_plus(k,5)))
   end
   disp(sprintf('-----------------------------------------------------------------------\n'))
   disp(sprintf('FITTING POINTS ALONG NEGATIVE SEMI-AXIS IN THE ROTATED NORMAL SPACE'))
   disp(sprintf('-----------------------------------------------------------------------'))
   disp(sprintf('Axis      abscissa       ordinate    Value of G      Curvatures'))
   for k = 1:nrv-1
    disp(sprintf('%3.d%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n',k,sorm_point_fitted_newton_results.U_prime_minus(k,2),sorm_point_fitted_newton_results.U_prime_minus(k,3),...
                                                                      sorm_point_fitted_newton_results.U_prime_minus(k,4),sorm_point_fitted_newton_results.U_prime_minus(k,5)))
   end
   disp(sprintf('-----------------------------------------------------------------------\n'))
   end
   diary off
   
   disp(['..............................................................................'])
   disp('  Do you want to see the parameters and the results now available ? ');
   disp('  1: YES');
   disp('  2: NO');
   sorm_all_results = input('  CHOOSE YES or NO ? ');
   disp(['..............................................................................'])
   disp(' ');
   
   switch sorm_all_results
      case 1  %----- YES -----------------------
   disp(' ');
   disp(['               #################################################################################'])
   disp(['               #     The following parameters are now available in your current workspace      #'])
   disp(['               #################################################################################'])
   disp(' ');
   disp('   sorm_point_fitted_newton_results.beta2_breitung     = Reliability index beta Breitung from SORM point-fitted analysis')
   disp('   sorm_point_fitted_newton_results.pf2_breitung       = Failure probability pf2_Breitung from SORM point-fitted analysis')
   disp('   sorm_point_fitted_newton_results.beta2_breitung_mod = Reliability index beta modified Breitung from SORM point-fitted analysis')
   disp('   sorm_point_fitted_newton_results.pf2_breitung_mod   = Failure probability pf2_modified_ Breitung from SORM point-fitted analysis')
   disp('   sorm_point_fitted_newton_results.beta2_tvedt_EI     = Reliability index beta2  according to Tvedt exact integral formula')
   disp('   sorm_point_fitted_newton_results.pf2_tvedt_EI       = Failure probability pf2  according to Tvedt exact integral formula')
   disp('   sorm_point_fitted_newton_results.U_prime_plus       = Coordinates, g values, and curvatures for the u_prime + points')
   disp('   sorm_point_fitted_newton_results.U_prime_minus      = Coordinates, g values, and curvatures for the u_prime - points')
   disp('   sorm_point_fitted_newton_results.kapa_plus          = Curvatures along positive semi-axis')
   disp('   sorm_point_fitted_newton_results.kapa_minus         = Curvatures along negative semi-axis')
   disp(' ');   
		case 2  %----------- NO ----------------------------
   disp(' ');      
   disp(' The analysis is over ')
   disp(' ');
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      
 end
         
      
	case 3 % ---- GRADIENT-BASED SORM ----------------------------------------------------------- 
      disp(' '); 
      %[res] = sormgradient(1,probdata,analysisopt,gfundata,femodel,randomfield);

   disp(' ');
   disp(' ');
   disp('  Will be available soon...  : ) ');
   disp(' ');

   
   

end
   
case 5 % ---- IMPORTANCE SAMPLING SIMULATION ----------------------------------------------
   
   % Clear screen and display message
   disp('Simulation analysis is running, please wait... (Ctrl+C breaks)')
   
   % Set point for importance sampling
   simpoint = analysisopt.sim_point;
   simpoint_previous = analysisopt.sim_point;
   switch lower(simpoint)
   case 'dspt'
      analysisopt.sim_point = formresults.dsptu;
   case 'origin'
      margsize = size(marg);
      analysisopt.sim_point = zeros(margsize(1),1);
   otherwise
      disp('ERROR: Invalid sampling point for simulation analysis.');
   end
   
   % Run simulation analysis
   [simulationresults] = simulation(1,probdata,analysisopt,gfundata,femodel,randomfield);
   results.simulation = simulationresults;
   analysisopt.sim_point = simpoint_previous;
   % Display results
   if test_ouputfile == 1
      diary(filename)
   end
   disp([' '])
   disp(' ');
   disp(['          ########################################################################################################'])
   disp(['          #                                 RESULTS FROM RUNNING SIMULATION ANALYSIS                             #'])
   disp(['          ########################################################################################################'])
   disp(' ');
   disp(['Failure Probability:  ',num2str(simulationresults.pf,'%0.5e'),''])
   disp([' '])
   disp(['Coefficient of variation of probability of failure:  ',num2str(simulationresults.cov_pf,'%0.5e'),''])
   disp([' '])
   disp(['Reliability index beta:  ',num2str(simulationresults.beta),''])
   disp([' '])
   disp(['Number of simulations:  ',num2str(simulationresults.num_sim'),''])
   diary off
   disp(['.................................................'])
   disp('The following parameters are now available in your current workspace:')
   disp('   simulationresults.pf = Failure probability from Monte Carlo simulation')
   disp('   simulationresults.cov_pf = Coefficient of variation for the failure probability from this simulation')
   disp('   simulationresults.beta = Reliability index beta from this simulation')
   disp('   simulationresults.num_sim = Number of simulations')
   disp(['.................................................'])
   disp([' '])
   
case 6 % ---- SYSTEMS ANALYSIS ----------------------------------------------
   
   % Display error message if the FERUMsystems package is not available
   if (exist('systems')~=2)
      disp('The FERUMsystems package must be on the Matlab path');
      disp('before a FERUM systems reliability analysis can be performed');
   end
   
   % Display message about 'being' in the same directory as the executable file
   disp(' ');
   disp('  The "scis_for.exe" file should exist in your current Matlab directory');
   disp('  (use the "cd" command or the "Set Path..." menu option to change current directory)');
   disp('  or you need to add the path to this file to your autoexec.bat file.');
   disp(' ');
   a = input('  Press 0 to exit, press 1 to continue: ');
   disp(' ');
   if (a==0)
      %break;
   end
   clear a;
   
   % Clear screen and display message
   disp('Systems analysis is running, please wait... (Ctrl+C breaks)')
   
   % Run the analysis
   [systemsresults] = systems(probdata,analysisopt,gfundata,femodel,randomfield,system);
   results.system = systemsresults;
   
   % Display results
   if test_ouputfile == 1
      diary(filename)
   end
   disp([' '])
   disp(' ');
   disp(['          ########################################################################################################'])
   disp(['          #                             RESULTS FROM RUNNING SYSTEM RELIABILITY ANALYSIS                         #'])
   disp(['          ########################################################################################################'])
   disp(' ');
   
   if size(systemsresults.pf,2) == 1,
    disp([' Failure Probability of the System is: ',num2str(systemsresults.pf,'%0.5e'),''])
    disp([' Corresponding Reliability Index is: ',num2str(systemsresults.beta),''])
   else
    disp([' Failure Probability of the General System exists between ',num2str(systemsresults.pf(1),'%0.5e'),' and ',num2str(systemsresults.pf(2),'%0.5e'),''])
    disp([' Corresponding Reliability Index exists between           ',num2str(systemsresults.beta(1)),' and ',num2str(systemsresults.beta(2)),''])
   end
   diary off
   disp(['.................................................'])
   disp('The following parameters are now available in your current workspace:')
   disp('(index i and j denotes limit-state function number)')
   disp('   systemsresults.Nscis(i,j)   = Actual Number of iterations of SCIS algorithm for P(C_i*C_j)')
   disp('   systemsresults.covscis(i,j) = Actual Coefficient of Variation of SCIS algorithm for P(C_i*C_j)')
   disp(['.................................................'])
   disp('   systemsresults.iter[i]       = Number of iterations')
   disp('   systemsresults.beta_form[i]  = Reliability index beta from FORM analysis')
   disp('   systemsresults.pf1[i]        = Failure probability pf1')
   disp('   systemsresults.dsptu[:,i]    = Design point u_star')
   disp('   systemsresults.alpha[:,i]    = Alpha vector')
   disp('   systemsresults.dsptx[:,i]    = Design point in original space')
   disp('   systemsresults.imptg[:,i]    = Importance vector gamma')
   disp('   systemsresults.gfcn[i]       = Recorded  values of the limit-state function')
   disp('   systemsresults.stpsz[i]      = Recorded step size values')
   disp(['.................................................'])
   disp([' '])
   
 

  
case 7 % ---- IVERSE FORM ANALYSIS ----------------------------------------------
      
% Clear screen and display message
   disp('FORM analysis is running, please wait... (Ctrl+C breaks)')
   
   % Run the analysis
   t = clock;
   [inverse_formresults] = inverse_form(1,probdata,analysisopt,gfundata,femodel,randomfield);
   results.inverse_form = inverse_formresults;
   time = etime(clock,t); 
   % Display results
   %clc
   if test_ouputfile == 1
      diary(filename)
   end
   disp(' ');
   disp(['          ########################################################################################################'])
   disp(['          #                          RESULTS FROM RUNNING INVERSE FORM RELIABILITY ANALYSIS                      #'])
   disp(['          ########################################################################################################'])
   disp(' ');
   disp(['Number of iterations: ',num2str(inverse_formresults.iter),''])
   disp([' '])
   disp(['Time to complete the analysis:    ',num2str(time),''])
   disp(' ');
   disp(['Reliability index beta1 target: ',num2str(inverse_formresults.beta1),''])
   disp([' '])
   disp(['Value of the deterministic parameter: ',num2str(inverse_formresults.theta,'%0.5e'),''])
   disp([' '])
    
   diary off
      
   disp(['.................................................'])
   disp('The following parameters are now available in your current workspace:')
   disp('   inverse_formresults.iter             = Number of iterations')
   disp('   inverse_formresults.beta1            = Reliability index beta1 target')
   disp('   inverse_formresults.theta            = Value of the deterministic parameter')
   disp('   inverse_formresults.theta_recorded   = Recorded  values of the deterministic parameter')
   disp('   inverse_formresults.u_recorded       = Recorded  values of vector u')
   disp('   inverse_formresults.norm_u_recorded  = Recorded value of the norm of vector u ')
   disp('   inverse_formresults.gfunc_values     = Recorded_G_function_values')
   disp(['.................................................'])
   disp([' '])
   
   
otherwise % -----------------------------------------------------------------------------------
   disp(' ');
   disp('  You entered an invalid choice.');





   disp(' ');
   
end % -----------------------------------------------------------------------------------------

