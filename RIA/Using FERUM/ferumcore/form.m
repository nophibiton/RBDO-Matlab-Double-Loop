function formresults = form(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

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

marg      = probdata.marg;
R         = probdata.correlation;

i_max     = analysisopt.ig_max;
e1        = analysisopt.e1;
e2        = analysisopt.e2;
step_code = analysisopt.step_code;
grad_flag = lower(analysisopt.grad_flag);

% Determine number of random variables
marg_dim = size(marg,1);

% Compute corrected correlation coefficients
Ro = mod_corr( probdata, R );

% Cholesky decomposition
Lo = (chol(Ro))';
iLo = inv(Lo);

% Compute starting point for the algorithm
x = ( marg(:,4) );     
u = x_to_u(x,probdata,iLo);

% Set parameters for the iterative loop
i = 1;         % Initialize counter.
conv_flag = 0; % Convergence is achieved when this flag is set to 1.

% Perform iterative loop to find the design point
while conv_flag == 0;
%    disp('..................................')
%    disp('Now carrying out iteration number:'),disp(i)
   
   % Transformation from u to x space
   x = u_to_x(u,probdata,Lo);
   J_u_x = jacobian(x,u,probdata,Lo,iLo);
   J_x_u = inv(J_u_x);
   
   % Evaluate limit-state function and its gradient
   [ G, grad_g ] = gfun(lsf,x,grad_flag,probdata,gfundata,femodel,randomfield);
   grad_G = (grad_g * J_x_u)';
   Recorded_grad_G_values(:,i) = grad_G;
   Recorded_G_function_values(i) = G;
   
   % Set scale parameter Go and inform about struct. resp.
   if i == 1
      Go = G;
%       disp('Value of limit-state function in the first step:')
%       disp(G)
   end
   
   % Compute alpha vector
   alpha = -grad_G / norm(grad_G);
      
   % Check convergence
   if ( (abs(G/Go)<e1) & (norm(u-alpha'*u*alpha)<e2) ) | i==i_max
      conv_flag = 1;
   end
   
   % Take a step if convergence is not achieved
   if conv_flag == 0;
      
      % Determine search direction
      d = search_dir(G,grad_G,u);
      search_direction(:,i)=d;
      
      % Determine step size
      if step_code == 0
         step = step_size(lsf,G,grad_G,u,d,Lo,probdata,gfundata,femodel,randomfield,analysisopt,J_x_u);
         %step = step_size_previous(lsf,G,grad_G,alpha,u,d,Lo,probdata,gfundata,femodel,randomfield,analysisopt,J_x_u);

      else
         step = step_code;
      end
      Recorded_step_size_values(i) = step;
      
      % Determine new trial point
      u_new = u + step * d;
      
      % Prepare for a new round in the loop
      u = u_new;
      i = i + 1;
      
   end
   
end   

if i ~= i_max
   
   %  Post-processing 
   formresults.iter= i;                                                   % Number_of_iterations
   formresults.beta1= alpha' * u;                                          % Reliability_index_beta 
   formresults.pf1= ferum_cdf(1,-formresults.beta1,0,1);                   % Failure_probability_pf1   
   formresults.dsptu= u;                                                  % Design_point_u_star
   formresults.dsptx= x';                                                 % Design_point_in_original_space
   formresults.alpha= alpha;                                              % Alpha_vector
    
   D_prime = diag(diag(   sqrt( J_x_u * J_x_u' )    ));                   % (intermediate computation)
   
   formresults.imptg =(alpha'*J_u_x*D_prime/norm(alpha'*J_u_x*D_prime))'; % Importance_vector_gamma
   formresults.gfcn = Recorded_G_function_values;
   formresults.grad_gfunc = Recorded_grad_G_values;
   formresults.stpsz = Recorded_step_size_values;
   
   
   % ajout provisoire
   %formresults.search_direction = search_direction;
   %%%%%%%%%%%%%%
   
   % Post-processing : Reliability sensitivities
   % with respect to the distribution parameters
   [J_u_theta, beta_sensitivities, pf_sensitivities] = jacobian_u_theta(x',probdata,Lo,iLo,formresults.alpha,formresults.beta1);
      
   formresults.beta_sensi_thetaf = beta_sensitivities;  % Beta sensitivities with respect to the distribution parameters
   formresults.pf_sensi_thetaf = pf_sensitivities;      % Probability of failure sensitivities with respect to the distribution parameters
   % with respect to the limit-state function parameters
   switch lower(gfundata(lsf).parameter)
   case 'yes'
   %n_thetag = length(gfundata(lsf).thetag);
   [ beta_sensi_thetag, pf_sensi_thetag ] = sensi_wrt_thetag(lsf,formresults.dsptx,formresults.beta1,J_x_u,probdata,analysisopt,gfundata,femodel,randomfield);
   formresults.beta_sensi_thetag = beta_sensi_thetag;
   formresults.pf_sensi_thetag = pf_sensi_thetag;
   %disp(sprintf('SENSITIVITIES WITH RESPECT TO THE LIMIT-STATE FUNCTION PARAMETERS'))
   %disp(sprintf('----------------------------------------------------------------------------------------------'))
   %disp(sprintf(' par          d(beta)/d(par)        d(beta)/d(par)  '))
   %for k = 1:n_thetag
   % disp(sprintf('%3.d%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n',k,beta_sensi_thetag(k),pf_sensi_thetag(k)))
   %end
   %disp(sprintf('----------------------------------------------------------------------------------------------\n'))  
   case 'no '
   otherwise
   %disp('Unknown type of limit-state function parameter input.')
   end

   
else
   disp('Maximum number of iterations was reached before convergence.');            
end
