function inverse_formresults = inverse_form(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

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
e3        = analysisopt.e3;

step_code = analysisopt.step_code;
grad_flag = lower(analysisopt.grad_flag);

beta_target = analysisopt.beta_target;

%theta = gfundata(1).deterministic_parameter_start_value;

% Determine number of random variables
marg_dim = size(marg,1);

% Compute corrected correlation coefficients
Ro = mod_corr( probdata, R );

% Cholesky decomposition
Lo = (chol(Ro))';
iLo = inv(Lo);

% Compute starting point for the algorithm
theta = gfundata(1).deterministic_parameter_start_value;
x = ( marg(:,4) );     
u = x_to_u(x,probdata,iLo);

% Set parameters for the iterative loop
i = 1;         % Initialize counter.
conv_flag = 0; % Convergence is achieved when this flag is set to 1.


% Perform iterative loop to find the design point
while conv_flag == 0;
   disp('..................................')
   disp('Now carrying out iteration number:'),disp(i)
   Recorded_theta_values(i) = theta;
   Recorded_u_values(i,:) = u';
   Recorded_norm_u_values(i) = norm(u);
   % Transformation from u to x space
   x = u_to_x(u,probdata,Lo);
   J_u_x = jacobian(x,u,probdata,Lo,iLo);
   J_x_u = inv(J_u_x);
   
   % Evaluate limit-state function and its gradient
   [ G, grad_g, dgdtheta ] = inverse_gfun(lsf,x,theta,grad_flag,probdata,gfundata,femodel,randomfield);
   grad_G = (grad_g * J_x_u)';
   Recorded_G_function_values(i) = G;
   
   % Set scale parameter Go and inform about struct. resp.
   if i == 1
      Go = G;
      disp('Value of limit-state function in the first step:')
      disp(G)
   end
   
   % Compute alpha vector
   alpha = - grad_G / norm(grad_G);
   
   %abs(G/Go)
   %disp(G)
   %beta_target
   %theta
   %norm(u)
  
   % Check convergence
   if i ~= 1
   %((norm(u - u_prev))^2 + (abs(theta - theta_prev))^2)^0.5/((norm(u))^2 + (abs(theta))^2)^0.5
   %abs(theta - theta_prev)
   %aaa = norm(u-u_prev)^2/(norm(u)^2 + (abs(theta))^2)
   %bbb = (abs(theta-theta_prev))^2/(norm(u)^2 + (abs(theta))^2)
   %(aaa + bbb)^0.5
   
   %norm(u-u_prev)^2/(norm(u)^2)
   %(abs(theta - theta_prev))/(abs(theta))
   %(abs(theta - theta_prev))
   %aa = (norm(u) - norm(u_prev))^2 +  (abs(theta-theta_prev))^2
   %bb =(norm(u))^2 + (abs(theta))^2
   %aa/bb
   %(aa/bb)^0.5
   %abs(G/Go)
   %abs(norm(u)- beta_target)
   %abs(theta - theta_prev)/abs(theta_prev)
   
   if ( (abs(G/Go)<e1) & (abs(norm(u)- beta_target) < e3) & ( abs(theta - theta_prev) < e3) ) | i==i_max
      %( (abs(G/Go)<e1) & (norm(u-alpha'*u*alpha)<e2)& ( abs(theta - theta_prev)/(abs(theta)) < e3) ) | i==i_max
%( (abs(G/Go)<e1) &  (norm(u)- beta_target < e3) & ( abs(theta - theta_prev)/(abs(theta)) < e3 ) | i==i_max

      %(((norm(u-u_prev))^2 + (abs(theta-theta_prev))^2)^0.5/((norm(u))^2 + (abs(theta))^2)^0.5 < e3)  | i==i_max

         %( (abs(G/Go)<e1) & (norm(u-alpha'*u*alpha)<e2)& (((norm(u - u_prev))^2- (abs(theta - theta_prev))^2)/((norm(u))^2- (abs(theta))^2) < e3) ) | i==i_max

         %(abs(theta - theta_prev) < e3)  | i==i_max

      %( (abs(G/Go)<e1) &  (norm(u)- beta_target < e3) ) | i==i_max
      %( (abs(G/Go)<e1) & (norm(u-alpha'*u*alpha)<e2)& ((norm(u - u_prev))^2- (abs(theta - theta_prev))^2)/((norm(u))^2- (abs(theta))^2) < e3) ) | i==i_max
      conv_flag = 1;
   end
end

   
   % Take a step if convergence is not achieved
   if conv_flag == 0;
      
      % Determine search direction
      [d_u, d_theta] = inverse_search_dir(G,grad_G,u,beta_target,dgdtheta);
      d = [d_u ; d_theta];
      search_direction(:,i)= d;
      
      % Determine step size
      if step_code == 0
         step = inverse_step_size(lsf,G,grad_G,u,theta,d_u,d_theta,Lo,probdata,gfundata,femodel,randomfield,analysisopt,J_x_u);
         

      else
         step = step_code;
      end
      Recorded_step_size_values(i) = step;
      
      % Determine new trial point
      u_new = u + step * d_u;
      theta_new = theta + step * d_theta;
      
      % Prepare for a new round in the loop
      u_prev = u;
      u = u_new;
      theta_prev = theta;
      theta = theta_new;
      i = i + 1;
      
   end
   
end   

if i ~= i_max
   
   inverse_formresults.iter = i;                                 % Number_of_iterations
   inverse_formresults.beta1 = beta_target;                      % Reliability index beta1 target
   inverse_formresults.theta = theta;                            % Value of the deterministic parameter
   inverse_formresults.theta_recorded = Recorded_theta_values;   % Recorded  values of the deterministic parameter
   inverse_formresults.u_recorded = Recorded_u_values;           % Recorded  values of vector u
   inverse_formresults.norm_u_recorded = Recorded_norm_u_values; % Recorded value of the norm of vector u 
   inverse_formresults.gfunc_values = Recorded_G_function_values;% Recorded G function values
        
else
   disp('Maximum number of iterations was reached before convergence.');            
end
