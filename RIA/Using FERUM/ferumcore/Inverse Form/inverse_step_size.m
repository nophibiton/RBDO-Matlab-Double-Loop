function step_size = inverse_step_size(lsf,G,grad_G,u,theta,d_u,d_theta,Lo,probdata,gfundata,femodel,randomfield,analysisopt,J_x_u)

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

grad_flag = lower(analysisopt.grad_flag);
il_max = analysisopt.il_max;

e3 = analysisopt.e3;
beta_target = analysisopt.beta_target;

% Parameter c of the merit function
c = ( beta_target*norm(u) / e3 ) * 2 + 10;

% Define the merit function
merit_u = 0.5 * (norm(u))^2 + c * abs(G);

trial_step_size = 1;

trial_u = u + trial_step_size * d_u;
trial_theta = theta + trial_step_size*d_theta;
trial_x = u_to_x(trial_u,probdata,Lo); 


[ trial_G, dummy1, dummy2] = inverse_gfun(lsf,trial_x,trial_theta,'no ',probdata,gfundata,femodel,randomfield);
%merit_u_new = 0.5 * (norm(trial_u))^2 + (( beta_target*norm(trial_u) / e3 ) * 2 + 10)* abs(trial_G);
merit_u_new = 0.5 * (norm(trial_u))^2 + c* abs(trial_G);

i = 1;

while ( merit_u_new > merit_u ) & ( i < il_max + 1 )
   
   trial_step_size = trial_step_size * 0.5;
   
   trial_u = u + trial_step_size * d_u;
   trial_theta = theta + trial_step_size*d_theta;

   [ trial_x ] = u_to_x(trial_u,probdata,Lo); 

	[ trial_G, dummy1, dummy2] = inverse_gfun(lsf,trial_x,trial_theta,'no ',probdata,gfundata,femodel,randomfield);
	%merit_u_new = 0.5 * (norm(trial_u))^2 + (( beta_target*norm(trial_u) / e3 ) * 2 + 10)* abs(trial_G);
    merit_u_new = 0.5 * (norm(trial_u))^2 + c* abs(trial_G);
   
   i = i + 1;
   
   if ( i == il_max  ) & ( merit_u_new > merit_u )
      factor = 2^(il_max); 
      disp(['The step size has been reduced by a factor of 1/',num2str(factor),' before continuing.'])
   end
   
end

step_size = trial_step_size;
