function step_size = step_size(lsf,G,grad_G,u,d,Lo,probdata,gfundata,femodel,randomfield,analysisopt,J_x_u)

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


c = ( norm(u) / norm(grad_G) ) * 2 + 10;

merit = 0.5 * (norm(u))^2 + c * abs(G);

trial_step_size = 1;

trial_u = u + trial_step_size * d;

trial_x = u_to_x(trial_u,probdata,Lo); 

[ trial_G, dummy ] = gfun(lsf,trial_x,'no ',probdata,gfundata,femodel,randomfield);

%[ trial_G, trial_grad_g ] = gfun(lsf,trial_x,grad_flag,probdata,gfundata,femodel,randomfield);
%trial_grad_G = (trial_grad_g * J_x_u)';

merit_new = 0.5 * (norm(trial_u))^2 + c * abs(trial_G);

%merit_new = 0.5 * (norm(trial_u))^2 + (( norm(trial_u) / norm(trial_grad_G) ) * 2 + 10)* abs(trial_G);


i = 1;

while ( merit_new > merit ) & ( i < il_max + 1 )
   
   trial_step_size = trial_step_size * 0.5;
   
   trial_u = u + trial_step_size * d;

   [ trial_x ] = u_to_x(trial_u,probdata,Lo); 

   [ trial_G, dummy ] = gfun(lsf,trial_x,'no ',probdata,gfundata,femodel,randomfield);

	%[ trial_G, trial_grad_g ] = gfun(lsf,trial_x,grad_flag,probdata,gfundata,femodel,randomfield);
	%trial_grad_G = (trial_grad_g * J_x_u)';


	merit_new = 0.5 * (norm(trial_u))^2 + c * abs(trial_G);

	%merit_new = 0.5 * (norm(trial_u))^2 + (( norm(trial_u) / norm(trial_grad_G) ) * 2 + 10)* abs(trial_G);
   
   i = i + 1;
   
   if ( i == il_max ) & ( merit_new > merit )
      factor = 2^(il_max); 
      disp(['The step size has been reduced by a factor of 1/',num2str(factor),' before continuing.'])
   end
   
end

step_size = trial_step_size;
