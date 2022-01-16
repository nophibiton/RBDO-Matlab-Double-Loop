function [ G, grad_g, dgdtheta] = inverse_gfun(lsf,x,theta,grad_flag,probdata,gfundata,femodel,randomfield)


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


switch lower(gfundata(lsf).evaluator)
   
   
case 'basic'
   G = eval(gfundata(lsf).expression);
   if grad_flag == 'no '
      grad_g = 0;
      dgdtheta = 0;
   elseif grad_flag == 'ffd'
      marg = probdata.marg;
      original_x = x;
      % Computes the derivative of the limit-state function with respect to parameter theta
      h_theta = 0.1*theta/2000;
      theta = theta + h_theta;
      G_a_theta_step_ahead = eval(gfundata(lsf).expression);
      dgdtheta = (G_a_theta_step_ahead - G)/h_theta;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for j = 1 : length(x)
         h = marg(j,3)/200;
         x = original_x;
         x(j) = x(j) + h;
         G_a_step_ahead = eval(gfundata(lsf).expression);
         grad_g(j) = (G_a_step_ahead - G)/h;
      end
   elseif grad_flag == 'ddm'
      for i = 1 : length(gfundata(lsf).dgdq)
         grad_g(i) = eval(gfundata(lsf).dgdq{i});
      end
      % Computes the derivative of the limit-state function with respect to parameter theta
      dgdtheta = eval(gfundata(lsf).dgdthetadeter{1});
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   else
      disp('ERROR: Invalid method for gradient computations');
   end
   
   
case 'ferumlinearfecode'
   % Set the right parameter value at each call of FERUMlinerfecode
   gfundata(lsf).lim = theta;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if grad_flag == 'no '
      [ G, dummy ] = FERUMlinearfecode(lsf,x,'no ',gfundata,femodel,randomfield);
      grad_g = 0;
      dgdtheta = 0;
   elseif grad_flag == 'ffd'
      marg = probdata.marg;
      [ G, dummy ] = FERUMlinearfecode(lsf,x,'no ',gfundata,femodel,randomfield);
      % The derivative of the limit-state function with respect to the parameter theta is always equal to 1
      % when using ferumlinearfecode limit-state function type
      dgdtheta = 1;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for j = 1 : length(x)
         h = marg(j,3)/200;
         perturbed_x = x;
         perturbed_x(j) = perturbed_x(j) + h;
         [ G_a_step_ahead, dummy ] = FERUMlinearfecode(lsf,perturbed_x,'no ',gfundata,femodel,randomfield);
         grad_g(j) = (G_a_step_ahead - G)/h;
      end
   elseif grad_flag == 'ddm'
      [ G, grad_g ] = FERUMlinearfecode(lsf,x,'yes',gfundata,femodel,randomfield);
      dgdtheta = 1;
   else
      disp('ERROR: Invalid method for gradient computations');
   end
   
   
case 'ferumnonlinearfecode'
   
   
case 'fedeas'
   marg     = probdata.marg;
   resp     = gfundata(lsf).resp;
   lim      = gfundata(lsf).lim;
   id       = femodel.id;
   Model    = femodel.Model;
   ElData   = femodel.ElData;
   Load     = femodel.Load;
   SolStrat = femodel.SolStrat;
   if grad_flag == 'no '
      [ G, dummy ] = gfunfedeas(x,'no ',id,resp,lim,Model,ElData,Load,SolStrat);
      grad_g = 0;
   elseif grad_flag == 'ffd'
      [ G, grad_g ] = gfunfedeas(x,'no ',id,resp,lim,Model,ElData,Load,SolStrat);
      for j = 1 : length(x)
         h = marg(j,3)/200;
         perturbed_x = x;
         perturbed_x(j) = perturbed_x(j) + h;
         [ G_a_step_ahead, dummy ] = gfunfedeas(perturbed_x,'no ',id,resp,lim,Model,ElData,Load,SolStrat);
         grad_g(j) = (G_a_step_ahead - G)/h;
      end
   elseif grad_flag == 'ddm'   
      [ G, grad_g ] = gfunfedeas(x,'yes',id,resp,lim,Model,ElData,Load,SolStrat);
   else
      disp('ERROR: Invalid method for gradient computations');
   end
   
   
otherwise
   disp('Unknown type of limit-state function evaluator.')
   
   
end


