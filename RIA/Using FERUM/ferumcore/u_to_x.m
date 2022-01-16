function  [ x ] = u_to_x(u,probdata,Lo)

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

marg = probdata.marg;
parameter = probdata.parameter;

marg_dim = size(marg);     nrv = marg_dim(1);

z = Lo * u;
for i = 1 : nrv
   switch marg(i,1)
      case 1        % Normal distribution
         x(i) = z(i) * parameter(i,2) + parameter(i,1);
      case 2        % Lognormal distribution
         xi = parameter(i,4);
         lambda = parameter(i,3);
         x(i) = exp(  z(i) * xi + lambda  );
      case 3  % Gamma marginal distribution
         lambda = parameter(i,3);
         k = parameter(i,4);
         mean = parameter(i,1);
         normal_val = ferum_cdf(1,z(i),0,1);
         %x(i) = fzero('zero_gamma',mean,[],k,lambda,normal_val);
         %x(i) = fzero('zero_gamma',mean,optimset('fsolve'),k,lambda,normal_val);
         x(i) = fzero('zero_gamma',mean,optimset('fzero'),k,lambda,normal_val);
      case 4  % Shifted Exponential marginal distribution
          lambda = parameter(i,3);
          x_zero = parameter(i,4);
          x(i)= x_zero + 1 / lambda * ( log( 1 / (1-ferum_cdf(1,z(i),0,1)) ) );
      case 5  % Shifted Rayleigh marginal distribution
          a = parameter(i,3);
          x_zero = parameter(i,4);
          x(i)= x_zero + a * ( 2*log( 1 / (1-ferum_cdf(1,z(i),0,1)) ) )^0.5;
      case 6  % Uniform marginal distribution
         a = parameter(i,3);
         b = parameter(i,4);
         x(i) = (b-a) * ferum_cdf(1,z(i),0,1) + a;
      case 7  % (Reserved for Beta marginal distribution)
         
         % extract interval bounds from marg
         a = parameter(i,5);
         b = parameter(i,6);
         q = parameter(i,3);
         r = parameter(i,4);
         
         x01(i) = betainv(ferum_cdf(1,z(i),0,1),q,r);
         % transform x01 from [0,1] to the interval [a,b]
         x(i) = a + x01(i)*(b-a);
         
      case 8 % Chi-square marginal distribution
         lambda = 0.5;
         nu = parameter(i,3);
         k = nu/2 ;
         mean = parameter(i,1);
         normal_val = ferum_cdf(1,z(i),0,1);
         %x(i) = fzero('zero_gamma',mean,[],k,lambda,normal_val);
         %x(i) = fzero('zero_gamma',mean,optimset('fsolve'),k,lambda,normal_val);
         x(i) = fzero('zero_gamma',mean,optimset('fzero'),k,lambda,normal_val);
      case 11 % Type I Largest Value marginal distribution
         a_n = parameter(i,4);
         u_n = parameter(i,3);
         x(i) = u_n - (1/a_n)*log(log(1/ferum_cdf(1,z(i),0,1))); 
      case 12 % Type I Smallest Value marginal distribution
         a_1 = parameter(i,4);
         u_1 = parameter(i,3);
         x(i) = u_1 + (1/a_1)*log(log(1/(1-ferum_cdf(1,z(i),0,1))));
      case 13 % Type II Largest Value marginal distribution
         u_n = parameter(i,3);
         k = parameter(i,4);
         x(i)=u_n * ( ( log(1 / ferum_cdf(1,z(i),0,1) ) ) ^ (-1/k) );
      case 14 % Type III Smallest Value marginal distribution
         u_1 = parameter(i,3);
         k = parameter(i,4);
         epsilon = parameter(i,5);
         x(i) = epsilon + (u_1-epsilon)*(log(1/(1-ferum_cdf(1,z(i),0,1)))^(1/k));
      
      case 15 % Gumbel marginal distribution
         a_n = parameter(i,4);
         u_n = parameter(i,3);
         x(i) = u_n - (1/a_n)*log(log(1/ferum_cdf(1,z(i),0,1))); 
      case 16 % Weibull marginal distribution
         u_1 = parameter(i,3);
         k = parameter(i,4);
         x(i) = u_1*(log(1/(1-ferum_cdf(1,z(i),0,1)))^(1/k));
      case 18 % (Reserved for Laplace marginal distribution)
      case 19 % (Reserved for Pareto marginal distribution)
         
      otherwise
   end
end

x = x';