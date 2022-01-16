function [ u ] = x_to_u(x,probdata,iLo);
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

marg_dim = size(marg);      nrv = marg_dim(1);

for i = 1 : nrv
   switch marg(i,1)
      case 1   % Normal marginal distribution
         z(i) =  ( x(i) - parameter(i,1) ) / parameter(i,2) ;
      case 2   % Lognormal marginal distribution
         xi = parameter(i,4);
         lambda = parameter(i,3);
         z(i) =  ( log(x(i)) - lambda ) / xi ;
      case 3  % Gamma marginal distribution
         lambda = parameter(i,3);
         k = parameter(i,4);
         z(i)=inv_norm_cdf(gammainc(lambda*x(i),k));
      case 4  % Shifted exponential marginal distribution
          lambda = parameter(i,3);
          x_zero = parameter(i,4);
          z(i)=inv_norm_cdf(1-exp(-lambda*(x(i)-x_zero)));
       case 5  % Shifted Rayleigh marginal distribution
          a = parameter(i,3) ;
          x_zero = parameter(i,4);
          z(i)=inv_norm_cdf(1-exp(-0.5*((x(i)-x_zero)/a)^2));
      case 6  % Uniform marginal distribution
         a = parameter(i,3);
         b = parameter(i,4);
         z(i) = inv_norm_cdf( (x(i) - a) / (b-a) );
      case 7  % (Reserved for Beta marginal distribution)
         
         a = parameter(i,5);
         b = parameter(i,6);
         q = parameter(i,3);
         r = parameter(i,4);
      
         % reduce x to the interval [0,1]
         x01(i) = (x(i)-a)/(b-a);
         z(i) = inv_norm_cdf(betacdf(x01(i),q,r));
         
      case 8 % Chi-square marginal distribution
         lambda = 0.5;
         nu = parameter(i,3);
         k = nu/2;
         z(i)=inv_norm_cdf(gammainc(lambda*x(i),k));
         
      case 11 % Type I Largest Value marginal distribution
         a_n = parameter(i,4);
         u_n = parameter(i,3);
         z(i) = inv_norm_cdf(exp(-exp(-a_n*(x(i)-u_n)))); 
      case 12 % Type I Smallest Value marginal distribution
         a_1 = parameter(i,4);
         u_1 = parameter(i,3);
         z(i) = inv_norm_cdf(1-exp(-exp(a_1*(x(i)-u_1))));
      case 13 % Type II Largest Value marginal distribution
         u_n = parameter(i,3);
         k = parameter(i,4);
         z(i) = inv_norm_cdf(exp(-(u_n/x(i))^k));
      case 14 % Type III Smallest Value marginal distribution
         u_1 = parameter(i,3);
         k = parameter(i,4);
         epsilon = parameter(i,5);
         z(i) = inv_norm_cdf(1-exp(-((x(i)-epsilon)/(u_1-epsilon))^k)); 
      
      case 15 % Gumbel marginal distribution
         a_n = parameter(i,4);
         u_n = parameter(i,3);
         z(i) = inv_norm_cdf(exp(-exp(-a_n*(x(i)-u_n))));
      case 16 % Weibull marginal distribution
         u_1 = parameter(i,3);
         k = parameter(i,4);
         z(i) = inv_norm_cdf(1-exp(-(x(i)/u_1)^k));
      case 18 % (Reserved for Laplace marginal distribution)
      case 19 % (Reserved for Pareto marginal distribution)

      otherwise
   end
end 


u = iLo * z';
