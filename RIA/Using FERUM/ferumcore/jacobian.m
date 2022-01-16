function  [ J_u_x ] = jacobian(x,u,probdata,Lo,iLo)

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
marg_dim = size(probdata.marg);       nrv = marg_dim(1);

z = Lo * u;

J_z_x = zeros(nrv); 

for i = 1 : nrv
   
   switch marg(i,1)
      case 1   % Normal distribution
         J_z_x(i,i) = 1/parameter(i,2);
      case 2   % Lognormal distribution
         xi = sqrt( log( 1 + ( parameter(i,2) / parameter(i,1) )^2 ) );
         J_z_x(i,i) = 1/( xi * x(i) );
      case 3  % (Gamma marginal distribution)
         lambda = parameter(i,3) ;
         k = parameter(i,4);
         pdf1 = ((lambda*(lambda*x(i))^(k-1))/gamma(k))*exp(-lambda*x(i));
         pdf2 = ferum_pdf(1,z(i),0,1);
         J_z_x(i,i)= pdf1/pdf2;
      case 4  % Shifted exponential marginal distribution
          lambda = parameter(i,3);
          x_zero = parameter(i,4);
          pdf1=lambda * exp(-lambda*(x(i)-x_zero));
          pdf2 = ferum_pdf(1,z(i),0,1);
          J_z_x(i,i) = pdf1/pdf2;
       case 5  % Shifted Rayleigh marginal distribution
          a = parameter(i,3) ;
          x_zero = parameter(i,4);
          pdf1=(x(i)-x_zero)/a^2 * exp(-0.5*((x(i)-x_zero)/a)^2);
          pdf2 = ferum_pdf(1,z(i),0,1);
          J_z_x(i,i) = pdf1/pdf2;
      case 6  % Uniform marginal distribution
         a = parameter(i,3);
         b = parameter(i,4);
         pdf1 = 1 / (b-a);
         pdf2 = ferum_pdf(1,z(i),0,1);
         J_z_x(i,i) = pdf1/pdf2;
      case 7  % ( Reserved Beta marginal distribution )
         % extract interval bounds from marg
         a = parameter(i,5);
         b = parameter(i,6);
         q = parameter(i,3);
         r = parameter(i,4);
         
         % reduce x to the interval [0,1]
         x01(i) = (x(i)-a)/(b-a);
        
         % computes the pdf of x
         pdf1 = betapdf(x01(i),q,r)/(b-a);
         % Here is another expression of the pdf of beta distribution
         % pdf1 = (((x-a)^(q-1))*(b-x)^(r-1))/((gamma(q)*gamma(r)/gamma(q+r))*(b-a)^(q+r-1));
         % computes the pdf of z
         pdf2 = ferum_pdf(1,z(i),0,1);
         J_z_x(i,i) = pdf1/pdf2;
         
      case 11 % Type I Largest Value marginal distribution
         a_n = parameter(i,4);
         u_n = parameter(i,3);
			pdf1 = a_n*exp(-a_n*(x(i)-u_n)-exp(-a_n*(x(i)-u_n)));
         pdf2 = ferum_pdf(1,z(i),0,1);
         J_z_x(i,i) = pdf1/pdf2;
      case 12 % Type I Smallest Value marginal distribution
         a_1 = parameter(i,4);
         u_1 = parameter(i,3);
			pdf1 = a_1*exp(a_1*(x(i)-u_1)-exp(a_1*(x(i)-u_1)));
         pdf2 = ferum_pdf(1,z(i),0,1);
         J_z_x(i,i) = pdf1/pdf2;
      case 13 % Type II Largest Value marginal distribution
         %cov = marg(i,3) / marg(i,2);
         %f = inline('(gamma(1-2/k)-(gamma(1-1/k))^2)^0.5/(gamma(1-1/k))')
         %guess_zero = 3+1/cov;
         %k = fzero('k_type2largest',guess_zero,[],cov);
         %u_n = marg(i,2)/gamma(1-1/k);
         u_n = parameter(i,3);
         k = parameter(i,4);         
         pdf1=(k/u_n)* (u_n/x(i))^(k+1) * exp(-(u_n/x(i))^k);
         pdf2 = ferum_pdf(1,z(i),0,1);
         J_z_x(i,i) = pdf1/pdf2;
      case 14 % Type III Smallest Value marginal distribution
         u_1 = parameter(i,3);
         k = parameter(i,4);
         epsilon = parameter(i,5);
         pdf1 = (k/(u_1-epsilon))* (((x(i)-epsilon)/(u_1-epsilon))^(k-1)) ...
                * exp(-((x(i)-epsilon)/(u_1-epsilon))^k);
         pdf2 = ferum_pdf(1,z(i),0,1);
         J_z_x(i,i) = pdf1/pdf2;
      case 8 % Chi squared distribution
         lambda = 0.5 ;
         nu = parameter(i,3);
         k = nu/2;
         pdf1 = ((lambda*(lambda*x(i))^(k-1))/gamma(k))*exp(-lambda*x(i));
         pdf2 = ferum_pdf(1,z(i),0,1);
         J_z_x(i,i)= pdf1/pdf2;
      case 15 % Gumbel marginal distribution
         a_n = parameter(i,4);
         u_n = parameter(i,3);
			pdf1 = a_n*exp(-a_n*(x(i)-u_n)-exp(-a_n*(x(i)-u_n)));
         pdf2 = ferum_pdf(1,z(i),0,1);
         J_z_x(i,i) = pdf1/pdf2;
      case 16 % Weibull marginal distribution
		 u_1 = parameter(i,3);
         k = parameter(i,4);
         pdf1 = (k/u_1)* ((x(i)/u_1)^(k-1))* exp(-(x(i)/u_1)^k);
         pdf2 = ferum_pdf(1,z(i),0,1);
         J_z_x(i,i) = pdf1/pdf2;
      otherwise
   end
      
end

J_u_x = iLo * J_z_x;