function f = betacorr(rho0,rho_target,margi,margj)

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

f = rho_target - integral(rho0,margi,margj);
return


function rho = integral(rho0,margi,margj)
% computes the double integral
% integration limits (should be -infty,infty on both axes...)
zmin = -4;
zmax = -zmin;
% determinant of the jacobian of the transformation between [z1max,z1min]x[z2max,z2min] and [-1,1]x[-1,1]
detJ = (zmax-zmin)^2/4;
% get integration points and weight in [-1,1], nIP is the number of integration pts
[xIP,wIP] = Gauss(16);
% transform integration points coordinates from [-1,1] to [zmax,zmin]
z1 = zmin.*ones(size(xIP)) + (zmax-zmin).*(xIP+ones(size(xIP)))./2;
z2 = z1;
switch margi(1)
case 1 % normal
   x1 = z1*margi(3)+margi(2);
case 2 % lognormal
   delta = margi(3)/margi(2);
   zeta = sqrt(log(1+delta^2));
   lambda = log(margi(2)/sqrt(1+delta^2));
   x1 = exp(z1*zeta+lambda);
case 6 % uniform   
case 7 % beta
   % solve non linear set of equations for shape parameters m,n (collected in the vector par)
   par_guess = [1 1];
   par = fsolve('beta_parameters_eqn_set',par_guess,[],'',margi(5),margi(6),margi(2),margi(3));
   % marginally transform the z's in the beta distributed x's
   x1 = betainv(normcdf(z1,0,1),par(1),par(2));
   % inverse beta cdf in matlab gives quantiles in [0,1], transform them to [a,b]
   x1 = x1.*(margi(6)-margi(5))+margi(5);
otherwise
end
switch margj(1)
case 1 % normal
   x2 = z2*margj(3)+margj(2);
case 2 % lognormal
   delta = margj(3)/margj(2);
   zeta = sqrt(log(1+delta^2));
   lambda = log(margj(2)/sqrt(1+delta^2));
   x2 = exp(z2*zeta+lambda);
case 6 % uniform
case 7 % beta
   % solve non linear set of equations for shape parameters m,n (collected in the vector par)
   par_guess = [1 1];
   par = fsolve('beta_parameters_eqn_set',par_guess,[],'',margi(5),margi(6),margi(2),margi(3));
   % marginally transform the z's in the beta distributed x's
   x2 = betainv(normcdf(z2,0,1),par(1),par(2));
   % inverse beta cdf in matlab gives quantiles in [0,1], transform them to [a,b]
   x2 = x2.*(margj(6)-margj(5))+margj(5);
otherwise
end

% compute integral
rho = 0;
for i=1:16
   for j=1:16
      phi2(i,j) = 1/(6.28*sqrt(1-rho0^2))*exp(-1/(2*(1-rho0^2))*(z1(i)^2-2*rho0*z1(i)*z2(j)+z2(j)^2));
      rho = rho + (x1(i)-margi(2))/margi(3)*(x2(j)-margj(2))/margj(3)*phi2(i,j)*detJ*wIP(i)*wIP(j);
   end
end
return