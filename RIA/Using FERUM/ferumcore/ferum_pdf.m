function p = ferum_pdf(type,x,mean,stdv)

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

switch type
case 1  % (Normal marginal distribution)
   p = 1/( sqrt(2*pi) * stdv ) * exp(-1/2*((x-mean)/stdv)^2);
case 2  % (Lognormal marginal distribution)
case 3  % (Reserved for Gamma marginal distribution)
case 4  % (Reserved for Shifted exponential marginal distribution)
case 5  % (Reserved for Shifted Rayleigh marginal distribution)
case 6  % (Reserved for Uniform marginal distribution)
case 7  % (Reserved for Beta marginal distribution)
case 11 % (Reserved for Type I Largest Value marginal distribution)
case 12 % (Reserved for Type I Smallest Value marginal distribution)
case 13 % (Reserved for Type II Largest Value marginal distribution)
case 14 % (Reserved for Type III Smallest Value marginal distribution)
case 15 % (Reserved for Chi-square marginal distribution)
case 16 % (Reserved for Gumbel marginal distribution)
case 17 % (Reserved for Weibull marginal distribution)
case 18 % (Reserved for Laplace marginal distribution)
case 19 % (Reserved for Pareto marginal distribution)
otherwise
end
