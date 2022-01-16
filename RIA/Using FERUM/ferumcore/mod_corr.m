function [ Ro ] = mod_corr( probdata, R )

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

for i = 1 : nrv       % Loop through all elements
   for j = 1 : nrv    % of the correlation matrix.
      if i == j       % Diagonal terms in correlation matrix:
         Ro(i,j) = 1.0;
      else            % Off-diagonal terms in correlation matrix:
         
         if ( R(i,i) ~= 1.0 | R(j,j) ~= 1.0  )
            Ro(i,j) = R(i,j);
         else
            
            
            
            if marg(i,1) == 1    % Normal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  Ro(i,j) = R(i,j);
               case  2,    % Lognormal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  xi_j = sqrt( log( 1 + cov_j^2 ) );
                  Ro(i,j) = R(i,j) * cov_j/xi_j;
               case  3,    % Gamma marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor = 1.001-0.007*cov_j+0.018*cov_j^2;
                  Ro(i,j) = R(i,j) * F_factor;
               case  4,    % Shifted Exponential marginal distribution
                  F_factor =  1.107 ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  5,    % Shifted Rayleigh marginal distribution
                  F_factor =  1.014 ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  6,    % Uniform marginal distribution
                  F_factor =  1.023;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  7,    % Beta marginal distribution
                  % solve non linear integral equation for rho0
                  %rho0_guess = .5;
                  %Ro(i,j) = fzero('betacorr',rho0_guess,[],'',R(i,j),marg(i,:),marg(j,:));
                  ba = parameter(j,6) - parameter(j,5);
			      u1 = ( parameter(j,1) - parameter(j,5) ) / ba; 
			      s1 = parameter(j,2) / ba;
			      F_factor = (1.026 + 0.001*R(i,j) - 0.178*u1 + 0.268*s1 - 0.001*R(i,j)*R(i,j)+ 0.178*u1*u1 - 0.679*s1*s1 - 0.003*s1*R(i,j));
                  Ro(i,j) = R(i,j) * F_factor ;
               case 8,    % Chi-square marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor = 1.001-0.007*cov_j+0.018*cov_j^2;
                  Ro(i,j) = R(i,j) * F_factor; 
               case 11,    % Type I Largest Value marginal distribution
                  F_factor =  1.031 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  F_factor =  1.031 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor = 1.030 + 0.238*cov_j + 0.364*cov_j^2;
                  Ro(i,j) = R(i,j) * F_factor;
               case 14,    % Type III Smallest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor = 1.031 - 0.195*cov_j + 0.328*cov_j^2;
                  Ro(i,j) = R(i,j) * F_factor;
               
               case 15,    % Gumbel marginal distribution
                  F_factor =  1.031 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor = 1.031 - 0.195*cov_j + 0.328*cov_j^2;
                  Ro(i,j) = R(i,j) * F_factor;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 2    % Lognormal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  xi_i = sqrt( log( 1 + cov_i^2 ) );
                  Ro(i,j) = R(i,j) * cov_i/xi_i;
               case  2,    % Lognormal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  xi_i = sqrt( log( 1 + cov_i^2 ) );
                  xi_j = sqrt( log( 1 + cov_j^2 ) );
                  Ro(i,j) = 1/(xi_i*xi_j)*log(1+R(i,j)*cov_i*cov_j);
               case  3,    % Gamma marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.001 + 0.033*R(i,j) + 0.004*cov_i - 0.016*cov_j + 0.002*R(i,j)^2 + 0.223*cov_i^2 ...
                              + 0.130*cov_j^2 - 0.104*R(i,j)*cov_i + 0.029*cov_i*cov_j - 0.119*R(i,j)*cov_j; 
                  Ro(i,j) = R(i,j)*F_factor;
               case  4,    % Shifted Exponential marginal distribution
                   cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.098 + 0.003*R(i,j) + 0.019*cov_i + 0.025*R(i,j)^2 + 0.303*cov_i^2 - 0.437*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  5,    % Shifted Rayleigh marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.011 + 0.001*R(i,j) + 0.014*cov_i + 0.004*R(i,j)^2 + 0.231*cov_i^2 - 0.130*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  6,    % Uniform marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor = 1.019 + 0.014*cov_i + 0.010 * R(i,j)^2 + 0.249 * cov_i^2 ;
                  Ro(i,j) = R(i,j) * F_factor;
               case  7,    % Beta marginal distribution
                  % solve non linear integral equation for rho0
                  %rho0_guess = .5;
                  %Ro(i,j) = fzero('betacorr',rho0_guess,[],'',R(i,j),marg(i,:),marg(j,:));
                  ba = parameter(j,6) - parameter(j,5);
			      u1 = ( parameter(j,1) - parameter(j,5) ) / ba; 
			      s1 = parameter(j,2) / ba;
			      uu = u1*u1;
			      ss = s1*s1;
			      xu = R(i,j) * u1;
			      xs = R(i,j) * s1;
			      us = u1*s1;
                  cov_i = parameter(i,2)/parameter(i,1);
			      up = u1*cov_i;
			      sp = s1 * cov_i;
			      F_factor = 0.979 + 0.053*R(i,j) + 0.181*u1 + 0.293*s1 + 0.002*cov_i- 0.004*R(i,j)*R(i,j)...
                      - 0.181*uu + 5.796*ss + 0.277*cov_i*cov_i - 0.107*xu	- 0.619*xs - 0.190*R(i,j)*cov_i...
                      - 3.976*us - 0.097*up + 0.133*sp - 14.402*ss*s1 - 0.069*cov_i*cov_i*cov_i	+ 0.031*R(i,j)*xs...
                      + 0.015*R(i,j)*R(i,j)*cov_i + 3.976*uu*s1 + 0.097*uu*cov_i - 0.430*ss*cov_i + 0.113*sp*cov_i ...
                      + 1.239*R(i,j)*us + 0.380*R(i,j)*up;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 11,    % Type I Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.029 + 0.001*R(i,j) + 0.014*cov_i + 0.004*R(i,j)^2 + 0.233*cov_i^2 - 0.197*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.029 - 0.001*R(i,j) + 0.014*cov_i + 0.004*R(i,j)^2 + 0.233*cov_i^2 + 0.197*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.026 + 0.082*R(i,j) - 0.019*cov_i + 0.222*cov_j + 0.018*R(i,j)^2 + 0.288*cov_i^2 + 0.379*cov_j^2 - 0.441*R(i,j)*cov_i + 0.126*cov_i*cov_j - 0.277*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 14,    % Type III Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.031 + 0.052*R(i,j) + 0.011*cov_i - 0.210*cov_j + 0.002*R(i,j)^2 + 0.220*cov_i^2 + 0.350*cov_j^2 + 0.005*R(i,j)*cov_i + 0.009*cov_i*cov_j - 0.174*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 8,    % Chi-square marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.001 + 0.033*R(i,j) + 0.004*cov_i - 0.016*cov_j + 0.002*R(i,j)^2 + 0.223*cov_i^2 ...
                              + 0.130*cov_j^2 - 0.104*R(i,j)*cov_i + 0.029*cov_i*cov_j - 0.119*R(i,j)*cov_j; 
                  Ro(i,j) = R(i,j)*F_factor;
               case 15,    % Gumbel marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.029 + 0.001*R(i,j) + 0.014*cov_i + 0.004*R(i,j)^2 + 0.233*cov_i^2 - 0.197*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.031 + 0.052*R(i,j) + 0.011*cov_i - 0.210*cov_j + 0.002*R(i,j)^2 + 0.220*cov_i^2 + 0.350*cov_j^2 + 0.005*R(i,j)*cov_i + 0.009*cov_i*cov_j - 0.174*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
               
            elseif marg(i,1) == 3    % Gamma marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor = 1.001-0.007*cov_i+0.018*cov_i^2;
                  Ro(i,j) = R(i,j) * F_factor;
               case  2,    % Lognormal distribution      
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.001 + 0.033*R(i,j) + 0.004*cov_j - 0.016*cov_i + 0.002*R(i,j)^2 + 0.223*cov_j^2 ...
                              + 0.130*cov_i^2 - 0.104*R(i,j)*cov_j + 0.029*cov_j*cov_i - 0.119*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  3,    % Gamma marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.002 + 0.022*R(i,j) - 0.012*(cov_j+cov_i) + 0.001*R(i,j)^2 + 0.125*(cov_j^2+cov_i^2)...
                     - 0.077*R(i,j)*(cov_j+cov_i) + 0.014*cov_i*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  4,    % Shifted Exponential marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.104 + 0.003*R(i,j) - 0.008*cov_i + 0.014*R(i,j)^2 + 0.173*cov_i^2 - 0.296*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  5,    % Shifted Rayleigh marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.014 + 0.001*R(i,j) - 0.007*cov_i + 0.002*R(i,j)^2 + 0.126*cov_i^2 - 0.090*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  6,    % Uniform marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.023 - 0.007*cov_i + 0.002*R(i,j)^2 + 0.127*cov_i^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  7,    % Beta marginal distribution
                  ba = parameter(j,6) - parameter(j,5);
			      u1 = ( parameter(j,1) - parameter(j,5) ) / ba; 
			      s1 = parameter(j,2) / ba;
			
                  %cov_j = parameter(j,2)/parameter(j,1);
                  cov_i = parameter(i,2)/parameter(i,1);
                  
                  uu=u1*u1;
			      ss=s1*s1;
			      xu=R(i,j)*u1;
			      xs=R(i,j)*s1;
			      dus=u1*s1;
			      up=u1*cov_i;
			      sp=s1*cov_i;
			      xp=R(i,j)*cov_i;
            
			if R(i,j) > 0.0
				F_factor = 0.931+0.050*R(i,j)+0.366*u1+0.549*s1+0.181*cov_i...
                    -0.055*R(i,j)*R(i,j)-0.515*uu+4.804*ss-0.484*cov_i*cov_i-0.064*xu...
                    -0.637*xs+0.032*xp-4.059*us-0.319*up-0.211*sp...
                    +0.052*R(i,j)*R(i,j)*R(i,j)+0.227*uu*u2-10.220*ss*s2+0.559*cov_i*cov_i*cov_i-0.042*R(i,j)*xu...
                    +0.223*R(i,j)*xs-0.172*R(i,j)*xp+0.028*R(i,j)*uu+0.695*R(i,j)*ss+0.126*R(i,j)*cov_i*cov_i...
                    +3.845*uu*s2+0.019*uu*cov_i-1.244*us*s1+0.008*up*cov_i-2.075*ss*cov_i...
                    +0.167*sp*cov_i+0.666*R(i,j)*us+0.386*R(i,j)*up-0.517*R(i,j)*sp+2.125*us*cov_i;
			else
				F_factor = 1.025+0.050*R(i,j)-0.029*u1+0.047*s1-0.136*cov_i...
                    +0.069*R(i,j)*R(i,j)+0.178*uu+6.281*ss+0.548*cov_i*cov_i-0.027*xu...
                    -0.686*xs+0.046*xp-3.513*us+0.231*up+0.299*sp...
                    +0.063*R(i,j)*R(i,j)*R(i,j)-0.226*uu*u1-17.507*ss*s2-0.366*cov_i*cov_i*cov_i+0.051*R(i,j)*xu...
                    -0.246*R(i,j)*xs+0.186*R(i,j)*xp-0.001*R(i,j)*uu+0.984*R(i,j)*ss+0.121*R(i,j)*cov_i*cov_i...
                    +3.700*uu*s1+0.081*uu*cov_i+1.356*us*s1+0.002*up*cov_i+1.654*ss*cov_i...
                    -0.135*sp*cov_i+0.619*R(i,j)*us+0.410*R(i,j)*up-0.686*R(i,j)*sp-2.205*us*cov_i;
            end
        
			Ro(i,j) = F_factor*R(i,j);
               case 11,    % Type I Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.031 + 0.001*R(i,j)- 0.007*cov_i  + 0.003*R(i,j)^2 + 0.131*cov_i^2 - 0.132*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.031 - 0.001*R(i,j)- 0.007*cov_i  + 0.003*R(i,j)^2 + 0.131*cov_i^2 + 0.132*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.029 + 0.056*R(i,j)- 0.030*cov_i + 0.225*cov_j + 0.012*R(i,j)^2 + 0.174*cov_i^2 ...
                              + 0.379*cov_j^2 - 0.313*R(i,j)*cov_i + 0.075*cov_i*cov_j - 0.182*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case 14,    % Type III Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.032 + 0.034*R(i,j)- 0.007*cov_i - 0.202*cov_j + 0.121*cov_i^2 ...
                              + 0.339*cov_j^2 - 0.006*R(i,j)*cov_i + 0.003*cov_i*cov_j - 0.111*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case 8,    % Chi-square marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.002 + 0.022*R(i,j) - 0.012*(cov_j+cov_i) + 0.001*R(i,j)^2 + 0.125*(cov_j^2+cov_i^2)...
                     - 0.077*R(i,j)*(cov_j+cov_i) + 0.014*cov_i*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case 15,    % Gumbel marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.031 + 0.001*R(i,j)- 0.007*cov_i  + 0.003*R(i,j)^2 + 0.131*cov_i^2 - 0.132*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.032 + 0.034*R(i,j)- 0.007*cov_i - 0.202*cov_j + 0.121*cov_i^2 ...
                              + 0.339*cov_j^2 - 0.006*R(i,j)*cov_i + 0.003*cov_i*cov_j - 0.111*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 4    % Shifted Exponential marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  F_factor =  1.107 ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  2,    % Lognormal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.098 + 0.003*R(i,j) + 0.019*cov_j + 0.025*R(i,j)^2 + 0.303*cov_j^2 - 0.437*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  3,    % Gamma marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.104 + 0.003*R(i,j) - 0.008*cov_j + 0.014*R(i,j)^2 + 0.173*cov_j^2 - 0.296*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  4,    % Shifted Exponential marginal distribution
                  F_factor =  1.229 - 0.367*R(i,j) + 0.153*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  5,    % Shifted Rayleigh marginal distribution
                  F_factor =  1.123 - 0.100*R(i,j) + 0.021*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  6,    % Uniform marginal distribution
                  F_factor =  1.133 + 0.029*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  7,    % Beta marginal distribution
                  ba = parameter(j,6) - parameter(j,5);
			      u1 = ( parameter(j,1) - parameter(j,5) )/ ba; 
			      s1 = parameter(j,2) / ba;
                   
			     uu = u1*u1;
			     ss = s1*s1;
			     xu = R(i,j)*u1;
			     xs = R(i,j)*s1;
			     us = u1*s1;
			     F_factor = 1.082-0.004*R(i,j)+0.204*u1+0.432*s1-0.001*R(i,j)*R(i,j)...
                            - 0.204*uu+7.728*ss+0.008*xu-1.699*xs-5.338*us...
                            - 19.741*ss*s2+0.135*R(i,j)*xs+5.338*uu*s1+3.397*R(i,j)*us;
			      Ro(i,j) = F_factor*R(i,j);
               case 11,    % Type I Largest Value marginal distribution
                  F_factor =  1.142 - 0.154*R(i,j) + 0.031*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  F_factor =  1.142 + 0.154*R(i,j) + 0.031*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case 13,    % Type II Largest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.109 - 0.152*R(i,j) + 0.361*cov_j + 0.130*R(i,j)^2 + 0.455*cov_j^2 - 0.728*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 14,    % Type III Smallest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.147 + 0.145*R(i,j) - 0.271*cov_j + 0.010*R(i,j)^2 + 0.459*cov_j^2 - 0.467*R(i,j)*cov_j ;                  
                  Ro(i,j) = R(i,j) * F_factor ;
               case 8,    % Chi-square marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.104 + 0.003*R(i,j) - 0.008*cov_j + 0.014*R(i,j)^2 + 0.173*cov_j^2 - 0.296*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 15,    % Gumbel marginal distribution
                  F_factor =  1.142 - 0.154*R(i,j) + 0.031*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.147 + 0.145*R(i,j) - 0.271*cov_j + 0.010*R(i,j)^2 + 0.459*cov_j^2 - 0.467*R(i,j)*cov_j ;                  
                  Ro(i,j) = R(i,j) * F_factor ;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 5    % Shifted Rayleigh marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  F_factor = 1.014;
                  Ro(i,j) =  R(i,j)*F_factor;
               case  2,    % Lognormal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.011 + 0.001*R(i,j) + 0.014*cov_j + 0.004*R(i,j)^2 + 0.231*cov_j^2 - 0.130*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  3,    % Gamma marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.014 + 0.001*R(i,j) - 0.007*cov_j + 0.002*R(i,j)^2 + 0.126*cov_j^2 - 0.090*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  4,    % Shifted Exponential marginal distribution
                 F_factor =  1.123 - 0.100*R(i,j) + 0.021*R(i,j)^2 ;
                 Ro(i,j) = R(i,j) * F_factor ;
              case  5,    % Shifted Rayleigh marginal distribution
                 F_factor =  1.028 - 0.029*R(i,j);
                  Ro(i,j) = R(i,j) * F_factor ;
               case  6,    % Uniform marginal distribution
                  F_factor =  1.038 - 0.008*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  7,    % Beta marginal distribution
                  ba = parameter(j,6) - parameter(j,5);
			      u1 = ( parameter(j,1) - parameter(j,5) )/ ba; 
			      s1 = parameter(j,2) / ba;
                  
                  F_factor  = 1.037-0.042*R(i,j)-0.182*u1+0.369*s1-0.001*R(i,j)*R(i,j)+...
                              0.182*u1*u1-1.150*s1*s1+0.084*R(i,j)*u1;
			       Ro(i,j) = F_factor *R(i,j);
               case 11,    % Type I Largest Value marginal distribution
                  F_factor =  1.046 - 0.045*R(i,j) + 0.006*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  F_factor =  1.046 + 0.045*R(i,j) + 0.006*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.036 - 0.038*R(i,j) + 0.266*cov_j + 0.028*R(i,j)^2 + 0.383*cov_j^2 - 0.229*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 14,    % Type III Smallest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.047 + 0.042*R(i,j) - 0.212*cov_j + 0.353*cov_j^2 - 0.136*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 8,    % Chi-square marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.014 + 0.001*R(i,j) - 0.007*cov_j + 0.002*R(i,j)^2 + 0.126*cov_j^2 - 0.090*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 15,    % Gumbel marginal distribution
                  F_factor =  1.046 - 0.045*R(i,j) + 0.006*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.047 + 0.042*R(i,j) - 0.212*cov_j + 0.353*cov_j^2 - 0.136*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 6    % Uniform marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  Ro(i,j) = 1.023 * R(i,j);
               case  2,    % Lognormal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor = 1.019 + 0.014*cov_j + 0.010 * R(i,j)^2 + 0.249 * cov_j^2;
                  Ro(i,j) = R(i,j) * F_factor;
               case  3,    % Gamma marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor = 1.023 - 0.007*cov_j + 0.002 * R(i,j)^2 + 0.127 * cov_j^2;
                  Ro(i,j) = R(i,j) * F_factor;
               case  4,    % Shifted Exponential marginal distribution
                  F_factor = 1.133 + 0.029 * R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor;
               case  5,    % Shifted Rayleigh marginal distribution
                  F_factor =  1.038 - 0.008*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  6,    % Uniform marginal distribution
                  Ro(i,j) = R(i,j) * ( 1.047 - 0.047 * R(i,j)^2 );
               case  7,    % Beta marginal distribution
                   ba = parameter(j,6) - parameter(j,5);
			       u1 = ( parameter(j,1) - parameter(j,5) )/ ba; 
			       s1 = parameter(j,2) / ba;
                                     
			       F_factor = 1.040+0.015*R(i,j)-0.176*u1+0.432*s1-0.008*R(i,j)*R(i,j)+...
                              0.176*u1*u1-1.286*s1*s1-0.137*R(i,j)*s1;
			       Ro(i,j) = F_factor * R(i,j);
               case 11,    % Type I Largest Value marginal distribution
                  F_factor =  1.055 + 0.015*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  F_factor =  1.055 + 0.015*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.033 + 0.305*cov_j + 0.074*R(i,j)^2 + 0.405*cov_j^2 ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case 14,    % Type III Smallest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.061 - 0.237*cov_j - 0.005*R(i,j)^2 + 0.379*cov_j^2 ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case 8,    % Chi-square marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor = 1.023 - 0.007*cov_j + 0.002 * R(i,j)^2 + 0.127 * cov_j^2;
                  Ro(i,j) = R(i,j) * F_factor;
               case 15,    % Gumbel marginal distribution
                  F_factor =  1.055 + 0.015*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.061 - 0.237*cov_j - 0.005*R(i,j)^2 + 0.379*cov_j^2 ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 7    % Beta marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  ba = parameter(i,6) - parameter(i,5);
			      u1 = ( parameter(i,1) - parameter(i,5) ) / ba; 
			      s1 = parameter(i,2) / ba;
			      F_factor = (1.026 + 0.001*R(i,j) - 0.178*u1 + 0.268*s1 - 0.001*R(i,j)*R(i,j)+ 0.178*u1*u1 - 0.679*s1*s1 - 0.003*s1*R(i,j));
                  Ro(i,j) = R(i,j) * F_factor ;
                  
               case  2,    % Lognormal distribution
                  ba = parameter(i,6) - parameter(i,5);
			      u1 = ( parameter(i,1) - parameter(i,5) ) / ba; 
			      s1 = parameter(i,2) / ba;
			      uu = u1*u1;
			      ss = s1*s1;
			      xu = R(i,j)* u1;
			      xs = R(i,j)* s1;
			      us = u1*s1;
                  cov_j = parameter(j,2)/parameter(j,1);
			      up = u1*cov_j;
			      sp = s1 * cov_j;
			      F_factor = 0.979 + 0.053*R(i,j) + 0.181*u1 + 0.293*s1 + 0.002*cov_j- 0.004*R(i,j)*R(i,j)...
                      - 0.181*uu + 5.796*ss + 0.277*cov_j*cov_j - 0.107*xu	- 0.619*xs - 0.190*R(i,j)*cov_j...
                      - 3.976*us - 0.097*up + 0.133*sp - 14.402*ss*s1 - 0.069*cov_j*cov_j*cov_j	+ 0.031*R(i,j)*xs...
                      + 0.015*R(i,j)*R(i,j)*cov_j + 3.976*uu*s1 + 0.097*uu*cov_j - 0.430*ss*cov_j + 0.113*sp*cov_j ...
                      + 1.239*R(i,j)*us + 0.380*R(i,j)*up;
                  Ro(i,j) = R(i,j) * F_factor ;
                   
               case  3,    % Gamma marginal distribution
                   %double ba2 = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			%double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba2; 
			%double s2 = rv2Ptr->getStdv() / ba2;
			      ba = parameter(i,6) - parameter(i,5);
			      u1 = ( parameter(i,1) - parameter(i,5) ) / ba; 
			      s1 = parameter(i,2) / ba;
			
                  cov_j = parameter(j,2)/parameter(j,1);
                  %cov_i = parameter(i,2)/parameter(i,1);
                  
                  uu=u1*u1;
			      ss=s1*s1;
			      xu=R(i,j)*u1;
			      xs=R(i,j)*s1;
			      dus=u1*s1;
			      up=u1*cov_j;
			      sp=s1*cov_j;
			      xp=R(i,j)*cov_j;
            
			if R(i,j) > 0.0
				F_factor = 0.931+0.050*R(i,j)+0.366*u1+0.549*s1+0.181*cov_j...
                    -0.055*R(i,j)*R(i,j)-0.515*uu+4.804*ss-0.484*cov_j*cov_j-0.064*xu...
                    -0.637*xs+0.032*xp-4.059*us-0.319*up-0.211*sp...
                    +0.052*R(i,j)*R(i,j)*R(i,j)+0.227*uu*u2-10.220*ss*s2+0.559*cov_j*cov_j*cov_j-0.042*R(i,j)*xu...
                    +0.223*R(i,j)*xs-0.172*R(i,j)*xp+0.028*R(i,j)*uu+0.695*R(i,j)*ss+0.126*R(i,j)*cov_j*cov_j...
                    +3.845*uu*s2+0.019*uu*cov_j-1.244*us*s1+0.008*up*cov_j-2.075*ss*cov_j...
                    +0.167*sp*cov_j+0.666*R(i,j)*us+0.386*R(i,j)*up-0.517*R(i,j)*sp+2.125*us*cov_j;
			else
				F_factor = 1.025+0.050*R(i,j)-0.029*u1+0.047*s1-0.136*cov_j...
                    +0.069*R(i,j)*R(i,j)+0.178*uu+6.281*ss+0.548*cov_j*cov_j-0.027*xu...
                    -0.686*xs+0.046*xp-3.513*us+0.231*up+0.299*sp...
                    +0.063*R(i,j)*R(i,j)*R(i,j)-0.226*uu*u1-17.507*ss*s2-0.366*cov_j*cov_j*cov_j+0.051*R(i,j)*xu...
                    -0.246*R(i,j)*xs+0.186*R(i,j)*xp-0.001*R(i,j)*uu+0.984*R(i,j)*ss+0.121*R(i,j)*cov_j*cov_j...
                    +3.700*uu*s1+0.081*uu*cov_j+1.356*us*s1+0.002*up*cov_j+1.654*ss*cov_j...
                    -0.135*sp*cov_j+0.619*R(i,j)*us+0.410*R(i,j)*up-0.686*R(i,j)*sp-2.205*us*cov_j;
            end
        
			Ro(i,j) = F_factor*R(i,j);
                   
               case  4,    % Shifted Exponential marginal distribution
                  ba = parameter(i,6) - parameter(i,5);
			      u1 = ( parameter(i,1) - parameter(i,5) )/ ba; 
			      s1 = parameter(i,2) / ba;
                   
			     uu = u1*u1;
			     ss = s1*s1;
			     xu = R(i,j)*u1;
			     xs = R(i,j)*s1;
			     us = u1*s1;
			     F_factor = 1.082-0.004*R(i,j)+0.204*u1+0.432*s1-0.001*R(i,j)*R(i,j)...
                            - 0.204*uu+7.728*ss+0.008*xu-1.699*xs-5.338*us...
                            - 19.741*ss*s2+0.135*R(i,j)*xs+5.338*uu*s1+3.397*R(i,j)*us;
			      Ro(i,j) = F_factor*R(i,j);
               case  5,    % Shifted Rayleigh marginal distribution
                  ba = parameter(i,6) - parameter(i,5);
			      u1 = ( parameter(i,1) - parameter(i,5) )/ ba; 
			      s1 = parameter(i,2) / ba;
                  
                  F_factor  = 1.037-0.042*R(i,j)-0.182*u1+0.369*s1-0.001*R(i,j)*R(i,j)+...
                              0.182*u1*u1-1.150*s1*s1+0.084*R(i,j)*u1;
			       Ro(i,j) = F_factor *R(i,j);
               case  6,    % Uniform marginal distribution
                   ba = parameter(i,6) - parameter(i,5);
			       u1 = ( parameter(i,1) - parameter(i,5) )/ ba; 
			       s1 = parameter(i,2) / ba;
                                     
			       F_factor = 1.040+0.015*R(i,j)-0.176*u1+0.432*s1-0.008*R(i,j)*R(i,j)+...
                              0.176*u1*u1-1.286*s1*s1-0.137*R(i,j)*s1;
			       Ro(i,j) = F_factor * R(i,j);
                   
               case  7,    % Beta marginal distribution
                  ba1 = parameter(i,6) - parameter(i,5);
			      u1 = ( parameter(i,1) - parameter(i,5) )/ ba1; 
			      s1 = parameter(i,2) / ba1;
                  
                  ba2 = parameter(j,6) - parameter(j,5);
			      u2 = ( parameter(j,1) - parameter(j,5) )/ ba2; 
			      s2 = parameter(j,2) / ba2;
                  
			          
            cov_j = parameter(j,2)/parameter(j,1);
            cov_i = parameter(i,2)/parameter(i,1);
            
			o=u1;
			p=s1;
			q=u2;
			r=s2;
			u=o+q;
			s=p+r;
			oo=o*o;
			pp=p*cov_i;
			qq=q*q;
			rr=r*r;
			oq=o*q;
			pr=p*r;
            
			us=oo+cov_j*cov_j;
			ss=cov_i*cov_i+rr;
			uc=oo*o+cov_j*cov_j*q;
			sc=cov_i*cov_i*cov_i+rr*r;
			x=R(i,j);
			xx=R(i,j)*R(i,j);
			
			if R(i,j) > 0.0
                F_factor =1.030-0.050*x-0.056*u+0.094*s+0.009*xx-0.084*us+2.583*ss...
                    +0.100*x*u+0.010*x*s-0.485*u*s+0.270*oq-0.688*pr-0.011*xx*x...
                    +0.024*uc-10.786*sc+0.013*xx*u-0.035*xx*s+0.001*x*us-0.069*x*ss...
                    +1.174*us*s+0.004*oq*u+0.227*ss*u+2.783*pr*s+0.058*x*s*u...
                    -0.260*x*oq-0.352*x*pr-1.609*oq*s+0.194*pr*u;
            else
                F_factor=0.939-0.023*x+0.147*u+0.668*s+0.035*xx-0.008*us+3.146*ss...
                    +0.103*x*u-0.126*x*s-1.866*u*s-0.268*oq-0.304*pr+0.011*xx*x...
                    -0.024*uc-10.836*sc-0.013*xx*u-0.035*xx*s-0.001*x*us+0.069*x*ss...
                    +1.175*us*s-0.005*oq*u-0.270*ss*u+2.781*pr*s+0.058*x*u*s...
                    -0.259*x*oq+0.352*x*pr+1.608*oq*s-0.189*pr*u;
            end
            Ro(i,j) = F_factor* R(i,j);
                   
               case 11,    % Type I Largest Value marginal distribution
                   ba = parameter(i,6) - parameter(i,5);
			       u1 = ( parameter(i,1) - parameter(i,5) )/ ba; 
			       s1 = parameter(i,2) / ba;
                                      
			       F_factor = 1.055-0.066*R(i,j)-0.194*u1+0.391*s1+0.003*R(i,j)*R(i,j)...
                             +0.194*u1*u1-1.134*s1*s1+0.130*R(i,j)*u1+0.003*R(i,j)*s1;
			       Ro(i,j) = F_factor * R(i,j);
                   
               case 12,    % Type I Smallest Value marginal distribution
                   ba = parameter(i,6) - parameter(i,5);
			       u1 = ( parameter(i,1) - parameter(i,5) )/ ba; 
			       s1 = parameter(i,2) / ba;
			       F_factor = 1.055+0.066*R(i,j)-0.194*u1+0.391*s1+0.003*R(i,j)*R(i,j)...
                              +0.194*u1*u1-1.134*s1*s1-0.130*R(i,j)*u1-0.003*R(i,j)*s1;
			       Ro(i,j) = F_factor * R(i,j);
                   
               case 13,    % Type II Largest Value marginal distribution
                   ba = parameter(i,6) - parameter(i,5);
			       u1 = ( parameter(i,1) - parameter(i,5) )/ ba; 
			       s1 = parameter(i,2) / ba;
			       uu = u1*u1;
			       ss = s1*s1;
			       xu = R(i,j)*u1;
			       xs = R(i,j)*s1;
			       us = u1*s1;
                   cov_j = parameter(j,2)/parameter(j,1);
			       uq = u1*cov_j;
			       sq = s1*cov_j;
			       xq=R(i,j)*cov_j;
			       F_factor = 1.005 + 0.091*R(i,j) + 0.285*u1+ 0.260*s1+ 0.199*cov_j...
                             - 0.023*R(i,j)*R(i,j) - 0.285*uu + 8.180*ss + 0.543*cov_j*cov_j - 0.181*xu...
                             - 1.744*xs - 0.336*xq - 5.450*us - 0.265*uq + 0.514*sq...
                             -19.661*ss*s1- 0.178*cov_j*cov_j*cov_j...
                             + 0.244*R(i,j)*xs + 0.066*R(i,j)*R(i,j)*cov_j - 0.001*R(i,j)*ss...
                             + 5.450*uu*s1+ 0.265*uu*cov_j - 0.986*ss*cov_j...
                             + 0.133*sq*cov_j + 3.488*R(i,j)*us + 0.671*R(i,j)*uq;
			       Ro(i,j) = F_factor * R(i,j);
               case 14,    % Type III Smallest Value marginal distribution
                   ba = parameter(i,6) - parameter(i,5);
			       u1 = ( parameter(i,1) - parameter(i,5) )/ ba; 
			       s1 = parameter(i,2) / ba;
                   cov_j = parameter(j,2)/parameter(j,1);
                   
			       F_factor = 1.054+0.002*R(i,j)-0.176*u1+0.366*s1-0.201*cov_j...
                       -0.002*R(i,j)*R(i,j)+0.176*u1*u1-1.098*s1*s1+0.340*cov_j*cov_j...
                       -0.004*R(i,j)*u1-0.029*s1*cov_j;
			       Ro(i,j) = F_factor* R(i,j);
               case 8,    % Chi-square marginal distribution
                        %double ba2 = rv2Ptr->getParameter4() - rv2Ptr->getParameter3();
			%double u2 = ( rv2Ptr->getMean() - rv2Ptr->getParameter3() ) / ba2; 
			%double s2 = rv2Ptr->getStdv() / ba2;
			      ba = parameter(i,6) - parameter(i,5);
			      u1 = ( parameter(i,1) - parameter(i,5) ) / ba; 
			      s1 = parameter(i,2) / ba;
			
                  cov_j = parameter(j,2)/parameter(j,1);
                  %cov_i = parameter(i,2)/parameter(i,1);
                  
                  uu=u1*u1;
			      ss=s1*s1;
			      xu=R(i,j)*u1;
			      xs=R(i,j)*s1;
			      dus=u1*s1;
			      up=u1*cov_j;
			      sp=s1*cov_j;
			      xp=R(i,j)*cov_j;
            
			if R(i,j) > 0.0
				F_factor = 0.931+0.050*R(i,j)+0.366*u1+0.549*s1+0.181*cov_j...
                    -0.055*R(i,j)*R(i,j)-0.515*uu+4.804*ss-0.484*cov_j*cov_j-0.064*xu...
                    -0.637*xs+0.032*xp-4.059*us-0.319*up-0.211*sp...
                    +0.052*R(i,j)*R(i,j)*R(i,j)+0.227*uu*u2-10.220*ss*s2+0.559*cov_j*cov_j*cov_j-0.042*R(i,j)*xu...
                    +0.223*R(i,j)*xs-0.172*R(i,j)*xp+0.028*R(i,j)*uu+0.695*R(i,j)*ss+0.126*R(i,j)*cov_j*cov_j...
                    +3.845*uu*s2+0.019*uu*cov_j-1.244*us*s1+0.008*up*cov_j-2.075*ss*cov_j...
                    +0.167*sp*cov_j+0.666*R(i,j)*us+0.386*R(i,j)*up-0.517*R(i,j)*sp+2.125*us*cov_j;
			else
				F_factor = 1.025+0.050*R(i,j)-0.029*u1+0.047*s1-0.136*cov_j...
                    +0.069*R(i,j)*R(i,j)+0.178*uu+6.281*ss+0.548*cov_j*cov_j-0.027*xu...
                    -0.686*xs+0.046*xp-3.513*us+0.231*up+0.299*sp...
                    +0.063*R(i,j)*R(i,j)*R(i,j)-0.226*uu*u1-17.507*ss*s2-0.366*cov_j*cov_j*cov_j+0.051*R(i,j)*xu...
                    -0.246*R(i,j)*xs+0.186*R(i,j)*xp-0.001*R(i,j)*uu+0.984*R(i,j)*ss+0.121*R(i,j)*cov_j*cov_j...
                    +3.700*uu*s1+0.081*uu*cov_j+1.356*us*s1+0.002*up*cov_j+1.654*ss*cov_j...
                    -0.135*sp*cov_j+0.619*R(i,j)*us+0.410*R(i,j)*up-0.686*R(i,j)*sp-2.205*us*cov_j;
            end
        
			Ro(i,j) = F_factor*R(i,j);
               case 15,    % Gumbel marginal distribution
                   ba = parameter(i,6) - parameter(i,5);
			       u1 = ( parameter(i,1) - parameter(i,5) )/ ba; 
			       s1 = parameter(i,2) / ba;
                                      
			       F_factor = 1.055-0.066*R(i,j)-0.194*u1+0.391*s1+0.003*R(i,j)*R(i,j)...
                             +0.194*u1*u1-1.134*s1*s1+0.130*R(i,j)*u1+0.003*R(i,j)*s1;
			       Ro(i,j) = F_factor * R(i,j);
               case 16,    % Weibull marginal distribution
                   ba = parameter(i,6) - parameter(i,5);
			       u1 = ( parameter(i,1) - parameter(i,5) )/ ba; 
			       s1 = parameter(i,2) / ba;
                   cov_j = parameter(j,2)/parameter(j,1);
                   
			       F_factor = 1.054+0.002*R(i,j)-0.176*u1+0.366*s1-0.201*cov_j...
                       -0.002*R(i,j)*R(i,j)+0.176*u1*u1-1.098*s1*s1+0.340*cov_j*cov_j...
                       -0.004*R(i,j)*u1-0.029*s1*cov_j;
			       Ro(i,j) = F_factor* R(i,j);
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 11    % Type I Largest Value marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  F_factor =  1.031 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  2,    % Lognormal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.029 + 0.001*R(i,j) + 0.014*cov_j + 0.004*R(i,j)^2 + 0.233*cov_j^2 - 0.197*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  3,    % Gamma marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.031 + 0.001*R(i,j)- 0.007*cov_j  + 0.003*R(i,j)^2 + 0.131*cov_j^2 - 0.132*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  4,    % Shifted Exponential marginal distribution
                  F_factor =  1.142 - 0.154*R(i,j) + 0.031*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  5,    % Shifted Rayleigh marginal distribution
                  F_factor =  1.046 - 0.045*R(i,j) + 0.006*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  6,    % Uniform marginal distribution
                  F_factor =  1.055 + 0.015*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  7,    % Beta marginal distribution
                   ba = parameter(j,6) - parameter(j,5);
			       u1 = ( parameter(j,1) - parameter(j,5) )/ ba; 
			       s1 = parameter(j,2) / ba;
                                      
			       F_factor = 1.055-0.066*R(i,j)-0.194*u1+0.391*s1+0.003*R(i,j)*R(i,j)...
                             +0.194*u1*u1-1.134*s1*s1+0.130*R(i,j)*u1+0.003*R(i,j)*s1;
			       Ro(i,j) = F_factor * R(i,j);
               case 11,    % Type I Largest Value marginal distribution
                  F_factor =  1.064 - 0.069*R(i,j) + 0.005*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  F_factor =  1.064 + 0.069*R(i,j) + 0.005*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.056 - 0.060*R(i,j)+ 0.263*cov_j  + 0.020*R(i,j)^2 + 0.383*cov_j^2 - 0.332*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 14,    % Type III Smallest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.064 + 0.065*R(i,j)- 0.210*cov_j  + 0.003*R(i,j)^2 + 0.356*cov_j^2 - 0.211*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 8,    % Chi-square marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.031 + 0.001*R(i,j)- 0.007*cov_j  + 0.003*R(i,j)^2 + 0.131*cov_j^2 - 0.132*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 15,    % Gumbel marginal distribution
                  F_factor =  1.064 - 0.069*R(i,j) + 0.005*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.064 + 0.065*R(i,j)- 0.210*cov_j  + 0.003*R(i,j)^2 + 0.356*cov_j^2 - 0.211*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 12    % Type I Smallest Value marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  F_factor =  1.031;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  2,    % Lognormal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.029 - 0.001*R(i,j) + 0.014*cov_j + 0.004*R(i,j)^2 + 0.233*cov_j^2 + 0.197*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  3,    % Gamma marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.031 - 0.001*R(i,j)- 0.007*cov_j  + 0.003*R(i,j)^2 + 0.131*cov_j^2 + 0.132*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  4,    % Shifted Exponential marginal distribution
                  F_factor =  1.142 + 0.154*R(i,j) + 0.031*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  5,    % Shifted Rayleigh marginal distribution
                  F_factor =  1.046 + 0.045*R(i,j) + 0.006*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  6,    % Uniform marginal distribution
                  F_factor =  1.055 + 0.015*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  7,    % Beta marginal distribution
                   ba = parameter(j,6) - parameter(j,5);
			       u1 = ( parameter(j,1) - parameter(j,5) )/ ba; 
			       s1 = parameter(j,2) / ba;
			       F_factor = 1.055+0.066*R(i,j)-0.194*u1+0.391*s1+0.003*R(i,j)*R(i,j)...
                              +0.194*u1*u1-1.134*s1*s1-0.130*R(i,j)*u1-0.003*R(i,j)*s1;
			       Ro(i,j) = F_factor * R(i,j);
               case 11,    % Type I Largest Value marginal distribution
                  F_factor =  1.064 + 0.069*R(i,j) + 0.005*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  F_factor =  1.064 - 0.069*R(i,j) + 0.005*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.056 + 0.060*R(i,j) + 0.263*cov_j  + 0.020*R(i,j)^2 + 0.383*cov_j^2 + 0.332*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 14,    % Type III Smallest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.064 - 0.065*R(i,j) - 0.210*cov_j  + 0.003*R(i,j)^2 + 0.356*cov_j^2 + 0.211*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 8,    % Chi-square marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.031 - 0.001*R(i,j)- 0.007*cov_j  + 0.003*R(i,j)^2 + 0.131*cov_j^2 + 0.132*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 15,    % Gumbel marginal distribution
                  F_factor =  1.064 + 0.069*R(i,j) + 0.005*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.064 - 0.065*R(i,j) - 0.210*cov_j  + 0.003*R(i,j)^2 + 0.356*cov_j^2 + 0.211*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 13    % Type II Largest Value marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor = 1.030 + 0.238*cov_i + 0.364*cov_i^2;
                  Ro(i,j) = R(i,j) * F_factor;
               case  2,    % Lognormal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.026 + 0.082*R(i,j) - 0.019*cov_j + 0.222*cov_i + 0.018*R(i,j)^2 + 0.288*cov_j^2 ...
                              + 0.379*cov_i^2 - 0.441*R(i,j)*cov_j + 0.126*cov_i*cov_j - 0.277*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  3,    % Gamma marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.029 + 0.056*R(i,j)- 0.030*cov_j + 0.225*cov_i + 0.012*R(i,j)^2 + 0.174*cov_j^2 ...
                              + 0.379*cov_i^2 - 0.313*R(i,j)*cov_j + 0.075*cov_i*cov_j - 0.182*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  4,    % Shifted Exponential marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.109 - 0.152*R(i,j) + 0.361*cov_i + 0.130*R(i,j)^2 + 0.455*cov_i^2 - 0.728*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  5,    % Shifted Rayleigh marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.036 - 0.038*R(i,j) + 0.266*cov_i + 0.028*R(i,j)^2 + 0.383*cov_i^2 - 0.229*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  6,    % Uniform marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.033 + 0.305*cov_i + 0.074*R(i,j)^2 + 0.405*cov_i^2 ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  7,    % Beta marginal distribution
                   ba = parameter(j,6) - parameter(j,5);
			       u1 = ( parameter(j,1) - parameter(j,5) )/ ba; 
			       s1 = parameter(j,2) / ba;
			       uu = u1*u1;
			       ss = s1*s1;
			       xu = R(i,j)*u1;
			       xs = R(i,j)*s1;
			       us = u1*s1;
                   cov_i = parameter(i,2)/parameter(i,1);
			       uq = u1*cov_i;
			       sq = s1*cov_i;
			       xq = R(i,j)*cov_i;
			       F_factor = 1.005 + 0.091*R(i,j) + 0.285*u1+ 0.260*s1+ 0.199*cov_i...
                             - 0.023*R(i,j)*R(i,j) - 0.285*uu + 8.180*ss + 0.543*cov_i*cov_i - 0.181*xu...
                             - 1.744*xs - 0.336*xq - 5.450*us - 0.265*uq + 0.514*sq...
                             -19.661*ss*s1- 0.178*cov_i*cov_i*cov_i...
                             + 0.244*R(i,j)*xs + 0.066*R(i,j)*R(i,j)*cov_i - 0.001*R(i,j)*ss...
                             + 5.450*uu*s1+ 0.265*uu*cov_i - 0.986*ss*cov_i...
                             + 0.133*sq*cov_i + 3.488*R(i,j)*us + 0.671*R(i,j)*uq;
			       Ro(i,j) = F_factor * R(i,j);
               case 11,    % Type I Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.056 - 0.060*R(i,j)+ 0.263*cov_i  + 0.020*R(i,j)^2 + 0.383*cov_i^2 - 0.332*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.056 + 0.060*R(i,j)+ 0.263*cov_i  + 0.020*R(i,j)^2 + 0.383*cov_i^2 + 0.332*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.086 + 0.054*R(i,j)+ 0.104*(cov_i+cov_j)...
                     			- 0.055*R(i,j)^2 + 0.662*(cov_i^2+cov_j^2)- 0.570*R(i,j)*(cov_i+cov_j) + 0.203*cov_i*cov_j ...
                     			- 0.020*R(i,j)^3 - 0.218*(cov_i^3+cov_j^3)- 0.371*R(i,j)*(cov_i^2+cov_j^2) + 0.257*(cov_i+cov_j)*R(i,j)^2 ...
                              + 0.141*cov_i*cov_j*(cov_i+cov_j);
                  Ro(i,j) = R(i,j) * F_factor ;
               case 14,    % Type III Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.065 + 0.146*R(i,j)+ 0.241*cov_i - 0.259*cov_j...
                    				+ 0.013*R(i,j)^2 + 0.372*cov_i^2 + 0.435*cov_j^2 ...
                  				+ 0.005*R(i,j)*cov_i + 0.034*cov_i*cov_j - 0.481*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 8,    % Chi-square marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.029 + 0.056*R(i,j)- 0.030*cov_j + 0.225*cov_i + 0.012*R(i,j)^2 + 0.174*cov_j^2 ...
                              + 0.379*cov_i^2 - 0.313*R(i,j)*cov_j + 0.075*cov_i*cov_j - 0.182*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case 15,    % Gumbel marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.056 - 0.060*R(i,j)+ 0.263*cov_i  + 0.020*R(i,j)^2 + 0.383*cov_i^2 - 0.332*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.065 + 0.146*R(i,j)+ 0.241*cov_i - 0.259*cov_j...
                    				+ 0.013*R(i,j)^2 + 0.372*cov_i^2 + 0.435*cov_j^2 ...
                  				+ 0.005*R(i,j)*cov_i + 0.034*cov_i*cov_j - 0.481*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 14    % Type III Smallest Value marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor = 1.031 - 0.195*cov_i + 0.328*cov_i^2;
                  Ro(i,j) = R(i,j) * F_factor;
               case  2,    % Lognormal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.031 + 0.052*R(i,j) + 0.011*cov_j - 0.210*cov_i + 0.002*R(i,j)^2 + 0.220*cov_j^2 + 0.350*cov_i^2 + 0.005*R(i,j)*cov_j + 0.009*cov_i*cov_j - 0.174*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  3,    % Gamma marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.032 + 0.034*R(i,j)- 0.007*cov_j - 0.202*cov_i + 0.121*cov_j^2 ...
                              + 0.339*cov_i^2 - 0.006*R(i,j)*cov_j + 0.003*cov_i*cov_j - 0.111*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  4,    % Shifted Exponential marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.147 + 0.145*R(i,j) - 0.271*cov_i + 0.010*R(i,j)^2 + 0.459*cov_i^2 - 0.467*R(i,j)*cov_i ;                  
                  Ro(i,j) = R(i,j) * F_factor ;
               case  5,    % Shifted Rayleigh marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.047 + 0.042*R(i,j) - 0.212*cov_i + 0.353*cov_i^2 - 0.136*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  6,    % Uniform marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.061 - 0.237*cov_i - 0.005*R(i,j)^2 + 0.379*cov_i^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  7,    % Beta marginal distribution
                   ba = parameter(j,6) - parameter(j,5);
			       u1 = ( parameter(j,1) - parameter(j,5) )/ ba; 
			       s1 = parameter(j,2) / ba;
                   cov_i = parameter(i,2)/parameter(i,1);
                   
			       F_factor = 1.054+0.002*R(i,j)-0.176*u1+0.366*s1-0.201*cov_i...
                       -0.002*R(i,j)*R(i,j)+0.176*u1*u1-1.098*s1*s1+0.340*cov_i*cov_i...
                       -0.004*R(i,j)*u1-0.029*s1*cov_i;
			       Ro(i,j) = F_factor* R(i,j);
               case 11,    % Type I Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.064 + 0.065*R(i,j)- 0.210*cov_i  + 0.003*R(i,j)^2 + 0.356*cov_i^2 - 0.211*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.064 - 0.065*R(i,j) - 0.210*cov_i  + 0.003*R(i,j)^2 + 0.356*cov_i^2 + 0.211*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.065 + 0.146*R(i,j)+ 0.241*cov_j - 0.259*cov_i...
                    				+ 0.013*R(i,j)^2 + 0.372*cov_j^2 + 0.435*cov_i^2 ...
                  				+ 0.005*R(i,j)*cov_j + 0.034*cov_i*cov_j - 0.481*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 14,    % Type III Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.063 - 0.004*R(i,j)- 0.200*(cov_i+cov_j)...
                    				- 0.001*R(i,j)^2 + 0.337*(cov_i^2+cov_j^2)...
                  				+ 0.007*R(i,j)*(cov_i+cov_j) - 0.007*cov_i*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 8,    % Chi-square marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.032 + 0.034*R(i,j)- 0.007*cov_j - 0.202*cov_i + 0.121*cov_j^2 ...
                              + 0.339*cov_i^2 - 0.006*R(i,j)*cov_j + 0.003*cov_i*cov_j - 0.111*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 15,    % Gumbel marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.064 + 0.065*R(i,j)- 0.210*cov_i  + 0.003*R(i,j)^2 + 0.356*cov_i^2 - 0.211*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.063 - 0.004*R(i,j)- 0.200*(cov_i+cov_j)...
                    				- 0.001*R(i,j)^2 + 0.337*(cov_i^2+cov_j^2)...
                  				+ 0.007*R(i,j)*(cov_i+cov_j) - 0.007*cov_i*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 8    % Chi-square marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor = 1.001-0.007*cov_i+0.018*cov_i^2;
                  Ro(i,j) = R(i,j) * F_factor;
               case  2,    % Lognormal distribution      
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.001 + 0.033*R(i,j) + 0.004*cov_j - 0.016*cov_i + 0.002*R(i,j)^2 + 0.223*cov_j^2 ...
                              + 0.130*cov_i^2 - 0.104*R(i,j)*cov_j + 0.029*cov_j*cov_i - 0.119*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  3,    % Gamma marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.002 + 0.022*R(i,j) - 0.012*(cov_j+cov_i) + 0.001*R(i,j)^2 + 0.125*(cov_j^2+cov_i^2)...
                     - 0.077*R(i,j)*(cov_j+cov_i) + 0.014*cov_i*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  4,    % Shifted Exponential marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.104 + 0.003*R(i,j) - 0.008*cov_i + 0.014*R(i,j)^2 + 0.173*cov_i^2 - 0.296*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  5,    % Shifted Rayleigh marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.014 + 0.001*R(i,j) - 0.007*cov_i + 0.002*R(i,j)^2 + 0.126*cov_i^2 - 0.090*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  6,    % Uniform marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.023 - 0.007*cov_i + 0.002*R(i,j)^2 + 0.127*cov_i^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  7,    % Beta marginal distribution
                  ba = parameter(j,6) - parameter(j,5);
			      u1 = ( parameter(j,1) - parameter(j,5) ) / ba; 
			      s1 = parameter(j,2) / ba;
			
                  %cov_j = parameter(j,2)/parameter(j,1);
                  cov_i = parameter(i,2)/parameter(i,1);
                  
                  uu=u1*u1;
			      ss=s1*s1;
			      xu=R(i,j)*u1;
			      xs=R(i,j)*s1;
			      dus=u1*s1;
			      up=u1*cov_i;
			      sp=s1*cov_i;
			      xp=R(i,j)*cov_i;
            
			if R(i,j) > 0.0
				F_factor = 0.931+0.050*R(i,j)+0.366*u1+0.549*s1+0.181*cov_i...
                    -0.055*R(i,j)*R(i,j)-0.515*uu+4.804*ss-0.484*cov_i*cov_i-0.064*xu...
                    -0.637*xs+0.032*xp-4.059*us-0.319*up-0.211*sp...
                    +0.052*R(i,j)*R(i,j)*R(i,j)+0.227*uu*u2-10.220*ss*s2+0.559*cov_i*cov_i*cov_i-0.042*R(i,j)*xu...
                    +0.223*R(i,j)*xs-0.172*R(i,j)*xp+0.028*R(i,j)*uu+0.695*R(i,j)*ss+0.126*R(i,j)*cov_i*cov_i...
                    +3.845*uu*s2+0.019*uu*cov_i-1.244*us*s1+0.008*up*cov_i-2.075*ss*cov_i...
                    +0.167*sp*cov_i+0.666*R(i,j)*us+0.386*R(i,j)*up-0.517*R(i,j)*sp+2.125*us*cov_i;
			else
				F_factor = 1.025+0.050*R(i,j)-0.029*u1+0.047*s1-0.136*cov_i...
                    +0.069*R(i,j)*R(i,j)+0.178*uu+6.281*ss+0.548*cov_i*cov_i-0.027*xu...
                    -0.686*xs+0.046*xp-3.513*us+0.231*up+0.299*sp...
                    +0.063*R(i,j)*R(i,j)*R(i,j)-0.226*uu*u1-17.507*ss*s2-0.366*cov_i*cov_i*cov_i+0.051*R(i,j)*xu...
                    -0.246*R(i,j)*xs+0.186*R(i,j)*xp-0.001*R(i,j)*uu+0.984*R(i,j)*ss+0.121*R(i,j)*cov_i*cov_i...
                    +3.700*uu*s1+0.081*uu*cov_i+1.356*us*s1+0.002*up*cov_i+1.654*ss*cov_i...
                    -0.135*sp*cov_i+0.619*R(i,j)*us+0.410*R(i,j)*up-0.686*R(i,j)*sp-2.205*us*cov_i;
            end
        
			Ro(i,j) = F_factor*R(i,j);
               case 11,    % Type I Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.031 + 0.001*R(i,j)- 0.007*cov_i  + 0.003*R(i,j)^2 + 0.131*cov_i^2 - 0.132*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.031 - 0.001*R(i,j)- 0.007*cov_i  + 0.003*R(i,j)^2 + 0.131*cov_i^2 + 0.132*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.029 + 0.056*R(i,j)- 0.030*cov_i + 0.225*cov_j + 0.012*R(i,j)^2 + 0.174*cov_i^2 ...
                              + 0.379*cov_j^2 - 0.313*R(i,j)*cov_i + 0.075*cov_i*cov_j - 0.182*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case 14,    % Type III Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.032 + 0.034*R(i,j)- 0.007*cov_i - 0.202*cov_j + 0.121*cov_i^2 ...
                              + 0.339*cov_j^2 - 0.006*R(i,j)*cov_i + 0.003*cov_i*cov_j - 0.111*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case 8,    % Chi-square marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.002 + 0.022*R(i,j) - 0.012*(cov_j+cov_i) + 0.001*R(i,j)^2 + 0.125*(cov_j^2+cov_i^2)...
                     - 0.077*R(i,j)*(cov_j+cov_i) + 0.014*cov_i*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case 15,    % Gumbel marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.031 + 0.001*R(i,j)- 0.007*cov_i  + 0.003*R(i,j)^2 + 0.131*cov_i^2 - 0.132*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.032 + 0.034*R(i,j)- 0.007*cov_i - 0.202*cov_j + 0.121*cov_i^2 ...
                              + 0.339*cov_j^2 - 0.006*R(i,j)*cov_i + 0.003*cov_i*cov_j - 0.111*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 15    % Gumbel marginal distribution
                switch ( marg(j,1) )
               case  1,    % Normal distribution
                  F_factor =  1.031 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  2,    % Lognormal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.029 + 0.001*R(i,j) + 0.014*cov_j + 0.004*R(i,j)^2 + 0.233*cov_j^2 - 0.197*R(i,j)*cov_j ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  3,    % Gamma marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.031 + 0.001*R(i,j)- 0.007*cov_j  + 0.003*R(i,j)^2 + 0.131*cov_j^2 - 0.132*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  4,    % Shifted Exponential marginal distribution
                  F_factor =  1.142 - 0.154*R(i,j) + 0.031*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  5,    % Shifted Rayleigh marginal distribution
                  F_factor =  1.046 - 0.045*R(i,j) + 0.006*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  6,    % Uniform marginal distribution
                  F_factor =  1.055 + 0.015*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  7,    % Beta marginal distribution
                  ba = parameter(j,6) - parameter(j,5);
			      u1 = ( parameter(j,1) - parameter(j,5) )/ ba; 
			      s1 = parameter(j,2) / ba;
                                      
			      F_factor = 1.055-0.066*R(i,j)-0.194*u1+0.391*s1+0.003*R(i,j)*R(i,j)...
                             +0.194*u1*u1-1.134*s1*s1+0.130*R(i,j)*u1+0.003*R(i,j)*s1;
			      Ro(i,j) = F_factor * R(i,j);
               case 11,    % Type I Largest Value marginal distribution
                  F_factor =  1.064 - 0.069*R(i,j) + 0.005*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  F_factor =  1.064 + 0.069*R(i,j) + 0.005*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.056 - 0.060*R(i,j)+ 0.263*cov_j  + 0.020*R(i,j)^2 + 0.383*cov_j^2 - 0.332*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 14,    % Type III Smallest Value marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.064 + 0.065*R(i,j)- 0.210*cov_j  + 0.003*R(i,j)^2 + 0.356*cov_j^2 - 0.211*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 8,    % Chi-square marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.031 + 0.001*R(i,j)- 0.007*cov_j  + 0.003*R(i,j)^2 + 0.131*cov_j^2 - 0.132*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 15,    % Gumbel marginal distribution
                  F_factor =  1.064 - 0.069*R(i,j) + 0.005*R(i,j)^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.064 + 0.065*R(i,j)- 0.210*cov_j  + 0.003*R(i,j)^2 + 0.356*cov_j^2 - 0.211*R(i,j)*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 16    % Weibull marginal distribution
               switch ( marg(j,1) )
              case  1,    % Normal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor = 1.031 - 0.195*cov_i + 0.328*cov_i^2;
                  Ro(i,j) = R(i,j) * F_factor;
               case  2,    % Lognormal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.031 + 0.052*R(i,j) + 0.011*cov_j - 0.210*cov_i + 0.002*R(i,j)^2 + 0.220*cov_j^2 + 0.350*cov_i^2 + 0.005*R(i,j)*cov_j + 0.009*cov_i*cov_j - 0.174*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  3,    % Gamma marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.032 + 0.034*R(i,j)- 0.007*cov_j - 0.202*cov_i + 0.121*cov_j^2 ...
                              + 0.339*cov_i^2 - 0.006*R(i,j)*cov_j + 0.003*cov_i*cov_j - 0.111*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ; 
               case  4,    % Shifted Exponential marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.147 + 0.145*R(i,j) - 0.271*cov_i + 0.010*R(i,j)^2 + 0.459*cov_i^2 - 0.467*R(i,j)*cov_i ;                  
                  Ro(i,j) = R(i,j) * F_factor ;
               case  5,    % Shifted Rayleigh marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.047 + 0.042*R(i,j) - 0.212*cov_i + 0.353*cov_i^2 - 0.136*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  6,    % Uniform marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.061 - 0.237*cov_i - 0.005*R(i,j)^2 + 0.379*cov_i^2 ;
                  Ro(i,j) = R(i,j) * F_factor ;
               case  7,    % Beta marginal distribution
                  ba = parameter(j,6) - parameter(j,5);
			      u1 = ( parameter(j,1) - parameter(j,5) )/ ba; 
			      s1 = parameter(j,2) / ba;
                  cov_i = parameter(i,2)/parameter(i,1);
                   
			      F_factor = 1.054+0.002*R(i,j)-0.176*u1+0.366*s1-0.201*cov_i...
                       -0.002*R(i,j)*R(i,j)+0.176*u1*u1-1.098*s1*s1+0.340*cov_i*cov_i...
                       -0.004*R(i,j)*u1-0.029*s1*cov_i;
			      Ro(i,j) = F_factor* R(i,j);
              case 8,    % Chi-square marginal distribution
                   cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.032 + 0.034*R(i,j)- 0.007*cov_j - 0.202*cov_i + 0.121*cov_j^2 ...
                              + 0.339*cov_i^2 - 0.006*R(i,j)*cov_j + 0.003*cov_i*cov_j - 0.111*R(i,j)*cov_i ;
                  Ro(i,j) = R(i,j) * F_factor ;      
              case 11,    % Type I Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.064 + 0.065*R(i,j)- 0.210*cov_i  + 0.003*R(i,j)^2 + 0.356*cov_i^2 - 0.211*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 12,    % Type I Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.064 - 0.065*R(i,j) - 0.210*cov_i  + 0.003*R(i,j)^2 + 0.356*cov_i^2 + 0.211*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 13,    % Type II Largest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.065 + 0.146*R(i,j)+ 0.241*cov_j - 0.259*cov_i...
                    				+ 0.013*R(i,j)^2 + 0.372*cov_j^2 + 0.435*cov_i^2 ...
                  				+ 0.005*R(i,j)*cov_j + 0.034*cov_i*cov_j - 0.481*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 14,    % Type III Smallest Value marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.063 - 0.004*R(i,j)- 0.200*(cov_i+cov_j)...
                    				- 0.001*R(i,j)^2 + 0.337*(cov_i^2+cov_j^2)...
                  				+ 0.007*R(i,j)*(cov_i+cov_j) - 0.007*cov_i*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
                               
               case 15,    % Gumbel marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  F_factor =  1.064 + 0.065*R(i,j)- 0.210*cov_i  + 0.003*R(i,j)^2 + 0.356*cov_i^2 - 0.211*R(i,j)*cov_i;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 16,    % Weibull marginal distribution
                  cov_i = parameter(i,2)/parameter(i,1);
                  cov_j = parameter(j,2)/parameter(j,1);
                  F_factor =  1.063 - 0.004*R(i,j)- 0.200*(cov_i+cov_j)...
                    				- 0.001*R(i,j)^2 + 0.337*(cov_i^2+cov_j^2)...
                  				+ 0.007*R(i,j)*(cov_i+cov_j) - 0.007*cov_i*cov_j;
                  Ro(i,j) = R(i,j) * F_factor ;
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
               
            elseif marg(i,1) == 18    % Laplace marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
               case  2,    % Lognormal distribution
               case  3,    % Gamma marginal distribution
               case  4,    % Shifted Exponential marginal distribution
               case  5,    % Shifted Rayleigh marginal distribution
               case  6,    % Uniform marginal distribution
               case  7,    % Beta marginal distribution
               case 8,    % Chi-square marginal distribution
               case 11,    % Type I Largest Value marginal distribution
               case 12,    % Type I Smallest Value marginal distribution
               case 13,    % Type II Largest Value marginal distribution
               case 14,    % Type III Smallest Value marginal distribution
               
               case 15,    % Gumbel marginal distribution
               case 16,    % Weibull marginal distribution
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
               
            elseif marg(i,1) == 19    % Pareto marginal distribution
               switch ( marg(j,1) )
               case  1,    % Normal distribution
               case  2,    % Lognormal distribution
               case  3,    % Gamma marginal distribution
               case  4,    % Shifted Exponential marginal distribution
               case  5,    % Shifted Rayleigh marginal distribution
               case  6,    % Uniform marginal distribution
               case  7,    % Beta marginal distribution
               case 8,    % Chi-square marginal distribution
               case 11,    % Type I Largest Value marginal distribution
               case 12,    % Type I Smallest Value marginal distribution
               case 13,    % Type II Largest Value marginal distribution
               case 14,    % Type III Smallest Value marginal distribution
               
               case 15,    % Gumbel marginal distribution
               case 16,    % Weibull marginal distribution
               case 18,    % Laplace marginal distribution
               case 19,    % Pareto marginal distribution
               otherwise,
                  Ro(i,j) = R(i,j);
               end
            end
            
         end
      end
   end
end

