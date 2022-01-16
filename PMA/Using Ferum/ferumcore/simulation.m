function simulationresults = simulation(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

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

% Exctract model data
marg = probdata.marg;
R = probdata.correlation;

point = analysisopt.sim_point;
stdv1 = analysisopt.stdv_sim;
num_sim = analysisopt.num_sim;
target_cov = analysisopt.target_cov;

% Find number of random variables
nrv = length(point);

% Modify correlation matrix and perform Cholesky decomposition
Ro = mod_corr( probdata, R );
Lo = (chol(Ro))';

% Establish covariance matrix, its Cholesky decomposition, and its inverse
covariance = stdv1^2 * eye(nrv);
chol_covariance = chol(covariance);
inv_covariance = inv(covariance);

% Initializations
sum_q = 0;
sum_q_squared = 0;

% Pre-compute some factors to minimize computations inside simulation loop
factor1 = 1 / ( (2*pi)^(nrv/2) );
factor2 = 1 / ( (2*pi)^(nrv/2) * sqrt(det(covariance)) );


k = 1;
cov_of_q_bar(k) = 1.0;

while( k<=num_sim  & cov_of_q_bar(length(cov_of_q_bar))>target_cov ) 
   
   % Do the simulation (create array of random numbers)
   u = point + chol_covariance * randn(nrv,1);
   
   % Transform into original space
   [ x ] = u_to_x(u,probdata,Lo);
   
   % Evaluate limit-state function
   [ g, dummy ] = gfun(lsf,x,'no ',probdata,gfundata,femodel,randomfield);
   
   % Collect result of sampling
   if g < 0
      I = 1;
   else
      I = 0;
   end
   
   % Compute values of joint distributions at the u-point
   phi = factor1 * exp( -0.5 * u' * u );
   h   = factor2 * exp( -0.5 * (u-point)' * inv_covariance * (u-point) );
   
   % Update sums
   q = I * phi / h;
   sum_q = sum_q + q;
   sum_q_squared = sum_q_squared + q^2;
   
   % Compute coefficient of variation (of pf) for each simulation
   if sum_q > 0
      q_bar = 1/k * sum_q;
      variance_of_q_bar = 1/k * ( 1/k * sum_q_squared - (1/k*sum_q)^2);
      cov_of_q_bar(k) = sqrt(variance_of_q_bar) / q_bar;
      if cov_of_q_bar(k) == 0
         cov_of_q_bar(k) = 1.0;
      end
   else
      cov_of_q_bar(k) = 1.0;
   end
   
   k = k + 1;
end

if sum_q > 0
   % Compute probability of failure and reliability index
   pf = q_bar;
   cov = cov_of_q_bar(k-1);
   beta = -inv_norm_cdf(pf);
   
   % Plot how the coefficient of variation varies with number of simulations
   figure
   plot( (1:1:length(cov_of_q_bar)),cov_of_q_bar,'b-');
   axis([0 length(cov_of_q_bar) 0 1]);
   title('C.O.V. of probability of failure (\delta_p_f)')
   xlabel('Number of simulations');
   ylabel('Coefficient of variation');
else
   pf = 0;
   cov = 0;
   beta = 0;
end

simulationresults.pf     = pf;
simulationresults.cov_pf = cov;
simulationresults.beta   = beta;
simulationresults.num_sim = k-1;