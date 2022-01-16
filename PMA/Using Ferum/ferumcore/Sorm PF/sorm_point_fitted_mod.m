function  sorm_point_fitted_newton_results = sorm_point_fitted_mod(lsf,formresults,probdata,analysisopt,gfundata,femodel,randomfield)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the improved version of the point fitting SORM    %
% method suggested by Araya and Der Kiureghian (1988)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%[formresults] = form(lsf,probdata,analysisopt,gfundata,femodel,randomfield);


% Useful available parameters after Form analyis 
	beta = formresults.beta1;
	pf1 = formresults.pf1;
	alpha = formresults.alpha;
	iter = formresults.iter;
	%design point
	dsptx = formresults.dsptx;
	dsptu = formresults.dsptu;

   
% Rotation matrix obtained by Gram-Schmidt scheme   
R1 = gram_schmidt(alpha');  % Rotation matrix computation

marg = probdata.marg;
nrv = size(marg,1);

% These parameters are necessary to compute G and Grad_G
	R         = probdata.correlation;
	grad_flag = lower(analysisopt.grad_flag);
	% Determine number of random variables
	marg_dim = size(marg,1);
	% Compute corrected correlation coefficients
	Ro = mod_corr( probdata, R );
	% Cholesky decomposition
	Lo = (chol(Ro))';
	iLo = inv(Lo);

% Determination of the coefficient k (limit-state surface intersection with path Pi)
if abs(beta)<1
   k = 1 / abs(beta);
elseif abs(beta)>=1 & abs(beta) <= 3
   k = 1;
else
   k = 3 / abs(beta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial trials points of ordinates +beta %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_prime_1 = [ [ -k*beta*eye(nrv-1) ; beta*ones(1,nrv-1) ] [ k*beta*eye(nrv-1) ; beta*ones(1,nrv-1) ] ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determination of the fitting points in the rotated space %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the fitting points on the negative side of axes and then on the positive side of axes
for i=1:2*(nrv-1)
   	num = i;
   	[u_prime_i,G_u] = fitting_point(Ro,Lo,iLo,k,R1,U_prime_1,num,lsf,formresults,probdata,analysisopt,gfundata,femodel,randomfield);
      U_prime_final(:,i) = u_prime_i;
      g_final(i) = G_u;
   end

 U_prime_final_negative = U_prime_final(:,[1:nrv-1]);
 U_prime_final_positive = U_prime_final(:,[nrv:2*(nrv-1)]);

% Computes the curvatures a_i_+/-
for i=1:nrv-1
a_curvatures_minus(i) = 2*(U_prime_final_negative(nrv,i)-beta)/(U_prime_final_negative(i,i))^2;
a_curvatures_plus(i) = 2*(U_prime_final_positive(nrv,i)-beta)/(U_prime_final_positive(i,i))^2;
end
a_curvatures_plus;
a_curvatures_minus;

kappa_plus_minus(1,:)  = a_curvatures_plus;
kappa_plus_minus(2,:) = a_curvatures_minus;

% Curvatures ai
%lala = (1./((1/2 * ( (1./(1+beta*a_curvatures_plus).^0.5) + (1./(1+beta*a_curvatures_minus).^0.5) )).^2)-1)./beta
%lala/2
%lala_1 = (1./(1+beta*a_curvatures_plus).^0.5)
%lala_2 = (1./(1+beta*a_curvatures_minus).^0.5)
%0.5*(1/beta).*((4./(lala_1 + lala_2).^2)-1)

% Coordinates of fitting points in rotated space 
% along with the value of lsf at the fitting points and the curvature a_i_+/- in a same matrix

% Along minus axis
U_prime_minus = zeros(nrv-1,5);

U_prime_minus(:,1)= round([1:1:nrv-1]');
for i = 1:nrv-1
   U_prime_minus(i,2) = U_prime_final(i,i);
end
U_prime_minus(:,3)= U_prime_final_negative(nrv,:)';
U_prime_minus(:,4)= g_final([1:(nrv-1)])';
U_prime_minus(:,5)= a_curvatures_minus';
U_prime_minus;

%Along plus axis
U_prime_plus = zeros(nrv-1,5);
U_prime_plus(:,1)= [1:1:(nrv-1)]';
for i = nrv:(2*(nrv-1))
   U_prime_plus(i-nrv+1,2) = U_prime_final(i-nrv+1,i);
end
U_prime_plus(:,3)= U_prime_final_positive(nrv,:)';
U_prime_plus(:,4)= g_final([nrv:2*(nrv-1)])';
U_prime_plus(:,5)= a_curvatures_plus';
U_prime_plus;

% Probability of failure & beta Breitung
pf2_breitung = ferum_cdf(1,-beta,0,1) * prod( 1/2 * ( (1./(1+beta*a_curvatures_plus).^0.5) + (1./(1+beta*a_curvatures_minus).^0.5) ) );
beta2_breitung = -inv_norm_cdf(pf2_breitung);

% Probability of failure & beta Breitung modified by Hohenbichler / Rackwitz
kk = ferum_pdf(1,beta,0,1) / ferum_cdf(1,-beta,0,1);
pf2_breitung_mod = ferum_cdf(1,-beta,0,1) * prod( 1/2 * ( (1./(1+kk*a_curvatures_plus).^0.5) + (1./(1+kk*a_curvatures_minus).^0.5) ) );
beta2_breitung_mod = -inv_norm_cdf(pf2_breitung_mod);


us = saddle_point(beta,kappa_plus_minus);
pf2_tvedt_EI = pf2_Tvedt_EI(us,beta,kappa_plus_minus);
beta2_tvedt_EI = -inv_norm_cdf(pf2_tvedt_EI);


% Point-fitted SORM Post-processing
   
   	% FORM results
      %  Post-processing 
   	  sorm_point_fitted_newton_results.iter  = iter;                                % Number_of_iterations in FORM analysis
   	  sorm_point_fitted_newton_results.beta1 = beta;                                % Reliability_index_beta 
   	  sorm_point_fitted_newton_results.pf1   = pf1;                                 % Failure_probability_pf1   
   	  
      % SORM results
      %sorm_point_fitted_newton_results.max_iter_fitting = max_iter_fitting;		% Number_of_iterations for each fitting points
      
      sorm_point_fitted_newton_results.R1 = R1;								        % Rotation Matrix (G-S scheme) 
      
	  sorm_point_fitted_newton_results.beta2_breitung = beta2_breitung;             %Beta Breitung
	  sorm_point_fitted_newton_results.pf2_breitung = pf2_breitung;                 %Probability of failure Breitung

	  sorm_point_fitted_newton_results.beta2_breitung_mod = beta2_breitung_mod;     %Beta Breitung modified by Hohenbichler / Rackwitz
	  sorm_point_fitted_newton_results.pf2_breitung_mod = pf2_breitung_mod;         %Failure Probability Breitung modified by Hohenbichler / Rackwitz
      
      sorm_point_fitted_newton_results.beta2_tvedt_EI = beta2_tvedt_EI;             %Beta Tvedt Exact integral
      sorm_point_fitted_newton_results.pf2_tvedt_EI = pf2_tvedt_EI;                 %Failure Probability Exact integral
      
      sorm_point_fitted_newton_results.U_prime_plus = U_prime_plus;
      sorm_point_fitted_newton_results.U_prime_minus = U_prime_minus;
      
      sorm_point_fitted_newton_results.kapa_plus = a_curvatures_plus;               % Curvature along positive semi-axis
      sorm_point_fitted_newton_results.kapa_minus = a_curvatures_minus;             % Curvature along negative semi-axis
      
      
      %sorm_point_fitted_newton_results.kappa_plus_minus = kappa_plus_minus;
