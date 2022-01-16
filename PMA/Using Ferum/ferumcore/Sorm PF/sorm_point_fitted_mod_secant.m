function  sorm_point_fitted_secant_results = sorm_point_fitted_mod_secant(lsf,formresults,probdata,analysisopt,gfundata,femodel,randomfield)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the improved version of the point fitting SORM    %
% method suggested by Araya and Der Kiureghian (1988)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximum number of iteration for each fitting points
threshold = 10^-7;
stop_flag = 0;

%[formresults] = form(lsf,probdata,analysisopt,gfundata,femodel,randomfield);


% Useful available parameters after Form analyis 
    iter = formresults.iter;	
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
%U_prime_final = zeros(nrv,2*(nrv-1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determination of the fitting points in the rotated space %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Points on the negative sides of axes
for i=1:2*(nrv-1)
    stop_flag = 0;
   if i <= nrv-1         % find the point on Negative axis
   	counter = i;
   	sign = -1;
	else                 % find the point on Positive axis
   	counter = i-nrv+1;
   	sign = 1;
   end
   
   u_prime_i = U_prime_1(:,i);
   u = R1'*u_prime_i;
   x = u_to_x(u,probdata,Lo);
   
   g1 = gfun(lsf,x,'no ',probdata,gfundata,femodel,randomfield);
   eta1 = u_prime_i(nrv);
   
   if g1 > 0  % trying point in the safe domain
      
      u_prime_i(nrv)= (1+0.5*k^2)*beta;
      u_prime_i(counter)= sign*((k*beta)^2-(u_prime_i(nrv)-beta)^2)^0.5;
      
      u = R1'*u_prime_i;
      x = u_to_x(u,probdata,Lo);
      
      eta2 = u_prime_i(nrv);
      g2 = gfun(lsf,x,'no ',probdata,gfundata,femodel,randomfield);
      
      
      eta3 = eta1 - g1*(eta1-eta2)/(g1-g2);
      
      u_prime_i(nrv)= eta3;
      u_prime_i(counter)= sign*((k*beta)^2-(u_prime_i(nrv)-beta)^2)^0.5;
      u = R1'*u_prime_i;
      x = u_to_x(u,probdata,Lo);

      g3 = gfun(lsf,x,'no ',probdata,gfundata,femodel,randomfield);
       
      % initialize function while parameters
      g_1 = g2;
	  g_2 = g3;
	  n_1 = eta2;
	  n_2 = eta3;
                       
      while stop_flag == 0; 
        n_new = n_1 - g_1 * (n_1-n_2) / (g_1-g_2);
        u_prime_i(nrv)= n_new;
      	u_prime_i(counter)= sign*((k*beta)^2-(u_prime_i(nrv)-beta)^2)^0.5;
      	u = R1'*u_prime_i;
        x = u_to_x(u,probdata,Lo);
        
        g_new = gfun(lsf,x,'no ',probdata,gfundata,femodel,randomfield);
   
   		g_1 = g_2;
		g_2 = g_new;
		n_1 = n_2;
		n_2 = n_new;
   		
        % Test to break loop in function while
        if abs(g_new) < threshold
            stop_flag = 1;
        end
        
	  end 
      
      U_prime_final(:,i) = u_prime_i;
      g_final(i)= g_new;
      
   else % trying point in the failure domain
      
      u_prime_i(nrv)= (1-0.5*k^2)*beta;
      u = R1'*u_prime_i;
      x = u_to_x(u,probdata,Lo);
      g2 = gfun(lsf,x,'no ',probdata,gfundata,femodel,randomfield);
      eta2 = u_prime_i(nrv);
      
      eta3 = eta1 - g1*(eta1-eta2)/(g1-g2);
      
      u_prime_i(nrv)= eta3;
      
      u = R1'*u_prime_i;
      x = u_to_x(u,probdata,Lo);
      g3 = gfun(lsf,x,'no ',probdata,gfundata,femodel,randomfield);
       
      % initialize function while parameters
      g_1 = g2;
	  g_2 = g3;
	  n_1 = eta2;
	  n_2 = eta3;
                 
      while stop_flag == 0;  
         n_new = n_1 - g_1 * (n_1-n_2) / (g_1-g_2);
         
         u_prime_i(nrv)= n_new;
      	
         u = R1'*u_prime_i;
         x = u_to_x(u,probdata,Lo);
   		 g_new = gfun(lsf,x,'no ',probdata,gfundata,femodel,randomfield);
   
   		 g_1 = g_2;
		 g_2 = g_new;
		 n_1 = n_2;
		 n_2 = n_new;
   		            
         if abs(g_new) < threshold
           stop_flag = 1;
         end
         
		end 
      
      u_prime_i;
      U_prime_final(:,i)=u_prime_i;
      g_final(i)= g_new;
     
   end
end

%U_prime_final;
U_prime_final_negative = U_prime_final(:,[1:nrv-1]);
U_prime_final_positive = U_prime_final(:,[nrv:2*(nrv-1)]);

% Curvatures a_i_+/-
for i=1:nrv-1
a_curvatures_minus(i) = 2*(U_prime_final_negative(nrv,i)-beta)/(U_prime_final_negative(i,i))^2;
a_curvatures_plus(i) = 2*(U_prime_final_positive(nrv,i)-beta)/(U_prime_final_positive(i,i))^2;
end
a_curvatures_plus;
a_curvatures_minus;

kappa_plus_minus(1,:)  = a_curvatures_plus;
kappa_plus_minus(2,:) = a_curvatures_minus;


% Coordinates of fitting points in rotated space
% Along minus axis
U_prime_minus = zeros(nrv-1,5);
U_prime_minus(:,1)= [1:1:nrv-1]';
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
      sorm_point_fitted_secant_results.iter = iter;                                  % Number_of_iterations in FORM analysis
   	  sorm_point_fitted_secant_results.beta_form = beta;                             % Reliability_index_beta 
   	  sorm_point_fitted_secant_results.pf1_form = pf1;                               % Failure_probability_pf1   
   	  %sorm_point_fitted_secant_results.dsptu = dsptu;                               % Design_point_u_star
      %sorm_point_fitted_secant_results.alpha = alpha;                               % Alpha_vector
      %sorm_point_fitted_secant_results.dsptx = dsptx;                               % Design_point_in_original_space 
   	
      % SORM results
      sorm_point_fitted_secant_results.R1 = R1;								         % Rotation Matrix (G-S scheme) 
      sorm_point_fitted_secant_results.beta2_breitung = beta2_breitung;              %Beta Breitung
	  sorm_point_fitted_secant_results.pf2_breitung = pf2_breitung;                  %Probability of failure Breitung
	  sorm_point_fitted_secant_results.beta2_breitung_mod = beta2_breitung_mod;      %Beta Breitung modified by Hohenbichler / Rackwitz
	  sorm_point_fitted_secant_results.pf2_breitung_mod = pf2_breitung_mod;          %Probability of failure Breitung modified by Hohenbichler / Rackwitz
      
      sorm_point_fitted_secant_results.beta2_tvedt_EI = beta2_tvedt_EI;
      sorm_point_fitted_secant_results.pf2_tvedt_EI = pf2_tvedt_EI;
      
      
      sorm_point_fitted_secant_results.U_prime_plus = U_prime_plus;
      sorm_point_fitted_secant_results.U_prime_minus = U_prime_minus;
      
      sorm_point_fitted_secant_results.kapa_plus = a_curvatures_plus;
      sorm_point_fitted_secant_results.kapa_minus = a_curvatures_minus;
            
      %sorm_point_fitted_secant_results.kappa_plus_minus = kappa_plus_minus;




