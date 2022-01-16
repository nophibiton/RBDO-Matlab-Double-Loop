function sorm_curfit_results = sorm_curvature_fitting(lsf,formresults,probdata,analysisopt,gfundata,femodel,randomfield)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function compute an approximation of the probability of failure % 
% using the curvature-fitting method                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_complex = sqrt(-1);
   %[formresults] = form(lsf,probdata,analysisopt,gfundata,femodel,randomfield);
   
   grad_flag = lower(analysisopt.grad_flag);

   
   iter = formresults.iter;
   beta = formresults.beta1;
   pf1 = formresults.pf1;
   dsptu = formresults.dsptu;
   alpha = formresults.alpha;
   dsptx = formresults.dsptx';
   imptg = formresults.imptg;
   gfcn = formresults.gfcn;
   stpsz = formresults.stpsz;
   
   
   marg      = probdata.marg;
   R         = probdata.correlation;
      
   %Compute the gradient of G
    % Determine number of random variables
	marg_dim = size(marg,1);

	% Compute corrected correlation coefficients
	Ro = mod_corr( probdata, R );

	% Cholesky decomposition
	Lo = (chol(Ro))';
   	iLo = inv(Lo);
   
   	J_u_x = jacobian(dsptx,dsptu,probdata,Lo,iLo);
   	J_x_u = inv(J_u_x);
   	[ G, grad_g ] = gfun(lsf,dsptx,grad_flag,probdata,gfundata,femodel,randomfield);
   	grad_G = (grad_g * J_x_u)';   % Gradient of G evaluate at the design point

   
   % Compute the rotation matrix
   R1 = gram_schmidt(alpha');
   
   % Compute the hessian matrix....
   hess_G = hessian_G(1,dsptx,dsptu,G,probdata,gfundata,femodel,randomfield);
   
   A = R1*hess_G*R1' / norm(grad_G);

   [eigenvectors,D] = eig(A([1:(marg_dim-1)],[1:(marg_dim-1)]));
   
   for i=1:marg_dim-1
      kappa(i) = D(i,i);
   end
  
   %[V,D] = EIG(X) produces a diagonal matrix D of eigenvalues and a
   %full matrix V whose columns are the corresponding eigenvectors so
   %that X*V = V*D.

% Breitung
pf2_breitung = ferum_cdf(1,-beta,0,1) * prod(1./(1+beta*kappa).^0.5);
beta2_breitung = -inv_norm_cdf(pf2_breitung);

% Breitung modified by Hohenbichler / Rackwitz
kk = ferum_pdf(1,beta,0,1) / ferum_cdf(1,-beta,0,1);
pf2_breitung_mod = ferum_cdf(1,-beta,0,1) * prod(1./(1+kk*kappa).^0.5);
beta2_breitung_mod = -inv_norm_cdf(pf2_breitung_mod);

% Tvedt three-term approximation
A1 = ferum_cdf(1,-beta,0,1) * prod(1./(1+beta*kappa).^0.5);
A2 = (beta*ferum_cdf(1,-beta,0,1) - ferum_pdf(1,beta,0,1))*(prod(1./(1+beta*kappa).^0.5)-prod(1./(1+(beta+1)*kappa).^0.5));
A3 = (beta+1)*(beta*ferum_cdf(1,-beta,0,1) - ferum_pdf(1,beta,0,1))*(prod(1./(1+beta*kappa).^0.5)-real(prod(1./(1+(beta+i_complex)*kappa).^0.5)));
pf2_tvedt = A1 + A2 + A3;
beta_tvedt = -inv_norm_cdf(pf2_tvedt);

% Tvedt single integral formula using 10 point from Gauss-Laguerre quadrature
%num_gauss = 4;
%[s,w] = GaussLaguerrePoints(num_gauss);
%sum = 0;
%for k = 1:num_gauss
%    bb = beta^2 + 2*s(k);
%    sum = sum + w(k)*real(prod((1 + (bb^0.5 + i_complex)*kappa).^(-0.5)))*bb^(-0.5);
%end
%pf2_tvedt_single_integral = ferum_pdf(1,beta,0,1)*sum;
%beta2_tvedt_single_integral = -inv_norm_cdf(pf2_tvedt_single_integral);

% Tvedt exact integral formula using saddle-point method and trapezoidale rule
us = saddle_point(beta,kappa);
pf2_tvedt_EI = pf2_Tvedt_EI(us,beta,kappa);
beta2_tvedt_EI = -inv_norm_cdf(pf2_tvedt_EI);  

% SORM Post-processing
   
   	% FORM results
   	sorm_curfit_results.iter = iter;                                % Number_of_iterations in FORM analysis
   	sorm_curfit_results.beta1 = beta;                                % Reliability_index_beta FORM analysis
   	sorm_curfit_results.pf1 = pf1;                  				% Failure_probability_pf1   
        
    % SORM results
    sorm_curfit_results.R1 = R1;									 % Rotation matrix 
    sorm_curfit_results.hess_G = hess_G;                             % Hessian of G evaluated at the design point
    sorm_curfit_results.kappa = kappa';                              % Curvatures in the (n-1)(n-1) space
    sorm_curfit_results.eig_vectors=eigenvectors;                    % Matrix of eigenvectors
    sorm_curfit_results.beta2_breitung = beta2_breitung;             % Reliability index beta2  according to Breitung formula 
   	sorm_curfit_results.pf2_breitung = pf2_breitung;                 % Failure probability pf2  according to Breitung formula
    sorm_curfit_results.beta2_breitung_mod = beta2_breitung_mod;     % Reliability index beta2  according to improved  Breitung formula 
   	sorm_curfit_results.pf2_breitung_mod = pf2_breitung_mod;         % Failure probability pf2  according to improved Breitung formula
    sorm_curfit_results.beta2_tvedt_EI = beta2_tvedt_EI;             % Reliability index beta2  according to Tvedt exact integral formula 
   	sorm_curfit_results.pf2_tvedt_EI = pf2_tvedt_EI;                 % Failure probability pf2  according to Tvedt exact integral formula

    
  
