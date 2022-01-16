function hess_G = hessian_G(lsf,dsptx,dsptu,G,probdata,gfundata,femodel,randomfield)



marg = probdata.marg;
marg_dim = size(marg,1);
R        = probdata.correlation;

% Compute corrected correlation coefficients
Ro = mod_corr( probdata, R );
% Cholesky decomposition
Lo = (chol(Ro))';
iLo = inv(Lo);


allG_a_step_ahead = zeros(marg_dim);
hess_G = zeros(marg_dim);

h = 1/2000;

for i = 1:marg_dim

   perturbed_u_a_step_ahead = dsptu;
   perturbed_u_a_step_back = dsptu;

   perturbed_u_a_step_ahead(i) = perturbed_u_a_step_ahead(i) + h;
   perturbed_u_a_step_back(i) = perturbed_u_a_step_back(i) - h;
   
   % ----------- Transformation from u to x space ---------------%
   
   perturbed_x_a_step_ahead = u_to_x(perturbed_u_a_step_ahead,probdata,Lo);
   perturbed_x_a_step_back = u_to_x(perturbed_u_a_step_back,probdata,Lo);

   [ G_a_step_ahead, dummy ] = gfun(lsf,perturbed_x_a_step_ahead,'no ',probdata,gfundata,femodel,randomfield);
   [ G_a_step_back, dummy ] = gfun(lsf,perturbed_x_a_step_back,'no ',probdata,gfundata,femodel,randomfield);
   
   allG_a_step_ahead(i,i) = G_a_step_ahead;
   
   hess_G(i,i) = (G_a_step_ahead-2*G+G_a_step_back)/h^2;
   
   perturbed_u_a_step_ahead_i = perturbed_u_a_step_ahead;
   
   for j = 1:(i-1)
      
      perturbed_u_a_step_ahead = perturbed_u_a_step_ahead_i;
      
      perturbed_u_a_step_ahead(j) = perturbed_u_a_step_ahead(j) + h;
      
      % ----------- Transformation from u to x space ---------------%
      perturbed_x_a_step_ahead = u_to_x(perturbed_u_a_step_ahead,probdata,Lo);
      
      [ G_a_step_ahead, dummy ] = gfun(lsf,perturbed_x_a_step_ahead,'no ',probdata,gfundata,femodel,randomfield);
      
      allG_a_step_ahead(i,j) = G_a_step_ahead;
   end      
end

for i = 1:marg_dim
   for j = 1:(i-1)
      hess_G(i,j) = ( (allG_a_step_ahead(i,j)-allG_a_step_ahead(j,j)) - (allG_a_step_ahead(i,i)- G) ) / h^2;
      hess_G(j,i) = hess_G(i,j);
   end
end


   
