function  [ beta_sensitivities_thetag, pf_sensitivities_thetag ] = sensi_wrt_thetag(lsf,x,beta,J_x_u,probdata,analysisopt,gfundata,femodel,randomfield)

grad_flag = lower(analysisopt.grad_flag);


switch lower(gfundata(lsf).evaluator)
   
   
case 'basic'

if grad_flag == 'ddm'
      for i = 1 : length(gfundata(lsf).dgthetag)
  		dgthetag(i) = eval(gfundata(lsf).dgthetag{i}) 
		end
elseif grad_flag == 'ffd' 
[ G, dummy ] = gfun(lsf,x','no ',probdata,gfundata,femodel,randomfield);
   for i = 1 : length(gfundata(lsf).thetag)
      h = (gfundata(lsf).thetag(i)*0.1) /200;
      thetag_old(i) = gfundata(lsf).thetag(i);
      gfundata(lsf).thetag(i) = gfundata(lsf).thetag(i) + h;
      
      [ G_a_step_ahead, dummy ] = gfun(lsf,x','no ',probdata,gfundata,femodel,randomfield);
      dgthetag(i) = (G_a_step_ahead - G)/h;
      
      gfundata(lsf).thetag(i)= thetag_old(i);
   end
end

[ G, grad_g ] = gfun(lsf,x',grad_flag,probdata,gfundata,femodel,randomfield);
normgradG = norm((grad_g * J_x_u)');

beta_sensitivities_thetag = (1/normgradG)*dgthetag;
pf_sensitivities_thetag = -ferum_pdf(1,beta,0,1)*beta_sensitivities_thetag;

case 'ferumlinearfecode'
dgthetag = 1;

[ G, grad_g ] = gfun(lsf,x',grad_flag,probdata,gfundata,femodel,randomfield);
normgradG = norm((grad_g * J_x_u)');

beta_sensitivities_thetag = (1/normgradG)*dgthetag;
pf_sensitivities_thetag = -ferum_pdf(1,beta,0,1)*beta_sensitivities_thetag;

case 'ferumnonlinearfecode'
   
   
case 'fedeas'
dgthetag = 1;

[ G, grad_g ] = gfun(lsf,x',grad_flag,probdata,gfundata,femodel,randomfield);
normgradG = norm((grad_g * J_x_u)');

beta_sensitivities_thetag = (1/normgradG)*dgthetag;
pf_sensitivities_thetag = -ferum_pdf(1,beta,0,1)*beta_sensitivities_thetag;  
otherwise
   disp('Unknown type of limit-state function evaluator.')
   
   
end   
