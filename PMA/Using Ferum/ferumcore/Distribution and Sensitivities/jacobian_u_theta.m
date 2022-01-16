function  [ J_u_theta, beta_sensitivities, pf_sensitivities ] = jacobian_u_theta(xstar,probdata,Lo,iLo,alpha,beta)


% This function computes the partial derivatives of the transformation with respect to the distribution parameters
% and then computes the jacobian of the probability distribution

% 
marg = probdata.marg;
parameter = probdata.parameter;

marg_dim = size(marg,1);

% Detect the type of the distribution and then give the number of distribution parameters
%parameter_dim = 0; sensitivity = zeros(marg_dim,1)
%for i=1:marg_dim
%   if marg(i,1)== 1 | marg(i,1)== 2 | marg(i,1)== 3 | marg(i,1)== 4 | marg(i,1)== 5 | marg(i,1)== 6 | ...
%      marg(i,1)== 11 | marg(i,1)== 12 | marg(i,1)== 13 | marg(i,1)== 15 | marg(i,1)== 17
%      parameter_dim = parameter_dim + 4; % variable mean stdv p1 and p2
%   elseif marg(i,1)== 7 | marg(i,1)== 14
%      parameter_dim = parameter_dim + 5;  % variable mean stdv p1, p2 and p3
%      sensitivity_flag(i) = 1;
%   end
%end


% Number of distribution parameters: for each random variable mean stdv p1 and p2

parameter_dim = 6 * marg_dim;
J_z_theta=zeros(marg_dim,parameter_dim);

for i=1:marg_dim
   
   switch marg(i,1)
      
   case 1   % (Normal marginal distribution)
      %  z(i) =  ( x(i) - marg(i,2) ) / marg(i,3) ;
      mean = parameter(i,1);
      stdv = parameter(i,2);
      
      dzdmean = -1/stdv;
      dzdstdv = -(xstar(i) - mean )/(stdv)^2;
      
      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzdmean;
      J_z_theta(i,6*(i-1)+4) = dzdstdv;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;
      
      %partial_Z(i,2*i-1) = -1/marg(i,3);
      %partial_Z(i,2*i) = -( xstar - marg(i,2) ) / marg(i,3)^2 ;
      
   case 2   % Lognormal marginal distribution
      %    z(i) =  ( log(x(i)) - lambda ) / zeta ;
      mean = parameter(i,1);
      stdv = parameter(i,2);
      lambda = parameter(i,3);
      zeta = parameter(i,4);
      
     
      dzdlambda = -1/zeta;
      dzdzeta = -(log(xstar(i))-lambda)/zeta^2;
      
      dzdmean = dzdlambda*(1/mean + (stdv^2)/(mean^3+(stdv^2)*mean)) + dzdzeta*(-stdv^2)/((mean^3+(stdv^2)*mean)*(log(1+(stdv/mean)^2))^0.5);
      dzdstdv = dzdlambda*(-stdv/(stdv^2+mean^2)) + dzdzeta*stdv/((mean^2+stdv^2)*(log(1+(stdv/mean)^2))^0.5);

	  J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzdlambda;
      J_z_theta(i,6*(i-1)+4) = dzdzeta;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;


      %partial_Z(i,2*i-1)=-1/zeta;
      %partial_Z(i,2*i)=-(log(xstar)-lambda)/zeta^2;
      
   case 3  % Gamma marginal distribution
      %    z(i)= inv_norm_cdf(gammainc(lambda*xstar(i),k));
      mean = parameter(i,1);
      stdv = parameter(i,2); 
      lambda = parameter(i,3);
      k = parameter(i,4);
      
      ccc = inv_norm_cdf(gammainc(lambda*xstar(i),k));
      dzdlambda = ((xstar(i)*(lambda*xstar(i))^(k-1))/gamma(k))*exp(-lambda*xstar(i))/ferum_pdf(1,ccc,0,1);
      % Computed dzdk with a forward finite difference scheme
      h = 1/2000;
      cc = (gammainc(lambda*xstar(i),k + h) - gammainc(lambda*xstar(i),k))/h;
      dzdk = cc/ferum_pdf(1,ccc,0,1);   
      
      dzdmean = dzdk*(2*mean/stdv^2) + dzdlambda*(1/stdv^2) ;
      dzdstdv = dzdk*(-2*(mean^2)/stdv^3) + dzdlambda*(-2*mean/stdv^3) ;
      
      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzdlambda;
      J_z_theta(i,6*(i-1)+4) = dzdk;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;


      %partial_Z(i,2*i-1)=1/0.00001*(norminv(gammainc((lambda+0.00001)*xstar(i),k),0,1)-norminv(gammainc(lambda*xstar(i),k),0,1));
      %partial_Z(i,2*i)=1/0.00001*(norminv(gammainc(lambda*xstar(i),k+0.00001),0,1)-norminv(gammainc(lambda*xstar(i),k),0,1));

      
   case 4  % Shifted exponential marginal distribution
      %    z(i)=norminv(1-exp(-lambda*(x(i)-x_zero)),0,1);
      mean = parameter(i,1);
      stdv = parameter(i,2);
      lambda = parameter(i,3);
      x_zero = parameter(i,4);
      
      dzdlambda = (xstar(i)-x_zero)*exp(-lambda*(xstar(i)-x_zero))*1/ferum_pdf(1,inv_norm_cdf(1-exp(-lambda*(xstar(i)-x_zero))),0,1);
      dzdx_zero = -lambda*exp(-lambda*(xstar(i)-x_zero))*1/ferum_pdf(1,inv_norm_cdf(1-exp(-lambda*(xstar(i)-x_zero))),0,1);

		dzdmean = dzdlambda*0 + dzdx_zero*1;
      dzdstdv = dzdlambda*(-1/stdv^2) - dzdx_zero*1;
      
      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzdlambda;
      J_z_theta(i,6*(i-1)+4) = dzdx_zero;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;


		%partial_Z(i,2*i-1)=(xstar(i)-x_zero)*exp(-lambda*(xstar(i)-x_zero))/normpdf(norminv(1-exp(-lambda*(xstar(i)-x_zero)),0,1),0,1);
      %partial_Z(i,2*i)=-lambda*exp(-lambda*(xstar(i)-x_zero))/normpdf(norminv(1-exp(-lambda*(xstar(i)-x_zero)),0,1),0,1);
      
   case 5  % Shifted Rayleigh marginal distribution
      %    z(i)=norminv(1-exp(-0.5*((xstar(i)-x_zero)/a)^2),0,1);
      mean = parameter(i,1);
      stdv = parameter(i,2);
      a = parameter(i,3);
      x_zero = parameter(i,4);
      
      
      dzda = -((xstar(i)-x_zero)^2/a^3)*exp(-0.5*((xstar(i)-x_zero)/a)^2)*1/ferum_pdf(1,inv_norm_cdf(1-exp(-0.5*((xstar(i)-x_zero)/a)^2)),0,1);
      dzdx_zero = -((xstar(i)-x_zero)/a^2)*exp(-0.5*((xstar(i)-x_zero)/a)^2)*1/ferum_pdf(1,inv_norm_cdf(1-exp(-0.5*((xstar(i)-x_zero)/a)^2)),0,1);
      
      dzdmean = dzda*0 + dzdx_zero*1;
      dzdstdv = dzda*(1/(2-pi/2)^0.5) + dzdx_zero*(-(pi/(4-pi))^0.5);
      
      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzda;
      J_z_theta(i,6*(i-1)+4) = dzdx_zero;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;


      
      %partial_Z(i,2*i-1)=-(xstar(i)-x_zero)^2/a^3*exp(-0.5*((xstar(i)-x_zero)/a)^2)/normpdf(norminv(1-exp(-0.5*((xstar(i)-x_zero)/a)^2),0,1),0,1);
      %partial_Z(i,2*i)=-(xstar(i)-x_zero)/a^2*exp(-0.5*((xstar(i)-x_zero)/a)^2)/normpdf(norminv(1-exp(-0.5*((xstar(i)-x_zero)/a)^2),0,1),0,1);
      
   case 6  % Uniform marginal distribution
      %    z(i) = inv_norm_cdf( (x(i) - a) / (b-a) );
      mean = parameter(i,1);
      stdv = parameter(i,2);
      a = parameter(i,3);
      b = parameter(i,4);
      
      dzda = ((xstar(i)-b)/(b-a)^2)*1/ferum_pdf(1,inv_norm_cdf((xstar(i)-a)/(b-a)),0,1);
      dzdb = ((xstar(i)-a)/(b-a)^2)*(-1/ferum_pdf(1,inv_norm_cdf((xstar(i)-a)/(b-a)),0,1));
      
      dzdmean = dzda*1 + dzdb*1;
      dzdstdv = -dzda*3^0.5 + dzdb*3^0.5;

      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzda;
      J_z_theta(i,6*(i-1)+4) = dzdb;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;
      
      %partial_Z(i,2*i-1)=((xstar-b)/(b-a)^2)*1/ferum_pdf(1,inv_norm_cdf((xstar-a)/(b-a)),0,1);
      %partial_Z(i,2*i)=((xstar-a)/(b-a)^2)*(-1/ferum_pdf(1,inv_norm_cdf((xstar-a)/(b-a)),0,1));
      
   case 7  % (Reserved for Beta marginal distribution)
         mean = parameter(i,1);
         stdv = parameter(i,2);
         q = parameter(i,3);
         r = parameter(i,4);
         a = parameter(i,5);
         b = parameter(i,6);
         
        x0 = (xstar(i)-a)/(b-a);
        ccc = inv_norm_cdf(betacdf(x0,q,r));
        
        %ccc = inv_norm_cdf(gammainc(lambda*xstar(i),k));
        h = 1/2000;
        
        x0a = (xstar(i)-a-h)/(b-a-h);
        cca = (betacdf(x0a,q,r) - betacdf(x0,q,r))/h;
        dzda = cca/ferum_pdf(1,ccc,0,1);  
        
        x0b = (xstar(i)-a)/(b + h -a);
        ccb = (betacdf(x0b,q,r) - betacdf(x0,q,r))/h;
        dzdb = ccb/ferum_pdf(1,ccc,0,1);
        
        ccq = (betacdf(x0,q+h,r) - betacdf(x0,q,r))/h;
        dzdq = ccq/ferum_pdf(1,ccc,0,1);
        
        ccr = (betacdf(x0,q,r+h) - betacdf(x0,q,r))/h;
        dzdr = ccr/ferum_pdf(1,ccc,0,1);
        
         dzdmean = 0; % still to be implemented
         dzdstdv = 0; % still to be implemented
                 
         J_z_theta(i,6*(i-1)+1) = dzdmean;
         J_z_theta(i,6*(i-1)+2) = dzdstdv;
         J_z_theta(i,6*(i-1)+3) = dzdq;
         J_z_theta(i,6*(i-1)+4) = dzdr;
         J_z_theta(i,6*(i-1)+5) = dzda;
         J_z_theta(i,6*(i-1)+6) = dzdb;
      
   case 11 % Type I Largest Value marginal distribution
      %    z(i) = inv_norm_cdf(exp(-exp(-a_n*(x(i)-u_n))))
      mean = parameter(i,1);
      stdv = parameter(i,2);
      u_n = parameter(i,3);
      a_n = parameter(i,4);
      
      ccc = exp(-exp(-a_n*(xstar(i)-u_n)));
      
      dzdu_n = -a_n*exp(-a_n*(xstar(i)-u_n))*ccc/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
      dzda_n = (xstar(i)-u_n)*exp(-a_n*(xstar(i)-u_n))*ccc/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
      
      dzdmean = dzdu_n*1 + dzda_n*0;
      dzdstdv = dzdu_n*(-0.5772*(6^0.5)/pi) + dzda_n*(-pi/((6^0.5)*stdv^2));
      
      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzdu_n;
      J_z_theta(i,6*(i-1)+4) = dzda_n;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;
      
   case 12 % Type I Smallest Value marginal distribution
      %    z(i) = inv_norm_cdf(1-exp(-exp(a_1*(x(i)-u_1))))
      mean = parameter(i,1);
      stdv = parameter(i,2);
      u_1 = parameter(i,3);
      a_1 = parameter(i,4);
      
      ccc = 1 - exp(-exp(a_1*(xstar(i)-u_1)));
      
      dzdu_1 = -a_1*exp(a_1*(xstar(i)-u_1))*exp(-exp(a_1*(xstar(i)-u_1)))/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
      dzda_1 = (xstar(i)-u_1)*exp(a_1*(xstar(i)-u_1))*exp(-exp(a_1*(xstar(i)-u_1)))/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
            
      dzdmean = dzdu_1*1 + dzda_1*0;
      dzdstdv = dzdu_1*0.5772*(6^0.5)/pi + dzda_1*(-pi/((6^0.5)*stdv^2));
      
      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzdu_1;
      J_z_theta(i,6*(i-1)+4) = dzda_1;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;
      
   case 13 % Type II Largest Value marginal distribution
      %    z(i) = inv_norm_cdf(exp(-(u_n/x(i))^k));
      mean = parameter(i,1);
      stdv = parameter(i,2);
      
      u_n = parameter(i,3);
      k = parameter(i,4);
            
      ccc = exp(-(u_n/xstar(i))^k);
      
      dzdu_n = (-(u_n/xstar(i))^k)*(k/u_n)*exp(-(u_n/xstar(i))^k)/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
      dzdk = (-(u_n/xstar(i))^k)*log(u_n/xstar(i))*exp(-(u_n/xstar(i))^k)/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
      
      
      r = 1./k;
      
      
      g1 = gamma(1.-r);
      g2 = gamma(1.-r*2);
      
      rp =1./(k + 0.001);
      rm=1./ (k - 0.001);
      
      gp1 = gamma(1.-rp);
      gm1 = gamma(1.-rm);
      
      gp2 = gamma(1.-rp*2.);
      gm2 = gamma(1.-rm*2.);
      
      du1=g1;
      du2=(gp1-gm1)*u_n/0.002;
      ds1 = sqrt(g2-g1*g1);
      ds2=(sqrt(gp2-gp1*gp1)-sqrt(gm2-gm1*gm1))*u_n/0.002;
      
      det=du1*ds2-du2*ds1;
      
      dzdmean = (ds2*dzdu_n-ds1*dzdk)/det;
      dzdstdv = (-du2*dzdu_n +du1*dzdk)/det;

      
      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzdu_n;
      J_z_theta(i,6*(i-1)+4) = dzdk;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;
      
      % partial_Z(i,2*i-1)=(-k/u)*(u/xstar(i))^k * exp(-(u/xstar(i))^k)/normpdf(norminv(exp(-(u/xstar(i))^k),0,1),0,1);
      % partial_Z(i,2*i)=(u/xstar(i))^k * log(xstar(i)/u) * exp(-(u/xstar(i))^k)/normpdf(norminv(exp(-(u/xstar(i))^k),0,1),0,1);

	case 14 % (Reserved for Type III Smallest Value marginal distribution)
   		  % z(i) = inv_norm_cdf(1-exp(-((x(i)-epsilon)/(u_1-epsilon))^k));
      mean = parameter(i,1);
      stdv = parameter(i,2);
      u_1 = parameter(i,3);
      k = parameter(i,4);
      epsilon = parameter(i,5);
      
      ccc = 1-exp(-((xstar(i)-epsilon)/(u_1-epsilon))^k);
      cst1 = ((xstar(i) - epsilon)/(u_1 - epsilon))^k ;
      cst2 = exp(-cst1);
      cst3 = (xstar(i)-epsilon)/(u_1 - epsilon)^2 - 1/(u_1 - epsilon);
      
      dzdu_1 = - cst1*cst2*(k/(u_1-epsilon))/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
      dzdk = cst1*cst2*log((xstar(i) - epsilon)/(u_1 - epsilon))/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
      dzdepsilon = cst1*cst2*k*cst3*(u_1 - epsilon)/((xstar(i)-epsilon)*ferum_pdf(1,inv_norm_cdf(ccc),0,1));
      
      
      %ux=xx/par(1)
      %uxk=ux**par(2)
      %e=uxk*dexp(-uxk)/pz
      %s(3)=-par(2)*e/par(1)
      %s(4)=dlog(ux)*e
      
      test = dzdu_1 + dzdepsilon
      r=1./k;
      
      %g1=dgama(1.+r)
      %g2=dgama(1.+r*2.)
      g1 = gamma(1.+r);
      g2 = gamma(1.+r*2);
      
      
      
      
      %dr=0.001d0*par(2)
      dr = 0.001*k;
      
      %rp=1.d0/(par(2)+dr)
      rp=1./(k+dr);
      %rm=1.d0/(par(2)-dr)
      rm=1./(k-dr);
      
      %gp1=dgama(1.+rp)
      %gm1=dgama(1.+rm)
      %gp2=dgama(1.+rp*2.)
      %gm2=dgama(1.+rm*2.)
      gp1 = gamma(1.+rp);
      gm1 = gamma(1.+rm);
      
      gp2 = gamma(1.+rp*2.);
      gm2 = gamma(1.+rm*2.);
      
      %du1=g1
      du1=g1;
      %du2=(gp1-gm1)*par(1)/2.d0/dr
      du2=(gp1-gm1)*u_1/2./dr;
      %ds1=dsqrt(g2-g1*g1)
      ds1 = sqrt(g2-g1*g1);
      %ds2=(dsqrt(gp2-gp1*gp1)-dsqrt(gm2-gm1*gm1))*par(1)/2.d0/dr
      ds2=(sqrt(gp2-gp1*gp1)-sqrt(gm2-gm1*gm1))*u_1/2./dr;
      
      %det=du1*ds2-du2*ds1
      det=du1*ds2-du2*ds1;
      
      %s(1)=(ds2*s(3)-ds1*s(4))/det
      %s(2)=(-du2*s(3)+du1*s(4))/det
      
      
      
      dzdmean = (ds2*dzdu_1-ds1*dzdk)/det;  
      dzdstdv = (-du2*dzdu_1+du1*dzdk)/det;	
      
      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzdu_1;
      J_z_theta(i,6*(i-1)+4) = dzdk;
      J_z_theta(i,6*(i-1)+5) = dzdepsilon;
      J_z_theta(i,6*(i-1)+6) = 0;
      
   case 8 % (Reserved for Chi-square marginal distribution)
       %    z(i)= inv_norm_cdf(gammainc(0..5*xstar(i),nu/2));
      mean = parameter(i,1);
      stdv = parameter(i,2); 
      nu = parameter(i,3);
      lambda = 0.5;
      %k = nu /2;
      
      ccc = inv_norm_cdf(gammainc(lambda*xstar(i),nu/2));
      dzdlambda = 0;
      %dzdlambda = ((xstar(i)*(lambda*xstar(i))^(k-1))/gamma(k))*exp(-lambda*xstar(i))/ferum_pdf(1,ccc,0,1);
      % Computed dzdk with a forward finite difference scheme
      h = 1/500;
      cc = (gammainc(lambda*xstar(i),(nu + h)/2) - gammainc(lambda*xstar(i),nu/2))/h;
      dzdnu = cc/ferum_pdf(1,ccc,0,1);   
      
      dzdmean = dzdnu*(4*mean/stdv^2);
      dzdstdv = dzdnu*(-4*(mean^2)/stdv^3);
      
      %dzdmean = dzdnu*(4*mean/stdv^2) + dzdlambda*(1/stdv^2) ;
      %dzdstdv = dzdnu*(-4*(mean^2)/stdv^3) + dzdlambda*(-2*mean/stdv^3) ;
      
      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzdnu;
      J_z_theta(i,6*(i-1)+4) = dzdlambda;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;
      
   case 15 % Gumbel marginal distribution
      mean = parameter(i,1);
      stdv = parameter(i,2);
      u_n = parameter(i,3);
      a_n = parameter(i,4);
      
      ccc = exp(-exp(-a_n*(xstar(i)-u_n)));
      
      dzdu_n = -a_n*exp(-a_n*(xstar(i)-u_n))*ccc/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
      dzda_n = (xstar(i)-u_n)*exp(-a_n*(xstar(i)-u_n))*ccc/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
      
      dzdmean = dzdu_n*1 + dzda_n*0;
      dzdstdv = dzdu_n*(-0.5772*(6^0.5)/pi) + dzda_n*(-pi/((6^0.5)*stdv^2));
      
      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzdu_n;
      J_z_theta(i,6*(i-1)+4) = dzda_n;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;
      
   case 16 % (Reserved for Weibull marginal distribution)
      % z(i) = inv_norm_cdf(1-exp(-(x(i)/u_1)^k));
      mean = parameter(i,1);
      stdv = parameter(i,2);
      u_1 = parameter(i,3);
      k = parameter(i,4);
      
      ccc = 1 - exp(-(xstar(i)/u_1)^k);
      
      dzdu_1 = -((xstar(i)/u_1)^k)*(k/u_1)*exp(-(xstar(i)/u_1)^k)/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
      dzdk = ((xstar(i)/u_1)^k)*log(xstar(i)/u_1)*exp(-(xstar(i)/u_1)^k)/ferum_pdf(1,inv_norm_cdf(ccc),0,1);
      
      r=1./k;
      
      g1 = gamma(1.+r);
      g2 = gamma(1.+r*2);
      
      dr = 0.001*k;
      
      rp = 1./(k+dr);
      rm = 1./(k-dr);
      
      gp1 = gamma(1.+rp);
      gm1 = gamma(1.+rm);
      
      gp2 = gamma(1.+rp*2.);
      gm2 = gamma(1.+rm*2.);
      
      du1=g1;
      du2=(gp1-gm1)*u_1/2./dr;

      ds1 = sqrt(g2-g1*g1);
      ds2=(sqrt(gp2-gp1*gp1)-sqrt(gm2-gm1*gm1))*u_1/2./dr;
      
      det=du1*ds2-du2*ds1;
  
      dzdmean = (ds2*dzdu_1-ds1*dzdk)/det;  
      dzdstdv = (-du2*dzdu_1+du1*dzdk)/det;	
      
      J_z_theta(i,6*(i-1)+1) = dzdmean;
      J_z_theta(i,6*(i-1)+2) = dzdstdv;
      J_z_theta(i,6*(i-1)+3) = dzdu_1;
      J_z_theta(i,6*(i-1)+4) = dzdk;
      J_z_theta(i,6*(i-1)+5) = 0;
      J_z_theta(i,6*(i-1)+6) = 0;
      
   case 18 % (Reserved for Laplace marginal distribution)
      
   case 19 % (Reserved for Pareto marginal distribution)
      
      
   otherwise
      
   end
   
end

J_u_theta = iLo * J_z_theta;

sensitivities = (alpha'*J_u_theta)';

for j = 1 : marg_dim
   beta_sensitivities(j,1) = sensitivities(6*(j-1) + 1);
   beta_sensitivities(j,2) = sensitivities(6*(j-1) + 2);
   beta_sensitivities(j,3) = sensitivities(6*(j-1) + 3);
   beta_sensitivities(j,4) = sensitivities(6*(j-1) + 4);
   beta_sensitivities(j,5) = sensitivities(6*(j-1) + 5);
   beta_sensitivities(j,6) = sensitivities(6*(j-1) + 6);
end

pf_sensitivities = -ferum_pdf(1,beta,0,1)*beta_sensitivities;
