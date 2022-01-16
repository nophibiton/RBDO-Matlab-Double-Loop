function [ parameter ]  = distribution_parameter(marg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the parameters (resp. the mean and stdv) corresponding to the marginal distributions knowing %
% the mean and std. dev. ( the parameters)of the distribution. 
%
% The distribution library is available in Ferum's new user's guide. 
%
% It also computes the jacobian of the standard parameters (mean and sd) wrt these other parameters, in order to be used 
% in the sensitivity analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

marg_dim = size(marg,1);

J_teta_tetaprime=zeros(2*marg_dim,2*marg_dim);
parameter = zeros(marg_dim,6);

for i=1 : marg_dim
   switch marg(i,1)
   case 1 % Normal distribution
      if marg(i,9)== 1
      	mean = marg(i,2);
      	stdv = marg(i,3);
      elseif marg(i,9)== 0
         mean = marg(i,2);
         stdv = marg(i,3);
      end
            
      parameter(i,1) = mean;
      parameter(i,2) = stdv;
      parameter(i,3) = mean;
      parameter(i,4) = stdv;
            
   case 2 % Lognormal distribution
      if marg(i,9)== 1
         lambda = marg(i,5);
         zeta = marg(i,6);
         mean = exp(lambda + (zeta^2)/2);
         stdv = exp(lambda + (zeta^2)/2)*(exp(zeta^2)-1)^0.5;
      elseif marg(i,9)== 0  
      	mean = marg(i,2);
      	stdv = marg(i,3);
      	cov = marg(i,3) / marg(i,2);
      	zeta =(log(1+cov^2))^0.5;
      	lambda = log(marg(i,2))-0.5*zeta^2;
      end
            
      parameter(i,1) = mean;
      parameter(i,2) = stdv;
      parameter(i,3) = lambda;
      parameter(i,4) = zeta;
            
      %J_teta_tetaprime(2*i-1,2*i-1)=marg(i,2);
      %J_teta_tetaprime(2*i-1,2*i)=marg(i,2)*zeta;
      %J_teta_tetaprime(2*i,2*i-1)=marg(i,3);
      %J_teta_tetaprime(2*i,2*i)=marg(i,3)*zeta+marg(i,2)^2*zeta*exp(zeta^2)/marg(i,3);
      
   case 3 % Gamma marginal distribution
      if marg(i,9)== 1
         lambda = marg(i,5);
         k = marg(i,6);
         mean = k/lambda;
         stdv = (k^0.5)/lambda;
      elseif marg(i,9)== 0  
      	mean = marg(i,2);
      	stdv = marg(i,3);
      	lambda = marg(i,2) / marg(i,3)^2 ;
      	k = marg(i,2)^2 / marg(i,3)^2;
      end
          
      parameter(i,1)= mean;
      parameter(i,2)= stdv;
      parameter(i,3)= lambda;
      parameter(i,4)= k;
      
      %J_teta_tetaprime(2*i-1,2*i-1)= -k / lambda^2 ;
      %J_teta_tetaprime(2*i-1,2*i)= 1 / lambda ;
      %J_teta_tetaprime(2*i,2*i-1)= -k^0.5 / lambda^2 ;
      %J_teta_tetaprime(2*i,2*i)= 1 / (2*lambda*k^0.5) ;
        
   case 4 % Shifted Exponential marginal distribution      
      if marg(i,9)== 1
         lambda = marg(i,5);
         x_zero = marg(i,6);
         mean = x_zero + 1/lambda;
         stdv = 1/lambda;
      elseif marg(i,9)== 0
         mean = marg(i,2);
      	stdv = marg(i,3);
      	x_zero= marg(i,2) - marg(i,3);
      	lambda= 1 / marg(i,3) ;
      end
           
      parameter(i,1) = marg(i,2);
      parameter(i,2) = marg(i,3);
      parameter(i,3) = lambda;
      parameter(i,4) = x_zero;
      
      %J_teta_tetaprime(2*i-1,2*i-1)= -1 / lambda^2 ;
      %J_teta_tetaprime(2*i-1,2*i)= 1 ;
      %J_teta_tetaprime(2*i,2*i-1)= -1 / lambda^2 ;
      %J_teta_tetaprime(2*i,2*i)= 0 ;
      
   case 5 % Shifted Rayleigh marginal distribution
      if marg(i,9)== 1
         a = marg(i,5);
         x_zero = marg(i,6);
		 mean = x_zero + a*(pi/2)^0.5;
      	 stdv = a*(2-pi/2)^0.5;          
      elseif marg(i,9)== 0
         mean = marg(i,2);
      	 stdv = marg(i,3);      
      	 a = marg(i,3) / (2-pi/2)^0.5 ;
      	 x_zero = marg(i,2) - (pi/(4-pi))^0.5 * marg(i,3);
      end
     
      parameter(i,1) = mean;
      parameter(i,2) = stdv;
      parameter(i,3) = a;
      parameter(i,4) = x_zero;
      
      %J_teta_tetaprime(2*i-1,2*i-1)= (pi/2)^0.5 ;
      %J_teta_tetaprime(2*i-1,2*i)= 1  ;
      %J_teta_tetaprime(2*i,2*i-1)= (2-pi/2)^0.5 ;
      %J_teta_tetaprime(2*i,2*i)= 0 ;
      
      
   case 6 % Uniform marginal distribution
      if marg(i,9)== 1
         a = marg(i,5);
         b = marg(i,6);
			mean = (a+b)/2 ;
      	stdv = (b-a)/(2*(3)^0.5);          
      elseif marg(i,9)== 0
         mean = marg(i,2);
      	stdv = marg(i,3);
      	a = marg(i,2) - 3^0.5 * marg(i,3);
      	b = marg(i,2) + 3^0.5 * marg(i,3);
      end
           
      parameter(i,1) = mean;
      parameter(i,2) = stdv;
      parameter(i,3)= a;
      parameter(i,4)= b;
           
      %J_teta_tetaprime(2*i-1,2*i-1)=0.5;
      %J_teta_tetaprime(2*i-1,2*i)=0.5;
      %J_teta_tetaprime(2*i,2*i-1)=-1/2/3^0.5;
      %J_teta_tetaprime(2*i,2*i)= 1/2/3^0.5;
      
   case 7 % (Reserved for Beta marginal distribution)
      if marg(i,9)== 1
         q = marg(i,5);
         r = marg(i,6);
         a = marg(i,7);
         b = marg(i,8);
		   mean = a + q*(b-a)/(q+r);
         stdv = ((b-a)/(q+r))*(q*r/(q+r+1))^0.5;          
      elseif marg(i,9)== 0
         mean = marg(i,2);
      	stdv = marg(i,3);
         a = marg(i,7);
         b = marg(i,8);
         % solve non linear set of equations for shape parameters m,n (collected in the vector par)
         par_guess = [1 1];
         par = fsolve('betaparam',par_guess,optimset('fsolve'),a,b,mean,stdv);
         r = par(1);
         q = par(2);
              	
      end
           
      parameter(i,1) = mean;
      parameter(i,2) = stdv;
      parameter(i,3) = q;
      parameter(i,4) = r;
      parameter(i,5) = a;
      parameter(i,6) = b;
      
   case 8 % Chi-square Value marginal distribution
      if marg(i,9)== 1
         lambda = 0.5;
         nu = marg(i,5);
         mean = nu/(2*lambda);
         stdv = ((nu/2)^0.5)/lambda;
      elseif marg(i,9)== 0  
      	mean = marg(i,2);
      	stdv = marg(i,3);
        lambda = 0.5;
        mean_test = lambda*stdv^2
        if mean/mean_test < 0.95 | mean/mean_test > 1.05
            error('ERROR WHEN USING CHI-SQUARE DISTRIBUTION MEAN AND ST.DEV SHOULD BE GIVEN SUCH THAT MEAN = 0.5*ST.DEV.^2')
        end
      	nu = 2*(marg(i,2)^2 / marg(i,3)^2);
      end
     
      parameter(i,1)= mean;
      parameter(i,2)= stdv;
      parameter(i,3)= nu;
      %parameter(i,4)= k;
      
      %J_teta_tetaprime(2*i-1,2*i-1)= -k / lambda^2 ;
      %J_teta_tetaprime(2*i-1,2*i)= 1 / lambda ;
      %J_teta_tetaprime(2*i,2*i-1)= -k^0.5 / lambda^2 ;
      %J_teta_tetaprime(2*i,2*i)= 1 / (2*lambda*k^0.5) ;
      
   case 11 % Type I Largest Value marginal distribution
      if marg(i,9)== 1
         u_n = marg(i,5);
         a_n = marg(i,6);
			mean = u_n + 0.5772/a_n ;
      	stdv = pi/(a_n*6^0.5);          
      elseif marg(i,9)== 0
         mean = marg(i,2);
      	stdv = marg(i,3);
      	a_n = pi/(marg(i,3)*sqrt(6));
      	u_n = marg(i,2) - (0.5772*marg(i,3)*sqrt(6))/pi;
      	%p3 = ;
      	%P4 = ;
      end

         
      parameter(i,1) = mean;
      parameter(i,2) = stdv;
      parameter(i,3) = u_n;
      parameter(i,4) = a_n;
      %parameter(i,5) = p3;
      %parameter(i,6) = P4;
      
      
   case 12 % Type I Smallest Value marginal distribution
      if marg(i,9)== 1
         u_1 = marg(i,5);
         a_1 = marg(i,6);
			mean = u_1 - 0.5772/a_1 ;
      	stdv = pi/(a_1*6^0.5);          
      elseif marg(i,9)== 0
         mean = marg(i,2);
      	stdv = marg(i,3);
      	a_1 = pi/(marg(i,3)*sqrt(6));
      	u_1 = marg(i,2) + (0.5772*marg(i,3)*sqrt(6))/pi;
      	%p3 = ;
      	%P4 = ;
      end
 
      parameter(i,1) = mean;
      parameter(i,2) = stdv;
      parameter(i,3) = u_1;
      parameter(i,4) = a_1;
      %parameter(i,5) = p3;
      %parameter(i,6) = P4;
      
   case 13 % Type II Largest Value marginal distribution
      if marg(i,9)== 1
         u_n = marg(i,5);
         k = marg(i,6);
		 mean = u_n*gamma(1-1/k); 
      	 stdv = u_n*(gamma(1-2/k)-(gamma(1-1/k))^2)^0.5;          
      elseif marg(i,9)== 0
         mean = marg(i,2);
      	 stdv = marg(i,3);
                           
         %par = fsolve('typIIlargestpar',par_guess,optimset('fsolve'),mean,stdv);
         mean_test = mean  ;
         par_guess = [2.000001 10e+07];
         par = fzero('typIIlargestpar',par_guess,optimset('fzero'),mean,stdv);
         k = par;
         u_n = mean/gamma(1-1/k);
         %u_n = par(1);
         %k = par(2);
      	
      end

      parameter(i,1) = mean;
      parameter(i,2) = stdv;
      parameter(i,3) = u_n;
      parameter(i,4) = k;

      %J_teta_tetaprime(2*i-1,2*i-1)= gamma(1-1/k) ;
      %J_teta_tetaprime(2*i-1,2*i)= u/0.001*(gamma(1-1/(k+0.001))-gamma(1-1/k))  ;
      %J_teta_tetaprime(2*i,2*i-1)= (gamma(1-2/k)-(gamma(1-1/k))^2)^0.5 ;
      %J_teta_tetaprime(2*i,2*i)= u/0.001*((gamma(1-2/(k+0.001))-(gamma(1-1/(k+0.001)))^2)^0.5-(gamma(1-2/k)-(gamma(1-1/k))^2)^0.5) ;
      
    case 14 % Type III Smallest Value marginal distribution
       if marg(i,9)== 1
         u_1 = marg(i,5);
         k = marg(i,6);
         epsilon = marg(i,7);
         mean = epsilon + (u_1-epsilon)*gamma(1+1/k); 
      	 stdv = (u_1-epsilon)*(gamma(1+2/k)-gamma(1+1/k)^2)^0.5;        
         
      elseif marg(i,9)== 0
         mean = marg(i,2);
      	 stdv = marg(i,3);
         epsilon = marg(i,7);
         
         meaneps = mean - epsilon;
         
         par_guess = [0.1 10e+07];
         par = fzero('typIIIsmallestpar',par_guess,optimset('fzero'),meaneps,stdv);
         k = par;
         u_1 = meaneps/gamma(1+1/k)+epsilon;
         %par = fsolve('typIIIsmallestpar',par_guess,optimset('fsolve'),mean,stdv,epsilon);
         %u_1 = par(1);
         %k = par(2);

         
      end

      parameter(i,1) = mean;
      parameter(i,2) = stdv;
      parameter(i,3) = u_1;
      parameter(i,4) = k;
      parameter(i,5) = epsilon;

    
         
    case 15 % Gumbel marginal distribution
      if marg(i,9)== 1
         a_n = marg(i,5);
         u_n = marg(i,6);
		 mean = u_n + 0.5772/a_n ;
      	 stdv = pi/(a_n*6^0.5);          
      elseif marg(i,9)== 0
         mean = marg(i,2);
      	 stdv = marg(i,3);
      	 a_n = pi/(marg(i,3)*sqrt(6));
      	 u_n = marg(i,2) - (0.5772*marg(i,3)*sqrt(6))/pi;
      	 
      end

         
      parameter(i,1) = mean;
      parameter(i,2) = stdv;
      parameter(i,3) = u_n;
      parameter(i,4) = a_n;
      %parameter(i,5) = p3;
      %parameter(i,6) = P4;

         
    case 16 % Weibull marginal distribution
       if marg(i,9)== 1
         u_1 = marg(i,5);
         k = marg(i,6);
         mean = u_1*gamma(1+1/k);
      	stdv = u_1*(gamma(1+2/k)-gamma(1+1/k)^2)^0.5;          
      elseif marg(i,9)== 0
         mean = marg(i,2);
         stdv = marg(i,3);
         epsilon = 0;
      	 meaneps = mean - epsilon;
         
         par_guess = [0.1 10e+07];
         par = fzero('typIIIsmallestpar',par_guess,optimset('fzero'),meaneps,stdv);
         k = par;
         u_1 = meaneps/gamma(1+1/k)+epsilon;
         
              	
      end

      parameter(i,1) = mean;
      parameter(i,2) = stdv;
      parameter(i,3) = u_1;
      parameter(i,4) = k;
           
    case 18 % (Reserved for Laplace marginal distribution)
    case 19 % (Reserved for Pareto marginal distribution)
      
      
   otherwise   
      
      parameter(i,1)=marg(i,2);
      parameter(i,2)=marg(i,3);
      %J_teta_tetaprime(2*i-1,2*i-1)=1;
      %J_teta_tetaprime(2*i-1,2*i)=0;
      %J_teta_tetaprime(2*i,2*i-1)=0;
      %J_teta_tetaprime(2*i,2*i)=1;
      
   end
end

