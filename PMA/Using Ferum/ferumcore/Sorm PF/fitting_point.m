function [u_prime_i , G_u] = fitting_point(Ro,Lo,iLo,k,R1,U_prime_1,num,lsf,formresults,probdata,analysisopt,gfundata,femodel,randomfield)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the fitting points :                                     %
% - the parameter k which define the path is needed                               %
% - the matrix U_prime_1 defining the starting points on each axis is also needed %
% - R1 is the rotation matrix                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


threshold = 10^-6;
stop_flag = 0;

%Useful available parameters after Form analyis : formresults 
beta = formresults.beta1;
alpha = formresults.alpha;
% Probdata parameters
marg = probdata.marg;
% These parameters are necessary to compute G and Grad_G
grad_flag = lower(analysisopt.grad_flag);
%Determine number of random variables
marg_dim = size(marg,1);
nrv = size(marg,1);
   
% Determine the side of the axis the fitting point belongs to with respect to the matrix U_prime_1 construction
i = num;

if num <= nrv-1      % Negative axis
   counter = i;
   sign = -1;
else                 % Positive axis
   counter = i-nrv+1;
   sign = 1;
end


vect=zeros(nrv,1);             % Useful to compute the scalar product in the computation of b
u_prime_i = U_prime_1(:,i);    % Define the coordinates of the fitting point in the rotated space
a = u_prime_i(counter);
b = u_prime_i(nrv);
         
u = R1'*u_prime_i;				 % Rotation to the standard normal space
x = u_to_x(u,probdata,Lo);
J_u_x = jacobian(x,u,probdata,Lo,iLo);
J_x_u = inv(J_u_x);
   
[ G, grad_g ] = gfun(lsf,x,grad_flag,probdata,gfundata,femodel,randomfield);
grad_G = R1*(grad_g * J_x_u)';

% Case where G is negative at the starting points
if G < 0															
    a = sign*k*beta;
    dadb = 0;
      
    vect(counter) = dadb;
   	vect(nrv) = 1;

    while stop_flag == 0;  
    %for j = 1:4
      b = b - G/((grad_G)'*vect);
      a = sign*k*beta;
      
      % new point coordinates
      u_prime_i(counter)= a;
	  u_prime_i(nrv)= b;
   	        
      u = R1'*u_prime_i;
   	  x = u_to_x(u,probdata,Lo);
   	  J_u_x = jacobian(x,u,probdata,Lo,iLo);
   	  J_x_u = inv(J_u_x);
   
   	  [ G, grad_g ] = gfun(lsf,x,grad_flag,probdata,gfundata,femodel,randomfield);
   	  grad_G = R1*(grad_g * J_x_u)';
      
      
      if abs(G) < threshold
        stop_flag = 1;
      end
      
   	  %dadb = 0;
   	  %vect(i) = dadb;
   	  %vect(nrv) = 1;
      
   	  end
      u_prime_i;
     G_u = G;  

% Case where G is positive at the starting points      
else  % G is positive
          
      a = sign * ((k*beta)^2-(b-beta)^2)^0.5;
      dadb = sign * (-(b-beta)/((k*beta)^2-(b-beta)^2)^0.5);
      
      vect(counter) = dadb;
   	  vect(nrv) = 1;

    while stop_flag == 0;
    %  for j = 1:4
   	
   	   b = b - G/((grad_G)'*vect);
       a = sign * ((k*beta)^2-(b-beta)^2)^0.5;
      
       % New point coordinates
   	   u_prime_i(counter)= a;
	   u_prime_i(nrv)= b;
   	
       u = R1'*u_prime_i;
       x = u_to_x(u,probdata,Lo);
   	   J_u_x = jacobian(x,u,probdata,Lo,iLo);
   	   J_x_u = inv(J_u_x);
   
   	   [ G, grad_g ] = gfun(lsf,x,grad_flag,probdata,gfundata,femodel,randomfield);
   	   grad_G = R1*(grad_g * J_x_u)';
    
       if abs(G) < threshold
           stop_flag = 1;
       end
   
   	   dadb = sign * (-(b-beta)/((k*beta)^2-(b-beta)^2)^0.5);
   	   vect(counter) = dadb;
       	      
   	end
	u_prime_i;
	G_u = G;
end
