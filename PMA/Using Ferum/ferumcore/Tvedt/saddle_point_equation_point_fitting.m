function dpsidu = saddle_point_equation_point_fitting(u,beta,kapa_plus,kapa_minus)


m = - beta;
 
    
dpsidu = -m + u - 1/u + sum((kapa_plus./(2*(1-kapa_plus*u).^1.5) + kapa_minus./(2*(1-kapa_minus*u).^1.5))/(1./(1-kapa_plus*u).^0.5 + 1./(1-kapa_minus*u).^0.5) );






