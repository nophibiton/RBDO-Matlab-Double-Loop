function f = tvedt_parabol_integrand(u,beta,kapa)


i_complex = i;
m = -beta;

size_kapa = size(kapa);

if size_kapa(1) == 1
    f = (-1/(i_complex*pi*u))*exp(- m*u)*exp(u^2/2)*prod(1./(1-kapa.*u).^0.5);
elseif size_kapa(1)== 2
    a_curvatures_plus = kapa(1,:);
    a_curvatures_minus = kapa(2,:);
    f = (-1/(i_complex*pi*u))*exp(- m*u)*exp(u^2/2)*prod( 1/2 * ( (1./(1-a_curvatures_plus*u).^0.5) + (1./(1-a_curvatures_minus*u).^0.5)) );
else
    ERROR('ERROR IN THE INPUT CURVATURES')
end
