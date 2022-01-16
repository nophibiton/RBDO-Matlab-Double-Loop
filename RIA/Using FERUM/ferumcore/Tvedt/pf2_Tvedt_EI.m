function pf2_Tvedt_EI = pf2_Tvedt_EI(u,beta,kapa)

i_complex = i;
h0 = 1;
h = 1;
nu = i_complex;

p = 0;

j = 1;
while abs(tvedt_parabol_integrand(u + j*nu*h,beta,kapa)) >= 10e-12
%for j = 1:12
%abs(tvedt_parabol_integrand(u + j*nu*h,beta,kapa))
p = p + tvedt_parabol_integrand(u + j*nu*h,beta,kapa);
j = j + 1;
end


pf2_Tvedt_EI = real(nu*h*(tvedt_parabol_integrand(u,beta,kapa)/2 + p));