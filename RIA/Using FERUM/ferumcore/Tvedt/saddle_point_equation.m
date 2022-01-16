function dpsidu = saddle_point_equation(u,beta,kapa)


m = - beta;

size_kapa = size(kapa);

if size_kapa(1) == 1
    dpsidu = -m + u - 1/u + sum((kapa./2)./(1-kapa*u));
elseif size_kapa(1)== 2
    a_curvatures_plus = kapa(1,:);
    a_curvatures_minus = kapa(2,:);
        
    %num = a_curvatures_plus./(2*(1-a_curvatures_plus*u).^1.5) + a_curvatures_minus./(2*(1-a_curvatures_minus*u).^1.5);
    %denom = 1./(1-a_curvatures_plus*u).^0.5 + 1./(1-a_curvatures_minus*u).^0.5;
    sum((a_curvatures_plus./(2*(1-a_curvatures_plus*u).^1.5) + a_curvatures_minus./(2*(1-a_curvatures_minus*u).^1.5))/(1./(1-a_curvatures_plus*u).^0.5 + 1./(1-a_curvatures_minus*u).^0.5) )    
    f = -m + u - 1/u + sum((a_curvatures_plus./(2*(1-a_curvatures_plus*u).^1.5) + a_curvatures_minus./(2*(1-a_curvatures_minus*u).^1.5))/(1./(1-a_curvatures_plus*u).^0.5 + 1./(1-a_curvatures_minus*u).^0.5) );
else
    ERROR('ERROR IN THE INPUT CURVATURES')
end






% Psi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sum = 0;
%for i=1:kapa_dim
%    sum = sum + (log(1-k(j)*u))^(-0.5)
%end

%psi = -m*u + 0.5*u*u + sum((log(1-kapa.*u)).^(-0.5)) -log(u);

% dpsi/du
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sum1 = 0;
%for i=1:kapa_dim
%    sum1 = sum1 + (k(j)/2)/(1-k(j)*u)
%end




% d2psi/du
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sum2 = 0;
%for i=1:kapa_dim
%    sum2 = sum2 + 2*(k(j)^2)*(2-2k(j)*u)^(-2)
%end


%d2psidu2 = 1 + 1/u^2 + sum(2*(kapa.^2).*(2-2.*kapa.*u).^(-2));