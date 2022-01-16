function us = saddle_point(beta,kapa)


size_kapa = size(kapa);

if size_kapa(1) == 1
    vect = sort(kapa);
    counter = 1;
    counter1 = 1;

    for j=1:length(vect)
        if sign(vect(j))== -1
            minus_kapa(counter)=vect(j);
            counter = counter + 1 ;
        else
            plus_kapa(counter1)=vect(j);
            counter1 = counter1 + 1 ;
        end
    end
    
    
    if exist('minus_kapa') == 0
        lower_bound = -10^15;
    else
        minus_kapa = 1./minus_kapa;
        lower_bound = max(minus_kapa)+ 10^(-6);
    end

    upper_bound = -10^(-15);  
    
 guess = [lower_bound, upper_bound];

 us  = fzero('saddle_point_equation',guess,optimset('fzero'),beta,kapa);  
    
elseif size_kapa(1)== 2
    
    vect = [kapa(1,:), kapa(2,:)];
    vect = sort(vect);
    
    counter = 1;
    counter1 = 1;

    for j=1:length(vect)
        if sign(vect(j))== -1
            minus_kapa(counter)=vect(j);
            counter = counter + 1 ;
        else
            plus_kapa(counter1)=vect(j);
            counter1 = counter1 + 1 ;
        end
    end
    
     if exist('minus_kapa') == 0
        lower_bound = -10^15;
    else
        minus_kapa = 1./minus_kapa;
        lower_bound = max(minus_kapa)+ 10^(-6);
    end

    upper_bound = -10^(-15);
    
    guess = [lower_bound, upper_bound];
    
    kapa_plus = kapa(1,:);
    kapa_minus = kapa(2,:);

    us  = fzero('saddle_point_equation_point_fitting',guess,optimset('fzero'),beta,kapa_plus,kapa_minus);  
end



