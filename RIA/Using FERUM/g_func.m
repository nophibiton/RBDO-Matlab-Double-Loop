function val = g_func(x1,x2,d1,d2)

% limit state function
val = (1/5)*d1*d2*x2^2-x1;
end