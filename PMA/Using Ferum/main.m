clc,clear; format compact; format shortG;

fun = @(x) x(1)^2+x(2)^2;

lb = [0,0];
ub = [15,15];

A = [];
b = [];
Aeq = [];
beq = [];

x0 = [7.5,7.5];

nonlcon = @calc_const;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)