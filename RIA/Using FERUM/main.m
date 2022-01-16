clc,clear; format compact; format shortG;

fun = @calc_obj;       % objective function
nonlcon = @calc_const; % constraint function

lb = [0,0];
ub = [15,15];

A = [];
b = [];
Aeq = [];
beq = [];

x0 = [7.5,7.5];   % intial guess

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)