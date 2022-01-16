function [c,ceq] = calc_const(d)

%inequality nonlinear constraint
ferum_inputfile;
c = formresults.pf1 - 1/100;

% equality nonlinear constraint
ceq = [];
end