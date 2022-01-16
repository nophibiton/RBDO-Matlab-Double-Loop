function [c,ceq] = calc_const(d)

%inequality nonlinear constraint
ferum_inputfile_inverse;
c = -inverse_formresults.theta;

% equality nonlinear constraint
ceq = [];
end