function [ zero_gamma ] = zero_gamma(x,k,lambda,normal_val)

% used for the gamma distribution. 
% Notice: the definition of gammainc in matlab


zero_gamma = gammainc(lambda*x,k)- normal_val;
             