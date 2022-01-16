function [ k_type2largest ] = k_type2largest(x,cov)

% used for the type II largest value distribution. 

k_type2largest =(gamma(1-2/x)-(gamma(1-1/x))^2)^0.5/(gamma(1-1/x))-cov;