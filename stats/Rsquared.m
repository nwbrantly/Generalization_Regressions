function [rsq]= Rsquared(data,data_hat)

% [rsq]= (data,data_hat)
% This function calculates thecorrelation coefficient in a stride by stride
% base 
%  R^2= 1 - sum(data - data_hat)^2./sum(data - mean(data))^2

% Input:
%       data    matrix of observed data 
%       data_hat    matrix of reconstructed/predicted data
%
% Output:
%       ursqr   r^2 vallue 
% Written by: DMMO  April 22, 2022
% Last modified:

rsq= 1 - sum((data- data_hat).^2)./sum((data- mean(data)).^2);

end 

