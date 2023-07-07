function [ursqr]=rsqr_uncentered(data,data_rec)
% [r_sqr]=rsqr_uncentered(data,data_rec)
% This function calculates the uncetered correlation coefficient using "Cluster"
% method.  
% Input:
%       data    matrix of observed data  (e.g., data=[mus pert_dir])
%       data_rec    matrix of reconstructed/predicted data (e.g., data_rec
%       = [mus pert_dir]
%
% Output:
%       ursqr   matrix with uncentered correlation coefficients
%
% called functions:
%       std_mean0.m (optional)
%       
% 
% this function is called by:
%       funur.m
%
% Written by: GTO May 24th, 2006
% Last modified:
%
%

% Shift dimensions for regression purposes data magnitudes in rows and
% data channels in columns
warning off;
data=data';
data_rec=data_rec';

% Like it is done in Cluster 3.0
% dim_data=size(data);
% for i=1:dim_data(2)
%     % Calculate standard deviation assuming mean of 0
%     datastd0=std_mean0(data(:,i));
%     datarecstd0=std_mean0(data_rec(:,i));
% 
%     X=[data(:,i)/datastd0 data_rec(:,i)/datarecstd0];
%     n=length(X);
%     ur(i)=(sum(prod(X,2)))/n;
%     ursqr(i)=ur(i)^2;
% end

% Zar book
dim_data=size(data);
for i=1:dim_data(2)
    
    X=[data(:,i) data_rec(:,i)];
    n=length(X);
    ur(i)=(sum(prod(X,2)))^2/(sum(data(:,i).^2)*sum(data_rec(:,i).^2));
   
%        x=data(:,i)';
%     y=data_rec(:,i)';
%     [r,p]=corrcoef(x,y);
%     if length(x)==1
%         pvalue_m(i)=p(1);
%     else
%     pvalue_m(i)=p(1,2);
%     end
%     if pvalue_m(i)<0.05
%         rsqure_m(i)=r(1,2)^2;
%     else
%         rsqure_m(i)=NaN;
%         ur(i)=NaN;
%     end
    
end

ursqr=ur';
