function [ d ] = bhatt_weibull( lambda_k_1,lambda_k_2 )
%BHATT_WEIBULL -- Computes the Bhattacharryya inspired distance measure
%                 between two sets of Weibull parameters
%
% Usage
%  d = BHATT_WEIBULL(lambda_k_1,lambda_k_2)
%
% Inputs
%  lambda_k_1, lambda_2: Two Nx2 matrices of Weibull parameters. N is the 
%                        number of distributions to be compared
%                             
% Outputs
%  d:                    Distance measure

l1=lambda_k_1(:,1);
l2=lambda_k_2(:,1);
k1=lambda_k_1(:,2);
k2=lambda_k_2(:,2);
k=(k1+k2)/2; %k1=k; k2=k;
%l=(l1+l2)/2; l1=l; l2=l;
%d=-sum(log(4)+log((l1.^k1).*(l2.^k2).*k1.*k2)/2-log(k1+k2)-log(l1.^k1+l2.^k2));
%d=-sum(log(2)+k.*log(l1.*l2)/2-log(l1.^k+l2.^k));
d=-sum(log(4)+(log(k1)+log(k2)+k.*log(l1)+k.*log(l2))/2-log(l1.^k+l2.^k)-log(k1+k2));
%d=2*((l1.*l2).^(k/2))./(l1.^k+l2.^k);
end

