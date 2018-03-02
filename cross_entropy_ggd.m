function [ d ] = cross_entropy_ggd( alpha_beta_1, alpha_beta_2 )
%CROSS_ENTROPY_GGD -- Computes the cross entropy distance measure between 
%                     two sets of GGD parameters
%
% Usage
%  d = CROSS_ENTROPY_GGD(alpha_beta_1, alpha_beta_2)
%
% Inputs
%  alpha_beta_1, alpha_beta_2: Two Nx2 matrices of GGF parameters. N is the
%                              number of distributions to be compared
%                             
% Outputs
%  d:                          Distance measure         
assert(size(alpha_beta_1,2)==size(alpha_beta_2,2),'alpha-beta-Vectors must have same length');
d=0;
alpha_1=alpha_beta_1(:,1); alpha_2=alpha_beta_2(:,1);
beta_1=alpha_beta_1(:,2); beta_2=alpha_beta_2(:,2); one=ones(size(alpha_1));
d=log(2*alpha_2.*gamma(one./beta_2)./(beta_2)) ... 
    +((alpha_1./alpha_2).^beta_2).*gamma((beta_2+one)./beta_1)./gamma(one./beta_1);
d=sum(d);


end

