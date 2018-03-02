function [ S ] = normalize_scat( S,mu )
%NORMALIZE_SCAT -- Converts the regular Scattering Transform of a signal to
%                  its Normalized Scattering Transform
%
% Usage
%  Sn = NORMALIZE_SCAT(S)
%
% Inputs
%  S:  A scattering cell array as produced by scat
%
% Outputs
%  Sn: S, normalized

assert(length(S)>2)

p2_cmp=S{3}.meta.theta(1,:)*1000+S{3}.meta.j(1,:)*100+S{3}.meta.theta(2,:)*10+S{3}.meta.j(2,:);
if length(S)>3
    for k=1:length(S{4}.signal)
        Theta2=S{4}.meta.theta(2,k);
        J2=S{4}.meta.j(2,k);
        Theta1=S{4}.meta.theta(3,k);
        J1=S{4}.meta.j(3,k);
        p2=Theta2*1000+J2*100+Theta1*10+J1;
        S{4}.signal{k}=S{4}.signal{k}./(S{3}.signal{p2_cmp==p2}+eps);
    end
end

p_cmp=S{2}.meta.theta(1,:)*10+S{2}.meta.j(1,:);
for k=1:length(S{3}.signal)
    Theta=S{3}.meta.theta(2,k);
    J=S{3}.meta.j(2,k);
    p=Theta*10+J;
    S{3}.signal{k}=S{3}.signal{k}./(S{2}.signal{p_cmp==p}+eps);
end

for k=1:length(S{2}.signal)
    S{2}.signal{k}=S{2}.signal{k}/mu;
end