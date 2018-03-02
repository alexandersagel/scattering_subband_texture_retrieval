function [ alpha_beta ] = fwt2alphabeta( Img, decLevel )
%FWT2ALPHABETA -- Performs a Subband Decomposition and calculates the GGD
%parameters of for each subband
%
% Usage:
%  alpha_beta = FWT2ALPHABETA(Img, decLevel)
%
% Inputs
%  Img:       Input image
%  decLevel:  Depth of the Subband Decomposition
%
% Outputs
% alpha_beta: (decLevel*3)x2 matrix of the parameters alpha and beta for
%             each subband

assert(size(Img,1)==size(Img,2), 'Assuming quadratic images here')
nSide=size(Img,1);
if nargin>1 && decLevel>0
    h0=[.48296291314483155 .83651630373770821 .224143868041921833 -.12940952255095486]; %Daubechies D2
    %h0=[0.47046721 1.14111692 0.650365 -0.19093442 -0.12083221 0.0498175]; h0=h0/norm(h0);
    Img = subbandDecomposition(Img, h0,h0(end:-1:1),decLevel);
else
    decLevel=round(log2(nSide/16));
end
alpha_beta=zeros(decLevel*3,2);
abi=0;
for k=decLevel:-1:1
    X=Img(1:nSide/(2^k),nSide/(2^k)+1:nSide/(2^(k-1))); X=X(:);
    abi=abi+1; L=length(X);
    alpha_beta(abi,2)=ggd_beta(X); alpha_beta(abi,1)=norm(X,alpha_beta(abi,2))*(alpha_beta(abi,2)/L)^(1/alpha_beta(abi,2));
    X=Img(nSide/(2^k)+1:nSide/(2^(k-1)),1:nSide/(2^k)); X=X(:);
    abi=abi+1;
    alpha_beta(abi,2)=ggd_beta(X); alpha_beta(abi,1)=norm(X,alpha_beta(abi,2))*(alpha_beta(abi,2)/L)^(1/alpha_beta(abi,2));
    X=Img(nSide/(2^k)+1:nSide/(2^(k-1)),nSide/(2^k)+1:nSide/(2^(k-1))); X=X(:);
    abi=abi+1;
    alpha_beta(abi,2)=ggd_beta(X); alpha_beta(abi,1)=norm(X,alpha_beta(abi,2))*(alpha_beta(abi,2)/L)^(1/alpha_beta(abi,2));
end

end

