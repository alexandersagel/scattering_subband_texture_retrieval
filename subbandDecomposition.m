function [ I_dec ] = subbandDecomposition(I, h0,g0,n_layers )
%SUBBANDDECOMPOSITION -- Decomposes an image into a multiresolution
%                        representation
%
% Usage
%  I_dec = subbandDecomposition(I, h0,g0,n_layers)
% 
% Inputs
%  I:        Image to be decomposed
%  h0:       Lowpass filter
%  g0:       Complementary filter to h0
%  n_layers: Decomposition depth
% 
% Outputs:
%  I_dec:    Decomposed Image, contains n_layers*3+1 subbands

assert(mod(size(I,1),2^(n_layers-1))==0,'Dimensions wrong');
assert(mod(size(I,2),2^(n_layers-1))==0,'Dimensions wrong');
assert(length(g0)==length(h0),'Filters must have same length');


h0=h0(:)'; g0=g0(:)';
%g1=h0.*cos((0:length(h0)-1)*pi);
h1=g0.*cos((1:length(g0))*pi);

n=length(h0)-1;
I=[I(floor(n/2):-1:1,floor(n/2):-1:1) I(floor(n/2):-1:1,:) I(floor(n/2):-1:1,end:-1:end-ceil(n/2)+1); ...
   I(:, floor(n/2):-1:1) I I(:,end:-1:end-ceil(n/2)+1); ...
   I(end:-1:end-ceil(n/2)+1,floor(n/2):-1:1) I(end:-1:end-ceil(n/2)+1,:)  I(end:-1:end-ceil(n/2)+1,end:-1:end-ceil(n/2)+1) ];


HH=conv2(h1,h1,I,'valid'); HH=HH(1:2:end,1:2:end);
LH=conv2(h0,h1,I,'valid'); LH=LH(1:2:end,1:2:end);
HL=conv2(h1,h0,I,'valid'); HL=HL(1:2:end,1:2:end);
LL=conv2(h0,h0,I,'valid'); LL=LL(1:2:end,1:2:end);

% HH=imfilter(I,h1,'symmetric'); HH=HH(:,1:2:end); HH=imfilter(HH,h1','symmetric'); HH=HH(1:2:end,:);
% LH=imfilter(I,h0,'symmetric'); LH=LH(:,1:2:end); LH=imfilter(LH,h1','symmetric'); LH=LH(1:2:end,:);
% HL=imfilter(I,h1,'symmetric'); HL=HL(:,1:2:end); HL=imfilter(HL,h0','symmetric'); HL=HL(1:2:end,:);
% LL=imfilter(I,h0,'symmetric'); LL=LL(:,1:2:end); LL=imfilter(LL,h0','symmetric'); LL=LL(1:2:end,:);
if n_layers<2
else
    LL=subbandDecomposition(LL, h0,g0,n_layers-1);
end

I_dec=[LL HL; LH, HH];
end

