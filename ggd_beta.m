function [ beta ] = ggd_beta( X )
%GGD_BETA -- Computes the beta parameter of a GGD from a set of Samples
%
% Usage
%  beta = GGD_BETA(X)
%
% Inputs
%  X:  Set of realizations of a random variable
%
% Outputs
%  beta
     X=X(:);
     m1=mean(abs(X));
     m2=mean(abs(X).^2);
     beta_small=0.1;
     beta_large=6;
     beta=m1/sqrt(m2);
     if beta>0.3
         beta=0.16/(0.88-beta);
     end
     if beta<0 || beta>4
         beta=4;
     end
     gX=g(X,beta); k=0;
     while abs(gX)>1e-13
         beta_old=beta;
         beta=beta-gX/gTick(X,beta);
         %if beta out of scale, fall back to rough guess
         if beta < beta_small
             beta=beta_small;
             break;
         elseif beta > beta_large
             beta=beta_large;
             break;
         end
         if k>1000 || isnan(beta)
             beta=beta_old;
             break;
         end
         gX=g(X,beta); k=k+1;
     end
end

function [ result ] = g(X, beta)
     X=X(:);
     result=1+psi(1/beta)/beta-sum(abs(X).^beta.*log(abs(X)))/sum(abs(X).^beta) ...
         +log(beta*mean(abs(X).^beta))/beta;
end

function [ result ] = gTick(X, beta)
     X=X(:);
     result=-psi(1/beta)/beta^2-psi(1,1/beta)/beta^3+1/beta^2 ...
         -sum(abs(X).^beta.*(log(abs(X)).^2))/sum(abs(X).^beta)+sum(abs(X).^beta.*log(abs(X)))^2/sum(abs(X).^beta)^2 ...
         +sum(abs(X).^beta.*log(abs(X)))/(beta*sum(abs(X).^beta))-log(beta*mean(abs(X).^beta))/beta^2;
end