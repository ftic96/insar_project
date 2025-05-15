%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate the parameters of a sample for the Pearson 3
% distribution, using the L-moments method. The Hosking and Wallis (1997) 
% version of the distribution is chosen. 
%
% Given the location (mu), scale (sigma) and shape (Gamma) parameters, we can
% estimate the three other parameters used in this version :
%   -alpha : 4 / sigma^2
%   -beta :0.5 * sigma * abs(Gamma)
%   -xi : mu - 2*sigma/Gamma
%
% If Gamma  > 0, the range of x is : xi <= x < Inf 
% If Gamma  = 0, the range of x is : -Inf < x < Inf 
% If Gamma  < 0, the range of x is : -Inf < x <= xi 
%
% Input :
%    -x : vector of values
%   
% Output : 
%   -alpha, beta, xi : parameters of the distribution
%   -Gamma : shape parameter. Depending on its value, the skewness is
%       positive (Gamma > 0) or negative (Gamma < 0). If Gamma = 0, the
%       distribution is normal where the mean is alpha and the standard
%       deviation is beta
%
% If Gamma = 0 (normal distribution case), alpha correspond to the mean 
% (mu), and beta to the standard deviation (sigma)
%
% Source : Hosking, J., & Wallis, J. (1997). Regional Frequency Analysis:
% An Approach Based on L-Moments. Cambridge: Cambridge University Press. 
% doi:10.1017/CBO9780511529443
%
% Guillaume Talbot, INRS-ETE 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha,beta,xi,Gamma]=pearson3_fit(x)

%Remove the missing values
x(isnan(x))=[];

%Number of non-missing values in the vector
n=length(x);
     
%Sort the vector x
x=sort(x);

j=1:n;
if size(j)~=size(x) %If j and x have different orientations, we transpose j
    j=j';
end

%Compute the L-moments : l1, l2, t, t3 and t4
b0=sum(x)/n;
b1j=(j-1)./(n-1).*x;
b2j=(j-1).*(j-2)./(n-1)./(n-2).*x;
% b3j=(j-1).*(j-2).*(j-3)./(n-1)./(n-2)./(n-3).*x; %Not necessary for this distribution

b1=sum(b1j)/n;
b2=sum(b2j)/n;
% b3=sum(b3j)/n; %Not necessary for this distribution

l1=b0;
l2=2*b1-b0;
l3=6*b2-6*b1+b0;
% l4=20*b3-30*b2+12*b1-b0;  %Not necessary for this distribution

% t=l2/l1; %Not necessary for this distribution
t3=l3/l2;
% t4=l4/l2; %Not necessary for this distribution

%Estimate the alpha parameter
if abs(t3)>=(1/3) && abs(t3)<1
    z=1-abs(t3);
    alpha=(0.36067.*z-0.59567.*z.^2 +0.25361.*z.^3)./(1-2.78861.*z+2.56096.*z.^2-0.77045.*z.^3);
else
    z=3.*pi.*t3.^2;
    alpha=(1+0.2906.*z)./(z+0.1882.*z.^2+0.0442.*z.^3);
end

%Estimate the Gamma,beta and xi parameters
Gamma=2./sqrt(alpha).*sign(t3);
mu=l1;
sigma=l2.*sqrt(pi).*sqrt(alpha).*gamma(alpha)./gamma(alpha+.5);
beta=0.5.*sigma.*abs(Gamma);
xi=mu-2.*sigma./Gamma;

%If alpha is bigger than ~171, beta is infinite (overflow) and Gamma is
%very small. So, we choose to select the case where Gamma = 0 and the
%distribution becomes Normal, with alpha = mu and beta = sigma.
if isnan(beta) 
   alpha=l1;
   beta=sqrt(pi).*l2;
   xi=nan;
   Gamma=0;    
end