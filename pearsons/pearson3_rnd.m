%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute a random sample of (m * n) values, using the the Pearson 3
% distribution. The Hosking and Wallis (1997) version of the distribution
% is chosen.  
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
%    -alpha, beta, xi : parameters of the distribution
%    -Gamma : shape parameter. Depending on its value, the skewness is
%       positive (Gamma > 0) or negative (Gamma < 0). If Gamma = 0, the
%       distribution is normal where the mean is alpha and the standard
%       deviation is beta
%    -m, n : number of rows (m) and colums (n) of the output matrix
%
% Output
%   -x : matrix of the random sample
%
% Source : Hosking, J., & Wallis, J. (1997). Regional Frequency Analysis:
% An Approach Based on L-Moments. Cambridge: Cambridge University Press. 
% doi:10.1017/CBO9780511529443
%
% Guillaume Talbot, INRS-ETE 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x]=pearson3_rnd(alpha,beta,xi,Gamma,m,n)

%% Generate a random matrix (m * n) with values between 0 and 1
p=rand(m,n);

%% If Gamma = 0, the Normal distribution is used
if Gamma==0 
    norm_inv=@(P,mu,sigma) sigma.*(-sqrt(2).*erfcinv(2.*P))+mu;
    x=norm_inv(p,alpha,beta);    
    return
end

%% Compute the inverse cdf of the matrix p to generate the random sample
if Gamma>0 %Case of positive skewness    
    X=gammaincinv(p,alpha.*ones(size(p)));
    x=(X.*beta)+xi;  
else %Case of negative skewness    
    X=gammaincinv(1-p,alpha.*ones(size(p)));    
    x=(-X.*beta)+xi;
end
