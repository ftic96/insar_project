clear

alpha=10;
beta=0.3;
xi=5;
Gamma=1;

%Generate a random sample using the Pearson 3 distribution
X=pearson3_rnd(alpha,beta,xi,Gamma,100000,1);

%Empirical histogram of the random sample
subplot(3,1,1)
hist(X,5:0.1:13)
title('Random sample','fontsize',14)
xlabel('x')
ylabel('Bin Count')

%Compute the corresponding probability density function 
x=5:0.01:13;
subplot(3,1,2)
plot(x,pearson3_pdf(x,alpha,beta,xi,Gamma));
title('Probability density function','fontsize',14)
xlabel('x')
ylabel('PDF')

%Compute the corresponding cumulative distribution function 
subplot(3,1,3)
plot(x,pearson3_cdf(x,alpha,beta,xi,Gamma));
title('Cumulative distribution function','fontsize',14)
xlabel('x')
ylabel('CDF')

%Estimate the parameters of the distribution for the random sample
[alpha_est,beta_est,xi_est,Gamma_est]=pearson3_fit(X);

%Compare the "real" parameters to the estimated parameters
Parameters_comparison=[alpha alpha_est ; beta beta_est ; xi xi_est] 

%Compute the theoretical 25th, 50th, 75th and 90th percentiles, using the Pearson 3 distribution 
p=[0.25 0.5 0.75 0.9];
theoretical_percentiles=pearson3_inv(p,alpha_est,beta_est,xi_est,Gamma_est);

%Compute the empirical 25th, 50th and 75th percentiles
empirical_percentiles=quantile(X,p);

%Compare the empirical and theoretical percentiles
Percentiles_comparison=[empirical_percentiles' theoretical_percentiles']