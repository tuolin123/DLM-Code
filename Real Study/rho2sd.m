function dif = rho2sd(rho, sigma, D)
cov = discrete_covariance(sigma, D);
dif = rho - cov(2);
end