%least-squares/gmm-type objective function
function f = gmmmoments(X,theta,W,moments)
%bins=number of bins
%X=data made up of integrals with different powers in the integrand
%W=weighting matrix
%theta=vector of polynomial coefficients

fitted=X*theta;%fitted polynomial values

%empirical values of moment conditions
m = zeros(moments,1);
for i=1:moments
    m(i) = (1/size(X,1))*sum(fitted.^i) - 1/(i+1);
end
 
f = m'*W*m;