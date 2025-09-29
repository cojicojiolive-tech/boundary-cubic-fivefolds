----------------------------------------
-- 0. Ring and polynomial
----------------------------------------
R = QQ[x0,x1,x2,x3,x4,x5,x6];

sigma1=2;
sigma2=3;
rho=4;
kappa=5;
mu=6;
lambda=7;




F = x2^3 + sigma1*x3^3 + sigma2*x4^3 + rho*x2*x3*x4
      + x1*x2*x5 + x0*x3*x6 + x0*x5^2 + x1^2*x6
      + kappa*x2^2*x3 + mu*x3^2*x4 + lambda*x2^2*x4;

m = ideal vars R;


----------------------------------------
-- 2. Projective singular-locus ideal (via J only)
----------------------------------------
J       = ideal jacobian matrix {{F}};         -- gradient ideal
IprojJ  = saturate(J, m);                      -- projectivized singular locus in P^6


minsJ  = minimalPrimes IprojJ;

print netList minsJ; 




scan(minsJ, P -> (
  print concatenate(
    "codim=", toString codim P,
    ", degree=", toString degree P,
    ", ideal=", toString P
  )
));