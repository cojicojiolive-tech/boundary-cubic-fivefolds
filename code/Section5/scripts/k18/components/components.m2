----------------------------------------
-- 0. Ring and polynomial
----------------------------------------
R = QQ[x0,x1,x2,x3,x4,x5,x6];

sigma1=2;
sigma2=3;
rho=4;
alpha=5;
beta=6;
gamma=7;



F = x2^3 + sigma1*x3^3 + sigma2*x4^3 + rho*x2*x3*x4
      + x0*x2*x5 + x1*x3*x5 + alpha*x0*x4*x5
      + x0*x2*x6 + beta*x1*x3*x6 + gamma*x1*x4*x6;

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