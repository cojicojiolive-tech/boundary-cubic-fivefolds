----------------------------------------
-- 0. Ring and polynomial
----------------------------------------
R = QQ[x0,x1,x2,x3,x4,x5,x6];
rho=2;
sigma=3;
F = x3^3 + x0*x4^2 + x0*x5^2 + x0*x6^2 + x1^2*x4
      + rho*x1*x2*x5 + sigma*x2^2*x6;
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