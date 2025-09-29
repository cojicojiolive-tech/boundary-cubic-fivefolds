----------------------------------------
-- 0. Ring and polynomial
----------------------------------------
R = QQ[x0,x1,x2,x3,x4,x5,x6];

sigma=2;
tau=3;
rho=4;

F = x2^3 + sigma*x3^3 + tau*x4^3 + x5^3 + x2*x3*x4
      + x0^2*x6 + rho*x1^2*x6;

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