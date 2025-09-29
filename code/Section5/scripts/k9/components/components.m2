----------------------------------------
-- 0. Ring and polynomial
----------------------------------------
R = QQ[x0,x1,x2,x3,x4,x5,x6];

tau=2;
rho=3;


F = x2^2*x3 + tau*x2*x3^2 + x0*x4^2 + rho*x0*x5^2 + x1*x4^2 + x1*x5^2
      + x0*x2*x6 + x1*x3*x6;

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