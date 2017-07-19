%% See page 364 of 
% Higham, D. J. (2008). Modeling and simulating chemical reactions. 
    % SIAM review, 50(2), 347-368.

function yprime = mm_rre(t,y)
k1 = 1e6; k2 = 1e-4; k3 = 0.1;
% MM_RRE Michaelis-Menten Reaction Rate Equation
yprime = zeros(4,1);
yprime(1) = -k1*y(1)*y(2) + k2*y(3);
yprime(2) = -k1*y(1)*y(2) + (k2+k3)*y(3);
yprime(3) = k1*y(1)*y(2) - (k2+k3)*y(3);
yprime(4) = k3*y(3);
end