%% See page 364 of 
% Higham, D. J. (2008). Modeling and simulating chemical reactions. 
    % SIAM review, 50(2), 347-368.

% ODE15s solution the Reaction Rate Equation for
% the Michaelis-Menten system.
%
% Parameters from Chapter 7 of
% Stochastic Modelling for Systems Biology,
% by Darren J. Wilkinson, Chapman & Hall/CRC, 2006.
%
% Downloadable from
% http://www.maths.strath.ac.uk/˜aas96106/algfiles.html
% along with an extended version that produces graphical output.
tspan = [0 50]; yzero = [5e-7; 2e-7; 0; 0];
options = odeset('AbsTol',1e-8);

[t,y] = ode15s(@mm_rre,tspan,yzero,options);
% Recordor plot (t,y) at this stage
size(t)
size(y)
plot(t,y(:,1),'b',t,y(:,4),'r')
%end