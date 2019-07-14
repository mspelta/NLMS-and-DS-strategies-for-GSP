% ================================

S = 2;
P_up = 0.1;

kappa_1 = sqrt(2)*erfinv((1-P_up)^(1/S))
kappa_2 = (2/S) * gammaincinv( P_up, 0.5*S,'upper')

% ---------------------------------

pause

% -- CW -----------------------
updateProbability_1 = 1 - (erf(kappa_1/sqrt(2)))^S
% -- l2-norm ------------------
updateProbability_2 = gammainc(0.5*S*kappa_2,0.5*S,'upper')

% ================================