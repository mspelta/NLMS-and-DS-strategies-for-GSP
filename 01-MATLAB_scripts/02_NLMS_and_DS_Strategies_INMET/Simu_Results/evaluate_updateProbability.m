function [ updateProbability ] = evaluate_updateProbability(dataselection_algorithm, kappa, D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    S = sum(diag(D)); % D is the sampling matrix

    if(dataselection_algorithm == 1)
        % -- Component-wise -----------
        updateProbability = 1 - (erf(kappa/sqrt(2)))^S;
    else
        % -- l2-norm ------------------
        updateProbability = gammainc(0.5*S*kappa,0.5*S,'upper');
    end

end

