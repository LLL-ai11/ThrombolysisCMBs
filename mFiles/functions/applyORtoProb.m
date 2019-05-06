function p_out = applyORtoProb(p_in, oddsRatio)
%applyORtoProb: Summary of this function goes here
%   Detailed explanation goes here

odds_old = p_in / (1 - p_in);
odds_new = odds_old * oddsRatio;
p_out = odds_new / (1 + odds_new);
end

