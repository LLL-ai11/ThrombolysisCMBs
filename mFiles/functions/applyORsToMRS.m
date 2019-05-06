function [mRS_out, status] = applyORsToMRS(mRS_in, OddsRatio_disability, OddsRatio_mortality, cut)
%applyORsToProb: Applys odds ratio for disability / mortality to a given mRS distribution
%   Input parameters
%   useDisability: remnant flag, indicating if disability and/or mortality is considered in the transformation of mRS scores
%   useMortality: see above
%   cut: mRS cutpoint, assumed to be 2 (mRS 0-2 vs 3-6) unless otherwise specified

try
    cut;
catch
    cut = 2;
end

if cut == 2
    pOld_0to2 = sum(mRS_in(1:3)); pOld_3To6 = 1 - pOld_0to2;
    pOld_3to5 = sum(mRS_in(4:6));
    pOld_alive = sum(mRS_in(1:6)); pOld_dead = 1 - pOld_alive;
    
    pWithin_0to2 = mRS_in(1:3) ./ pOld_0to2; pWithin_3to6 = mRS_in(4:7) ./ pOld_3To6; pWithin_3to5 = mRS_in(4:6) ./ pOld_3to5;
    oddsOld_3to6 = pOld_3To6 / pOld_0to2;
    oddsNew_3to6 = oddsOld_3to6 * OddsRatio_disability;
    pNew_3to6 = oddsNew_3to6 / (1 + oddsNew_3to6); pNew_0to2 = 1 - pNew_3to6;
    
    pWithin_alive = mRS_in(1:6) ./ pOld_alive;
    oddsOld_dead = pOld_dead / pOld_alive;
    oddsNew_dead = oddsOld_dead * OddsRatio_mortality;
    pNew_dead = oddsNew_dead / (1 + oddsNew_dead); pNew_alive = 1 - pNew_dead;
    
    pNew_3to5 = 1 - pNew_0to2 - pNew_dead;
    
    if OddsRatio_disability ~= 0 && OddsRatio_mortality ~= 0
        if pNew_3to5 < 0
            mRS_out = 0; status = 0;
            return;
        end
        mRS_out = [pNew_0to2 * pWithin_0to2, pNew_3to5 * pWithin_3to5, pNew_dead];
        status = 1;
    elseif OddsRatio_disability ~= 0
        mRS_out = [pNew_0to2 * pWithin_0to2, pNew_3to6 * pWithin_3to6];
        status = 1;
    elseif OddsRatio_mortality ~= 0
        mRS_out = [pNew_alive * pWithin_alive, pNew_dead];
        status = 1;
    else
        mRS_out = 0; status = 0;
    end
elseif cut == 1
    pOld_0to1 = sum(mRS_in(1:2)); pOld_2To6 = 1 - pOld_0to1;
    pOld_2to5 = sum(mRS_in(3:6));
    pOld_alive = sum(mRS_in(1:6)); pOld_dead = 1 - pOld_alive;
    
    pWithin_0to1 = mRS_in(1:2) ./ pOld_0to1; pWithin_2to6 = mRS_in(3:7) ./ pOld_2To6; pWithin_2to5 = mRS_in(3:6) ./ pOld_2to5;
    oddsOld_2to6 = pOld_2To6 / pOld_0to1;
    oddsNew_2to6 = oddsOld_2to6 * OddsRatio_disability;
    pNew_2to6 = oddsNew_2to6 / (1 + oddsNew_2to6); pNew_0to1 = 1 - pNew_2to6;
    
    pWithin_alive = mRS_in(1:6) ./ pOld_alive;
    oddsOld_dead = pOld_dead / pOld_alive;
    oddsNew_dead = oddsOld_dead * OddsRatio_mortality;
    pNew_dead = oddsNew_dead / (1 + oddsNew_dead); pNew_alive = 1 - pNew_dead;
    
    pNew_2to5 = 1 - pNew_0to1 - pNew_dead;
    
    if OddsRatio_disability ~= 0 && OddsRatio_mortality ~= 0
        if pNew_2to5 < 0
            mRS_out = 0; status = 0;
            return;
        end
        mRS_out = [pNew_0to1 * pWithin_0to1, pNew_2to5 * pWithin_2to5, pNew_dead];
        status = 1;
    elseif OddsRatio_disability ~= 0
        mRS_out = [pNew_0to1 * pWithin_0to1, pNew_2to6 * pWithin_2to6];
        status = 1;
    elseif OddsRatio_mortality ~= 0
        mRS_out = [pNew_alive * pWithin_alive, pNew_dead];
        status = 1;
    else
        mRS_out = 0; status = 0;
    end
end
mRS_out = mRS_out / sum(mRS_out);
end

