function [impactCMB_many_vs_no, impacttPA_yes_vs_no_noCMBs, impacttPA_yes_vs_no_manyCMBs] = CMB_tPA_selection_v2(useDisability, useMortality, ...
    mRS_control, mRS_tPA, ...
    sensiFlag, univSensiAnalysis, v, a)
%CMB_tPA_selection: Core of the multi-step algorithm
%   Input parameters
%   useDisability:      remnant flag, indicating if disability and/or mortality is considered in the transformation of mRS scores
%   useMortality:       see above
%   mRS_control:        estimated 90-day mRS-distributions for a patient with given age and NIHSS without tPA treatment
%   mRS_tPA:            estimated 90-day mRS-distributions for a patient with given age and NIHSS without tPA treatment
%   sensiFlag:          flag for probabilistic sensitivity analzsis (1 / 0)
%   univSensiAnalysis:  array containing flags for univariate sensitivity analyses
%   v                   verbose mode yes/no
%   a                   index in age array, needed for age-specific prevalence of >10 CMBs

%   Output
%   see code below


%% Initialize
accuracySolver = 0.0001;
qualityWeights = [0.975, 0.92, 0.795, 0.635, 0.375, 0.055, 0];  % Dijkland, Chaisinanunkul

% Odds ratios / effects
if v; fprintf('Initializing odds ratios / effect data ...'); end
% increased risk of poor outcome with 1-10 CMBs in treated patients, Charidimou et al.
OR_someCMBs_mRS_3to6 = 1;
OR_someCMBs_mRS_dead = 1;
% increased risk of poor outcome with >10 CMBs in treated patients, Charidimou et al.
OR_manyCMBs_mRS_3to6 = 3.99; OR_manyCMBs_mRS_3to6_95CI = [1.55, 10.22];
OR_manyCMBs_mRS_dead = 2.44; OR_manyCMBs_mRS_dead_95CI = [0.92, 6.49];

% increased risk to ICH with 1-10 CMBs, Charidimou et al.
OR_someCMBs_ICH = 1;
% increased risk to ICH with >10 CMBs, Charidimou et al.
OR_manyCMBs_ICH = 3.65; OR_manyCMBs_ICH_95CI = [1.17, 11.42];

% increased risk of poor outcome when sICH occurs, Strbian et al.
OR_ICH_mRS_3to6 = 2.68; OR_ICH_mRS_3to6_95CI = [2.35, 3.19];
OR_ICH_mRS_dead = 5.19;OR_ICH_mRS_dead_95CI = [3.55, 7.65];

% increased risk to ICH when given tPA vs controls, Emberson et al.
OR_tPA_ICH = 5.55;  OR_tPA_ICH_95CI = [4.01, 7.70];

% Vary parameters in sensitivity analyses.
if sensiFlag
    r = rand();
    m = log(OR_manyCMBs_mRS_3to6);     CI = log(OR_manyCMBs_mRS_3to6_95CI);
    if r > 0.5; SD = (CI(2) - m) / 1.96;
    else;       SD = (m - CI(1)) / 1.96; end
    OR_manyCMBs_mRS_3to6 = max(1, exp(norminv(r, m, SD)));
    m = log(OR_manyCMBs_mRS_dead);     CI = log(OR_manyCMBs_mRS_dead_95CI);
    if r > 0.5; SD = (CI(2) - m) / 1.96;
    else;       SD = (m - CI(1)) / 1.96; end
    OR_manyCMBs_mRS_dead = max(1, exp(norminv(r, m, SD)));
    m = log(OR_manyCMBs_ICH);     CI = log(OR_manyCMBs_ICH_95CI);
    if r > 0.5; SD = (CI(2) - m) / 1.96;
    else;       SD = (m - CI(1)) / 1.96; end
    OR_manyCMBs_ICH = max(1, exp(norminv(r, m, SD)));
       
    r = rand();
    m = log(OR_ICH_mRS_3to6);     CI = log(OR_ICH_mRS_3to6_95CI);
    if r > 0.5; SD = (CI(2) - m) / 1.96;
    else;       SD = (m - CI(1)) / 1.96; end
    OR_ICH_mRS_3to6 = max(1, exp(norminv(r, m, SD)));
    m = log(OR_ICH_mRS_dead);    CI = log(OR_ICH_mRS_dead_95CI);
    if r > 0.5; SD = (CI(2) - m) / 1.96;
    else;       SD = (m - CI(1)) / 1.96; end
    OR_ICH_mRS_dead = max(1, exp(norminv(r, m, SD)));
    
    r = rand();
    m = log(OR_tPA_ICH);    CI = log(OR_tPA_ICH_95CI);
    if r > 0.5; SD = (CI(2) - m) / 1.96;
    else;       SD = (m - CI(1)) / 1.96; end
    OR_tPA_ICH = max(1, exp(norminv(r, m, SD)));
    
    clear r mu SD
end
if v; fprintf(' completed.\n'); end

% Raw probabilities
if v; fprintf('Setting raw probabilities ...'); end
% age-adjusted prevalenca of >10 CMBs; Charidimou et al. & Poels et al.
pCMB_manyCMBs_array = [0.00372	0.00766	0.01385	0.02227	0.03294];  % corresponding to age groups 50:10:90
pCMB_manyCMBs_params_array =    [61.26, 16310; ...
                                180.9, 23390; ...
                                109, 7741; ...
                                65.07, 2844; ...
                                45.51, 1327];
pCMB_manyCMBs = pCMB_manyCMBs_array(a); pCMB_manyCMBs_params = pCMB_manyCMBs_params_array(a,:);        
pCMB_someCMBs = .25;  % dummy parameter

if sensiFlag
    r = rand(); pCMB_manyCMBs = betainv(r, pCMB_manyCMBs_params(1), pCMB_manyCMBs_params(2));
end
pCMB_noCMBs = 1 - pCMB_manyCMBs - pCMB_someCMBs;

pICH_ICH_tPA_avCMB = 0.07; pICH_ICH_tPA_avCMB_params = [69, 985 - 69];   % Strbian et al.
if sensiFlag
    r = rand(); pICH_ICH_tPA_avCMB = betainv(r, pICH_ICH_tPA_avCMB_params(1), pICH_ICH_tPA_avCMB_params(2));
end
pICH_ICH_control_avCMB = applyORtoProb(pICH_ICH_tPA_avCMB, 1/OR_tPA_ICH);

pICH_noICH_tPA_avCMB = 1 - pICH_ICH_tPA_avCMB;
pICH_noICH_control_avCMB = 1 - pICH_ICH_control_avCMB;
if v; fprintf(' completed.\n'); end

% Conditional probabilities (derived)
if v; fprintf('Preparing conditional probabilities ... \n'); end
if v; fprintf('... Probability of sICH, given CMB status, control group ...'); end
upper_p = 1; lower_p = 0; % initialize for solver
cnt = 1;
while true
    pICH_ICH_control_noCMB = (upper_p + lower_p) / 2;
    % switch according to sensitivity analysis flag
    switch univSensiAnalysis(1)
        case 0 % base case model
            pICH_ICH_control_someCMB = applyORtoProb(pICH_ICH_control_noCMB, OR_someCMBs_ICH);
            pICH_ICH_control_manyCMB = applyORtoProb(pICH_ICH_control_noCMB, OR_manyCMBs_ICH);
        case 1
            pICH_ICH_control_someCMB = applyORtoProb(pICH_ICH_control_noCMB, 1);
            pICH_ICH_control_manyCMB = applyORtoProb(pICH_ICH_control_noCMB, 1);
    end
    pICH_ICH_control_avCMB_expected = pICH_ICH_control_noCMB * pCMB_noCMBs + pICH_ICH_control_someCMB * pCMB_someCMBs + pICH_ICH_control_manyCMB * pCMB_manyCMBs;
    delta = pICH_ICH_control_avCMB_expected - pICH_ICH_control_avCMB;
    if abs(delta) > accuracySolver/100
        if delta > 0; upper_p = pICH_ICH_control_noCMB;
        else; lower_p = pICH_ICH_control_noCMB;
        end
        cnt = cnt + 1;
    else
        break;
    end
end
pICH_noICH_control_noCMB = 1- pICH_ICH_control_noCMB;
pICH_noICH_control_someCMB = 1- pICH_ICH_control_someCMB;
pICH_noICH_control_manyCMB = 1- pICH_ICH_control_manyCMB;
if v; fprintf(' completed in n = %d runs. Delta: %.6f.\n', cnt, delta); end
clear cnt pICH_ICH_control_avCMB_expected delta upper_p lower_p

if v; fprintf('... Probability of sICH, given CMB status, tPA group ...'); end
upper_p = 1; lower_p = 0; % initialize for solver

cnt = 1;
while true
    pICH_ICH_tPA_noCMB = (upper_p + lower_p) / 2;
    pICH_ICH_tPA_someCMB = applyORtoProb(pICH_ICH_tPA_noCMB, OR_someCMBs_ICH);
    pICH_ICH_tPA_manyCMB = applyORtoProb(pICH_ICH_tPA_noCMB, OR_manyCMBs_ICH);
    pICH_ICH_tPA_avCMB_expected = pICH_ICH_tPA_noCMB * pCMB_noCMBs + pICH_ICH_tPA_someCMB * pCMB_someCMBs + pICH_ICH_tPA_manyCMB * pCMB_manyCMBs;
    delta = pICH_ICH_tPA_avCMB_expected - pICH_ICH_tPA_avCMB;
    if abs(delta) > accuracySolver/100
        if delta > 0; upper_p = pICH_ICH_tPA_noCMB;
        else; lower_p = pICH_ICH_tPA_noCMB;
        end
        cnt = cnt + 1;
    else
        break;
    end
end
pICH_noICH_tPA_noCMB = 1- pICH_ICH_tPA_noCMB;
pICH_noICH_tPA_someCMB = 1- pICH_ICH_tPA_someCMB;
pICH_noICH_tPA_manyCMB = 1- pICH_ICH_tPA_manyCMB;
if v; fprintf(' completed in n = %d runs. Delta: %.6f.\n', cnt, delta); end
clear cnt pICH_ICH_tPA_avCMB_expected delta upper_p lower_p

if v; fprintf('... Probability of CMB, given sICH status, control group ...'); end
pCMB_noCMBs_control_noICH = pICH_noICH_control_noCMB * pCMB_noCMBs / pICH_noICH_control_avCMB;
pCMB_noCMBs_control_ICH = pICH_ICH_control_noCMB * pCMB_noCMBs / pICH_ICH_control_avCMB;
pCMB_someCMBs_control_noICH = pICH_noICH_control_someCMB * pCMB_someCMBs / pICH_noICH_control_avCMB;
pCMB_someCMBs_control_ICH = pICH_ICH_control_someCMB * pCMB_someCMBs / pICH_ICH_control_avCMB;
pCMB_manyCMBs_control_noICH = pICH_noICH_control_manyCMB * pCMB_manyCMBs / pICH_noICH_control_avCMB;
pCMB_manyCMBs_control_ICH = pICH_ICH_control_manyCMB * pCMB_manyCMBs / pICH_ICH_control_avCMB;
if v; fprintf(' completed.\n'); end

if v; fprintf('... Probability of CMB, given sICH status, tPA group ...'); end
pCMB_noCMBs_tPA_noICH = pICH_noICH_tPA_noCMB * pCMB_noCMBs / pICH_noICH_tPA_avCMB;
pCMB_noCMBs_tPA_ICH = pICH_ICH_tPA_noCMB * pCMB_noCMBs / pICH_ICH_tPA_avCMB;
pCMB_someCMBs_tPA_noICH = pICH_noICH_tPA_someCMB * pCMB_someCMBs / pICH_noICH_tPA_avCMB;
pCMB_someCMBs_tPA_ICH = pICH_ICH_tPA_someCMB * pCMB_someCMBs / pICH_ICH_tPA_avCMB;
pCMB_manyCMBs_tPA_noICH = pICH_noICH_tPA_manyCMB * pCMB_manyCMBs / pICH_noICH_tPA_avCMB;
pCMB_manyCMBs_tPA_ICH = pICH_ICH_tPA_manyCMB * pCMB_manyCMBs / pICH_ICH_tPA_avCMB;
if v; fprintf(' completed.\n'); end


%% Main 2 - tPA group
% Primary subgroups according to CMB-number
if v; fprintf('Estimating primary subgroups according to CMB-number in tPA patients ...'); end
cnt = 1;
upper_mRS_3to6 = 2; lower_mRS_3to6 = 0;  % initialize for solver
upper_mRS_dead = 2; lower_mRS_dead = 0;
while true
    try
        OR_noCMBs_mRS_3to6 = (upper_mRS_3to6 + lower_mRS_3to6) / 2;
        OR_noCMBs_mRS_dead = (upper_mRS_dead + lower_mRS_dead) / 2;
        mRS_tPA_noCMBs = applyORsToMRS(mRS_tPA, OR_noCMBs_mRS_3to6 * useDisability, OR_noCMBs_mRS_dead * useMortality);
        mRS_tPA_someCMBs = applyORsToMRS(mRS_tPA_noCMBs, OR_someCMBs_mRS_3to6 * useDisability, OR_someCMBs_mRS_dead * useMortality);
        mRS_tPA_manyCMBs = applyORsToMRS(mRS_tPA_noCMBs, OR_manyCMBs_mRS_3to6 * useDisability, OR_manyCMBs_mRS_dead * useMortality);
        
        mRS_0to2_expected = sum(mRS_tPA_noCMBs(1:3)) * pCMB_noCMBs + sum(mRS_tPA_someCMBs(1:3)) * pCMB_someCMBs + sum(mRS_tPA_manyCMBs(1:3)) * pCMB_manyCMBs;
        mRS_0to2_observed = sum(mRS_tPA(1:3));
        delta_mRS_0to2 = mRS_0to2_expected - mRS_0to2_observed;
        
        mRS_dead_expected = mRS_tPA_noCMBs(7) * pCMB_noCMBs + mRS_tPA_someCMBs(7) * pCMB_someCMBs + mRS_tPA_manyCMBs(7) * pCMB_manyCMBs;
        mRS_dead_observed = mRS_tPA(7);
        delta_dead = mRS_dead_expected - mRS_dead_observed;
    catch
        impactCMB_many_vs_no = NaN*[1 1 1]; impacttPA_yes_vs_no_noCMBs = NaN*[1 1 1]; impacttPA_yes_vs_no_manyCMBs = NaN*[1 1 1];
        % fprintf('      try catch. \n')
        return;
    end
    breakLoop = 1;
    if useDisability && abs(delta_mRS_0to2) > accuracySolver
        if delta_mRS_0to2 < 0
           upper_mRS_3to6 = OR_noCMBs_mRS_3to6; 
        else
           lower_mRS_3to6 = OR_noCMBs_mRS_3to6;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if useMortality && abs(delta_dead) > accuracySolver
        if delta_dead > 0
            upper_mRS_dead = OR_noCMBs_mRS_dead;
        else
            lower_mRS_dead = OR_noCMBs_mRS_dead;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if breakLoop
        break;
    end
end
if v; fprintf(' completed in n = %d runs. Delta: %.6f / %.6f.\n', cnt, delta_mRS_0to2, delta_dead); end
clear cnt mRS_0to2_expected mRS_0to2_observed delta_mRS_0to2 mRS_dead_expected mRS_dead_observed delta_dead
clear upper_mRS_3to6 lower_mRS_3to6 upper_mRS_dead lower_mRS_dead
clear breakLoop
clear OR_noCMBs_mRS_3to6 OR_noCMBs_mRS_dead

% Secondary subgroups according to sICH-status
if v; fprintf('Estimating secondary subgroups according to sICH-status in tPA patients ... \n'); end
if v; fprintf('... patients with no CMBs ...'); end
cnt = 1;
upper_mRS_3to6 = 2; lower_mRS_3to6 = 0;  % initialize for solver
upper_mRS_dead = 2; lower_mRS_dead = 0;
while true
    try
        OR_noICH_mRS_3to6 = (upper_mRS_3to6 + lower_mRS_3to6) / 2;
        OR_noICH_mRS_dead = (upper_mRS_dead + lower_mRS_dead) / 2;
        mRS_tPA_noCMBs_noICH = applyORsToMRS(mRS_tPA_noCMBs, OR_noICH_mRS_3to6 * useDisability, OR_noICH_mRS_dead * useMortality);
        mRS_tPA_noCMBs_ICH = applyORsToMRS(mRS_tPA_noCMBs_noICH, OR_ICH_mRS_3to6 * useDisability, OR_ICH_mRS_dead * useMortality);
        
        mRS_0to2_expected = sum(mRS_tPA_noCMBs_noICH(1:3)) * pICH_noICH_tPA_noCMB + sum(mRS_tPA_noCMBs_ICH(1:3)) * pICH_ICH_tPA_noCMB;
        mRS_0to2_observed = sum(mRS_tPA_noCMBs(1:3));
        delta_mRS_0to2 = mRS_0to2_expected - mRS_0to2_observed;
        
        mRS_dead_expected = mRS_tPA_noCMBs_noICH(7) * pICH_noICH_tPA_noCMB + mRS_tPA_noCMBs_ICH(7) * pICH_ICH_tPA_noCMB;
        mRS_dead_observed = mRS_tPA_noCMBs(7);
        delta_dead = mRS_dead_expected - mRS_dead_observed;
    catch
        impactCMB_many_vs_no = NaN*[1 1 1]; impacttPA_yes_vs_no_noCMBs = NaN*[1 1 1]; impacttPA_yes_vs_no_manyCMBs = NaN*[1 1 1];
        % fprintf('      try catch. \n')
        return;
    end
    breakLoop = 1;
    if useDisability && abs(delta_mRS_0to2) > accuracySolver
        if delta_mRS_0to2 < 0
           upper_mRS_3to6 = OR_noICH_mRS_3to6; 
        else
           lower_mRS_3to6 = OR_noICH_mRS_3to6;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if useMortality && abs(delta_dead) > accuracySolver
        if delta_dead > 0
            upper_mRS_dead = OR_noICH_mRS_dead;
        else
            lower_mRS_dead = OR_noICH_mRS_dead;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if breakLoop
        break;
    end
end
if v; fprintf(' completed in n = %d runs. Delta: %.6f / %.6f.\n', cnt, delta_mRS_0to2, delta_dead); end
clear cnt mRS_0to2_expected mRS_0to2_observed delta_mRS_0to2 mRS_dead_expected mRS_dead_observed delta_dead
clear upper_mRS_3to6 lower_mRS_3to6 upper_mRS_dead lower_mRS_dead
clear breakLoop
clear OR_noICH_mRS_3to6 OR_noICH_mRS_dead

if v; fprintf('... patients with some CMBs ...'); end
cnt = 1;
upper_mRS_3to6 = 2; lower_mRS_3to6 = 0;  % initialize for solver
upper_mRS_dead = 2; lower_mRS_dead = 0;
while true
    try
        OR_noICH_mRS_3to6 = (upper_mRS_3to6 + lower_mRS_3to6) / 2;
        OR_noICH_mRS_dead = (upper_mRS_dead + lower_mRS_dead) / 2;
        mRS_tPA_someCMBs_noICH = applyORsToMRS(mRS_tPA_someCMBs, OR_noICH_mRS_3to6 * useDisability, OR_noICH_mRS_dead * useMortality);
        mRS_tPA_someCMBs_ICH = applyORsToMRS(mRS_tPA_someCMBs_noICH, OR_ICH_mRS_3to6 * useDisability, OR_ICH_mRS_dead * useMortality);
        
        mRS_0to2_expected = sum(mRS_tPA_someCMBs_noICH(1:3)) * pICH_noICH_tPA_someCMB + sum(mRS_tPA_someCMBs_ICH(1:3)) * pICH_ICH_tPA_someCMB;
        mRS_0to2_observed = sum(mRS_tPA_someCMBs(1:3));
        delta_mRS_0to2 = mRS_0to2_expected - mRS_0to2_observed;
        
        mRS_dead_expected = mRS_tPA_someCMBs_noICH(7) * pICH_noICH_tPA_someCMB + mRS_tPA_someCMBs_ICH(7) * pICH_ICH_tPA_someCMB;
        mRS_dead_observed = mRS_tPA_someCMBs(7);
        delta_dead = mRS_dead_expected - mRS_dead_observed;
    catch
        impactCMB_many_vs_no = NaN*[1 1 1]; impacttPA_yes_vs_no_noCMBs = NaN*[1 1 1]; impacttPA_yes_vs_no_manyCMBs = NaN*[1 1 1];
        % fprintf('      try catch. \n')
        return;
    end
    breakLoop = 1;
    if useDisability && abs(delta_mRS_0to2) > accuracySolver
        if delta_mRS_0to2 < 0
           upper_mRS_3to6 = OR_noICH_mRS_3to6; 
        else
           lower_mRS_3to6 = OR_noICH_mRS_3to6;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if useMortality && abs(delta_dead) > accuracySolver
        if delta_dead > 0
            upper_mRS_dead = OR_noICH_mRS_dead;
        else
            lower_mRS_dead = OR_noICH_mRS_dead;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if breakLoop
        break;
    end
end
if v; fprintf(' completed in n = %d runs. Delta: %.6f / %.6f.\n', cnt, delta_mRS_0to2, delta_dead); end
clear cnt mRS_0to2_expected mRS_0to2_observed delta_mRS_0to2 mRS_dead_expected mRS_dead_observed delta_dead
clear upper_mRS_3to6 lower_mRS_3to6 upper_mRS_dead lower_mRS_dead
clear breakLoop
clear OR_noICH_mRS_3to6 OR_noICH_mRS_dead

if v; fprintf('... patients with many CMBs ...'); end
cnt = 1;
upper_mRS_3to6 = 2; lower_mRS_3to6 = 0;  % initialize for solver
upper_mRS_dead = 2; lower_mRS_dead = 0;
while true
    try
        OR_noICH_mRS_3to6 = (upper_mRS_3to6 + lower_mRS_3to6) / 2;
        OR_noICH_mRS_dead = (upper_mRS_dead + lower_mRS_dead) / 2;
        mRS_tPA_manyCMBs_noICH = applyORsToMRS(mRS_tPA_manyCMBs, OR_noICH_mRS_3to6 * useDisability, OR_noICH_mRS_dead * useMortality);
        mRS_tPA_manyCMBs_ICH = applyORsToMRS(mRS_tPA_manyCMBs_noICH, OR_ICH_mRS_3to6 * useDisability, OR_ICH_mRS_dead * useMortality);
        
        mRS_0to2_expected = sum(mRS_tPA_manyCMBs_noICH(1:3)) * pICH_noICH_tPA_manyCMB + sum(mRS_tPA_manyCMBs_ICH(1:3)) * pICH_ICH_tPA_manyCMB;
        mRS_0to2_observed = sum(mRS_tPA_manyCMBs(1:3));
        delta_mRS_0to2 = mRS_0to2_expected - mRS_0to2_observed;
        
        mRS_dead_expected = mRS_tPA_manyCMBs_noICH(7) * pICH_noICH_tPA_manyCMB + mRS_tPA_manyCMBs_ICH(7) * pICH_ICH_tPA_manyCMB;
        mRS_dead_observed = mRS_tPA_manyCMBs(7);
        delta_dead = mRS_dead_expected - mRS_dead_observed;
    catch
        impactCMB_many_vs_no = NaN*[1 1 1]; impacttPA_yes_vs_no_noCMBs = NaN*[1 1 1]; impacttPA_yes_vs_no_manyCMBs = NaN*[1 1 1];
        % fprintf('      try catch. \n')
        return;
    end
    breakLoop = 1;
    if useDisability && abs(delta_mRS_0to2) > accuracySolver
        if delta_mRS_0to2 < 0
            upper_mRS_3to6 = OR_noICH_mRS_3to6;
        else
            lower_mRS_3to6 = OR_noICH_mRS_3to6;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if useMortality && abs(delta_dead) > accuracySolver
        if delta_dead > 0
            upper_mRS_dead = OR_noICH_mRS_dead;
        else
            lower_mRS_dead = OR_noICH_mRS_dead;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if breakLoop
        break;
    end
end
if v; fprintf(' completed in n = %d runs. Delta: %.6f / %.6f.\n', cnt, delta_mRS_0to2, delta_dead); end
clear cnt mRS_0to2_expected mRS_0to2_observed delta_mRS_0to2 mRS_dead_expected mRS_dead_observed delta_dead
clear upper_mRS_3to6 lower_mRS_3to6 upper_mRS_dead lower_mRS_dead
clear breakLoop
clear OR_noICH_mRS_3to6 OR_noICH_mRS_dead


if v; fprintf('Calculating ORs for poor outcome with CMBs, adjusted for sICH ...'); end
% increased risk of poor outcome with 1-10 CMBs, adj. for ICH
OR_someCMBs_adjICH_mRS_3to6 = sum(mRS_tPA_someCMBs_noICH(4:7)) / sum(mRS_tPA_someCMBs_noICH(1:3)) / sum(mRS_tPA_noCMBs_noICH(4:7)) * sum(mRS_tPA_noCMBs_noICH(1:3));
OR_someCMBs_adjICH_mRS_dead = mRS_tPA_someCMBs_noICH(7) / sum(mRS_tPA_someCMBs_noICH(1:6)) / mRS_tPA_noCMBs_noICH(7) * sum(mRS_tPA_noCMBs_noICH(1:6));

% increased risk of poor outcome with >10 CMBs, adj. for ICH
OR_manyCMBs_adjICH_mRS_3to6 = sum(mRS_tPA_manyCMBs_noICH(4:7)) / sum(mRS_tPA_manyCMBs_noICH(1:3)) / sum(mRS_tPA_noCMBs_noICH(4:7)) * sum(mRS_tPA_noCMBs_noICH(1:3));
OR_manyCMBs_adjICH_mRS_dead = mRS_tPA_manyCMBs_noICH(7) / sum(mRS_tPA_manyCMBs_noICH(1:6)) / mRS_tPA_noCMBs_noICH(7) * sum(mRS_tPA_noCMBs_noICH(1:6));

if v; fprintf(' completed.\n'); end


%% Main 3 - control group
% Primary subgroups according to CMB-number
if v; fprintf('Estimating primary subgroups according to sICH-number in control patients ...'); end
cnt = 1;
upper_mRS_3to6 = 2; lower_mRS_3to6 = 0;  % initialize for solver
upper_mRS_dead = 2; lower_mRS_dead = 0;
while true
    try
        OR_noICH_mRS_3to6 = (upper_mRS_3to6 + lower_mRS_3to6) / 2;
        OR_noICH_mRS_dead = (upper_mRS_dead + lower_mRS_dead) / 2;
        mRS_control_noICH = applyORsToMRS(mRS_control, OR_noICH_mRS_3to6 * useDisability, OR_noICH_mRS_dead * useMortality);
        mRS_control_ICH = applyORsToMRS(mRS_control_noICH, OR_ICH_mRS_3to6 * useDisability, OR_ICH_mRS_dead * useMortality);
        
        mRS_0to2_expected = sum(mRS_control_noICH(1:3)) * pICH_noICH_control_avCMB + sum(mRS_control_ICH(1:3)) * pICH_ICH_control_avCMB;
        mRS_0to2_observed = sum(mRS_control(1:3));
        delta_mRS_0to2 = mRS_0to2_expected - mRS_0to2_observed;
        
        mRS_dead_expected = mRS_control_noICH(7) * pICH_noICH_control_avCMB + mRS_control_ICH(7) * pICH_ICH_control_avCMB;
        mRS_dead_observed = mRS_control(7);
        delta_dead = mRS_dead_expected - mRS_dead_observed;
    catch
        impactCMB_many_vs_no = NaN*[1 1 1]; impacttPA_yes_vs_no_noCMBs = NaN*[1 1 1]; impacttPA_yes_vs_no_manyCMBs = NaN*[1 1 1];
        % fprintf('      try catch. \n')
        return;
    end
    breakLoop = 1;
    if useDisability && abs(delta_mRS_0to2) > accuracySolver
        if delta_mRS_0to2 < 0
           upper_mRS_3to6 = OR_noICH_mRS_3to6; 
        else
           lower_mRS_3to6 = OR_noICH_mRS_3to6;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if useMortality && abs(delta_dead) > accuracySolver
        if delta_dead > 0
            upper_mRS_dead = OR_noICH_mRS_dead;
        else
            lower_mRS_dead = OR_noICH_mRS_dead;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if breakLoop
        break;
    end
end
if v; fprintf(' completed in n = %d runs. Delta: %.6f / %.6f.\n', cnt, delta_mRS_0to2, delta_dead); end
clear cnt mRS_0to2_expected mRS_0to2_observed delta_mRS_0to2 mRS_dead_expected mRS_dead_observed delta_dead
clear upper_mRS_3to6 lower_mRS_3to6 upper_mRS_dead lower_mRS_dead
clear breakLoop
clear OR_noICH_mRS_3to6 OR_noICH_mRS_dead

% Primary subgroups according to CMB-number
if v; fprintf('Estimating secondary subgroups according to CMB-number in control patients ...'); end
if v; fprintf('... patients with no sICH ...'); end
cnt = 1;
upper_mRS_3to6 = 2; lower_mRS_3to6 = 0;  % initialize for solver
upper_mRS_dead = 2; lower_mRS_dead = 0;
while true
    try
        OR_noCMBs_mRS_3to6 = (upper_mRS_3to6 + lower_mRS_3to6) / 2;
        OR_noCMBs_mRS_dead = (upper_mRS_dead + lower_mRS_dead) / 2;
        mRS_control_noICH_noCMBs = applyORsToMRS(mRS_control_noICH, OR_noCMBs_mRS_3to6 * useDisability, OR_noCMBs_mRS_dead * useMortality);
        mRS_control_noICH_someCMBs = applyORsToMRS(mRS_control_noICH_noCMBs, OR_someCMBs_adjICH_mRS_3to6 * useDisability, OR_someCMBs_adjICH_mRS_dead * useMortality);
        mRS_control_noICH_manyCMBs = applyORsToMRS(mRS_control_noICH_noCMBs, OR_manyCMBs_adjICH_mRS_3to6 * useDisability, OR_manyCMBs_adjICH_mRS_dead * useMortality);
        
        mRS_0to2_expected = sum(mRS_control_noICH_noCMBs(1:3)) * pCMB_noCMBs_control_noICH + sum(mRS_control_noICH_someCMBs(1:3)) * pCMB_someCMBs_control_noICH + sum(mRS_control_noICH_manyCMBs(1:3)) * pCMB_manyCMBs_control_noICH;
        mRS_0to2_observed = sum(mRS_control_noICH(1:3));
        delta_mRS_0to2 = mRS_0to2_expected - mRS_0to2_observed;
        
        mRS_dead_expected = mRS_control_noICH_noCMBs(7) * pCMB_noCMBs_control_noICH + mRS_control_noICH_someCMBs(7) * pCMB_someCMBs_control_noICH + mRS_control_noICH_manyCMBs(7) * pCMB_manyCMBs_control_noICH;
        mRS_dead_observed = mRS_control_noICH(7);
        delta_dead = mRS_dead_expected - mRS_dead_observed;
    catch
        impactCMB_many_vs_no = NaN*[1 1 1]; impacttPA_yes_vs_no_noCMBs = NaN*[1 1 1]; impacttPA_yes_vs_no_manyCMBs = NaN*[1 1 1];
        % fprintf('      try catch. \n')
        return;
    end
    breakLoop = 1;
    if useDisability && abs(delta_mRS_0to2) > accuracySolver
        if delta_mRS_0to2 < 0
           upper_mRS_3to6 = OR_noCMBs_mRS_3to6; 
        else
           lower_mRS_3to6 = OR_noCMBs_mRS_3to6;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if useMortality && abs(delta_dead) > accuracySolver
        if delta_dead > 0
            upper_mRS_dead = OR_noCMBs_mRS_dead;
        else
            lower_mRS_dead = OR_noCMBs_mRS_dead;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if breakLoop
        break;
    end
end
if v; fprintf(' completed in n = %d runs. Delta: %.6f / %.6f.\n', cnt, delta_mRS_0to2, delta_dead); end
clear cnt mRS_0to2_expected mRS_0to2_observed delta_mRS_0to2 mRS_dead_expected mRS_dead_observed delta_dead
clear upper_mRS_3to6 lower_mRS_3to6 upper_mRS_dead lower_mRS_dead
clear breakLoop
clear OR_noCMBs_mRS_3to6 OR_noCMBs_mRS_dead

if v; fprintf('... patients with sICH ...'); end
cnt = 1;
upper_mRS_3to6 = 2; lower_mRS_3to6 = 0;  % initialize for solver
upper_mRS_dead = 2; lower_mRS_dead = 0;
while true
    try
        OR_noCMBs_mRS_3to6 = (upper_mRS_3to6 + lower_mRS_3to6) / 2;
        OR_noCMBs_mRS_dead = (upper_mRS_dead + lower_mRS_dead) / 2;
        mRS_control_ICH_noCMBs = applyORsToMRS(mRS_control_ICH, OR_noCMBs_mRS_3to6 * useDisability, OR_noCMBs_mRS_dead * useMortality);
        mRS_control_ICH_someCMBs = applyORsToMRS(mRS_control_ICH_noCMBs, OR_someCMBs_adjICH_mRS_3to6 * useDisability, OR_someCMBs_adjICH_mRS_dead * useMortality);
        mRS_control_ICH_manyCMBs = applyORsToMRS(mRS_control_ICH_noCMBs, OR_manyCMBs_adjICH_mRS_3to6 * useDisability, OR_manyCMBs_adjICH_mRS_dead * useMortality);
        
        mRS_0to2_expected = sum(mRS_control_ICH_noCMBs(1:3)) * pCMB_noCMBs_control_ICH + sum(mRS_control_ICH_someCMBs(1:3)) * pCMB_someCMBs_control_ICH + sum(mRS_control_ICH_manyCMBs(1:3)) * pCMB_manyCMBs_control_ICH;
        mRS_0to2_observed = sum(mRS_control_ICH(1:3));
        delta_mRS_0to2 = mRS_0to2_expected - mRS_0to2_observed;
        
        mRS_dead_expected = mRS_control_ICH_noCMBs(7) * pCMB_noCMBs_control_ICH + mRS_control_ICH_someCMBs(7) * pCMB_someCMBs_control_ICH + mRS_control_ICH_manyCMBs(7) * pCMB_manyCMBs_control_ICH;
        mRS_dead_observed = mRS_control_ICH(7);
        delta_dead = mRS_dead_expected - mRS_dead_observed;
    catch
        impactCMB_many_vs_no = NaN*[1 1 1]; impacttPA_yes_vs_no_noCMBs = NaN*[1 1 1]; impacttPA_yes_vs_no_manyCMBs = NaN*[1 1 1];
        % fprintf('      try catch. \n')
        return;
    end
    breakLoop = 1;
    if useDisability && abs(delta_mRS_0to2) > accuracySolver
        if delta_mRS_0to2 < 0
           upper_mRS_3to6 = OR_noCMBs_mRS_3to6; 
        else
           lower_mRS_3to6 = OR_noCMBs_mRS_3to6;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if useMortality && abs(delta_dead) > accuracySolver
        if delta_dead > 0
            upper_mRS_dead = OR_noCMBs_mRS_dead;
        else
            lower_mRS_dead = OR_noCMBs_mRS_dead;
        end
        cnt = cnt + 1;
        breakLoop = 0;
    end
    if breakLoop
        break;
    end
end
if v; fprintf(' completed in n = %d runs. Delta: %.6f / %.6f.\n', cnt, delta_mRS_0to2, delta_dead); end
clear cnt mRS_0to2_expected mRS_0to2_observed delta_mRS_0to2 mRS_dead_expected mRS_dead_observed delta_dead
clear upper_mRS_3to6 lower_mRS_3to6 upper_mRS_dead lower_mRS_dead
clear breakLoop
clear OR_noCMBs_mRS_3to6 OR_noCMBs_mRS_dead

% Averaging accros sICH categories
mRS_control_noCMBs = mRS_control_ICH_noCMBs .* pICH_ICH_control_noCMB + mRS_control_noICH_noCMBs .* pICH_noICH_control_noCMB;
mRS_control_someCMBs = mRS_control_ICH_someCMBs .* pICH_ICH_control_someCMB + mRS_control_noICH_someCMBs .* pICH_noICH_control_someCMB;
mRS_control_manyCMBs = mRS_control_ICH_manyCMBs .* pICH_ICH_control_manyCMB + mRS_control_noICH_manyCMBs .* pICH_noICH_control_manyCMB;


%% Main 4 - Outcomes
if v; fprintf('Calculating ORs for good outcome / mortality according to CMB-status ...'); end
OR_tPA_avCMBs_mrs_0to2 = sum(mRS_tPA(1:3)) / sum(mRS_tPA(4:7)) / sum(mRS_control(1:3)) * sum(mRS_control(4:7));
OR_tPA_noCMBs_mrs_0to2 = sum(mRS_tPA_noCMBs(1:3)) / sum(mRS_tPA_noCMBs(4:7)) / sum(mRS_control_noCMBs(1:3)) * sum(mRS_control_noCMBs(4:7));
OR_tPA_someCMBs_mrs_0to2 = sum(mRS_tPA_someCMBs(1:3)) / sum(mRS_tPA_someCMBs(4:7)) / sum(mRS_control_someCMBs(1:3)) * sum(mRS_control_someCMBs(4:7));
OR_tPA_manyCMBs_mrs_0to2 = sum(mRS_tPA_manyCMBs(1:3)) / sum(mRS_tPA_manyCMBs(4:7)) / sum(mRS_control_manyCMBs(1:3)) * sum(mRS_control_manyCMBs(4:7));

OR_tPA_avCMBs_mrs_dead = mRS_tPA(7) / sum(mRS_tPA(1:7)) / mRS_control(7) * sum(mRS_control(1:7));
OR_tPA_noCMBs_mrs_dead = mRS_tPA_noCMBs(7) / sum(mRS_tPA_noCMBs(1:7)) / mRS_control_noCMBs(7) * sum(mRS_control_noCMBs(1:7));
OR_tPA_someCMBs_mrs_dead = mRS_tPA_someCMBs(7) / sum(mRS_tPA_someCMBs(1:7)) / mRS_control_someCMBs(7) * sum(mRS_control_someCMBs(1:7));
OR_tPA_manyCMBs_mrs_dead = mRS_tPA_manyCMBs(7) / sum(mRS_tPA_manyCMBs(1:7)) / mRS_control_manyCMBs(7) * sum(mRS_control_manyCMBs(1:7));

% common odds ratio, remnant,  not longer used
% OR_tPA_avCMBs_mrs_cOR = commonOR(mRS_control, mRS_tPA);
% OR_tPA_noCMBs_mrs_cOR = commonOR(mRS_control_noCMBs, mRS_tPA_noCMBs);
% OR_tPA_someCMBs_mrs_cOR = commonOR(mRS_control_someCMBs, mRS_tPA_someCMBs);
% OR_tPA_manyCMBs_mrs_cOR = commonOR(mRS_control_manyCMBs, mRS_tPA_manyCMBs);

quali_tPA_avCMBs = sum(mRS_tPA .* qualityWeights);
quali_tPA_noCMBs = sum(mRS_tPA_noCMBs .* qualityWeights);
quali_tPA_someCMBs = sum(mRS_tPA_someCMBs .* qualityWeights);
quali_tPA_manyCMBs = sum(mRS_tPA_manyCMBs .* qualityWeights);

quali_control_avCMBs = sum(mRS_control .* qualityWeights);
quali_control_noCMBs = sum(mRS_control_noCMBs .* qualityWeights);
quali_control_someCMBs = sum(mRS_control_someCMBs .* qualityWeights);
quali_control_manyCMBs = sum(mRS_control_manyCMBs .* qualityWeights);

impactCMB_many_vs_no_disability = (OR_tPA_manyCMBs_mrs_0to2 - OR_tPA_noCMBs_mrs_0to2) / OR_tPA_noCMBs_mrs_0to2;
impactCMB_many_vs_no_mortality = (OR_tPA_manyCMBs_mrs_dead - OR_tPA_noCMBs_mrs_dead) / OR_tPA_noCMBs_mrs_dead;
% impactCMB_many_vs_no_cOR = (OR_tPA_manyCMBs_mrs_cOR - OR_tPA_noCMBs_mrs_cOR) / OR_tPA_noCMBs_mrs_cOR;
impactCMB_many_vs_no_qualies = ( (quali_tPA_manyCMBs - quali_control_manyCMBs) - (quali_tPA_noCMBs - quali_control_noCMBs) ) / (quali_tPA_noCMBs - quali_control_noCMBs);
if useDisability == 1 && useMortality == 0
    impactCMB_many_vs_no = impactCMB_many_vs_no_disability;
elseif useDisability == 0 && useMortality == 1
    impactCMB_many_vs_no = impactCMB_many_vs_no_mortality;
elseif useDisability == 1 && useMortality == 1
    impactCMB_many_vs_no = [impactCMB_many_vs_no_disability, impactCMB_many_vs_no_mortality, impactCMB_many_vs_no_qualies];
end

impacttPA_yes_vs_no_manyCMBs_disability = OR_tPA_manyCMBs_mrs_0to2;
impacttPA_yes_vs_no_manyCMBs_mortality = OR_tPA_manyCMBs_mrs_dead;
% impacttPA_yes_vs_no_manyCMBs_cOR = OR_tPA_manyCMBs_mrs_cOR;
impacttPA_yes_vs_no_manyCMBs_qualies = quali_tPA_manyCMBs - quali_control_manyCMBs;
if useDisability == 1 && useMortality == 0
    impacttPA_yes_vs_no_manyCMBs = impacttPA_yes_vs_no_manyCMBs_disability;
elseif useDisability == 0 && useMortality == 1
    impacttPA_yes_vs_no_manyCMBs = impacttPA_yes_vs_no_manyCMBs_mortality;
elseif useDisability == 1 && useMortality == 1
    impacttPA_yes_vs_no_manyCMBs = [impacttPA_yes_vs_no_manyCMBs_disability, impacttPA_yes_vs_no_manyCMBs_mortality, impacttPA_yes_vs_no_manyCMBs_qualies];
end

impacttPA_yes_vs_no_noCMBs_disability = OR_tPA_noCMBs_mrs_0to2;
impacttPA_yes_vs_no_noCMBs_mortality = OR_tPA_noCMBs_mrs_dead;
% impacttPA_yes_vs_no_noCMBs_cOR = OR_tPA_noCMBs_mrs_cOR;
impacttPA_yes_vs_no_noCMBs_qualies = quali_tPA_noCMBs - quali_control_noCMBs;
if useDisability == 1 && useMortality == 0
    impacttPA_yes_vs_no_noCMBs = impacttPA_yes_vs_no_noCMBs_disability;
elseif useDisability == 0 && useMortality == 1
    impacttPA_yes_vs_no_noCMBs = impacttPA_yes_vs_no_noCMBs_mortality;
elseif useDisability == 1 && useMortality == 1
    impacttPA_yes_vs_no_noCMBs = [impacttPA_yes_vs_no_noCMBs_disability, impacttPA_yes_vs_no_noCMBs_mortality, impacttPA_yes_vs_no_noCMBs_qualies];
end
if v; fprintf(' completed.\n'); end

end