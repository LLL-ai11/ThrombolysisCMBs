function [mRS_out, status] = baseline_mRS(age, NIHSS, delay, tPA_given, useDisability, useMortality, sensiFlag, univSensiAnalysis)
%baseline_mRS: Estimates mRS distribution according to clinical parameters and treatment (IVT yes/no)
%   Imput parameters
%   age: age, in years
%   NIHSS: NIHSS score: 3 - 24
%   delay: treatment delay to start of IVT, in minutes
%   tPA_given: 1 / 0 (yes / no)
%   useDisability: remnant flag, indicating if disability and/or mortality is considered in the transformation of mRS scores
%   useMortality: see above
%   sensiFlag: flag for probabilistic sensitivity analzsis (1 / 0)
%   univSensiAnalysis: array containing flags for univariate sensitivity analyses


%% Initialize
oddsRatio_tPA_0to1 = 0.0097 * (delay/60)^2 - 0.2372 * (delay/60) + 2.1077;   % Emberson et al.
oddsRatio_tPA_0to1_95CI = [-0.0189 * (delay/60)^2 + 0.0269 * (delay/60) + 1.3517, ...
                           0.0557 * (delay/60)^2 - 0.647 * (delay/60) + 3.1807];

oddsRatio_tPA_death = 1.11; oddsRatio_tPA_death_95CI = [0.99, 1.25];  % Emberson et al.
                                           
oddsRatio_age = 0.001298*age^2 - 0.2416*age + 11.61; % Knoflach et al;
oddsRatio_age_95CI = [max(0.0001, -0.04636*age + 4.288), ...
                       43.01*exp(-0.05293*age)];
oddsRatio_age_95CI = [min(oddsRatio_age_95CI(1), oddsRatio_age), max(oddsRatio_age_95CI(2), oddsRatio_age)];  % adjustment due to imperfect fit
                       
if sensiFlag
    r = rand();
    m = log(oddsRatio_tPA_0to1);    CI = log(oddsRatio_tPA_0to1_95CI);
    if r > 0.5; SD = (CI(2) - m) / 1.96;
    else;       SD = (m - CI(1)) / 1.96; end
    oddsRatio_tPA_0to1 = exp(norminv(r, m, SD));
    m = log(oddsRatio_tPA_death);    CI = log(oddsRatio_tPA_death_95CI);
    if r > 0.5; SD = (CI(2) - m) / 1.96;
    else;       SD = (m - CI(1)) / 1.96; end
    oddsRatio_tPA_death = exp(norminv(r, m, SD));
    
    r = rand();
    m = log(oddsRatio_age);    CI = log(oddsRatio_age_95CI);
    if r > 0.5; SD = (CI(2) - m) / 1.96;
    else;       SD = (m - CI(1)) / 1.96; end
    oddsRatio_age = exp(norminv(r, m, SD));
    clear r m CI
end
if univSensiAnalysis(2)
    oddsRatio_tPA_0to1 = oddsRatio_tPA_0to1 * (-0.0028*age + 1.198);  % Emberson et al.
end
oddsRatio_tPA_2to6 = 1/oddsRatio_tPA_0to1;

%% Main
% obtain functional outcome (mRS distribution )according to NIHSS score vis interpolation
NIHSS_array = [2.5, 7.5, 13, 18.5, 24];  % Whitely et al.
mrs_control_NIHSS =     [0.268	0.321	0.237	0.090	0.034	0.028	0.022; ...
                        0.183	0.247	0.178	0.146	0.087	0.087	0.072; ...
                        0.072	0.145	0.123	0.165	0.194	0.146	0.156; ...
                        0.025	0.057	0.061	0.125	0.227	0.213	0.292; ...
                        0.002	0.022	0.029	0.070	0.163	0.272	0.441];
f = find(NIHSS_array > NIHSS);
upper = f(1); lower = f(1) - 1;
mrs_control_NIHSS = mrs_control_NIHSS(lower,:) + (mrs_control_NIHSS(upper,:) - mrs_control_NIHSS(lower,:)) * (NIHSS - NIHSS_array(lower)) / (NIHSS_array(upper) - NIHSS_array(lower));
mrs_control_NIHSS = mrs_control_NIHSS / sum(mrs_control_NIHSS);
clear f upper lower

% adjust for age, (inverse OR as it relates to good outcome instead of disability)
mRS = applyORsToMRS(mrs_control_NIHSS, 1/oddsRatio_age, 0);

% adjust for tPA-treatment
if tPA_given
    mRS_tPA = applyORsToMRS(mRS, oddsRatio_tPA_2to6 * useDisability, oddsRatio_tPA_death * useMortality, 1);
    mRS_out = mRS_tPA;
else
   mRS_out = mRS;
end   
status = 1;
end