function cOR = commonOR(mRS_control,mRS_tPA)
%commonOR Summary of this function goes here
%   Detailed explanation goes here
y = zeros(6, 1);
for cut = 1 : 6
    y(cut) = sum(mRS_tPA(1:cut)) / sum(mRS_tPA(cut+1:end)) / sum(mRS_control(1:cut)) * sum(mRS_control(cut+1:end));
end
fit1 = fitlm(ones(size(y)),y,'y~1');
cOR = fit1.Coefficients.Estimate(1);
end

