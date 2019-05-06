%% CMB_tPA_selection
%   c Ludwig Schlemm, 2019

%% Definitions
clearvars
age_array = 50:10:90; NIHSS_array = 5:5:20;
useDisability = 1; useMortality = 1;
nSensi = 1;
delay_array = 30:30:270; nDelay = length(delay_array);
univSensiAnalysis = [1, 1]; % first entry: no effect of >10 CMBs on sICH in patients without IVT
                            % second entry: modeling age-dependent treatment effect of IVT

%% Main
metaResults = cell(length(age_array), length(NIHSS_array));
for a = 1: length(age_array)
    fprintf('a: %d of %d\n', a, length(age_array))
    age = age_array(a);
    for N = 1: length(NIHSS_array)
        fprintf('   N: %d of %d\n', N, length(NIHSS_array))
        NIHSS = NIHSS_array(N);
        
        %Initialize
        results = cell(3,2);                        % 1: favourable outcome; 2: mortality; 3: combined  x 1: no, 2: many, 3: manyVSno
        for k = 1: 3
            for l = 1: 3
                results{k,l} = zeros(nSensi, nDelay);
            end
        end
        clear k l
        
        % Main
        for ii = 1: nDelay
            fprintf('      ii: %d of %d\n', ii, nDelay)
            delay = delay_array(ii);
            for jj = 1: nSensi
                %                 fprintf('         jj: %d of %d\n', jj, nSensi)
                mRS_control = baseline_mRS(age, NIHSS, delay, 0, useDisability, useMortality, jj>1, univSensiAnalysis);
                mRS_tPA = baseline_mRS(age, NIHSS, delay, 1, useDisability, useMortality, jj>1, univSensiAnalysis);
                [impactCMB_many_vs_no, impacttPA_yes_vs_no_noCMBs, impacttPA_yes_vs_no_manyCMBs] = CMB_tPA_selection_v2(useDisability, useMortality, ...
                    mRS_control, mRS_tPA, jj>1, univSensiAnalysis, 0, a);
                
                results{1,1}(jj, ii) = impacttPA_yes_vs_no_noCMBs(1);
                results{2,1}(jj, ii) = impacttPA_yes_vs_no_noCMBs(2);
                results{3,1}(jj, ii) = impacttPA_yes_vs_no_noCMBs(3);
                
                results{1,2}(jj, ii) = impacttPA_yes_vs_no_manyCMBs(1);
                results{2,2}(jj, ii) = impacttPA_yes_vs_no_manyCMBs(2);
                results{3,2}(jj, ii) = impacttPA_yes_vs_no_manyCMBs(3);
                
                results{1,3}(jj, ii) = impactCMB_many_vs_no(1);
                results{2,3}(jj, ii) = impactCMB_many_vs_no(2);
                results{3,3}(jj, ii) = impactCMB_many_vs_no(3);
            end
        end
        metaResults{a, N} = results;
    end
end
save(['..' filesep 'output' filesep 'output_s1_' num2str(univSensiAnalysis(1)) '_s2_' num2str(univSensiAnalysis(2)) '_' datestr(now, 'yyyy-mm-dd HH-MM_SS')], 'metaResults', 'age_array', 'NIHSS_array', 'delay_array', 'nSensi', 'univSensiAnalysis')


%% Plots - 1
pCI = .95;  % width of credible interval
panels = [2,1; 2, length(NIHSS_array)-1; length(age_array)-1, 1; length(age_array)-1, length(NIHSS_array)-1];
XTickArray = 60:90:240;
YTickArray1 = .5:.5:3;
YTickArray2 = 0:0.1:0.1;
styleArray = {'-', ':'};

% CI_colors = colormap('lines');
% CI_colors = CI_colors(1:3,:) .* [.8 .8 .8; 1 1 1; 1 1 1];
CI_colors = colormap('jet');
CI_colors = CI_colors([1,60, 60],:);

yLabels = {'OR mRS \leq 2', 'OR death', 'Net benefit'};
figure('Color', [1 1 1], 'Units', 'Centimeter', 'Position', [4 4 17 15]);
t = tight_subplot(3,size(panels, 1),[.01 .02], [.1 .06], [.1 .02]);

relDifffavOutcome = zeros(size(panels, 1), nSensi, length(delay_array));
relDiffMortality = zeros(size(panels, 1), nSensi, length(delay_array));
maxORfavOutcome = 0; minORfavOutcome = 100;
maxORmortality = 0; minORmortality = 100;

for p = 1: size(panels, 1)
    results = metaResults{panels(p, 1), panels(p, 2)};
    relDifffavOutcome(p, :, :) = ( results{1,2}(:,:) - results{1,1}(:,:) ) ./ results{1,1}(:,:);
    relDiffMortality(p, :, :) = ( results{2,2}(:,:) - results{2,1}(:,:) ) ./ results{2,1}(:,:);
    
    maxORfavOutcome = max([results{1,2}(1,:), maxORfavOutcome]); minORfavOutcome = min([results{1,2}(1,:), minORfavOutcome]);
    maxORmortality = max([results{2,2}(1,:), maxORmortality]); minORmortality = min([results{2,2}(1,:), minORmortality]);
    
    for s = 1 : 3      % s: disability, 2: mortality, 3: combined
        for ii = 2:-1:1 % no CMBs / many CMBs
            axes(t( (s-1)*size(panels, 1) + p ))
            x = delay_array;
            m = results{s,ii}(1,:);
            results_s = sort(results{s,ii});
            lower = prctile(results_s, 100*(1 - pCI)/2);
            upper = prctile(results_s, 100-100*(1 - pCI)/2);
            if nSensi > 1
                x2 = [x, fliplr(x)]; inBetween = [upper, fliplr(lower)];
                fill(x2, inBetween, CI_colors(ii, :), 'FaceAlpha', .2, 'EdgeColor', CI_colors(ii, :), 'EdgeAlpha', 0);
                hold on;
            end
            if s <= 2 && ii == 1
                plot([x(1) x(end)], [1 1], '-', 'LineWidth', 1, 'Color', .5* [1 1 1]);
            elseif s == 3 && ii == 1
                plot([x(1) x(end)], [0 0], '-', 'LineWidth', 1, 'Color', .5* [1 1 1]);
            end
            hold on;
            plot(x, m, 'Color', CI_colors(ii, :), 'LineWidth', 2, 'LineStyle', styleArray{ii});
            hold on
            if p == 4 && s == 2
                xLeg = [100 150];
                yLeg = [1.8 1.7];
                textLeg = {'\leq10 CMBs', ' >10 CMBs'};
                plot(xLeg, yLeg(ii)*[1 1], 'Color', CI_colors(ii, :), 'LineWidth', 4, 'LineStyle', styleArray{ii});
                text(xLeg(2)*1.05, yLeg(ii), textLeg(ii), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                    'FontSize', 8);
            end
            if p == 1
                ylabel(yLabels{s})
                if s < 3
                    set(gca, 'YTick', YTickArray1)
                else
                    set(gca, 'YTick', YTickArray2)
                end
            else
                set(gca, 'YTick', [])
            end
            xlim([min(x)-10,max(x)])
            ytickformat('%.2f')
            if s < 3
                ylim([0.8 1.9 ]);
            else
                ylim([-.08 .14]);
            end
            if s == 3
                set(gca, 'XTick', XTickArray);
                set(gca, 'XTickLabel', XTickArray)
                if p == 3
                    xl = xlabel('Treatment delay [min]');
                    xlp = get(xl, 'Position');
                    set(xl, 'Position', [0, xlp(2:3)])
                end
            else
                set(gca, 'XTick', XTickArray);
                set(gca, 'XTickLabels', []);
            end
            box off
            if s == 1 && ii == 2
                fprintf('Lower bound (disability) crosses 1 at delay = %d\n', delay_array(find(lower < 1, 1))) ;
            end
            if s == 3 && ii == 2
                beginnMeanNegative = delay_array(find(m < 0, 1));
                beginnUpperNegative = delay_array(find(upper < 0, 1));
                plot([beginnMeanNegative - (beginnMeanNegative>delay_array(1))*0.5*(delay_array(2) - delay_array(1)), delay_array(end)], -0.07*[1 1], 'Color', [154, 29, 17]/255, 'Linewidth', 4);
                plot([delay_array(1), beginnMeanNegative - (beginnMeanNegative>delay_array(1))*0.5*(delay_array(2) - delay_array(1))], -0.07*[1 1], 'Color', [169, 210, 137]/255, 'Linewidth', 4);
                fprintf('Upper bound (net effect) crosses 0 at delay = %d\n', beginnUpperNegative);
                fprintf('Mean (net effect) crosses 0 at delay = %d\n', beginnMeanNegative);
                clear beginnMeanNegative beginnUpperNegative
            end
        end
        if s == 1
            text(150, 1.95, sprintf('Age: %d, NIHSS: %d', age_array(panels(p, 1)), NIHSS_array(panels(p, 2))), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
        end
    end
end
fprintf('\nRel. Difference of fav. outcome: %.0f%% - %.0f%%\n', min(min(min(relDifffavOutcome(:,1,:))))*100, max(max(max(relDifffavOutcome(:,1,:))))*100)
fprintf('Range of OR fav outcome: %.2f - %.2f\n', minORfavOutcome, maxORfavOutcome)

fprintf('\nRel. Difference of mortality: %.0f%% - %.0f%%\n', min(min(min(relDiffMortality(:,1,:))))*100, max(max(max(relDiffMortality(:,1,:))))*100)
fprintf('Range of OR mortality: %.2f - %.2f\n', minORmortality, maxORmortality)
% export_fig(['..' filesep 'Figures' filesep 'Figure_1' filesep 'Figure_1_v1_s1_' num2str(univSensiAnalysis(1)) '_s2_' num2str(univSensiAnalysis(2))], '-nocrop', '-tif', '-r600')


%% Plots - 2
yLabelTime = 'Maximum time to obtain MRI [min]';
pCI = .95;
pCMB_manyCMBs_array = [0.00372	0.00766	0.01385	0.02227	0.03294];  % corresponding to age groups 50:10:90
pCMB_manyCMBsUpper_array = [0.00482	0.00884	0.01653	0.02788	0.04288];  % corresponding to age groups 50:10:90

CI_colors2 = colorgrad(10000, 'blue_up');
CI_colors2(1,:) = .9*[1 1 1];

excludeOldAges = 0;
% figure('Color', [1 1 1], 'Units', 'Centimeter', 'Position', [4 4 7.272 11]);  % 7, 5.23
% t = tight_subplot(length(age_array)-excludeOldAges,1,[.01 .12], [.25 .05], [.05 .04]);
figure('Color', [1 1 1], 'Units', 'Centimeter', 'Position', [4 4 9 11]);  % 7, 5.23
t = tight_subplot(length(age_array)-excludeOldAges,2,[.01 .12], [.25 .05], [.15 .04]);

cMax = 25;
if all(univSensiAnalysis == [0, 0])
   excelColumns = {'E', 'F'};
elseif all(univSensiAnalysis == [1, 0])
   excelColumns = {'G', 'H'};
elseif all(univSensiAnalysis == [0, 1])
   excelColumns = {'I', 'J'};
elseif all(univSensiAnalysis == [1, 1])
   excelColumns = {'K', 'L'};
end
for a = 1: length(age_array)-excludeOldAges
    plotImage = zeros(length(NIHSS_array)+1, length(delay_array)+1, 2);
    for N = 1: length(NIHSS_array)
        results = metaResults{a, N};
        slopeCombinedNoCMBs = ( results{3,1}(1,1) - results{3,1}(1,end) ) / ( delay_array(end) - delay_array(1) );
        x = delay_array;
        meanHarm = results{3,2}(1,:);
        results_s = sort(results{3,2});
        lowerHarm = prctile(results_s, 100*(1 - pCI)/2);
        upperHarm = prctile(results_s, 100-100*(1 - pCI)/2);
        meanTimeToMRI = -(pCMB_manyCMBs_array(a) .* meanHarm) ./ ((1 - pCMB_manyCMBs_array(a)) .* slopeCombinedNoCMBs);
%         upperTimeToMRI = -(pCMB_manyCMBsUpper_array(a) .* upperHarm) ./ ((1 - pCMB_manyCMBsUpper_array(a)) .* slopeCombinedNoCMBs);
        lowerTimeToMRI = -(pCMB_manyCMBsUpper_array(a) .* lowerHarm) ./ ((1 - pCMB_manyCMBsUpper_array(a)) .* slopeCombinedNoCMBs);
        plotImage(N, 1:end-1, 1) = max(0, meanTimeToMRI); plotImage(N, 1:end-1, 2) = max(0, lowerTimeToMRI);
        
        % display info for univ. sensi. analysis figure and save to excel
        displayAge = [60, 80]; displayNIHSS = [5, 15]; displayDelay = [30, 240];
        if any(age_array(a) == displayAge) && any(NIHSS_array(N) == displayNIHSS)
%            fprintf('   Info displ. sensi figure:\n') 
           for l = 1: 2
               pe = plotImage(N, delay_array == displayDelay(l) , 1);
               upBo = plotImage(N, delay_array == displayDelay(l) , 2);
%                fprintf('   Age:  %.0f, NIHSS: %.0f, delay: %.0f, Point estimate: %.2f \n', age_array(a), NIHSS_array(N), displayDelay(l), pe);
%                fprintf('   Age: %.0f, NIHSS: %.0f, delay: %.0f, Upper 95%% CI: %.2f \n', age_array(a), NIHSS_array(N), displayDelay(l), upBo);
               if 0
                   rowIndex = 8+ 4*(find(age_array(a) == displayAge)-1) + 2*(find(NIHSS_array(N) == displayNIHSS)-1) + l;
                   xlswrite(['..' filesep 'Figures' filesep 'Figure_3' filesep 'Figure_3.xlsx'], pe, 'Sheet1', [excelColumns{1}, num2str(rowIndex)])
                   xlswrite(['..' filesep 'Figures' filesep 'Figure_3' filesep 'Figure_3.xlsx'], upBo, 'Sheet1', [excelColumns{2}, num2str(rowIndex)])
               end
           
           end
           fprintf('\n')
           clear l pe upBo rowIndex
        end
        % end display info for univ. sensi. analysis figure and save to excel
    end
    fprintf('\nRange meanTimeToMRI for age: %.0f: %.2f - %.2f. \n', age_array(a), min(min(plotImage(:,:,1))), max(max(plotImage(:,:,1))));
    fprintf('Max upper bound timeToMRI for age: %.0f: %.2f. \n', age_array(a), max(max(plotImage(:,:,2))));

    plotImage(end, end, 1) = cMax;
    plotImage(end, end, 2) = cMax;
    
    for k = 1:2
        axes(t(2*((length(age_array)-excludeOldAges - a + 1)-1) + k))
%         axes(t( length(age_array)-excludeOldAges - a + 1))
        imagesc(plotImage(:,:,k))
        colormap(CI_colors2)
        hold on
        for i = 1:length(delay_array);  plot([i-.5,i-.5],[.5,length(NIHSS_array)+.5],'k-', 'LineWidth', .3, 'Color', 0*[1 1 1]);       end
        for i = 1:length(NIHSS_array);    plot([.5,length(delay_array)+.5], [i-.5,i-.5],'k-', 'LineWidth', .3, 'Color', 0* [1 1 1]);    end
        
        % markers on empty point estimates
        if k == 1
            for i = 1:length(delay_array)
                for j = 1:length(NIHSS_array)
                   if plotImage(j,i,1) == 0
                      plot(i,j, 'x', 'Color', .5*[1 1 1]); 
                      hold on;
                   end
                end
            end
            clear i j
        end
        
        set(gca, 'TickLength', [0 0])
        axis('image'); axis('xy')
        xlim([0,length(delay_array)]+.5); ylim([0, length(NIHSS_array)]+.5);
        if a == 1
            set(gca, 'XTick', 2: 2: length(delay_array))
            set(gca, 'XTickLabel', delay_array(get(gca, 'XTick')))
            xlabel('Treatment delay [min]')
        else
            set(gca, 'XTick', [])
        end
        if k == 1
            set(gca, 'YTick', 1: length(NIHSS_array))
            set(gca, 'YTickLabel', NIHSS_array)
            if a == 3; ylabel('NIHSS'); end
        else
            set(gca, 'YTick', [])
        end
        
        s = {'Point estimate', 'Upper bound of 95% CI'};
        if a == length(age_array)-excludeOldAges
           text(length(delay_array)/2+.5, length(NIHSS_array)+1.5, s{k}, 'HorizontalAlignment', 'center') 
        end
        
        if k == 1
           text(length(delay_array)+2.2, length(NIHSS_array)/2+.5, sprintf('Age\n%d', age_array(a)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        end
        
        if k == 1 && a == 1
            cb = colorbar('Location', 'South');
%             cb.Limits = [cMax/10000, cMax];
            cb.AxisLocation = 'out';
            cb.TickDirection = 'out';
%             cb.XTick = [cMax/10000 5:5:25];
%             cb.XTickLabel = [0 cb.XTick(2:end)];
            cb.Position = [.2, .06, .62/1, .03];
            BarPos = cb.Position;
            yl = ylabel(cb, yLabelTime);
            pos = yl.Position;
            set(yl, 'Position', [pos(1) 2.5 0])            
        end
    end
end
export_fig(['..' filesep 'Figures' filesep 'Figure_2' filesep 'Figure_2_v2_s1_' num2str(univSensiAnalysis(1)) '_s2_' num2str(univSensiAnalysis(2))], '-nocrop', '-tif', '-r600')



return;
%% temporary plots - mRS distributions
n_array = 3:20;
age_array_fig = [50 80];
XTickArrayMRS = 0.2:.2:.8;
figure('Color', [1 1 1], 'Units', 'Centimeter', 'Position', [4 1 18 9]);
t = tight_subplot(length(n_array),2,[.005 .01], [.15 .08], [.1 .02]);
c = colorgrad(7);
for ii = 1: length(n_array)
    NIHSS = n_array(ii);
    for a = 1: length(age_array_fig)
        age = age_array_fig(a);
        axes(t(2*(ii-1) + a))
        mRS_control = baseline_mRS(age, NIHSS, 30, 0, 1, 1, 0, [0, 0]);
        bh = barh([1, NaN], [mRS_control; [1 1 1 1 1 1 1]*NaN], 'stacked');
        for jj = 1: 7
            bh(jj).FaceColor = c(jj,:);
        end
        xlim([0, 1])
        ylim([0.5 1.5])
        set(gca, 'TickLength', [0, 0])
        box off
        if a == 1
            set(gca, 'YTick', 1); set(gca, 'YTickLabel', NIHSS);
        else
            set(gca, 'YTick', 1); set(gca, 'YTickLabel', []);
        end
        if a == 1 && ii == ceil(length(n_array)/2)
            ylabel('NIHSS')
        end
        if ii == 1
            text(.5, 1.8, sprintf('Age: %d years', age), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
        if ii < length(n_array)
            set(gca, 'XTick', []);
            set(gca, 'XTickLabel', [])
            set(gca, 'XColor', [1 1 1])
        else
            set(gca, 'XTick', XTickArrayMRS);
            set(gca, 'XTickLabel', XTickArrayMRS)
            xlabel('mRS [%]')
        end
    end
end
% export_fig(['..' filesep 'Figures' filesep 'Figure_mRS_distribution_v1_'], '-nocrop', '-tif', '-r300')
