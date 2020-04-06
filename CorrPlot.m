detect = sum(isnan(r), 4);
mean_r = mean(r, 4, 'omitnan');
mean_var = mean(var, 4, 'omitnan');

c = lines(7);
col = c([1, 6, 2, 4], :);
fntSize = 12;
rScale = [0, 0.07];
varScale = [0, 0.03];
detectScale = [60, 100];
corrScale = [0, 1];

%Localization bias
figure
subplot(3, 1, 1)
hold on
for i = 1:2:3
    plot(corrSet, squeeze(mean_r(i,1,:)), 'Color', col(i,:), 'LineWidth', 3)
    scatter(corrSet, squeeze(mean_r(i,1,:)), ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim(rScale)
    xlim(corrScale)
end
%plot([4, 4], [0 0.07], '--', 'Color', 'k')
set(gca,'FontSize', fntSize)
xlabel('Correlation between sources')
ylabel('Meters, m')

%Spreading area
subplot(3, 1, 2)
hold on
for i = 1:2:3
    plot(corrSet, squeeze(mean_var(i,1,:)), 'Color', col(i,:), 'LineWidth', 3)
    scatter(corrSet, squeeze(mean_var(i,1,:)), ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim(varScale)
    xlim(corrScale)
end
set(gca,'FontSize', fntSize)
%plot([4, 4], [0 0.03], '--', 'Color', 'k')
xlabel('Correlation between sources')
ylabel('Meters, m')

%Detection ratio
subplot(3, 1, 3)
hold on
for i = 1:2:3
    plot(corrSet,(1-(squeeze(detect(i,1,:))./Nmc))*100, 'Color', col(i,:), 'LineWidth', 3)
    scatter(corrSet, (1-(squeeze(detect(i,1,:))./Nmc))*100, ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim(detectScale)
    xlim(corrScale)
end
%plot([4, 4], [60 100], '--', 'Color', 'k')
set(gca,'FontSize', fntSize)
xlabel('Correlation between sources')
ylabel('Percents, %')

% Localization bias histogram for sources with 0.55 correlation
% ReciPSIICOS beamformer
% r(methodN, correlationN, ~, frac, simulationN)
hist(sqrt(squeeze(r(1, 6, :, :))),100) 