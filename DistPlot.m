detect = sum(isnan(r), 5);
mean_r = mean(r, 5, 'omitnan');
mean_var = mean(var, 5, 'omitnan');

c = lines(7);
col = c([1, 6, 2, 4], :);
fntSize = 12;
rScale = [0, 0.07];
varScale = [0, 0.03];
detectScale = [60, 100];
dScale = [dSet(1) - 0.005, dSet(length(dSet)) + 0.005];

%Localization bias for synchronous sources
figure
subplot(3, 2, 1)
hold on
for i = methodSet
    plot(dSet, squeeze(mean_r(i,1,:)), 'Color', col(i,:), 'LineWidth', 3)
    scatter(dSet, squeeze(mean_r(i,1,:)), ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim(rScale)
    xlim(dScale)
end
%plot([4, 4], [0 0.07], '--', 'Color', 'k')
set(gca,'FontSize', fntSize)
xlabel('Sources distance')
ylabel('Meters, m')

%Spreading area for synchronous sources
subplot(3,2,3)
hold on
for i = methodSet
    plot(dSet, squeeze(mean_var(i,1,:)), 'Color', col(i,:), 'LineWidth', 3)
    scatter(dSet, squeeze(mean_var(i,1,:)), ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim(varScale)
    xlim(dScale)
end
set(gca,'FontSize', fntSize)
%plot([4, 4], [0 0.03], '--', 'Color', 'k')
xlabel('Sources distance')
ylabel('Meters, m')

%Detection ratio for synchronous sources
subplot(3, 2, 5)
hold on
for i = methodSet
    plot(dSet,(1-(squeeze(detect(i,1,:))./Nmc))*100, 'Color', col(i,:), 'LineWidth', 3)
    scatter(dSet, (1-(squeeze(detect(i,1,:))./Nmc))*100, ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim(detectScale)
    xlim(dScale)
end
%plot([4, 4], [60 100], '--', 'Color', 'k')
set(gca,'FontSize', fntSize)
xlabel('Sources distance')
ylabel('Percents, %')


%Localization bias for asynchronous sources
subplot(3, 2, 2)
hold on
for i = methodSet
    plot(dSet, squeeze(mean_r(i,2,:)), 'Color', col(i,:), 'LineWidth', 3)
    scatter(dSet, squeeze(mean_r(i,2,:)), ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim(rScale)
    xlim(dScale)
end
%plot([4, 4], [0 0.05], '--', 'Color', 'k')
set(gca,'FontSize', fntSize)
xlabel('Sources distance')
ylabel('Meters, m')

%Spreading area for asynchronous sources
subplot(3, 2, 4)
hold on
for i = methodSet
    plot(dSet, squeeze(mean_var(i,2,:)), 'Color', col(i,:), 'LineWidth', 3)
    scatter(dSet, squeeze(mean_var(i,2,:)), ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim(varScale)
    xlim(dScale)
end
set(gca,'FontSize', fntSize)
%plot([4, 4], [0 0.03], '--', 'Color', 'k')
xlabel('Sources distance')
ylabel('Meters, m')

%Detection ratio for asynchronous sources
subplot(3, 2, 6)
hold on
for i = methodSet
    plot(dSet,(1-(squeeze(detect(i,2,:))./Nmc))*100, 'Color', col(i,:), 'LineWidth', 3)
    scatter(dSet, (1-(squeeze(detect(i,2,:))./Nmc))*100, ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim(detectScale)
    xlim(dScale)
end
%plot([4, 4], [60 100], '--', 'Color', 'k')
set(gca,'FontSize', fntSize)
xlabel('Sources distance')
ylabel('Percents, %')

% Localization bias histogram for 4cm distance between synchronous sources of
% ReciPSIICOS beamformer
% r(methodN, synchN, distanceN, ~, frac, simulationN)
hist(sqrt(squeeze(r(1, 1, 8, :, :))),100) 