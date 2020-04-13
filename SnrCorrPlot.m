%%%%%
% Aleksandra Kuznetsova,  Alexei Ossadtchi*, Grigoriy Mozgov
% *ossadtchi@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detect = sum(isnan(r), 5);
mean_r = mean(r, 5, 'omitnan');
mean_var = mean(var, 5, 'omitnan');

c = lines(7);
col = c([1, 6, 2, 4], :);
fntSize = 12;
rScale = [0, 0.12];
varScale = [0, 0.12];
detectScale = [60, 100];
YScale = [snrSet(1) - 0.25, snrSet(length(snrSet)) + 0.25];
dim = length(corrSet);
figure

%Localization bias
for n = 1:dim
    subplot(3, dim, n)
    hold on
    for i = methodSet
        plot(snrSet, squeeze(mean_r(i, :, n)), 'Color', col(i,:), 'LineWidth', 3)
        scatter(snrSet, squeeze(mean_r(i, :, n)), ...
        150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
        ylim(rScale)
        xlim(YScale)
    end
    %plot([4, 4], [0 0.07], '--', 'Color', 'k')
    set(gca,'FontSize', fntSize)
    xlabel('SNR')
    ylabel('Meters, m')
end

%Spreading area
for n = 1:dim
    subplot(3, dim, dim + n)
    hold on
    for i = methodSet
        plot(snrSet, squeeze(mean_var(i, :, n)), 'Color', col(i,:), 'LineWidth', 3)
        scatter(snrSet, squeeze(mean_var(i, :, n)), ...
        150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
        ylim(varScale)
        xlim(YScale)
    end
    set(gca,'FontSize', fntSize)
    %plot([4, 4], [0 0.03], '--', 'Color', 'k')
    xlabel('SNR')
    ylabel('Meters, m')
end

%Detection ratio
for n = 1:dim
    subplot(3, dim, 2*dim + n)
    hold on
    for i = methodSet
        plot(snrSet,(1-(squeeze(detect(i, :, n))./Nmc))*100, 'Color', col(i,:), 'LineWidth', 3)
        scatter(snrSet, (1-(squeeze(detect(i, :, n))./Nmc))*100, ...
        150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
        ylim(detectScale)
        xlim(YScale)
    end
    %plot([4, 4], [60 100], '--', 'Color', 'k')
    set(gca,'FontSize', fntSize)
    xlabel('SNR')
    ylabel('Percents, %')
end

% Localization bias histogram for sources with 0.55 correlation
% ReciPSIICOS beamformer
% r(methodN, correlationN, ~, frac, simulationN)
%hist(sqrt(squeeze(r(1, 6, :, :))),100) 