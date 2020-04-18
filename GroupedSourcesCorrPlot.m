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
rScale = [0, 0.07];
varScale = [0, 0.03];
detectScale = [20, 65];
YScale = [corrSet(length(corrSet)) - 0.1, corrSet(1) + 0.1];
dim = length(snrSet);
figure

%Localization bias
for n = 1:dim
    subplot(3, dim, n)
    hold on
    for i = methodSet
        lol(i) = plot(corrSet, squeeze(mean_r(i, n, :)), 'Color', col(i,:), 'LineWidth', 3);
        scatter(corrSet, squeeze(mean_r(i, n, :)), ...
        150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
        ylim(rScale)
        xlim(YScale)
    end
    set(gca,'FontSize', fntSize, 'Xdir', 'reverse')
    xlabel('Correlation')
    ylabel('Meters, m')
end

%Spreading area
for n = 1:dim
    subplot(3, dim, dim + n)
    hold on
    for i = methodSet
        plot(corrSet, squeeze(mean_var(i, n, :)), 'Color', col(i,:), 'LineWidth', 3)
        scatter(corrSet, squeeze(mean_var(i, n, :)), ...
        150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
        ylim(varScale)
        xlim(YScale)
    end
    set(gca,'FontSize', fntSize, 'Xdir', 'reverse')
    xlabel('Correlation')
    ylabel('Meters, m')
end

%Detection ratio
for n = 1:dim
    subplot(3, dim, 2*dim + n)
    hold on
    for i = methodSet
        plot(corrSet,(1-(squeeze(detect(i, n, :))./Nmc))*100, 'Color', col(i,:), 'LineWidth', 3)
        scatter(corrSet, (1-(squeeze(detect(i, n, :))./Nmc))*100, ...
        150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
        ylim(detectScale)
        xlim(YScale)
    end
    set(gca,'FontSize', fntSize, 'Xdir', 'reverse')
    xlabel('Correlation')
    ylabel('Percents, %')
end

%legend(lol, "ReciPSIICOS", "wReciPSIICOS", "LCMV")

% Localization bias histogram for sources with 0.55 correlation
% ReciPSIICOS beamformer
% r(methodN, correlationN, ~, frac, simulationN)
%hist(sqrt(squeeze(r(1, 6, :, :))),100) 