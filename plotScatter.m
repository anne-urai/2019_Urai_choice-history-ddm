
function [rho, tt, handles] = plotScatter(allresults, fld, legendWhere, doText)

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

% overall correlation
x = cat(1, allresults(:).(fld));
y = cat(1, allresults(:).criterionshift);
% show line
axis square;

% show lines to indicate origin
xlims = [min(x) max(x)];
ylims = [min(y) max(y)];
plot([0 0], ylims, 'color', [0.5 0.5 0.5], 'linewidth', 0.2);
plot(xlims, [0.5 0.5], 'color', [0.5 0.5 0.5], 'linewidth', 0.2); % if p(repeat), 0.5

for a = 1:length(allresults), % neutral last
    
    [rho, pval] = corr(allresults(a).(fld), allresults(a).criterionshift, ...
        'type', 'spearman', 'rows', 'complete');
    
    if pval < 0.05,
        
        % 18 SEPTEMBER 2018 - SWITCH TO PCA CORRELATION https://elifesciences.org/articles/00638
        b = deming(allresults(a).(fld)(~isnan(allresults(a).(fld))), allresults(a).criterionshift(~isnan(allresults(a).(fld))));
        %p = polyfit(allresults(a).(fld)(~isnan(allresults(a).(fld))), allresults(a).criterionshift(~isnan(allresults(a).(fld))), 1);
        xrangeextra = 0.15*range(allresults(a).(fld));
        xrange = linspace(min(allresults(a).(fld))- xrangeextra, ...
            max(allresults(a).(fld))+xrangeextra, 100);
        yrange = polyval(fliplr(b'), xrange);
        
        % NOW PLOT
        l = plot(xrange, yrange);
        l.Color = allresults(a).meancolor;
        l.LineWidth = 0.5;
        l.LineStyle = '-';
    end
    
    % PLOT ALL DATAPOINTS IN SPECIFIC COLOR
    s  = scatter(allresults(a).(fld), allresults(a).criterionshift, 7, allresults(a).scattercolor, allresults(a).marker);
	%set(s, 'markerfacecolor', allresults(a).scattercolor);
	handles{a} = s;
end

for a = 1:length(allresults), % neutral last
    % also add the group mean
    p = ploterr(nanmean(allresults(a).(fld)), nanmean(allresults(a).criterionshift), 2*nanstd(allresults(a).(fld)) ./ sqrt(length(allresults(a).(fld))), ...
        2*nanstd(allresults(a).criterionshift) ./ sqrt(length(allresults(a).criterionshift)), '.', 'abshhxy', 0);
    set(p(1), 'markersize', 0.1, 'color', allresults(a).meancolor); % tiny marker
    set(p(2), 'color', allresults(a).meancolor, 'linewidth', 1);
    set(p(3), 'color', allresults(a).meancolor, 'linewidth', 1);
end

axis tight; offsetAxes;

if doText,
    % PRINT THE CORRELATION COEFFICIENT
    txt = {sprintf('\\rho(%d) = %.3f', length(find(~isnan(y)))-2, rho) sprintf('p = %.3f', pval)};
    if pval < 0.001,
        txt = {sprintf('\\rho(%d) = %.3f', length(find(~isnan(y)))-2,rho) sprintf('p < 0.001')};
    end
    tt = text(min(get(gca, 'xlim')) + legendWhere*(range(get(gca, 'xlim'))), ...
        min(get(gca, 'ylim')) + 0.8*(range(get(gca, 'ylim'))), ...
        txt, 'fontsize', 5);
else
    tt = [];
end
set(gca, 'color', 'none');

end
