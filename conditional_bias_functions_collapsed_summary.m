function conditional_bias_functions_collapsed_summary

global colors
  thesecolors = {[0 0 0], [0.5 0.5 0.5],  colors(1, :), colors(2, :), mean(colors([1 2], :))};


load(sprintf('~/Data/serialHDDM/allds_cbfs_mediansplit.mat'));
periods    = {'fast', 'slow'};
models     = {'data', 'stimcoding_nohist', 'stimcoding_z_prevresp', 'stimcoding_dc_prevresp', 'stimcoding_dc_z_prevresp'};
modelnames = regexprep(regexprep(regexprep(models, 'stimcoding_', ''), '_', ' '), 'prevresp', '');
close all; 

for m = 2:length(models),
  subplot(4,8,m-1); hold on;

  for p = 1:2,
      bar(p, nanmean(allds.(periods{p})(:, m)), 'edgecolor', 'none', ...
        'facecolor', thesecolors{m}, 'basevalue', 0.5, 'barwidth', 0.3);
      errorbar(p, nanmean(allds.(periods{p})(:, 1)),  ...
      1.96* nanstd(allds.(periods{p})(:, 1)) ./ sqrt(size(allds.fast, 1)), ...
      'ko', 'capsize', 0, 'markerfacecolor', 'w');
  end

  if m == 2, ylabel('P(bias)'); end
  set(gca, 'xtick', 1:2, 'xticklabel', periods);
  title(sprintf('Model %s', modelnames{m}));
  axis tight; disp(get(gca, 'ylim'));
  ylim([0.5 0.55]);
  offsetAxes;
  if m > 2, set(gca, 'yticklabel', []); end
  set(gca, 'ycolor', 'k', 'xcolor', 'k');
end
end

