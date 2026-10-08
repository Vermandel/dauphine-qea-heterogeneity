% Replot Kaplan, Violante & Weidner (2014), Tables 6 and 7.
% Left: PSID estimates of MPC out of transitory income shocks (Table 6,
% baseline row; whiskers are reported bootstrap standard errors).
% Right: quarterly MPC out of a $500 transfer in their two-asset model
% (Table 7), by hand-to-mouth status and age.

figure('Color','w','Position',[100 100 920 390]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile
mpc = [0.243 0.301 0.127];
se  = [0.065 0.048 0.036];
b = bar(mpc,'FaceColor',[0.38 0.55 0.72]); hold on
errorbar(1:3,mpc,se,'k.','LineWidth',1.4,'CapSize',10);
set(gca,'XTick',1:3,'XTickLabel',{'Poor HtM','Wealthy HtM','Non-HtM'});
ylabel('MPC'); ylim([0 0.42]); grid on
title({'PSID estimates','transitory income shock'});

nexttile
age_mpc = [0.38 0.42 0.08; 0.30 0.42 0.01; 0.39 0.51 0.13];
bar(age_mpc,'grouped');
set(gca,'XTick',1:3,'XTickLabel',{'Age <= 40','Age 40--60','Age > 60'});
ylabel('Quarterly MPC'); ylim([0 0.60]); grid on
legend({'Poor HtM','Wealthy HtM','Non-HtM'},'Location','northwest');
title({'Two-asset life-cycle model','USD 500 unexpected transfer'});

title(tl,'Liquidity status, not net worth alone, predicts MPCs');
outdir=fullfile('..','figures');
if ~isfolder(outdir), mkdir(outdir); end
exportgraphics(gcf,fullfile(outdir,'kaplan2014_mpc.pdf'),'ContentType','vector');
