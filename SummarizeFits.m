clear all

dataDir = './AfterFitData/clean_contaminated';
outDir  = './AfterFitData/clean_contaminated/Summary';

% create a folder for the output files
if 7 ~= exist(outDir,'dir')
    mkdir(outDir);
end
files = dir(strcat(dataDir,'/*.mat'));


temp2color_idx = [60, 1;70, 2; 80, 3;90, 4;100, 5];

% Extract D, D_err, Temp, C0
D       = zeros(length(files),1);
D_err   = zeros(length(files),1);
TempC   = zeros(length(files),1);
Cs0      = zeros(length(files),1);
Cs0_err  = zeros(length(files),1);
C_max   = zeros(length(files),1);
dfname   = cell(length(files),1);
for k=1:length(files)
    fname = files(k).name;
    [filepath,name,extension] = fileparts(fname);
    mfile = matfile(fullfile(dataDir,fname));
    D(k)        = mfile.D_fit;
    D_err(k)    = mfile.D_err;
    TempC(k)    = mfile.tempC;
    Cs0(k)       = mfile.C0_fit;
    Cs0_err(k)   = mfile.C0_err;
%     C_t         = mfile.C_t;
    depth       = mfile.depth_um;
    idx         = find(depth(:,end) == mfile.thickness);
    C_max(k)   =  mfile.C_t(idx,1000);
    dfname{k}   = string(fname);
end

data = table(TempC,D,D_err,Cs0,Cs0_err,C_max,dfname);
data = sortrows(data,{'TempC'},{'ascend'});

c = jet(length(temp2color_idx));

mks = ['o','s','d','^','p','v','h','>'];

FigureVfb = figure;
lbls = cell(1,height(data));
chi2 = zeros(height(data),1);
chi2p = zeros(height(data),1);

maxTime = -1;

hold on
for i=1:height(data)
    filei = fullfile(dataDir,string(data.dfname(i)));
    mfile = matfile(filei);
    time_exp = mfile.time_exp;
    maxTime  = max(maxTime,max(time_exp));
    Vfb_exp  = mfile.Vsh;
    Vfb_std  = mfile.VfbStd;
    Vfb_fit  = mfile.ypred;
	[tbl,chi2,chi2p] = crosstab(Vfb_exp,Vfb_fit);
    delta    = mfile.delta;
    lower = Vfb_fit - delta;
    upper = Vfb_fit + delta;
    clidx = find(temp2color_idx(:,1) == data.TempC(i));
    lbls{i} = sprintf('T = %d °C Experimental',data.TempC(i));
    errorbar(time_exp./3600,Vfb_exp,Vfb_std,'Marker',mks(clidx),'Color',c(clidx,:),'LineWidth',1.5,'LineStyle','none');
    plot(time_exp./3600,Vfb_fit,'Color',c(clidx,:),'LineWidth',1.5,'HandleVisibility','off');
%     plot(time_exp./3600,[lower';upper'],'--','Color',[0.5,0.5,0.5],'LineWidth',1.5,'HandleVisibility','off')
end
hold off

xlabel('Time (hr)','FontSize',14,'FontWeight','bold')
ylabel('Flatband Shift (V)','FontSize',14,'FontWeight','bold')
leg = legend(lbls);
set(gca, 'FontSize', 14)
set(gca,'XMinorTick','on','YMinorTick','on')
legend('boxoff')
leg.Location = 'NorthEast';
leg.FontSize = 12;
xlim([0,maxTime/3600]);
ylim([-2.5,1]);
box on
ax = gca;
ax.LineWidth = 1.5;



% [tmp,ind]=sortrows(data.TempC');
% data=data(ind);
kB = 8.6173303E-5; % eV/K


logD        = log(data.D);
logD_err    = (data.D_err);
logCmax     = log(data.C_max);
weights     = logD_err.^(-2);
% weights     = weights./max(weights);
invT        = 1./(data.TempC+273.15);
[fitD,gof,fout] = fit(invT,logD,'poly1','Weights',weights)

ci_D        = confint(fitD);
pred_D      = predint(fitD,invT,0.95,'functional','on')

logD_err_f  = mean(min(ci_D(:,2)-fitD.p2));
slope_err   = mean(abs(ci_D(:,1)-fitD.p1));

s2 = var(log(data.D) - fitD(invT));
Cov=inv(fout.Jacobian'*fout.Jacobian)*s2;

D0          = exp(fitD.p2)
D0_err      = exp(logD_err_f)

EaD         = -fitD.p1*kB
EaD_err     = slope_err*kB

% Fit the maximum concentration in SiNx

fitC        = fit(invT,logCmax,'poly1')

ci_C        = confint(fitC);
pred_C      = predint(fitC,invT,0.95,'functional','on')

logC_err_f  = max(abs(ci_C(:,2)-fitC.p2));
slope_err   = max(abs(ci_C(:,1)-fitC.p1));

Cs0          = exp(fitC.p2)
Cs0_err      = exp(logC_err_f)

EaC         = -fitC.p1*kB
EaC_err     = slope_err*kB


% Arrhenius plot for D

FigureArrheniusD = figure;
hold on
errorbar(1000*invT,(data.D),(data.D_err),'bo','LineWidth',1.5);
plot(1000.*invT,exp(fitD(invT)),'b-','LineWidth',1.5);
plot(1000*invT,exp(pred_D),'b--','LineWidth',1.5);
hold off
set(gca, 'YScale', 'log')
% ylim([5E-15, 1E-9]);
xlabel('1000/T (K)','FontSize',14,'FontWeight','bold')
ylabel('D (cm^2/s)','FontSize',14,'FontWeight','bold')

d0_txt = sprintf('D_0 = %.2e ± %.2e cm^2/s\nE_a = %.2e ± %.2e eV',...
    D0,D0_err,EaD,EaD_err);

p1t_txt_D = text(0.025,0.97,d0_txt,'Units','normalized',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top');


% leg = legend('V_{fb} data','Fit','Prediction intervals','Location','Northeast');
% legend('boxoff')
% leg.FontSize = 12;

box on
ax = gca;
ax.LineWidth = 1.5;
set(gca, 'FontSize', 14)
set(gca,'XMinorTick','on','YMinorTick','on')

% Arrhenius plot for C

FigureArrheniusC = figure;
hold on
plot(1000*invT,(data.C_max),'o','LineWidth',1.5);
plot(1000.*invT,exp(fitC(invT)),'LineWidth',1.5);
plot(1000*invT,exp(pred_C),'r--','LineWidth',1.5);
hold off
set(gca, 'YScale', 'log')

xlabel('1000/T (K)','FontSize',14,'FontWeight','bold')
ylabel('C_{ss} (cm^{-3})','FontSize',14,'FontWeight','bold')

cs_txt = sprintf('C_{ss0} = %.2e ± %.2e cm^{-3}\nE_a = %.2e ± %.2e eV',...
    Cs0,Cs0_err,EaC,EaC_err);

p1t_txt_C = text(0.025,0.97,cs_txt,'Units','normalized',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top');


leg = legend('V_{fb} data','Fit','Prediction intervals','Location','Northeast');
legend('boxoff')
leg.FontSize = 12;

box on
ax = gca;
ax.LineWidth = 1.5;
set(gca, 'FontSize', 14)
set(gca,'XMinorTick','on','YMinorTick','on')


savefig(FigureArrheniusD,fullfile(outDir,'Fits_Summary_D.fig'))
print(FigureArrheniusD,fullfile(outDir,'Fits_Summary_D.png'),'-dpng','-r300')

savefig(FigureArrheniusC,fullfile(outDir,'Fits_Summary_C.fig'))
print(FigureArrheniusC,fullfile(outDir,'Fits_Summary_C.png'),'-dpng','-r300')

savefig(FigureVfb,fullfile(outDir,'Fits_Summary_Vfb.fig'))
print(FigureVfb,fullfile(outDir,'Fits_Summary_Vfb.png'),'-dpng','-r300')

save(fullfile(outDir,'Fits_Summary.mat'),'data','D0','D0_err','EaD','EaD_err',...
    'Cs0','Cs0_err','EaC','EaC_err','chi2','chi2p');
    