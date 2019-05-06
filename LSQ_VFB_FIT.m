clear all
% Set Constants
Kb = 8.6e-5; %eV/K
format long

global thickness;
global tempC;
global VBias;
global simulationTime;
global time_exp;
global fileDest;
global area_m;
global CN;
global Vsh;
global weights;
% global C0;

% Set Geometry
area_cm = pi*(.05)^2; %cm
area_m  = area_cm*1E-4;
% thickness = mean([66.85,76.96])*1E-3; % um
thickness = 80.26E-3; % um
% VBias = mean([6.685,7.696]); %V
VBias = 8.0; %V
E = VBias/thickness/100; %MV/cm

% Directory where the data is (if 3 files)
% inFileDirs = ["G:\My Drive\Research\PVRD1\Shared\Experimental Data\20180702",...
%               "G:\My Drive\Research\PVRD1\Shared\Experimental Data\20180724",...
%               "G:\My Drive\Research\PVRD1\Shared\Experimental Data\20180629"];

% Files to fit
% inFiles = ["A53_003_B_D14","A53_003_A_D5","A53_003_B_D13"];
% Temperature for each file
% tempCArray = [70,80,90];

% inFileDirs = ["G:\My Drive\Research\PVRD1\Shared\Experimental Data\20180831"];
% inFiles = ["A53_003A_D16D17D18D19_2_3_Stiched_OnlyCVfbData"];
% tempCArray = [60];

inFileDirs = [%"G:\My Drive\Research\PVRD1\Shared\Experimental Data\ZZcleanVScont\70C",...
              %"G:\My Drive\Research\PVRD1\Shared\Experimental Data\ZZcleanVScont\80C"];%,...
          %"G:\My Drive\#Shared\#Shared_Jonathan\Experimental Data_outdated_seeFenningMembers\ZZcleanVScont\80C"];
          "G:\My Drive\Research\PVRD1\DATA\CV\file_backup_20190416"];
              %"G:\My Drive\Research\PVRD1\Shared\Experimental Data\ZZcleanVScont\90C",...
              %"G:\My Drive\Research\PVRD1\Shared\Experimental Data\ZZcleanVScont\100C"];

inFiles = [%"comparison_A53_004C_D26D27D28D29_HR3_Na_3_D9D10D11D12_smoothed_11-14-2018",...
           %"CleanCont_80C_comparison_A53_004B_D5D6D7D8_HR3_Na_3_DD17D18D19D20"];%,...
           "AC1_D25D26D27D28_ANa1_D1D2D3D4_cont_minus_clean"];
           %"comparison_A53_004B_D1D2D3D4_HR3_Na_3_D13D14D15D16_11-12-2018_stitched_smoothed_11-14_usethis_smoothed_11-14-2018",...
           %"comparison_A53_004A_D15D16D17D18_20181001_A53_003C_D5D6D7D8"];
% tempCArray = [70,80,90,100];
tempCArray = [30];



%Load Experimental Vfb-Time Curve
fileDest = pwd;

FileSaveDir = "AfterFitData/clean_contaminated"; % A directory to save the output data
if 7 ~= exist(FileSaveDir,'dir')
    mkdir(FileSaveDir);
end

% Iterate over all the files
for i=1:length(tempCArray)

%%Set Migration Stresses
tempC = tempCArray(i); %C

dataFolder = inFileDirs(i);

FileName = inFiles(i); 
PinNum = 2; % The pin number ?


S   = load(fullfile(dataFolder,strcat(FileName,'.mat')));
% load(fullfile(dataFolder,strcat(FileName,'.mat')));
st  = 1; %Starting data array positon for fitting
% en  = length(S.Data(PinNum).VfbAve)-1; %Ending data array positon for fitting
en  = length(S.T12VfbPinAveShift)-1; %Ending data array positon for fitting

%en = 103;

% Estimate the flat band voltage shift to compare with the simulation
CN = mean(S.C(:,end)); %Average Accumulation Capacitance of Data (F) % Check this, we may need to average over all pins. This is the data of one pin only, the last pin of T2.
Vfb = S.T12VfbPinAveShift(st:en)'; %Average Flatband Voltage (V)
VfbStd = S.T1T2VfbAveStd(st:en); %Vfb Standard Deviation Error
Vsh = Vfb-Vfb(1); %Flatband Voltage Shift (V)
VfbStd = (VfbStd.^2+VfbStd(1)^2).^0.5;
tfb = S.tfb(st:en); % Time Data (s)
if isrow(tfb)
    tfb = tfb';
end
if isrow(Vsh)
    Vsh = Vsh';
end

% CN = mean(S.Data(PinNum).C(:,end)); %Average Accumulation Capacitance of Data (F)
% Vfb = S.Data(PinNum).VfbAve(st:en)'; %Average Flatband Voltage (V)
% VfbStd = S.Data(PinNum).VfbStd(st:en); %Vfb Standard Deviation Error
% Vsh = Vfb-Vfb(1); %Flatband Voltage Shift (V)
% VfbStd = (VfbStd.^2+VfbStd(1)^2).^0.5;
% tfb = S.Data(PinNum).tfb(st:en); % Time Data (s)
% if isrow(Vsh)
%     Vsh = Vsh';
% end

if isrow(VfbStd)
    VfbStd = VfbStd';
end

% FigureLoadData = figure;
% errorbar(tfb./3600,Vsh,VfbStd,'o','MarkerFaceColor',[1 1 1],'LineWidth',1.5)
% xlabel('Time (hr)','FontSize',14,'FontWeight','bold')
% ylabel('Flatband Shift (V)','FontSize',14,'FontWeight','bold')

% weights = (mean(VfbStd) - VfbStd).^(-2);
weights = (VfbStd).^(-2);
weights(isinf(weights)) = 1/min(weights);
%weights = (weights./max(weights))';
% weights = ones(length(VfbStd),1);

%%Set Migration Time
simulationTime = tfb(end);%/(3600); %hrs
tStep = (tfb(2)-tfb(1));%/3600; %Time Step During Experiment (hrs)

%Set simulated time array to start and end of experimental time
time_exp = (tfb(1):tStep:tfb(end))';

%%et Initial Kinetic Parameters Fot Fitting
%C0 = 5.401e16; %#/cm^3
%D0 = 4.099e-7; %cm^2/s
%Ea = .43;
%D = D0*exp(-Ea/(Kb*(T+273.15)));

C0  = 2.14E17; %Solubility (#/cm^3)
D   = 1.0e-17; %Diffusivity (cm^2/s)

% C0  = 2.14E17; %Solubility (#/cm^3)
% D   = 1.0e-13; %Diffusivity (cm^2/s)

x0 = log10([D,C0]);

fun = @(r)norm((getVFBSH_Matlab(10.^r)-Vsh).*weights);


% Minimize the problem to get initial values for 
% D and C0 using fminsearch
options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval},...
    'TolFun',1E-3,'TolX',1E-3);

x_simplex = fminsearch(fun,x0,options);

options = optimoptions('lsqnonlin','Display','iter-detailed',...
    'FunctionTolerance',1E-6,'Algorithm','levenberg-marquardt',...
    'StepTolerance', 1E-6,'FiniteDifferenceType','central',...
    'ScaleProblem','Jacobian','InitDamping',1E-2,...
    'PlotFcn',{'optimplotx','optimplotresnorm'});

[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@VFBSH_Residue,x_simplex,[],[],options);


% Get the output params from the x array
D_fit   = 10.^x(1)
C0_fit  = 10.^x(2)

ci = 10.^(nlparci(x,residual,'jacobian',jacobian));

D_err   = max(abs(ci(1,:)-D_fit));
C0_err  = max(abs(ci(2,:)-C0_fit));


[ypred,delta] = nlpredci('getVFBSH_Matlab_callback',time_exp,x,residual,'Jacobian',(jacobian),...
                         'SimOpt','on');
                     
                     
lower = ypred - delta;
upper = ypred + delta;

[C_t,VFB_fit,time_fd,depth_um] = FDNP_SiNxDevice1Dirichlet(D_fit,C0_fit,thickness,tempC,VBias,simulationTime,...
    CN,area_m);
                     
% Plot the results
FinalFig = figure;

hold on
errorbar(tfb./3600,Vsh,VfbStd,'o','MarkerFaceColor',[1 1 1],'LineWidth',1.5)
plot(time_exp./3600,ypred,'k','LineWidth',1.5)
plot(time_exp./3600,[lower';upper'],'r--','LineWidth',1.5)

hold off
xlabel('Time (hr)','FontSize',14,'FontWeight','bold')
ylabel('Flatband Shift (V)','FontSize',14,'FontWeight','bold')

d0_c0_txt = sprintf('D = %.2e ± %.2e cm^2/s\nC_0 = %.2e ± %.2e',...
    D_fit,D_err,C0_fit,C0_err);

p1t_txt = text(0.025,0.97,d0_c0_txt,'Units','normalized',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top');


leg = legend('Experiment','Model','Prediction bands','Location','Northeast');
legend('boxoff')
leg.FontSize = 12;

box on
ax = gca;
ax.LineWidth = 1.5;
set(gca, 'FontSize', 14)
set(gca,'XMinorTick','on','YMinorTick','on')

% fullFilename = fullfile('./FileSaveDir',strcat(FileName,'.fig'));



savefig(FinalFig,char(FileSaveDir+"/"+FileName+".fig"))
print(FinalFig,char(FileSaveDir+"/"+FileName+".png"),'-dpng','-r300')


SIMSFig = figure;
semilogy(depth_um,C_t(:,end),'o-')
xlabel('Depth (um)','FontSize',14,'FontWeight','bold')
ylabel('[Na^+] (cm^{-3})','FontSize',14,'FontWeight','bold')
title('Simulated SIMS Profile')
xlim([0,0.17])
ylim([1E12,1E20])
leg.FontSize = 12;
hold on
% plot([thickness thickness],ylim)
hold off

box on
ax = gca;
ax.LineWidth = 1.5;
set(gca, 'FontSize', 14)
set(gca,'XMinorTick','on','YMinorTick','on')

savefig(SIMSFig,char(FileSaveDir+"/"+FileName+"-SIMS.fig"))
print(SIMSFig,char(FileSaveDir+"/"+FileName+"-SIMS.png"),'-dpng','-r300')



save(char(FileSaveDir+"/"+FileName+".mat"),'D_fit','C0_fit','C_t',...
    'D_err','C0_err','VFB_fit','depth_um',...
    'time_exp','Vsh','VfbStd','VBias',...
    'thickness','tempC','area_m','CN','VBias','PinNum',...
    'time_fd','resnorm','residual','exitflag','output','lambda',...
    'jacobian','ypred','delta','x0','x_simplex','x');

end