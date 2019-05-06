% Make a movie of the VFB shift with time
% load the data
clear all
folder = './AfterFitData/saturation/';
file = 'A53_003A_D16D17D18D19_2_3_Stiched_OnlyCVfbData.mat';
fullFilename = fullfile(folder,file);

dataFile = matfile(fullFilename);

time_fd = dataFile.time_fd/3600;
time_fd = time_fd';
timeSteps = length(time_fd);
depth = dataFile.depth_um;
C_t = dataFile.C_t;
C = min(C_t(:,end))*ones(length(C_t(:,1)),1);
VFB_fit = dataFile.VFB_fit;
VFB_t = NaN(timeSteps,1);

FigureMov = figure;
outFilename = strcat('Simulated_VFB_',file);%,'-old-conditions');
movieFile  = fullfile(folder,strcat(outFilename,'.avi'));
subplot(2,1,1)

pl = semilogy(depth,C,'o-');

pl.XDataSource = 'depth';
pl.YDataSource = 'C';
xlim([0,max(depth)])
ylim([1E12,1E20])
hold on
plot([dataFile.thickness dataFile.thickness],[1E10,1E20],'r:')
plot([0,0.17],[max(C_t(:,end)) max(C_t(:,end))],'r:')
hold off
xlabel('Depth (um)','FontSize',11,'FontWeight','Bold')
ylabel('[Na^+] (cm^{-3})','FontSize',11,'FontWeight','Bold')
title('Simulated SIMS Profile')

p1t_txt = text(0.025,0.12,'Time = 0 (h)','Units','normalized','HorizontalAlignment','left');




box on
ax = gca;
ax.LineWidth = 1.5;
set(gca, 'FontSize', 14)
set(gca,'XMinorTick','on','YMinorTick','on')

subplot(2,1,2)
p2 = plot(time_fd,VFB_t,'o','LineWidth',1.5);
p2.XDataSource = 'time_fd';
p2.YDataSource = 'VFB_t';
xlabel('Time (hr)','FontSize',14,'FontWeight','bold')
ylabel('Flatband Shift (V)','FontSize',14,'FontWeight','bold')


box on
ax = gca;
ax.LineWidth = 1.5;
set(gca, 'FontSize', 14)
set(gca,'XMinorTick','on','YMinorTick','on')

% Create a video writer
v = VideoWriter(movieFile);
% Open the video writer
open(v);

for i=1:timeSteps
    C = C_t(:,i);
    VFB_t(i) = VFB_fit(i);
    refreshdata
    drawnow %limitrate

    time_txt = sprintf('Time = %g (hr)', time_fd(i));
    set(p1t_txt,'String',time_txt);

    % Make a video of the solution C(x,t)
    % Capture the plot as an image 
    frame = getframe(gcf); 
    writeVideo(v,frame);
end

drawnow
close(v);