%% Plot results
clear; close all; clc;

simyr  = [2001:2016];
addpath(genpath('../../_TOOLBOX/GENERAL'));
plotting = 0;
saveanim = 0;
metode   = 'closest';
dt       = 3;
shlon    = 0;
Ne       = 0;
TCDIR    = '../../_DATA/TC/';
tic;

load([num2str(mod(simyr(end),100),'%02d'),'_BASE1.mat'])
%%
for yi = 1:length(simyr)
    yr = simyr(yi);
    
    datefig  = [yr-1,11,01,0,0,0];
    Heraname = ['../../_DATA/ERA/025/H025_',num2str(yr)];
    Teraname = {['../../_DATA/ERA/025/Z025_',num2str(yr)];'mp2'};
    Deraname = ['../../_DATA/ERA/025/D025_',num2str(yr)];
    
    if yi == 1
        
        SWN_v  = [];
        ERA_v  = [];
        lon   = ncread([Heraname,'a_C.nc'],'longitude');
        lat   = ncread([Heraname,'a_C.nc'],'latitude');
        
        SH.x = zeros(size(Xp)); SH.x = SH.x(:);               % Maximum
        SP.x = zeros(size(SH.x));
        EH.x = zeros(length(lon)*length(lat),1);         
        EP.x = zeros(size(EH.x));
        SH.a = zeros(size(Xp)); SP.a = zeros(size(Xp));       % Average
        EH.a = zeros(length(lat),length(lon)); EP.a = EH.a;
        SH.n = zeros(size(SH.x)); SP.n = zeros(size(SH.x));   % Minimum
        EH.n = zeros(size(EH.x)); EP.n = zeros(size(EH.x));
        
%         BUOYl = [114.099, -21.6995]; %exmouth
        BUOYl = [...75.80996, 10.527197; ... % Kozhikode
            82.923165, 12.548534; ... % BD11
            ...73.108815, 8.144392; ... % AD09
            84.119604, 18.200825; ... % Vizag
            ...73.092898,16.897517; ... % Ratnagiri
            114.0986111,-21.699444; ... % Exmouth
            113.7441667,-23.174321; ... % Coral Bay
            110.547,-8.1362; ... % Baron Bay
            129.521667,28.451944; ... % 402
            127.6475, 26.258056; ... % 702
            127.965278, 26.24222; ... % 701, Okinawa
            124.102778, 24.365278; ... % 705
            125.23556, 24.860833]; % 706
        BUOYn = 'Okinawa';
    end
    
    % add tropical cyclone
%     dirname = pwd;
%     cd([TCDIR,num2str(yr),'a/']);
%     cfile   = dir('*.txt');
%     cd(dirname);
%     [NT,~]  = size(cfile);
%     for ii=1:NT
%         % read TC files
%         tcname = [[TCDIR,num2str(yr),'a/'] cfile(ii).name];
%         TC{2*yi-1,ii}  = readTC(tcname,round(dt));
%     end
%     
%     cd([TCDIR,num2str(yr),'b/']);
%     cfile   = dir('*.txt');
%     cd(dirname);
%     [NT,~]  = size(cfile); NT = NT;
%     for ii=1:NT
%         % read TC files
%         tcname = [[TCDIR,num2str(yr),'b/'] cfile(ii).name];
%         TC{2*yi,ii}  = readTC(tcname,round(dt));
%     end


%     %Collect
    [SH,SP,EH,EP,SWN_v,ERA_v] = minmaxav(yr,datefig,dt,metode,BUOYl,SWN_v,ERA_v,SH,SP,EH,EP,Xp,Yp,Heraname,Teraname,Deraname,[datenum(2017,1,1,0,0,0);1;0]);

    disp([num2str(yr),'-',num2str(toc)])
end
%% Reshaping
%  Maximum
SWNxH = reshape(SH.x,size(Xp));
ERAxH = reshape(EH.x,length(lat),length(lon));
SWNxP = reshape(SP.x,size(Xp));
ERAxP = reshape(EP.x,length(lat),length(lon));
%  Average
SWNaH = SH.a/length(SWN_v(:,1));
SWNaP = SP.a/length(SWN_v(:,1));
ERAaH = EH.a/length(ERA_v(:,1));
ERAaP = EP.a/length(ERA_v(:,1));
%  Minimum
SWNnH = reshape(SH.n,size(Xp));
ERAnH = reshape(EH.n,length(lat),length(lon));
SWNnP = reshape(SP.n,size(Xp));
ERAnP = reshape(EP.n,length(lat),length(lon));


[Nx,~] = meshgrid(lon,lat);

% shift longitude
if shlon==1
    Nx = Nx(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    Nx(:,1:floor(length(lon)/2))= Nx(:,1:floor(length(lon)/2))-360;
    ERAxH = ERAxH(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    ERAxP = ERAxP(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    ERAaH = ERAaH(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    ERAaP = ERAaP(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    ERAnH = ERAnH(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    ERAnP = ERAnP(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
end

% crop to fit domain
ERAxH = ERAxH((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));
ERAaH = ERAaH((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));
ERAnH = ERAnH((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));
ERAxP = ERAxP((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));
ERAaP = ERAaP((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));
ERAnP = ERAnP((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));

%% Plot Map
close all;
Nar  = round(size(Xp')/5);
Nn   = 128;
Amap = jet(Nn);
xx   = linspace(1,0,Nn);
xx   = [10.^xx,0];
yy   = linspace(xx(1),xx(end),Nn);
Amap = interp1(yy,Amap,xx);

fig1 = figure(1);
m_proj('miller','lat',double([min(Yp(:)), max(Yp(:))]), 'lon', double([min(Xp(:)), max(Xp(:))]));
set(gca,'color',[.9 .99 1],'FontSize',14);
m_pcolor(double(Xp),double(Yp),SWNaH-ERAaH); shading interp;
m_coast('patch',[0 0 0],'edgecolor','m');  %'none'
m_grid('linewi',2,'box','fancy','tickdir','in' ,'fontsize',14 );
hold on; m_plot(BUOYl(1),BUOYl(2),'or','MarkerFaceColor','g','MarkerSize',10); hold off;
hold off;
colormap(brmap(32,-1,1,0));
colorbar('location','Eastoutside','fontsize',12,'fontweight','bold');
caxis([-1 1]);
% title('Difference of average Tm_{01} [s] (SWAN-ERA)','FontSize',14)
title('Difference of average Hs [m] (SWAN-ERA)','FontSize',14)
set(fig1,'Position',[100 125 620 540])
saveas(fig1,'FIGURE/MAP-HSav.png')

fig2 = figure(2);
m_proj('miller','lat',double([min(Yp(:)), max(Yp(:))]), 'lon', double([min(Xp(:)), max(Xp(:))]));
set(gca,'color',[.9 .99 1],'FontSize',14);
m_pcolor(double(Xp),double(Yp),SWNxH-ERAxH); shading interp; %colormap(map);
m_coast('patch',[0 0 0],'edgecolor','m');  %'none'
m_grid('linewi',2,'box','fancy','tickdir','in','fontsize',14  );
hold on; m_plot(BUOYl(1),BUOYl(2),'or','MarkerFaceColor','g','MarkerSize',10); hold off;
% [Nx,Ny] = size(TC);
% for i =1:Nx*Ny
%     if ~isempty(TC{i})
%         ii = find(TC{i}.lon<=max(Xp(1,:))&TC{i}.lon>=min(Xp(1,:))&TC{i}.lat<=max(Yp(:,1))&TC{i}.lat>=min(Yp(:,1)));
%         if  ~isempty(ii)
%             hold on;
%             m_plot(TC{i}.lon,TC{i}.lat,'-o','MarkerSize',2,'Color',[0.2,0.2,0.2,0.2])
%         end
%     end
% end
hold off;
colormap(brmap(32,-5,5,0));
colorbar('location','Eastoutside','fontsize',12,'fontweight','bold');
caxis([-5 5]);
% title('Difference of maximum Tm_{01} [s] (SWAN-ERA)','FontSize',14)
title('Difference of maximum Hs [m] (SWAN-ERA)','FontSize',14)
% title(['Maximum Hs [m] (SWAN) during ',num2str(min(simyr)),'-',num2str(max(simyr))],'FontSize',14)
set(fig2,'Position',[800 125 620 540])
saveas(fig2,'FIGURE/MAP-HSmx02.png')

fig3 = figure(3);
m_proj('miller','lat',double([min(Yp(:)), max(Yp(:))]), 'lon', double([min(Xp(:)), max(Xp(:))]));
set(gca,'color',[.9 .99 1],'FontSize',14);
m_pcolor(double(Xp),double(Yp),SWNaP-ERAaP); shading interp;
m_coast('patch',[0 0 0],'edgecolor','m');  %'none'
m_grid('linewi',2,'box','fancy','tickdir','in' ,'fontsize',14 );
hold on; m_plot(BUOYl(1),BUOYl(2),'or','MarkerFaceColor','g','MarkerSize',10); hold off;
hold off;
colormap(brmap(32,-1,1,0));
colorbar('location','Eastoutside','fontsize',12,'fontweight','bold');
caxis([-1 1]);
title('Difference of average Tm_{01} [s] (SWAN-ERA)','FontSize',14)
set(fig3,'Position',[100 125 620 540])
saveas(fig3,'FIGURE/MAP-TMav_ori.fig')

fig4 = figure(4);
m_proj('miller','lat',double([min(Yp(:)), max(Yp(:))]), 'lon', double([min(Xp(:)), max(Xp(:))]));
set(gca,'color',[.9 .99 1],'FontSize',14);
m_pcolor(double(Xp),double(Yp),SWNxP-ERAxP); shading interp; %colormap(map);
m_coast('patch',[0 0 0],'edgecolor','m');  %'none'
m_grid('linewi',2,'box','fancy','tickdir','in','fontsize',14  );
hold on; m_plot(BUOYl(1),BUOYl(2),'or','MarkerFaceColor','g','MarkerSize',10); hold off;
hold off;
colormap(brmap(32,-5,5,0));
colorbar('location','Eastoutside','fontsize',12,'fontweight','bold');
caxis([-5 5]);
title('Difference of maximum Tm_{01} [s] (SWAN-ERA)','FontSize',14)
% title('Difference of maximum Hs [m] (SWAN-ERA)','FontSize',14)
% title(['Maximum Hs [m] (SWAN) during ',num2str(min(simyr)),'-',num2str(max(simyr))],'FontSize',14)
set(fig4,'Position',[800 125 620 540])
saveas(fig4,'FIGURE/MAP-TMmx_ori.fig')

%% Plot Time-Series
SWN_v  = double(SWN_v);
ERA_v  = double(ERA_v);
% BUOY     = dlmread('../../_DATA/BUOY/EXM52_2006-2011.txt');
BUOYv = [];
for ii = 1:length(simyr)
    buoyname = ['../../_DATA/JMA/MOORED/h701e.',num2str(mod(simyr(ii),1000),'%03d'),'.txt'];
    [BUOY]   = funreadJAPAN(buoyname);
    % BUOYv  = [round(BUOY(:,1)*48)/48+datenum(1900,1,0,0,0,0) BUOY(:,[2 4 12])];
    BUOYv    = [BUOYv; BUOY.date BUOY.Hs BUOY.Tave double(BUOY.Wdir)];
end
BUOYv        = BUOYv((BUOYv(:,2)<50 & BUOYv(:,3)<100),:);
[~,ie,is]    = intersect(ERA_v(:,1),SWN_v(:,1));
[~,ieb,ibe]  = intersect(ERA_v(:,1),BUOYv(:,1));
[~,isb,ibs]  = intersect(SWN_v(:,1),BUOYv(:,1));
[~,~,ies]    = intersect(ERA_v(ieb,1),SWN_v(:,1));

strdate = datenum(datefig);
% deviation analysis
Hvar_SE = var(SWN_v(is,4)-ERA_v(ie,4));
Hvar_SB = var(SWN_v(isb,4)-BUOYv(ibs,2));
Hvar_EB = var(ERA_v(ieb,4)-BUOYv(ibe,2));

figure(5);
plot(SWN_v(:,1),SWN_v(:,4),'b',...
    ERA_v(:,1),ERA_v(:,4),'r',...
    'LineWidth',2);
hold on; plot(BUOYv(:,1),BUOYv(:,2),'ok','MarkerSize',2)
hold on;
aa = corrcoef([BUOYv(ibe,2),SWN_v(ies,4),ERA_v(ieb,4)]);

text(min(SWN_v(:,1))+11,4,[num2str(aa(1,2),'corr(SWAN-BUOY) : %1.2f'),num2str(Hvar_SB,'| var : %1.2f')],'Color','b');
hold on;
text(min(SWN_v(:,1))+11,5,[num2str(aa(1,3),'corr(ERA-BUOY)  : %1.2f'),num2str(Hvar_EB,'| var : %1.2f')],'Color','r');
hold on;
text(min(SWN_v(:,1))+11,6,[num2str(aa(2,3),'corr(ERA-SWAN) : %1.2f'),num2str(Hvar_SE,'| var : %1.2f')],'Color','m');
hold off;
set(gca,'XTick',SWN_v(1:4*4*7:end,1),'XTickLabel',datestr(SWN_v(1:4*4*7:end,1),'dd-mm'));
ylabel('Wave Height [m]')
title('Time series of Significant Wave Height during 2008','FontSize',16)
axis([datenum(2007,1,1,0,0,0),datenum(2010,1,1,0,0,0),0,ceil(max(BUOYv(:,2)))+1])
h=legend('SWAN','ERA','BUOY'); set(h,'FontSize',14);
set(gcf,'Position',[1, 778, 1260, 330])
set(gca,'Position',[0.05 0.1 0.925 0.8],'FontSize',14)
grid on;

figure(6);
plot(SWN_v(:,1),SWN_v(:,5),'b',...
    ERA_v(:,1),ERA_v(:,5),'r',...
    'LineWidth',2);
hold on; plot(BUOYv(:,1),BUOYv(:,3),'ok','MarkerSize',2)
hold on;
aa = corrcoef([BUOYv(ibe,3),SWN_v(ies,5),ERA_v(ieb,5)]);
text(min(SWN_v(:,1))+11,11,num2str(aa(1,2),'corr(SWAN-BUOY) : %1.2f'),'Color','b');
hold on;
text(min(SWN_v(:,1))+11,12,num2str(aa(1,3),'corr(ERA-BUOY) : %1.2f'),'Color','r');
hold on;
text(min(SWN_v(:,1))+11,13,num2str(aa(2,3),'corr(ERA-SWAN) : %1.2f'),'Color','m');
hold off;
set(gca,'XTick',SWN_v(1:4*4*7:end,1),'XTickLabel',datestr(SWN_v(1:4*4*7:end,1),'dd-mm'));
ylabel('Period [s]')
title('Time series of Mean Wave Period (Tm_{01}) during 2008','FontSize',16)
axis([datenum(2007,11,1,0,0,0),datenum(2009,11,1,0,0,0),0,ceil(max(BUOYv(:,3)))+1])
h=legend('SWAN','ERA','BUOY'); set(h,'FontSize',14);
set(gcf,'Position',[1, 228, 1260, 330])
set(gca,'Position',[0.05 0.1 0.925 0.8],'FontSize',14)
grid on;

figure(7);
plot(SWN_v(:,1),SWN_v(:,6),'b',...
    ERA_v(:,1),ERA_v(:,6),'r',...
    'LineWidth',2);
hold on; plot(BUOYv(:,1),BUOYv(:,4),'ok','MarkerSize',2)
hold on;
aa = corrcoef([BUOYv(ibe,4),SWN_v(ies,6),ERA_v(ieb,6)]);
text(min(SWN_v(:,1))+11,20,num2str(aa(1,2),'corr(SWAN-BUOY) : %1.2f'),'Color','b');
hold on;
text(min(SWN_v(:,1))+11,40,num2str(aa(1,3),'corr(ERA-BUOY) : %1.2f'),'Color','r');
hold on;
text(min(SWN_v(:,1))+11,60,num2str(aa(2,3),'corr(ERA-SWAN) : %1.2f'),'Color','m');
hold off;
set(gca,'XTick',SWN_v(1:4*4*7:end,1),'XTickLabel',datestr(SWN_v(1:4*4*7:end,1),'dd-mm'));
ylabel('Direction [^\circ]')
title('Time series of Mean Wave Direction during 2008','FontSize',16)
axis([datenum(2007,11,1,0,0,0),datenum(2009,11,1,0,0,0),0,360])
h=legend('SWAN','ERA','BUOY'); set(h,'FontSize',14);
set(gcf,'Position',[1, 1028, 1260, 330])
set(gca,'Position',[0.05 0.1 0.925 0.8],'FontSize',14)
grid on;


%% qq-plot
figure(8);
subplot(1,2,1)
h=qqplot(BUOYv(:,2),SWN_v(:,4));
set(h,'MarkerSize',4)
set(gca,'FontSize',14)
hold on; plot([0 ceil(max(BUOYv(:,2)))+1],[0 ceil(max(BUOYv(:,2)))+1],'k')
xlabel([BUOYn,' Buoy']); ylabel('SWAN');
axis tight equal;
hold off;
box on;
title('Q-Q plot HS SWAN vs. BUOY','FontSize',14);
subplot(1,2,2)
h=qqplot(BUOYv(:,2),ERA_v(:,4));
set(h,'MarkerSize',4)
set(gca,'FontSize',14)
hold on; plot([0 ceil(max(BUOYv(:,2)))+1],[0 ceil(max(BUOYv(:,2)))+1],'k')
xlabel([BUOYn,' Buoy']); ylabel('ERA-Interim');
axis tight equal;
hold off;
box on;
title('Q-Q plot HS ERA-Interim vs. BUOY','FontSize',14);
set(gcf,'Position',[408   338   752   400])
saveas(gcf,['FIGURE/QQ-HS-',num2str(simyr(1)),'-',num2str(simyr(end)),'.png'])

figure(9);
subplot(1,2,1)
h=qqplot(BUOYv(:,3),SWN_v(:,5));
set(h,'MarkerSize',4)
set(gca,'FontSize',14)
hold on; plot([0 ceil(max(BUOYv(:,3)))+1],[0 ceil(max(BUOYv(:,3)))+1],'k')
xlabel([BUOYn,' Buoy']); ylabel('SWAN');
axis tight equal;
hold off;
box on;
title('Q-Q plot TM_{01} SWAN vs. BUOY','FontSize',14);
subplot(1,2,2)
h=qqplot(BUOYv(:,3),ERA_v(:,5));
set(h,'MarkerSize',4)
set(gca,'FontSize',14)
hold on; plot([0 ceil(max(BUOYv(:,3)))+1],[0 ceil(max(BUOYv(:,3)))+1],'k')
xlabel([BUOYn,' Buoy']); ylabel('ERA-Interim');
axis tight equal;
hold off;
box on;
title('Q-Q plot TM_{01} ERA-Interim vs. BUOY','FontSize',14);
set(gcf,'Position',[908   338   752   400])
saveas(gcf,['FIGURE/QQ-TM-',num2str(simyr(1)),'-',num2str(simyr(end)),'.png'])

%% Save
eval(['Dp = Depth',datestr(datenum(simyr(end)-1,11,1,0,0,0),'_yyyymmdd_HHMMSS'),'; clear(''Depth',datestr(datenum(simyr(end)-1,11,1,0,0,0),'_yyyymmdd_HHMMSS'),''')'])
eval('wb.ind = find(~isnan(Dp));'); % index of data point (which is regarded as ocean)
wb.SxH = SWNxH;
wb.SnH = SWNnH;
wb.ExH = ERAxH;
wb.EnH = ERAnH;
wb.SxP = SWNxP;
wb.SnP = SWNnP;
wb.ExP = ERAxP;
wb.EnP = ERAnP;

dlmwrite(strcat('TableSWN_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yymmdd'),'.txt'),SWN_v,'delimiter','\t','precision','%3.5f')
dlmwrite(strcat('TableERA_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yymmdd'),'.txt'),ERA_v,'delimiter','\t','precision','%3.5f')
save('MASK.mat','wb');