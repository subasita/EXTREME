%% Plot results
clear; close all; clc;

simyr  = [2008:2008];
addpath(genpath('../../_TOOLBOX/GENERAL'));
plotting = 0;
saveanim = 0;
metode   = 'closest';
dt       = 6;
shlon    = 0;
Ne       = 0;
tic;


load([num2str(mod(simyr(end),100),'%02d'),'_BASE1.mat'])

for yi = 1:length(simyr)
    yr = simyr(yi);
    
    datefig  = [yr-1,11,01,0,0,0];
    % datefig = [2014,06,15,00,00,00];
    BUOY     = dlmread('../../_DATA/BUOY/EXM52_2006-2011.txt');
    Heraname = ['../../_DATA/ERA/025/H025_',num2str(yr)];
    Teraname = ['../../_DATA/ERA/025/P025_',num2str(yr)];
    Deraname = ['../../_DATA/ERA/025/D025_',num2str(yr)];
    
    load([num2str(mod(yr,100),'%02d'),'_MAIN1.mat'])
    
    
    if yi == 1
        Nar  = round(size(Xp')/5);
        Nn   = 128;
        Amap = jet(Nn);
        xx   = linspace(1,0,Nn);
        xx   = [10.^xx,0];
        yy   = linspace(xx(1),xx(end),Nn);
        Amap = interp1(yy,Amap,xx);
        
        SWN_v  = [];
        ERA_v  = [];
        BUOYv  = [];
        SWNxH = zeros(size(Xp)); SWNxH = SWNxH(:);
        SWNxP = zeros(size(SWNxH));
        lon   = ncread([Heraname,'a_C.nc'],'longitude');
        lat   = ncread([Heraname,'a_C.nc'],'latitude');
        
        ERAxH = zeros(length(lon)*length(lat),1);
        ERAxP = zeros(size(ERAxH));
        SWNmH = zeros(size(Xp)); SWNmP = zeros(size(Xp));
        ERAmH = zeros(length(lat),length(lon)); ERAmP = ERAmH;
        BUOYl = [114.099, -21.6995];
    end
    % Collect
    i = 1;
    
    while exist(['Tm01_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS')],'var')
        if datenum(datefig+[0,0,0,i-1,0,0])<=datenum(datefig+[0,6,0,0,0,0])
            inf_ = 'a_C.nc';
        else
            inf_ = 'b_C.nc';
        end
        
        eval(['HS = Hsig_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Hsig_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        %     eval(['HSw = Hswell_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Hswell_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        %     eval(['Tz = Tm02_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Tm02_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        eval(['TP = TPsmoo_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''TPsmoo_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        eval(['Tm = Tm01_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Tm01_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        eval(['Dir = Dir_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Dir_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        %     eval(['Hu = Windv_x_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Windv_x_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        %     eval(['Hv = Windv_y_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Windv_y_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        
        % Make arrows
        Hu = -HS.^.25.*sind(Dir);
        Hv = -HS.^.25.*cosd(Dir);
        %     HS = sqrt(Hu.^2+Hv.^2);
        dx = round(length(Xp)/Nar(1));
        dy = round(length(Yp)/Nar(2));
        
        tt    = double(ncread([Heraname,inf_],'time'))/24+datenum(1900,1,1,0,0,0);
        
        Nt = find(round(BUOY(:,1)*48)/48+datenum(1900,1,0,0,0,0)==datenum(datefig+[0,0,0,(i-1),0,0]));
        NT = find(tt==datenum(datefig+[0,0,0,(i-1),0,0]));
        
        if (datenum(datefig+[0,0,0,(i-1),0,0])==datenum(yr,11,1,0,0,0)) && i ~= 1
            Nt = [];
        end
        
        if length(Nt)==1
            switch metode
                case 'closest'
                    % SWAN
                    Nx = find(sqrt((Xp(1,:)-BUOYl(1)).^2+(Yp(1,:)-BUOYl(2)).^2)==min(sqrt((Xp(1,:)-BUOYl(1)).^2+(Yp(1,:)-BUOYl(2)).^2)));
                    Ny = find(sqrt((Xp(:,1)-BUOYl(1)).^2+(Yp(:,1)-BUOYl(2)).^2)==min(sqrt((Xp(:,1)-BUOYl(1)).^2+(Yp(:,1)-BUOYl(2)).^2)));
                    %                 if isnan(HS(Ny(1),Nx(1))); Nx = Nx(1)+1; end
                    SWN_v = [SWN_v;datenum(datefig+[0,0,0,(i-1),0,0]),Xp(Ny(1),Nx(1)),Yp(Ny(1),Nx(1)),HS(Ny(1),Nx(1)),Tm(Ny(1),Nx(1)),TP(Ny(1),Nx(1)),Dir(Ny(1),Nx(1))];
                    
                    % ERA
                    Nx = find(sqrt((lon-BUOYl(1)).^2)==min(sqrt((lon-BUOYl(1)).^2)));
                    Ny = find(sqrt((lat-BUOYl(2)).^2)==min(sqrt((lat-BUOYl(2)).^2)));
                    if ~isempty(NT)
                        HS_ = ncread([Heraname,inf_],'swh',[Nx(1),Ny(1),NT],[1,1,1],[1,1,1]);
                        TM_ = ncread([Teraname,inf_],'mp1',[Nx(1),Ny(1),NT],[1,1,1],[1,1,1]);
                        Dr_ = ncread([Deraname,inf_],'mwd',[Nx(1),Ny(1),NT],[1,1,1],[1,1,1]);
                        ERA_v = [ERA_v;datenum(datefig+[0,0,0,(i-1),0,0]),lon(Nx(1)),lat(Ny(1)),HS_,TM_,Dr_];
                    end
                case 'interp'
                    % SWAN
                    HS_ = interp2(Xp,Yp,HS,BUOYl(1),BUOYl(2));
                    Tm_ = interp2(Xp,Yp,Tm,BUOYl(1),BUOYl(2));
                    Dr_ = interp2(Xp,Yp,Dir,BUOYl(1),BUOYl(2));
                    SWN_v = [SWN_v;datenum(datefig+[0,0,0,(i-1),0,0]),BUOYl(1),BUOYl(2),HS_,Tm_,Dr_];
                    
                    % ERA
                    if ~isempty(NT)
                        [Nx,Ny] = meshgrid(lon,lat);
                        HS_ = interp2(Nx,Ny,squeeze(ncread(Heraname,'swh',[1,1,NT],[Inf,Inf,1],[1,1,1]))',BUOYl(1),BUOYl(2));
                        TM_ = interp2(Nx,Ny,squeeze(ncread(Teraname,'mp1',[1,1,NT],[Inf,Inf,1],[1,1,1]))',BUOYl(1),BUOYl(2));
                        Dr_ = interp2(Nx,Ny,squeeze(ncread(Deraname,'mwd',[1,1,NT],[Inf,Inf,1],[1,1,1]))',BUOYl(1),BUOYl(2));
                        ERA_v = [ERA_v;datenum(datefig+[0,0,0,(i-1),0,0]),BUOYl(1),BUOYl(2),HS_,TM_,Dr_];
                    end
            end
            if ~isnan(BUOY(Nt,2))
                BUOYv = [BUOYv;datenum(datefig+[0,0,0,(i-1),0,0]),BUOYl(1),BUOYl(2),BUOY(Nt,2),BUOY(Nt,4),BUOY(Nt,8),BUOY(Nt,12)];
            end
            if ~isempty(NT)
                HS_   = squeeze(ncread([Heraname,inf_],'swh',[1,1,NT],[Inf,Inf,1],[1,1,1]))';
                TM_   = squeeze(ncread([Teraname,inf_],'mp1',[1,1,NT],[Inf,Inf,1],[1,1,1]))';
                Dr_   = squeeze(ncread([Deraname,inf_],'mwd',[1,1,NT],[Inf,Inf,1],[1,1,1]))';
            end
        else
            SWN_v = SWN_v;
            HS_   = zeros(size(ERAmH));
            TM_   = zeros(size(ERAmP));
        end
        
        if ~isempty(NT)
            HS_ = squeeze(ncread([Heraname,inf_],'swh',[1,1,NT],[Inf,Inf,1],[1,1,1]))';
            TM_ = squeeze(ncread([Teraname,inf_],'mp1',[1,1,NT],[Inf,Inf,1],[1,1,1]))';
            ERAmH = ERAmH + HS_;
            ERAmP = ERAmP + TM_;
        end
        SWNxH = max([SWNxH,double(HS(:))],[],2);
        SWNmH = SWNmH + double(HS);
        ERAxH = max([ERAxH,double(HS_(:))],[],2);
        SWNxP = max([SWNxP,double(Tm(:))],[],2);
        SWNmP = SWNmP + double(Tm);
        ERAxP = max([ERAxP,double(TM_(:))],[],2);
        Ne    = Ne + 1 ;
        
        if plotting==1
            fig1 = figure(1);
            m_proj('miller','lat',double([min(Yp(:)), max(Yp(:))]), 'lon', double([min(Xp(:)), max(Xp(:))]));
            set(gca,'color',[.9 .99 1]);
            m_pcolor(double(Xp),double(Yp),HS); shading interp; %colormap(map);
            hold on; m_quiver(Xp(1:dy:end,1:dx:end),Yp(1:dy:end,1:dx:end),...
                Hu(1:dy:end,1:dx:end),Hv(1:dy:end,1:dx:end),'k');
            hold on; m_plot(BUOY(:,1),BUOY(:,2),'ok','MarkerSize',4,'MarkerFaceColor','m','LineWidth',1.5);
            m_coast('patch',[0 0 0],'edgecolor','m');  %'none'
            
            m_grid('linewi',2,'box','fancy','tickdir','in' );
            colormap(Amap);
            colorbar('location','Eastoutside','fontsize',12,'fontweight','bold');
            caxis([0 15]);
            
            title(['HS of ',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'dd mmm yyyy [HH:MM]')],'FontSize',14)
            set(fig1,'renderer','zbuffer','Color','w','Position',[2 42 (length(Xp(1,:))/length(Yp(:,1)))*640 1.05*640]);
            
            if saveanim == 1
                filename = 'Animation00_HS.gif';
                drawnow
                frame = getframe(fig1);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256,'nodither');
                
                if i == 1
                    imwrite(imind,cm,filename,'gif', 'loopcount',inf,'DelayTime',0.25);
                else
                    imwrite(imind,cm,filename,'gif','writemode','append','DelayTime',0.25);
                end
            end
            hold off; pause(0.01)
            
        end
        i = i + dt;
    end
    disp([num2str(yr),'-',num2str(toc)])
end
%%
SWNxH = reshape(SWNxH,size(Xp));
ERAxH = reshape(ERAxH,length(lat),length(lon));
SWNxP = reshape(SWNxP,size(Xp));
ERAxP = reshape(ERAxP,length(lat),length(lon));
SWNmH = SWNmH/Ne;
SWNmP = SWNmP/Ne;
ERAmH = ERAmH/Ne;
ERAmP = ERAmP/Ne;

[Nx,Ny] = meshgrid(lon,lat);

% shift longitude
if shlon==1
    Nx = Nx(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    Nx(:,1:floor(length(lon)/2))= Nx(:,1:floor(length(lon)/2))-360;
    ERAxH = ERAxH(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    ERAxP = ERAxP(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    ERAmH = ERAmH(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    ERAmP = ERAmP(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
end

% crop
ERAxH = ERAxH((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));
ERAmH = ERAmH((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));
ERAxP = ERAxP((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));
ERAmP = ERAmP((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));

%%
fig1 = figure(1);
m_proj('miller','lat',double([min(Yp(:)), max(Yp(:))]), 'lon', double([min(Xp(:)), max(Xp(:))]));
set(gca,'color',[.9 .99 1]);
m_pcolor(double(Xp),double(Yp),SWNmH-ERAmH); shading interp;
m_coast('patch',[0 0 0],'edgecolor','m');  %'none'
m_grid('linewi',2,'box','fancy','tickdir','in' );
hold on; m_plot(BUOYl(1),BUOYl(2),'or','MarkerFaceColor','g'); hold off;
hold off;
colormap(brmap(32,-1,1,0));
colorbar('location','Eastoutside','fontsize',12,'fontweight','bold');
caxis([-1 1]);
% title('Difference of average Tm_{01} [s] (SWAN-ERA)','FontSize',14)
title('Difference of average Hs [m] (SWAN-ERA)','FontSize',14)
set(fig1,'Position',[100 125 620 540])
saveas(fig1,'FIGURE/difHSmean.fig')

fig2 = figure(2);
m_proj('miller','lat',double([min(Yp(:)), max(Yp(:))]), 'lon', double([min(Xp(:)), max(Xp(:))]));
set(gca,'color',[.9 .99 1]);
m_pcolor(double(Xp),double(Yp),SWNxH); shading interp; %colormap(map);
m_coast('patch',[0 0 0],'edgecolor','m');  %'none'
m_grid('linewi',2,'box','fancy','tickdir','in' );
hold on; m_plot(BUOYl(1),BUOYl(2),'or','MarkerFaceColor','g'); hold off;
hold off;
colormap(brmap(32,-5,5,0));
colorbar('location','Eastoutside','fontsize',12,'fontweight','bold');
caxis([0 15]);
% title('Difference of maximum Tm_{01} [s] (SWAN-ERA)','FontSize',14)
% title('Difference of maximum Hs [m] (SWAN-ERA)','FontSize',14)
title(['Maximum Hs [m] (SWAN) during ',num2str(min(simyr)),'-',num2str(max(simyr))],'FontSize',14)
set(fig2,'Position',[800 125 620 540])
saveas(fig2,'FIGURE/difHSmax.fig')

fig3 = figure(3);
m_proj('miller','lat',double([min(Yp(:)), max(Yp(:))]), 'lon', double([min(Xp(:)), max(Xp(:))]));
set(gca,'color',[.9 .99 1]);
m_pcolor(double(Xp),double(Yp),SWNmP-ERAmP); shading interp;
m_coast('patch',[0 0 0],'edgecolor','m');  %'none'
m_grid('linewi',2,'box','fancy','tickdir','in' );
hold on; m_plot(BUOYl(1),BUOYl(2),'or','MarkerFaceColor','g'); hold off;
hold off;
colormap(brmap(32,-1,1,0));
colorbar('location','Eastoutside','fontsize',12,'fontweight','bold');
caxis([-1 1]);
title('Difference of average Tm_{01} [s] (SWAN-ERA)','FontSize',14)
set(fig3,'Position',[100 125 620 540])
saveas(fig3,'FIGURE/difTMmean.fig')

fig4 = figure(4);
m_proj('miller','lat',double([min(Yp(:)), max(Yp(:))]), 'lon', double([min(Xp(:)), max(Xp(:))]));
set(gca,'color',[.9 .99 1]);
m_pcolor(double(Xp),double(Yp),SWNxP-ERAxP); shading interp; %colormap(map);
m_coast('patch',[0 0 0],'edgecolor','m');  %'none'
m_grid('linewi',2,'box','fancy','tickdir','in' );
hold on; m_plot(BUOYl(1),BUOYl(2),'or','MarkerFaceColor','g'); hold off;
hold off;
colormap(brmap(32,-5,5,0));
colorbar('location','Eastoutside','fontsize',12,'fontweight','bold');
caxis([-5 5]);
title('Difference of maximum Tm_{01} [s] (SWAN-ERA)','FontSize',14)
% title('Difference of maximum Hs [m] (SWAN-ERA)','FontSize',14)
% title(['Maximum Hs [m] (SWAN) during ',num2str(min(simyr)),'-',num2str(max(simyr))],'FontSize',14)
set(fig4,'Position',[800 125 620 540])
saveas(fig4,'FIGURE/difHSmax.fig')

%%
SWN_v = double(SWN_v);
ERA_v = double(ERA_v);
BUOYv = double(BUOYv);
strdate = datenum(datefig);
% deviation analysis
BUOY_HS = interp1(BUOYv(:,1),BUOYv(:,4),ERA_v(:,1));
BUOY_TM = interp1(BUOYv(:,1),BUOYv(:,5),ERA_v(:,1));
BUOY_Dr = interp1(BUOYv(:,1),BUOYv(:,6),ERA_v(:,1));
Hvar_SE = var(SWN_v(:,4)-ERA_v(:,4));
Hvar_SB = var(SWN_v(:,4)-BUOY_HS);
Hvar_EB = var(ERA_v(:,4)-BUOY_HS);

figure(5);
plot(SWN_v(:,1),SWN_v(:,4),'b',...
    ERA_v(:,1),ERA_v(:,4),'r',...
    ERA_v(:,1),BUOY_HS,'k',...
    'LineWidth',2);
hold on;
aa = corrcoef([BUOY_HS,SWN_v(:,4),ERA_v(:,4)]);

text(min(SWN_v(:,1))+11,4,[num2str(aa(1,2),'corr(SWAN) : %1.2f'),num2str(Hvar_SB,'| var : %1.2f')],'Color','b');
hold on;
text(min(SWN_v(:,1))+11,5,[num2str(aa(1,3),'corr(ERA)  : %1.2f'),num2str(Hvar_EB,'| var : %1.2f')],'Color','r');
hold on;
text(min(SWN_v(:,1))+11,6,[num2str(aa(2,3),'corr(ERA-SWAN) : %1.2f'),num2str(Hvar_SE,'| var : %1.2f')],'Color','m');
hold off;
set(gca,'XTick',SWN_v(1:2*4*7:end,1),'XTickLabel',datestr(SWN_v(1:2*4*7:end,1),'dd-mm'));
title('Significant Wave Height','FontSize',14)
axis([min(SWN_v(:,1)),max(SWN_v(:,1)),0,ceil(max(BUOYv(:,4)))+1])
legend('SWAN','ERA','BUOY')
set(gcf,'Position',[15,330,1230,330])
grid on;

figure(6);
plot(SWN_v(:,1),SWN_v(:,5),'b',...
    ERA_v(:,1),ERA_v(:,5),'r',...
    ERA_v(:,1),BUOY_TM,'k',...
    'LineWidth',2);
hold on;
aa = corrcoef([BUOY_TM,SWN_v(:,5),ERA_v(:,5)]);
text(min(SWN_v(:,1))+11,10,num2str(aa(1,2),'corr(SWAN) : %1.2f'),'Color','b');
hold on;
text(min(SWN_v(:,1))+11,11,num2str(aa(1,3),'corr(ERA) : %1.2f'),'Color','r');
hold on;
text(min(SWN_v(:,1))+11,12,num2str(aa(2,3),'corr(ERA-SWAN) : %1.2f'),'Color','m');
hold off;
set(gca,'XTick',SWN_v(1:2*4*7:end,1),'XTickLabel',datestr(SWN_v(1:2*4*7:end,1),'dd-mm'));
title('Period','FontSize',14)
axis([min(SWN_v(:,1)),max(SWN_v(:,1)),min([BUOYv(:,5);SWN_v(:,5)])-1,max(BUOYv(:,5))+1])
set(gcf,'Position',[15,80,1230,330])
grid on;

figure(7);
plot(SWN_v(:,1),SWN_v(:,6),'b',...
    ERA_v(:,1),ERA_v(:,6),'r',...
    ERA_v(:,1),BUOY_Dr,'k',...
    'LineWidth',2);
hold on;
aa = corrcoef([BUOY_Dr,SWN_v(:,6),ERA_v(:,6)]);
text(min(SWN_v(:,1))+11,10,num2str(aa(1,2),'corr(SWAN) : %1.2f'),'Color','b');
hold on;
text(min(SWN_v(:,1))+11,11,num2str(aa(1,3),'corr(ERA) : %1.2f'),'Color','r');
hold on;
text(min(SWN_v(:,1))+11,12,num2str(aa(2,3),'corr(ERA-SWAN) : %1.2f'),'Color','m');
hold off;
set(gca,'XTick',SWN_v(1:2*4*7:end,1),'XTickLabel',datestr(SWN_v(1:2*4*7:end,1),'dd-mm'));
title('Wave Direction','FontSize',14)
axis([min(SWN_v(:,1)),max(SWN_v(:,1)),0,360])
set(gcf,'Position',[15,80,1230,330])
grid on;


%% qq-plot
figure(8);
subplot(1,2,1)
h=qqplot(BUOYv(:,4),SWN_v(:,4));
set(h,'MarkerSize',4)
hold on; plot([0 ceil(max(BUOYv(:,4)))+1],[0 ceil(max(BUOYv(:,4)))+1],'k')
xlabel('buoy'); ylabel('swan');
axis tight equal;
subplot(1,2,2)
h=qqplot(BUOYv(:,4),ERA_v(:,4));
set(h,'MarkerSize',4)
hold on; plot([0 ceil(max(BUOYv(:,4)))+1],[0 ceil(max(BUOYv(:,4)))+1],'k')
xlabel('buoy'); ylabel('ERA');
axis tight equal;



% Save in Table
dlmwrite(strcat('TableSWN_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yymmdd'),'.txt'),SWN_v,'delimiter','\t','precision','%3.5f')
dlmwrite(strcat('TableERA_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yymmdd'),'.txt'),ERA_v,'delimiter','\t','precision','%3.5f')