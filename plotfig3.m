%% Plot results
clear; close all; clc;

simyr  = [2007,2008,2010]; 
addpath(genpath('../../_TOOLBOX/GENERAL'));
plotting = 0;
saveanim = 0;
metode   = 'closest';
dt       = 6;
shlon    = 0;
tic;
Ne       = 0;

for yi = 1:length(simyr)
yr = simyr(yi);
    
datefig  = [yr-1,11,01,0,0,0];
Ueraname = ['../../_DATA/ERA/025/U025_',num2str(yr)];
Veraname = ['../../_DATA/ERA/025/V025_',num2str(yr)];

load([num2str(mod(yr,100),'%02d'),'_BASE1.mat'])
load([num2str(mod(yr,100),'%02d'),'_WIND1.mat'])

i    = 1;

if yr == simyr(1)
Nar  = round(size(Xp')/5);
Nn   = 64;
Amap = jet(Nn);
xx   = linspace(1,0,Nn);
xx   = [10.^xx,0];
yy   = linspace(xx(1),xx(end),Nn);
Amap = interp1(yy,Amap,xx);

SWN_v  = [];
ERA_v  = [];
SWNxH = zeros(size(Xp)); SWNxH = SWNxH(:); 
lon   = ncread([Ueraname,'a_C.nc'],'longitude');
lat   = ncread([Ueraname,'a_C.nc'],'latitude');

ERAxH = zeros(length(lon)*length(lat),1);
SWNmH = zeros(size(Xp)); 
ERAmH = zeros(length(lat),length(lon)); 
end
%% Plot

while exist(['Windv_x_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS')],'var')
    if datenum(datefig+[0,0,0,i-1,0,0])<=datenum(datefig+[0,6,0,0,0,0])
       inf_ = 'a_C.nc'; 
    else
       inf_ = 'b_C.nc';
    end
    
%     eval(['HS = Hsig_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Hsig_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
%     eval(['HSw = Hswell_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Hswell_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
%     eval(['Tz = Tm02_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Tm02_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
%     eval(['TP = TPsmoo_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''TPsmoo_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
%     eval(['Tm = Tm01_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Tm01_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
%     eval(['Dir = Dir_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Dir_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
    eval(['Hu = Windv_x_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Windv_x_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
    eval(['Hv = Windv_y_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Windv_y_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
    
    % Make arrows
%     Hu = -HS.^.25.*sind(Dir);
%     Hv = -HS.^.25.*cosd(Dir);
    HS = sqrt(Hu.^2+Hv.^2);
    dx = round(length(Xp)/Nar(1));
    dy = round(length(Yp)/Nar(2));
    
    tt    = double(ncread([Veraname,inf_],'time'))/24+datenum(1900,1,1,0,0,0);

    NT = find(tt==datenum(datefig+[0,0,0,(i-1),0,0]));

  
    if ~isempty(NT) 
        Hu_ = squeeze(ncread([Ueraname,inf_],'u10',[1,1,NT],[Inf,Inf,1],[1,1,1]))';
        Hv_ = squeeze(ncread([Veraname,inf_],'v10',[1,1,NT],[Inf,Inf,1],[1,1,1]))';
        HS_ = sqrt(Hu_.^2+Hv_.^2);
        ERAmH = ERAmH + HS_;
    end
    
    SWNxH = max([SWNxH,double(HS(:))],[],2);
    SWNmH = SWNmH + double(HS);
    ERAxH = max([ERAxH,double(HS_(:))],[],2);
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

        title(['Tm_{01} of ',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'dd mmm yyyy [HH:MM]')],'FontSize',14)
        set(fig1,'renderer','zbuffer','Color','w','Position',[2 42 (length(Xp(1,:))/length(Yp(:,1)))*640 1.05*640]);
        
        if saveanim == 1
            filename = 'Animation00_TM.gif';
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
disp(num2str(toc))
end
%%
SWNxH = reshape(SWNxH,size(Xp));
ERAxH = reshape(ERAxH,length(lat),length(lon));
SWNmH = SWNmH/Ne;
ERAmH = ERAmH/Ne;

[Nx,Ny] = meshgrid(lon,lat);

% shift longitude
if shlon==1
    Nx = Nx(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    Nx(:,1:floor(length(lon)/2))= Nx(:,1:floor(length(lon)/2))-360;
    ERAxH = ERAxH(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
    ERAmH = ERAmH(:,[floor(length(lon)/2)+1:end,1:floor(length(lon)/2)]);
end

% crop
ERAxH = ERAxH((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));
ERAmH = ERAmH((lat>=min(Yp(:,1))&lat<=max(Yp(:,1))),(lon>=min(Xp(1,:))&lon<=max(Xp(1,:))));

%%
fig2 = figure(1);
m_proj('miller','lat',double([min(Yp(:)), max(Yp(:))]), 'lon', double([min(Xp(:)), max(Xp(:))]));
% m_proj('miller','lat',double([-50, 50]), 'lon', double([0, 360]));
%     m_proj('ortho','lat',-15','long',80');
set(gca,'color',[.9 .99 1]);
m_pcolor(double(Xp),double(Yp),SWNmH-ERAmH); shading interp;
m_coast('patch',[0 0 0],'edgecolor','m');  %'none'
m_grid('linewi',2,'box','fancy','tickdir','in' );
hold off;
colormap(brmap(32,-1,1,0));
colorbar('location','Eastoutside','fontsize',12,'fontweight','bold');
caxis([-1 1]);
% title('Difference of average Tm_{01} [s] (SWAN-ERA)','FontSize',14)
title('Difference of average Hs [m] (SWAN-ERA)','FontSize',14)
set(fig2,'Position',[100 125 620 540])

fig3 = figure(2);
m_proj('miller','lat',double([min(Yp(:)), max(Yp(:))]), 'lon', double([min(Xp(:)), max(Xp(:))]));
% m_proj('miller','lat',double([-30, 30]), 'lon', double([65, 120]));
%     m_proj('ortho','lat',-15','long',80');
set(gca,'color',[.9 .99 1]);
m_pcolor(double(Xp),double(Yp),SWNxH-ERAxH); shading interp; %colormap(map);
m_coast('patch',[0 0 0],'edgecolor','m');  %'none'
m_grid('linewi',2,'box','fancy','tickdir','in' );
hold off;
colormap(brmap(32,-10,10,0));
colorbar('location','Eastoutside','fontsize',12,'fontweight','bold');
caxis([-10 10]);
% title('Difference of maximum Tm_{01} [s] (SWAN-ERA)','FontSize',14)
title('Difference of maximum Hs [m] (SWAN-ERA)','FontSize',14)
set(fig3,'Position',[800 125 620 540])
