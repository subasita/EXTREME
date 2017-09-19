% read and plot table SWAN
clear; close all; clc;
simyr     = 2008;

tablename = [num2str(mod(simyr,100),'%02d'),'_EIO_1.tbl'];
buoyname  = '../../_DATA/BUOY/EXM52_2006-2011.txt';
Heraname  = ['../../_DATA/ERA/025/H025_',num2str(simyr)];
Teraname  = ['../../_DATA/ERA/025/P025_',num2str(simyr)];
plott     = 1;

varname   = {'Hsig';'Hswell';'Dir';'Tm01';'Tp';'Tm02';'Wu';'Wv'};
strdate   = datenum(simyr-1,11,1,0,0,0);
endate    = datenum(simyr,11,1,0,0,0);
dt        = datenum(0,0,0,1,0,0);
DT        = 6;

SWAN    = load(tablename);
BUOY    = dlmread(buoyname);
Nswan   = 1 ; %length(load('WaveStation2.txt'));
n_sta   = 1;

% LOCATION IN ERA
n_Loc   = [ 114.099, -21.6995];
n_Var   = length(SWAN(1,:));
Loc_x   = ncread([Heraname,'a_C.nc'],'longitude');
Loc_x   = find(abs(Loc_x-n_Loc(1))==min(abs(Loc_x-n_Loc(1)))); Loc_x = Loc_x(1);
Loc_y   = ncread([Heraname,'a_C.nc'],'latitude');
Loc_y   = find(abs(Loc_y-n_Loc(2))==min(abs(Loc_y-n_Loc(2)))); Loc_y = Loc_y(1);

% TIME
timeb   = BUOY(:,1)+datenum(1900,1,-1,0,0,0);
n_time  = length(SWAN(:,1))/Nswan;
timev   = strdate:dt:strdate+dt*(n_time-1);
timee   = double(ncread([Heraname,'a_C.nc'],'time'))/24+datenum(1900,1,1,0,0,0);
timee   = [timee(1:end-4);double(ncread([Heraname,'b_C.nc'],'time'))/24+datenum(1900,1,1,0,0,0)];
% HS
Hsig    = interp1(timeb,BUOY(:,2),timev);
Hsig    = smooth(Hsig,DT);
[fH,vH] = ecdf(Hsig);

% Tp
Tpeak   = interp1(timeb,BUOY(:,3),timev);
Tpeak   = smooth(Tpeak,DT);
[fTp,vTp] = ecdf(Tpeak);

% Tmean
Tmean   = interp1(timeb,BUOY(:,4),timev);
Tmean   = smooth(Tmean,DT);
[fT,vT] = ecdf(Tmean);

%%
close all;

for ii=1:n_Var
    if strcmp(varname{ii},'Hsig')
        HsSW  = SWAN(n_sta:Nswan:end,ii);
        HsEra = squeeze(ncread([Heraname,'a_C.nc'],'swh',[Loc_x,Loc_y,1],[1,1,Inf],[1,1,1]));
        HsEra = [HsEra(1:end-4,1); squeeze(ncread([Heraname,'b_C.nc'],'swh',[Loc_x,Loc_y,1],[1,1,Inf],[1,1,1]))];
        
        % timeseries
        figure(ii);
        plot(timev,HsSW,'b',timee,HsEra,'r',timeb,BUOY(:,2),'k','LineWidth',2)
        hold on;
        aa = corrcoef([Hsig(1:DT:end),HsSW(1:DT:end),HsEra(timee>=strdate&timee<=max(timev))]);
        bb = var(HsSW-Hsig);
        text(min(timev)+11,0.5,[num2str(aa(1,2),'corr(SWAN) : %1.2f'),num2str(bb,' | var = %1.2f')],'Color','b');
        hold on;
        bb = var(HsEra(timee>=strdate&timee<=max(timev))-Hsig(1:DT:end));
        text(min(timev)+11,2.75,[num2str(aa(1,3),'corr(ERA) : %1.2f'),num2str(bb,' | var = %1.2f')],'Color','r');
        hold on;
        bb = var(HsSW(1:DT:end)-HsEra(timee>=strdate&timee<=max(timev)));
        text(min(timev)+11,3.5,[num2str(aa(2,3),'corr(ERA-SWAN) : %1.2f'),num2str(bb,' | var = %1.2f')],'Color','m');
        hold on;
        patch([strdate,strdate,strdate,strdate],[0,0,ceil(max(BUOY(:,2))),ceil(max(BUOY(:,2)))],...
            zeros(4,1),'c','EdgeColor','none','FaceAlpha',0.25)
        hold off;
        set(gca,'XTick',timev(1:14*24:end),'XTickLabel',datestr(timev(1:14*24:end)','dd/mm'))
        grid on;
        set(gcf,'Position',[1 200 1300 400]);
        set(gca,'Position',[0.05 0.11 0.9 .8]);
        axis([min(timev) max(timev) 0 ceil(max(BUOY(:,2)))])
        legend('SWAN (0.25)','ERA (0.25)','BUOY')
        title('Significant Wave Height','FontSize',14)
        saveas(gcf,['FIGURE/BUOY_HsSeries_',num2str(simyr),'.fig'])
        pause(0.1)
        
        if plott==1
            % Q-Q Plot
            [fHs,vHs] = ecdf(HsSW(timev>=strdate)');

            figure(ii+10);
%             subplot(1,2,1)
%             plot(vHs,fHs,'b',vH,fH,'k','LineWidth',2)
%             h = legend('SWAN','Buoy'); set(h,'Location','SouthEast')
            subplot(1,2,1)
            h = qqplot(BUOY(timeb>=strdate&timeb<=endate,2),HsSW');
            set(h,'MarkerSize',3);
            hold on;
            plot([0.5 4],[0.5 4],'--k')
            hold off; 
            box on; axis([.5 4 .5 4])
            xlabel('Buoy'); ylabel('SWAN')
            title('HS : SWAN vs. Buoy')
            subplot(1,2,2)
            h = qqplot(BUOY(:,2),HsEra);
            set(h,'MarkerSize',3);
            hold on;
            plot([0.5 4],[0.5 4],'--k')
            hold off; 
            box on; axis([.5 4 .5 4])
            xlabel('Buoy'); ylabel('ERA')
            title('HS : ERA vs. Buoy')
            set(gcf,'Position',[403 310 660 356])
            saveas(gcf,['FIGURE/BUOY_HsQQ_',num2str(simyr),'.fig'])
            pause(0.1)
        end
    elseif strcmp(varname{ii},'Tm01')
        TmSW  = SWAN(n_sta:Nswan:end,ii);
        TmEra = squeeze(ncread([Teraname,'a_C.nc'],'mp1',[Loc_x,Loc_y,1],[1,1,Inf],[1,1,1]));
        TmEra = [TmEra(1:end-4,1);squeeze(ncread([Teraname,'b_C.nc'],'mp1',[Loc_x,Loc_y,1],[1,1,Inf],[1,1,1]))];
                
        % timeseries
        figure(ii);
        plot(timev,TmSW,'b',timee,TmEra,'r',timeb,BUOY(:,4),'k','LineWidth',2)
        hold on;
        aa = corrcoef([Tmean(1:DT:end),TmSW(1:DT:end),TmEra(timee>=strdate&timee<=max(timev))]);
        bb = var(TmSW-Tmean);
        text(min(timev)+11,1,[num2str(aa(1,2),'corr(SWAN) : %1.2f '),num2str(bb,' | var : %1.2f')],'Color','b');
        hold on;
        bb = var(TmEra(timee>=strdate&timee<=max(timev))-Tmean(1:DT:end));
        text(min(timev)+11,3.5,[num2str(aa(1,3),'corr(ERA) : %1.2f '),num2str(bb,' | var : %1.2f')],'Color','r');
        hold on;
        bb = var(TmSW(1:DT:end)-TmEra(timee>=strdate&timee<=max(timev)));
        text(min(timev)+11,2.25,[num2str(aa(2,3),'corr(ERA-SWAN) : %1.2f '),num2str(bb,' | var : %1.2f')],'Color',[1 0 1]);
        hold on;
        patch([strdate,strdate,strdate,strdate],...
            [0,0,ceil(max([BUOY(:,4);TmEra]))+1,ceil(max([BUOY(:,4);TmEra]))+1],...
            zeros(4,1),'c','EdgeColor','none','FaceAlpha',0.25)
        hold off;
        set(gca,'XTick',timev(1:14*24:end),'XTickLabel',datestr(timev(1:14*24:end)','dd/mm'))
        grid on;
        set(gcf,'Position',[1 200 1300 400]);
        set(gca,'Position',[0.05 0.11 0.9 .8]);
        axis([min(timev) max(timev) 0 ceil(max([BUOY(:,4);TmEra]))+1])
        h = legend('SWAN (0.25)','ERA (0.25)','BUOY');
        set(h,'Location','SouthEast')
        title('Spectral Mean Wave Period (Tm01)','FontSize',14)        
        saveas(gcf,['FIGURE/BUOY_Tm01Series_',num2str(simyr),'.fig'])
        pause(0.1)
        
        if plott==1
            % Q-Q Plot
            [fTm,vTm] = ecdf(TmSW(timev>=strdate)');

            figure(ii+10);
            title('Spectral Mean Wave Period (Tm01)')        
%             subplot(1,2,1)
%             plot(vTm,fTm,'b',vT,fT,'k','LineWidth',2)
%             h = legend('SWAN','Buoy'); set(h,'Location','SouthEast')
            subplot(1,2,1)
            h = qqplot(BUOY(timeb>=strdate&timeb<=endate,4),TmSW);
            set(h,'MarkerSize',3);
            hold on;
            plot([4 14],[4 14],'--k')
            hold off; 
            box on; axis([4 14 4 14])
            xlabel('Buoy'); ylabel('SWAN')
            title('Tm01 : SWAN vs. Buoy')
            subplot(1,2,2)
            h = qqplot(BUOY(:,4),TmEra);
            set(h,'MarkerSize',3);
            hold on;
            plot([4 14],[4 14],'--k')
            hold off; 
            box on; axis([4 14 4 14])
            xlabel('Buoy'); ylabel('ERA')
            title('Tm01 : ERA vs. Buoy')
            set(gcf,'Position',[403 310 660 356])
            saveas(gcf,['FIGURE/BUOY_Tm01QQ_',num2str(simyr),'.fig'])
            pause(0.1)
        end
    elseif strcmp(varname{ii},'Tp')
        TpSW  = SWAN(n_sta:Nswan:end,ii);
                
        % timeseries
        figure(ii);
        plot(timev,TpSW,'b',timev,Tpeak,'k','LineWidth',2)
        hold on;
        aa = corrcoef(Tpeak(timev>=strdate),TpSW(timev>=strdate));
        bb = var(TpSW(timev>=strdate)-Tpeak(timev>=strdate));
        text(min(timev)+11,1,[num2str(aa(1,2),'corr(SWAN) : %1.2f'),num2str(bb,'| var : %1.2f')],'Color','b');
        hold on;
        patch([strdate,strdate,strdate,strdate],...
            [0,0,ceil(max(BUOY(:,3)))+1,ceil(max(BUOY(:,3)))+1],...
            zeros(4,1),'c','EdgeColor','none','FaceAlpha',0.25)
        hold off;
        set(gca,'XTick',timev(1:14*24:end),'XTickLabel',datestr(timev(1:14*24:end)','dd/mm'))
        grid on;
        set(gcf,'Position',[1 200 1300 400]);
        set(gca,'Position',[0.05 0.11 0.9 .8]);
        axis([min(timev) max(timev) 0 ceil(max(BUOY(:,3)))+1])
        h = legend('SWAN (0.25)','BUOY');
        set(h,'Location','SouthEast')
        title('Spectral Peak Wave Period (Tp)','FontSize',14)        
        saveas(gcf,['FIGURE/BUOY_TpSeries_',num2str(simyr),'.fig'])
        pause(0.1)
        
        if plott==1
            % Q-Q Plot
            figure(ii+10);
            title('Spectral Peak Wave Period (Tp)')        
%             subplot(1,2,1)
%             plot(vTm,fTm,'b',vT,fT,'k','LineWidth',2)
%             h = legend('SWAN','Buoy'); set(h,'Location','SouthEast')
            h = qqplot(BUOY(timeb>=strdate&timeb<=endate,3),TpSW(timev>=strdate));
            set(h,'MarkerSize',3);
            hold on;
            plot([4 14],[4 14],'--k')
            hold off; 
            box on; axis([4 14 4 14])
            xlabel('Buoy'); ylabel('SWAN')
            title('Tpeak : SWAN vs. Buoy')
            saveas(gcf,['FIGURE/BUOY_TpQQ_',num2str(simyr),'.fig'])
            pause(0.1)
        end
    end
end
%%
% clear('kk')
Wu   = SWAN(:,end-1);
Wv   = SWAN(:,end);
Wm   = sqrt(Wu.^2+Wv.^2);
HsWm = 0.24*Wm.^2/9.813;

figure(90)
plot(timev,HsSW,'k',timev,HsWm,'b','LineWidth',2)
set(gca,'XTick',timev(1:24*14:end),'XTickLabel',datestr(timev(1:24*14:end)+1','dd/mm'))
axis([min(timev) max(timev) 0 6])
grid on;
legend('Hs','Hs from Ws')