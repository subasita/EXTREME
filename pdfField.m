%% Plot results
clear; close all; clc;
addpath(genpath('/Users/user/Documents/nunuk/_TOOLBOX/GENERAL/m_map'));

plotting = 0;
saveanim = 0;
dt       = 6;
stCal    = 10*24; % start calculation after 10 days
spanVar  = 36*4;
spanTm   = linspace(0,15,spanVar);
spanTp   = linspace(0,25,spanVar);
spanH    = linspace(0,10,spanVar);
spanD    = linspace(0,360,spanVar);

%% Plot
load('16_BASE1.mat')
aa = zeros(size(Xp)); aa = aa(:);
aa = aa*zeros(size(spanTm));
SWNpdf.HS = aa;
SWNpdf.TM = aa;
SWNpdf.TP = aa;
SWNpdf.Dir = aa;


simyr  = [2007,2008,2010,2011,2013:2016]; 

% a2 = '025';
for nn1 = 1 : length(simyr)
    %     for a1m = 1:2
    %         if a1m == 1; a1 = 11; else a1 = 05; end
    %
    %         if a1 == 11; a0 = simyr-1; a4 = 'a'; else a0 = simyr; a4 = 'b'; end
    %         if strcmp(a2,'025'); a5 = '_C'; else a5 = ''; end
    %         a3 = num2str(simyr);
    %
    datefig  = [simyr(nn1)-1,11,01,0,0,0];
    %         ERApath  = '/Users/user/Documents/nunuk/_DATA/ERA/'; %025/H025_2008a_C.nc';
    %         Heraname = [ERApath,a2,'/H',a2,'_',a3,a4,a5,'.nc'];
    %         Teraname = [ERApath,a2,'/T',a2,'_',a3,a4,a5,'.nc'];
    %         Deraname = [ERApath,a2,'/D',a2,'_',a3,a4,a5,'.nc'];
    
    tic;
    i = 1;
    
%     load([num2str(mod(simyr,100),'%02d'),'_MAIN1.mat'])
    load([num2str(mod(simyr(nn1),100),'%02d'),'_MAIN1.mat'])
    
    %         lon   = ncread(Heraname,'longitude');
    %         lat   = ncread(Heraname,'latitude');
    %         tt    = double(ncread(Heraname,'time'))/24+datenum(1900,1,1,0,0,0);
    
    %         ERApdf = zeros(length(lon)*length(lat),1);
    %         ERApdf = ERApdf*zeros(size(spanVar));
    
    while exist(['Tm01_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS')],'var')
        eval(['HS = Hsig_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Hsig_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        eval(['TP = TPsmoo_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''TPsmoo_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        eval(['TM = Tm01_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Tm01_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        eval(['Dir = Dir_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Dir_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        
        % start calculation after warmed-up
        if i>=stCal
               % input to PDF
                HS = HS(:); 
                TM = TM(:); 
                TP = TP(:); 
                Dir = Dir(:);
                
                for ii = 1:spanVar
                    if ii<spanVar
                        SWNpdf.HS(:,ii)  = SWNpdf.HS(:,ii)  + 1*(HS>=spanH(ii) &HS<spanH(ii+1));
                        SWNpdf.TM(:,ii)  = SWNpdf.TM(:,ii)  + 1*(TM>=spanTm(ii) &TM<spanTm(ii+1));
                        SWNpdf.TP(:,ii)  = SWNpdf.TP(:,ii)  + 1*(TP>=spanTp(ii) &TP<spanTp(ii+1));
                        SWNpdf.Dir(:,ii) = SWNpdf.Dir(:,ii) + 1*(Dir>=spanD(ii)&Dir<spanD(ii+1));                        
                    else
                        SWNpdf.HS(:,ii)  = SWNpdf.HS(:,ii)  + 1*(HS>=spanH(ii));
                        SWNpdf.TM(:,ii)  = SWNpdf.TM(:,ii)  + 1*(TM>=spanTm(ii));
                        SWNpdf.TP(:,ii)  = SWNpdf.TP(:,ii)  + 1*(TP>=spanTp(ii));
                        SWNpdf.Dir(:,ii) = SWNpdf.Dir(:,ii) + 1*(Dir>=spanD(ii));
                    end
                end
        end
        %   Pick data at Location
        if plotting==1
            HS = aa((Yp==-17.5),:);
            fig1 = figure(1);
            pcolor(spanH,Xp(1,:),HS); shading interp;
            
            if saveanim == 1
                filename = 'Animation00_TP.gif';
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
    %     end
end
%%
scanning = 'Lat';
Nn   = 2000;
Amap = parula(Nn);
xx   = linspace(1,0,Nn-1);
xx   = [50.^xx,0];
yy   = linspace(xx(1),xx(end),Nn);
Amap = interp1(yy,Amap,xx);
Amap(1,:) = [0.8,0.8,0.8];

if strcmp(scanning,'Lat')
    filename = 'PDF_TM_vs_Lon.gif';
    Ld = length(Yp(:,1));
    ii = length(Yp(Yp(:,1)>=-11,1));
else
    filename = 'PDF_TM_vs_Lat.gif';
    Ld = length(Xp(1,:));
    ii = length(Xp(Xp(1,:)<=87.5,1));
end
% cstdata  = load('D:\SMT01\BLTM\Exercises\Handin\BOWEN\WorldCoastline.dat');

if ii>1;
    if strcmp(scanning,'Lat')
        HS_ = SWNpdf.TM((Yp==Yp(ii,1)),:);
    else
        HS_ = SWNpdf.TM((Xp==Xp(1,ii)),:);
    end
    
    HS_ = HS_/max(sum(SWNpdf.HS,2))*100;
    
    fig1 = figure(1);
    subplot(1,2,1)
    eval(['pcolor(Xp,Yp,Depth',datestr(datenum(datefig),'_yyyymmdd_HHMMSS'),');'])
    shading interp;
    %     hold on; plot(cstdata(:,1),cstdata(:,2),'k');
    hold on;
    if strcmp(scanning,'Lon')
        plot(Xp(:,ii),Yp(:,ii),'r','LineWidth',2);
    elseif strcmp(scanning,'Lat')
        plot(Xp(ii,:),Yp(ii,:),'r','LineWidth',2);
    end
    title('Bathymetric Map of Site')
    
    hold off;
    axis equal;
    axis([min(Xp(:)) max(Xp(:)) min(Yp(:)) max(Yp(:))])
    xlabel('Lon (^\circ)');ylabel('Lat (^\circ)')
    
    subplot(1,2,2)
    if strcmp(scanning,'Lat')
        pcolor(spanH,Xp(1,:),HS_); shading interp;
        title(['Occurences of Tm_{01} at Lat: ',num2str(Yp(ii,1),'%1.2f'),'^\circ in %']);
        ylabel('Lon (^\circ)')
    elseif strcmp(scanning,'Lon')
        pcolor(spanH,Yp(:,1),HS_); shading interp;
        title(['Occurences of Tm_{01} at Lon: ',num2str(Xp(1,ii),'%1.2f'),'^\circ in %']);
        ylabel('Lat (^\circ)')
    end
    cb = colorbar; title(cb, '%'); colormap(Amap);
%     xlabel('Dir (^\circ)');
    xlabel('Tm_{01} (s)');
%     xlabel('Hs (m)');
    set(fig1,'renderer','zbuffer','Color','w','Position',[2 42 1220 510]);
    caxis([0 5])
    drawnow
    frame = getframe(fig1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256,'nodither');
    
    if ii == length(Xp(Xp(1,:)<=87.5,1))|| ii == length(Yp(Yp(:,1)>=-11,1))
        imwrite(imind,cm,filename,'gif', 'loopcount',inf,'DelayTime',0.25);
    else
        imwrite(imind,cm,filename,'gif','writemode','append','DelayTime',0.25);
    end
    hold off; pause(0.01) 
    
end
