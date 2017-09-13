%% Plot results
clear; close all; clc;
addpath(genpath('../../_TOOLBOX/GENERAL'));

plotting = 0;
saveanim = 0;
dt       = 6;
season   = 'allyear';
maxv     = 25;
spanVar  = maxv*5+1; % number of division in:
spanTm   = linspace(0,maxv,spanVar); % mean wave period
spanTp   = linspace(0,maxv,spanVar); % peak wave period
spanH    = linspace(0,maxv,spanVar); % significant wave height
spanD    = linspace(0,360,spanVar); % main wave direction 
spanP    = [10.0,50.0,90.0,95.0,99.0,99.9]; % in percent
savename = ['Dist_11Years_',season,'.mat'];

%% Plot
% Initial state
simyr  = 2006:2016;
load([num2str(mod(simyr(end),100),'%02d'),'_BASE1.mat'])
eval(['Dp = Depth',datestr(datenum(simyr(end)-1,11,1,0,0,0),'_yyyymmdd_HHMMSS'),'; clear(''Depth',datestr(datenum(simyr(end)-1,11,1,0,0,0),'_yyyymmdd_HHMMSS'),''')'])
eval('wb = zeros(size(Xp(~isnan(Dp))));');
wb = wb(:); wb = wb*zeros(size(spanTm));
SWNpdf.HS  = wb;
SWNpdf.TM  = wb;
SWNpdf.TP  = wb;
SWNpdf.Dir = wb;

% number of data point (which is regarded as ocean)
eval('wb = find(~isnan(Dp));');


% preparation for calculation of moment
EY_1 = zeros(size(wb))*zeros(1,4); % sum(y)
EY_2 = zeros(size(wb))*zeros(1,4); % sum(y^2)
EY_3 = zeros(size(wb))*zeros(1,4); % sum(y^3)
EY_4 = zeros(size(wb))*zeros(1,4); % sum(y^4)

for ni = 1 : length(simyr)
    tic;
    i = 1;
    datefig  = [simyr(ni)-1,11,01,0,0,0]; % date number
    if strcmp(season,'Summer')
        stCal = (datenum(simyr(ni),5,1,0,0,0)-datenum(datefig))*24;
    else
        stCal = 0;
    end

    % load SWAN result
    load([num2str(mod(simyr(ni),100),'%02d'),'_MAIN1.mat'])
    
    % filtering
    while exist(['Tm01_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS')],'var')
        eval(['HS  = Hsig_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Hsig_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        eval(['TP  = TPsmoo_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''TPsmoo_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        eval(['TM  = Tm01_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Tm01_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        eval(['Dir = Dir_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),'; clear(''Dir_',datestr(datenum(datefig+[0,0,0,i-1,0,0]),'yyyymmdd_HHMMSS'),''')'])
        
        % start calculation during specific year
%         if datenum(datefig+[0,0,0,i-1,0,0])<=datenum(simyr(ni),5,1,0,0,0) % winter NH
%         if datenum(datefig+[0,0,0,i-1,0,0])>=datenum(simyr(ni),5,1,0,0,0) % summer NH
        if i-1>=stCal; % all year  

            % input to PDF
            HS  = HS(wb);
            TM  = TM(wb);
            TP  = TP(wb);
            Dir = Dir(wb);
            
            % group into PDF
            for ii = 1:spanVar
                if ii<spanVar
                    SWNpdf.HS(:,ii)  = SWNpdf.HS(:,ii)  + 1*(HS>=spanH(ii) &HS<spanH(ii+1));
                    SWNpdf.TM(:,ii)  = SWNpdf.TM(:,ii)  + 1*(TM>=spanTm(ii) &TM<spanTm(ii+1));
                    SWNpdf.TP(:,ii)  = SWNpdf.TP(:,ii)  + 1*(TP>=spanTp(ii) &TP<spanTp(ii+1));
                    SWNpdf.Dir(:,ii) = SWNpdf.Dir(:,ii) + 1*(Dir>=spanD(ii)&Dir<spanD(ii+1));
                else % at boundary
                    SWNpdf.HS(:,ii)  = SWNpdf.HS(:,ii)  + 1*(HS>=spanH(ii));
                    SWNpdf.TM(:,ii)  = SWNpdf.TM(:,ii)  + 1*(TM>=spanTm(ii));
                    SWNpdf.TP(:,ii)  = SWNpdf.TP(:,ii)  + 1*(TP>=spanTp(ii));
                    SWNpdf.Dir(:,ii) = SWNpdf.Dir(:,ii) + 1*(Dir>=spanD(ii));
                end
            end
            
            % group of peaks for POT
            if i - 1 == stCal
                HSp = cell(size(HS));  % initialize the peak HS
                HS1 = HS;
            elseif i - 1 == stCal + dt
                HS2  = HS;
            elseif i -1 >= stCal + 2*dt
                HS3  = HS;
                
                % remarks every point which is a peak
                pi = find( HS3<=HS2 & HS2>=HS1 );
                for ip = 1 : length(pi)
                    HSp{pi(ip)} = [HSp{pi(ip)}; HS2(pi(ip))];
                end
                
                % back-substitution
                HS1  = HS2;
                HS2  = HS3;
            end
            
            % Calculate moment
            EY_1 = EY_1 + [HS, TM, TP, Dir]; % for moment 1
            EY_2 = EY_2 + [HS.^2, TM.^2, TP.^2, Dir.^2]; % for moment 2
            EY_3 = EY_3 + [HS.^3, TM.^3, TP.^3, Dir.^3]; % for moment 3
            EY_4 = EY_4 + [HS.^4, TM.^4, TP.^4, Dir.^4]; % for moment 4
            
        end
        
        % subtitute variables
        i = i + dt;
    end
    disp([num2str(simyr(ni)),'-',num2str(toc)])
end

%%  Calculate moment
N    =  max(sum(SWNpdf.HS,2));
HS2  = ones(size(Xp)); HS2 = HS2 * NaN;

%   Mean
Yavg    = EY_1./N;
HS2(wb) = Yavg(:,1); 

figure;
pcolor(Xp,Yp,HS2); shading interp; axis tight equal;
colormap(jet); colorbar;
caxis([0 3])
title('Average of the distribution of significant wave height (2006-2016)')

%   Variance
Yvar = (EY_2 - 2*Yavg.*EY_1 + N*Yavg.^2 )./N;
HS2(wb) = Yvar(:,1); 

figure;
pcolor(Xp,Yp,HS2); shading interp; axis tight equal;
colormap(jet); colorbar;
title('Variance of the distribution of significant wave height (2006-2016)')

%   Standard deviation
Ystd    = sqrt(Yvar);
HS2(wb) = Ystd(:,1); 

figure;
pcolor(Xp,Yp,HS2); shading interp; axis tight equal;
colormap(jet); colorbar;
title('Standard deviation of the distribution of significant wave height (2006-2016)')

%   Skewness
Yske    = (EY_3 - 3*Yavg.*EY_2 + 3*Yavg.^2.*EY_1 - N*Yavg.^3)./N;
Yske    = Yske./(Ystd.^3);
HS2(wb) = Yske(:,1); 

figure;
pcolor(Xp(1,:),Yp(:,1),HS2); shading interp; axis tight equal;
colormap(jet); colorbar; caxis([0 5])
title('Skewness of the distribution of significant wave height (2006-2016)')

%   Kurtosis
Ykur    = (EY_4 - 4*Yavg.*EY_3 + 6*Yavg.^2.*EY_2 - 4*Yavg.^3.*EY_1 + N*Yavg.^4)./N;
Ykur    = Ykur./(Yvar.^2);
HS2(wb) = Ykur(:,1); 

figure;
pcolor(Xp,Yp,HS2); shading interp; axis tight equal;
colormap(jet); colorbar; caxis([0 40])
title('Kurtosis of the distribution of significant wave height (2006-2016)')

%% Detect the location of percentage distribution
SWNdis.HS  = zeros(size(wb))*zeros(size(spanP));
SWNdis.TM  = zeros(size(wb))*zeros(size(spanP));
SWNdis.TP  = zeros(size(wb))*zeros(size(spanP));
SWNdis.Dir = zeros(size(wb))*zeros(size(spanP));

tic;
for ii = 1 : length(wb)
    if sum(SWNpdf.HS(ii,:))~=0
        SWNdis.HS(ii,:)  = pdf2perc(SWNpdf.HS(ii,:)./sum(SWNpdf.HS(ii,:)),spanH,spanP,'PDF');
    end
    if sum(SWNpdf.TM(ii,:))~=0
        SWNdis.TM(ii,:)  = pdf2perc(SWNpdf.TM(ii,:)./sum(SWNpdf.TM(ii,:)),spanH,spanP,'PDF');
    end
    if sum(SWNpdf.TP(ii,:))~=0
        SWNdis.TP(ii,:)  = pdf2perc(SWNpdf.TP(ii,:)./sum(SWNpdf.TP(ii,:)),spanH,spanP,'PDF');
    end
    if sum(SWNpdf.Dir(ii,:))~=0
        SWNdis.Dir(ii,:) = pdf2perc(SWNpdf.Dir(ii,:)./sum(SWNpdf.Dir(ii,:)),spanH,spanP,'PDF');
    end
end
disp(num2str(toc));



HS1 = ones(size(Xp)); HS1 = HS1 * NaN;
for ii = 1 : length(spanP)
    HS1(wb) = SWNdis.HS(:,ii);
    
    figure(200+ii);
    pcolor(Xp,Yp,HS1); shading interp; axis tight equal;
    colormap(jet(128)); caxis([0 5+ii]); colorbar;
    title([num2str(spanP(ii),'%2.1f'),'-^{th} Percentile'])
    saveas(gcf,['HS_',num2str(spanP(ii),'%2.1f'),'.fig'])
    pause(0.5)
end

%% plotting
scanning = 'Lon';
Nn   = 2000;
Amap = bone(Nn);
xx   = linspace(1,0,Nn-1);
xx   = [50.^xx,0];
yy   = linspace(xx(1),xx(end),Nn);
Amap = interp1(yy,Amap,xx);
Amap(1,:) = [0.8,0.8,0.8];

if strcmp(scanning,'Lat')
    filename = 'PDF_TM_vs_Lon.gif';
    Ld = length(Yp(:,1));
    ii = length(Yp(Yp(:,1)>=-15,1));
else
    filename = 'PDF_TM_vs_Lat.gif';
    Ld = length(Xp(1,:));
    ii = length(Xp(1,Xp(1,:)<=87.5));
end
% cstdata  = load('D:\SMT01\BLTM\Exercises\Handin\BOWEN\WorldCoastline.dat');
xx = Xp(wb); yy = Yp(wb);

if ii>1;
    if strcmp(scanning,'Lat')
        HS  = zeros(size(Xp(1,:)'))*zeros(size(spanP))*NaN;
        HS_ = zeros(size(Xp(1,:)'))*zeros(size(spanH))*NaN;
        ip  = xx(yy==Yp(ii,1)); ni = zeros(size(ip));
        for pi = 1:length(ip)
            ni(pi) = find(Xp(1,:)==ip(pi));
        end
        HS_(ni,:) = SWNpdf.HS((yy==Yp(ii,1)),:);
        HS(ni,:)  = SWNdis.HS((yy==Yp(ii,1)),:);
    else
        HS  = zeros(size(Yp(:,1)))*zeros(size(spanP))*NaN;
        HS_ = zeros(size(Yp(:,1)))*zeros(size(spanH))*NaN;
        ip  = yy(xx==Xp(1,ii)); ni = zeros(size(ip));
        for pi = 1:length(ip)
            ni(pi) = find(Yp(:,1)==ip(pi));
        end
        HS_(ni,:) = SWNpdf.HS((xx==Xp(1,ii)),:);
        HS(ni,:)  = SWNdis.HS((xx==Xp(1,ii)),:);
    end
    
    HS_ = HS_/max(sum(SWNpdf.HS,2))*100;
    
    fig1 = figure(1);
    subplot(1,5,1:3)
    pcolor(Xp,Yp,Dp)
    shading interp;
    %     hold on; plot(cstdata(:,1),cstdata(:,2),'k');
    hold on;
    if strcmp(scanning,'Lon')
        plot(Xp(:,ii),Yp(:,ii),'r','LineWidth',2);
    elseif strcmp(scanning,'Lat')
        plot(Xp(ii,:),Yp(ii,:),'r','LineWidth',2);
    end
    title('cross-sectional area (along red-line) in the domain')
    
    hold off;
    axis([min(Xp(:)) max(Xp(:)) min(Yp(:)) max(Yp(:))])
    xlabel('Lon (^\circ)');ylabel('Lat (^\circ)')
    
    subplot(1,5,[4 5])
    if strcmp(scanning,'Lat')
        imagesc(spanH,Xp(1,:),HS_); shading interp;
        hold on; pp = plot(HS,Xp(1,:),'LineWidth',2); hold off;
        set(pp,{'Color'},num2cell(jet(length(spanP)),2))
        legend(num2str(spanP'));
        set(gca,'YDir','normal')
        title(['Occurences of Hs at Lat: ',num2str(Yp(ii,1),'%1.2f'),'^\circ in %']);
        ylabel('Lon (^\circ)')
        axis([min(spanH),max(spanH),min(Xp(1,:)),max(Xp(1,:))])
    elseif strcmp(scanning,'Lon')
        imagesc(spanH,Yp(:,1),HS_); shading interp;
        hold on; pp = plot(HS,Yp(:,1),'LineWidth',2); hold off;
        set(pp,{'Color'},num2cell(jet(length(spanP)),2))
        legend(num2str(spanP'));
        set(gca,'YDir','normal')
        title(['Occurences of Hs at Lon: ',num2str(Xp(1,ii),'%1.2f'),'^\circ in %']);
        ylabel('Lat (^\circ)')
        axis([min(spanH),max(spanH),min(Yp(:,1)),max(Yp(:,1))])
    end
    cb = colorbar; title(cb, '%'); colormap(Amap);
    %     xlabel('Dir (^\circ)');
    %    xlabel('Tm_{01} (s)');
    xlabel('Hs (m)');
    set(fig1,'renderer','zbuffer','Color','w','Position',[2 42 1220 510]);
    
    caxis([0 10])
    drawnow
    frame = getframe(fig1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256,'nodither');
    imwrite(imind,cm,filename,'gif', 'loopcount',inf,'DelayTime',0.25);
    hold off; pause(0.01)
end

%% save important files
span.H   = spanH;
span.TM  = spanTm;
span.TP  = spanTp;
span.D   = spanD;
span.P   = spanP;

mom.Yavg = Yavg;
mom.Yvar = Yvar;
mom.Ystd = Ystd;
mom.Yske = Yske;
mom.Ykur = Ykur;
save(savename,'SWNpdf','SWNdis','HSp','Yavg','Yvar','Ystd','Yske','Ykur','Xp','Yp','Dp','wb','span');