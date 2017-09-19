%% Plot results
clear; close all; clc;
addpath(genpath('../../_TOOLBOX/GENERAL'));

season   = 'allyear';
plotting = 0;
saveanim = 0;
dt       = 6;
stCal    = 10*24; % start calculation after 10 days
maxv     = 25;
spanVar  = maxv*5+1; % number of division in:
spanTm   = linspace(0,maxv,spanVar); % mean wave period
spanTp   = linspace(0,maxv,spanVar); % peak wave period
spanH    = linspace(0,maxv,spanVar); % significant wave height
spanD    = linspace(0,360,spanVar); % main wave direction 
spanP    = [10.0,50.0,90.0,95.0,99.0,99.9]; % in percent
filename2load = ['Dist_11years_',season,'.mat'];

%% Plot
% Initial state
load(filename2load);

%% Plot of percentile distribution
HS1 = ones(size(Xp)); HS1 = HS1 * NaN;
for ii = 1 : length(spanP)
    HS1(wb) = SWNdis.HS(:,ii);
    
    figure(200+ii);
    pcolor(Xp,Yp,HS1); shading interp; axis tight equal;
    colormap(jet(128)); caxis([0 5+ii]); colorbar;
    title([num2str(spanP(ii),'%2.1f'),'-^{th} Percentile'])
    saveas(gcf,['FIGURE/HS_',num2str(spanP(ii),'%2.1f-'),season,'.fig'])
    pause(0.5)
end

%%  Calculate moment
N    =  max(sum(SWNpdf.HS,2));
HS2  = ones(size(Xp)); HS2 = HS2 * NaN;

%   Mean
HS2(wb) = mom.Yavg(:,1); 

figure(11);
pcolor(Xp,Yp,HS2); shading interp; axis tight equal;
colormap(jet); colorbar;
caxis([0 3])
title('Average of the distribution of significant wave height (2006-2016)')
saveas(gcf,['FIGURE/Moment_Average.fig'])

%   Variance
HS2(wb) = mom.Yvar(:,1); 

figure(12);
pcolor(Xp,Yp,HS2); shading interp; axis tight equal;
colormap(jet); colorbar;
title('Variance of the distribution of significant wave height (2006-2016)')
saveas(gcf,['FIGURE/Moment_Variance.fig'])

%   Standard deviation
HS2(wb) = mom.Ystd(:,1); 

figure(13);
pcolor(Xp,Yp,HS2); shading interp; axis tight equal;
colormap(jet); colorbar;
title('Standard deviation of the distribution of significant wave height (2006-2016)')
saveas(gcf,['FIGURE/Moment_StDev.fig'])

%   Skewness
HS2(wb) = mom.Yske(:,1); 

figure(14);
pcolor(Xp(1,:),Yp(:,1),HS2); shading interp; axis tight equal;
colormap(jet); colorbar; caxis([0 5])
title('Skewness of the distribution of significant wave height (2006-2016)')
saveas(gcf,['FIGURE/Moment_Skewness.fig'])

%   Kurtosis
HS2(wb) = mom.Ykur(:,1); 

figure(15);
pcolor(Xp,Yp,HS2); shading interp; axis tight equal;
colormap(jet); colorbar; caxis([0 40])
title('Kurtosis of the distribution of significant wave height (2006-2016)')
saveas(gcf,'FIGURE/Moment_Kurtosis.fig')

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
    filename = ['FIGURE/PDF_TM_vs_Lon-',season,'.gif'];
    Ld = length(Yp(:,1));
    ii = length(Yp(Yp(:,1)>=-20,1));
else
    filename = ['FIGURE/PDF_TM_vs_Lat-',season,'.gif'];
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
    xlabel('Lon (^\circ)'); ylabel('Lat (^\circ)')
    
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


%% plot at random position
clc; close all;
PS  = [90.0, 20; 87.5, 0.0; 85.0, -20.0; 122.0, -15.0; 110.0, -5; 132.0, 25; 115.0, 18.0];
dH2 = (spanH(2)-spanH(1))/2;
y   = [jet(50);hot(76)];
minx = [3.5,3,4,3,2,5,4.5];
for ii = 1:length(PS(:,1))
    x  = find(xx==PS(ii,1)&yy==(PS(ii,2)));
    
    % convert to CDF
    HS_ = (SWNpdf.HS(x,:))/max(sum(SWNpdf.HS,2));    
    HS1 = cumsum(HS_); % CDF
    
    
    % CDF of peaks
%     HS2   = sort(HSp{x,:});
%     [f,v] = ecdf(HS2);
    % fit to weibull
%     im  = fit(spanH',HS1','1 - exp(-((x)/a).^b)','StartPoint',[1,1],'Lower',[1,1],'Upper',[5,10])
    modelFun =  @(p,x) 1 - exp( -( (x) /p(1) ) .^ p(2));
    im   = nlinfit(spanH', HS1', modelFun,[5 3]);
    HS3 = 1-exp(-((spanH)/im(1)).^im(2));
    HS2 = [0 HS3(2:end)-HS3(1:end-1)]; 
    
    % fit to GDP to POT
%     im  = fit(double(v),double(f),'1-(1+b*(x-c)/a).^(-1/b)','StartPoint',[1,1,-1],'Lower',[0,0,-10],'Upper',[6,12,0])
%     HS3 = 1-(1+im.b*(spanH-im.c)/im.a).^(-1/im.b);
%     HS2 = [0 HS3(2:end)-HS3(1:end-1)];
    
    % Histogram HS
    figure(100+ii)
    for i = 1 : length(spanH)
        X = [spanH(i)-dH2,spanH(i)+dH2,spanH(i)+dH2,spanH(i)-dH2];
        Y = [ 0, 0, HS_(1,i)*100,HS_(1,i)*100]; 
        patch(X,Y,y(i,:));
        hold on;
    end
    
    text(5,10,['F(x) = ',num2str(im(2)/im(1),'%1.2f'),'(x/',num2str(im(1),'%1.2f'),')^{',num2str(im(2),'%1.2f'),'-1} exp(x/',num2str(im(1),'%1.2f'),')^{',num2str(im(2),'%1.2f'),'}']);
    hold on;
    plot(spanH,HS2*100,'m','LineWidth',2);
    hold off;
%     axis([minx(ii) minx(ii)+5 0 1]); %ceil(max(HS_(1,:)))*100]);
    axis([0 10 0 ceil(max(HS_(1,:)*10))*10]);
    xlabel('Significant Wave Height [m]'); ylabel('Percentage Occurrence [%]')
    box on;
    title(['Percentage of occurrence of significant wave height (Hs) during 2006-2016 ',num2str(ii,'[%01d]')]);
    saveas(gcf,['FIGURE/PDF location ',num2str(ii),'-',season,'.fig']);
end

