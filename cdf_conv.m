clear; close all; clc;
addpath(genpath('../../_TOOLBOX'))


load('SWANERA.mat')
Nb = length(BUOYl);

for ii = 1:Nb
    % ERA
    if ~isnan(ERA_v(1,ii+1))
        % Hs 
        [f,a] = ecdf(ERA_v(:,ii+1));                    % : CDF
        eval(['CDF.HSe',num2str(ii,'%03d'),' = [f a];'])
        [g,b] = findpeaks(ERA_v(:,ii+1));               % : POT
        [f,a] = ecdf(g);    
        eval(['POT.HSe',num2str(ii,'%03d'),' = {[f a],[g,b]};'])
        % Tz
        [f,a] = ecdf(ERA_v(:,Nb+ii+1));                 % : CDF
        eval(['CDF.TZe',num2str(ii,'%03d'),' = [f a];'])
        [g,b] = findpeaks(ERA_v(:,Nb+ii+1));            % : POT
        [f,a] = ecdf(g);
        eval(['POT.Tze',num2str(ii,'%03d'),' = {[f a],[g,b]};'])
        % Direction
        [f,a] = ecdf(ERA_v(:,2*Nb+ii+1));               % : CDF
        eval(['CDF.De',num2str(ii,'%03d'),' = [f a];'])
        [g,b] = findpeaks(ERA_v(:,2*Nb+ii+1));          % : POT
        [f,a] = ecdf(g);
        eval(['POT.De',num2str(ii,'%03d'),' = {[f a],[g,b]};'])
    end
    
    % SWAN
    if ~isnan(SWN_v(1,ii+1))
        % Hs
        [f,a] = ecdf(SWN_v(:,ii+1));
        eval(['CDF.HSs',num2str(ii,'%03d'),'= [f a];'])
        [g,b] = findpeaks(SWN_v(:,ii+1));               % : POT
        [f,a] = ecdf(g);    
        eval(['POT.HSs',num2str(ii,'%03d'),' = {[f a],[g,b]};'])
        % Tz
        [f,a] = ecdf(SWN_v(:,Nb+ii+1));
        eval(['CDF.TZs',num2str(ii,'%03d'),'= [f a];'])
        [g,b] = findpeaks(SWN_v(:,Nb+ii+1));            % : POT
        [f,a] = ecdf(g);
        eval(['POT.Tze',num2str(ii,'%03d'),' = {[f a],[g,b]};'])
        % Direction
        [f,a] = ecdf(SWN_v(:,2*Nb+ii+1));
        eval(['CDF.Ds',num2str(ii,'%03d'),'= [f a];'])
        [g,b] = findpeaks(SWN_v(:,2*Nb+ii+1));          % : POT
        [f,a] = ecdf(g);
        eval(['POT.De',num2str(ii,'%03d'),' = {[f a],[g,b]};'])
    end
        
end

%% 
clc; rmse = []; cic=[];
for i= 1:Nb%9%
    Nd      = length(SWN_v(:,i+1));
    
    % first moment
    av1(i)  = mean(SWN_v(:,i+1));  %     av2(i)  = sum(SWN_v(:,i+1))/Nd;

    
    % second moment
    var2(i) = (sum(SWN_v(:,i+1).^2)-2*av1(i)*sum(SWN_v(:,i+1))+Nd*av1(i).^2)/Nd;%     var1(i) = var(SWN_v(:,i+1));%     var3(i) = (sum(SWN_v(:,i+1).^2)-2*av2(i)*sum(SWN_v(:,i+1))+Nd*av2(i).^2)/Nd;
    SStot   = var2(i)*Nd;
    
    
    % third moment
    m3      = (sum(SWN_v(:,i+1).^3)-3*av1(i)*sum(SWN_v(:,i+1).^2)+ ...
        3*av1(i).^2*sum(SWN_v(:,i+1))-Nd*av1(i).^3)/Nd;
    sk2(i)  = m3./(var2(i).^(3/2));%     sk1(i)  = skewness(SWN_v(:,i+1));

    
    % fourth moment
    m4      = (sum(SWN_v(:,i+1).^4)-4*av1(i)*sum(SWN_v(:,i+1).^3)+ ...
        6*av1(i).^2*sum(SWN_v(:,i+1).^2)-4*av1(i).^3*sum(SWN_v(:,i+1))+ ...
        Nd*av1(i).^4)/Nd;
    ku2(i)  = m4./(var2(i).^2);%     ku1(i)  = kurtosis(SWN_v(:,i+1));


    % extract and plot 
    eval(['val1 = CDF.HSs',num2str(i,'%03d'),';']);  
    eval(['val2 = POT.HSs',num2str(i,'%03d'),'{2};']);
    
    
    % POT
    Rayleigh2 = '1-exp(-(x-b).^2/(2*a.^2))';
    GP        = '1-(1+a*x).^(-1/a)';
    th1     = linspace(0,1,1001);
    yi      = rand(100000,1);
    x2      = interp1(val1(:,1),val1(:,2),yi,'linear');
    
    for j=931%1:1001%
               
        th      = interp1(val1(:,1),val1(:,2),th1(j),'spline');
        x       = x2(x2>=th,1);
        [f,a]   = ecdf(x);
            
%         ft = fittype(Rayleigh2,...
%             'dependent',{'y'},'independent',{'x'},...
%             'coefficients',{'a'}); % ,'b'
        ft  = ('Generalized Extreme Value'); %Pareto');% 
%         ft = 'Lognormal';
        %     f1 = fit(val1(:,2),val1(:,1),ft,'StartPoint',[1,1],...
        %         'Lower',[1e-12,0],'Upper',[Inf,f(end)]);
%         f2 = fit(a,f,ft,'StartPoint',[1]);
        f2 = fitdist(x,ft);
        rmse(i) = sqrt(mean((f-cdf(f2,a)).^2));
        R2(i)   =  
        
        ci = paramci(f2);
        fc = th1(j)+f*(1-th1(j));
        mc = th1(j)+cdf(f2,a)*(1-th1(j));
        ml = th1(j)+cdf(makedist(ft,ci(1,1),ci(1,2),ci(1,3)),a)*(1-th1(j)); %
        mu = th1(j)+cdf(makedist(ft,ci(2,1),ci(2,2),ci(2,3)),a)*(1-th1(j)); %
        
        cic{i} = ci;      
        figure(i);
        plot(a,fc,'r',a,mc,'c','LineWidth',2);
        hold on;
        % confidence interval 95%
        plot(a,ml,'-.m',a,mu,'-.m')
        hold on; plot(val1(:,2),val1(:,1),'b');
        hold off
%         set(gca,'YTick',linspace(0,1,10),'YTickLabel',num2str(linspace(th1(j),1,10)','%1.4f'))
        title(sprintf('mean: %1.2f -- var: %1.2f -- skew: %1.2f -- kurt: %1.2f',av1(i),var2(i),sk2(i),ku2(i)))
    end
end
%%
figure(2017)
plot(sqrt(ku2),rmse,'or'); hold on;
text(sqrt(ku2),rmse,num2str([1:10]')); hold off;
