function show_MixRHLP_results(data,solution)
set(0,'defaultaxesfontsize',14);
[n, m] = size(data);
% fs = m/100
% t=1/fs:1/fs:m/fs;
t = 0:m-1;
% t=1:m;%

G = length(solution.param.alpha_g);
colors = {'r','b','g','m','c','k','y'};
% colors_cluster_means = {'m','c','k','b','r','g'};
colors_cluster_means = {[0.8 0 0],[0 0 0.8],[0 0.8 0],'m','c','k','y'};

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 560 scrsz(4)/2.15]);
plot(t, data');
title('Original time series')
ylabel('y')
xlabel('Time')

figure('Position',[scrsz(4)/2 scrsz(4)/2 560 scrsz(4)/2.15]);
%%%% cluster and means
for g=1:G
    cluster_g = data(solution.klas==g,:);
    plot(t,cluster_g','color',colors{g},'linewidth',0.001);    
    hold on
%     plot(t,solution.Ex_g(:,g),'color',colors{g},'linewidth',3)
    plot(t,solution.Ex_g(:,g),'color',colors_cluster_means{g},'linewidth',3)
end
title('Clustered and segmented time series')
ylabel('y')
xlabel('Time')

xlim([0 m-1])

% 
%    
% 
% %%%  means and logistic proportions
% for g=1:G
%     cluster_g = data(solution.klas==g,:); 
%    
%     figure(g+1)
%     subplot(211),plot(t,cluster_g','color',colors{g});    
%     hold on
% %     plot(t,solution.Ex_g(:,g),'color',colors{g},'linewidth',3)
%     plot(t,solution.Ex_g(:,g),'color',colors_cluster_means{g},'linewidth',3)
%     title(['Cluster ',int2str(g)])
%     ylabel('y');
%     %
%     %
%     subplot(212),plot(t,solution.param.pi_jgk(1:m,:,g),'linewidth',2);
%     ylabel('Logistic proportions')%['\pi_{jk}( w^g) , g = ',int2str(g)])
%     xlabel('Time')
%     set(gca,'ytick',[0:0.2:1])
% end

for g=1:G
    cluster_g = data(solution.klas==g,:); 
    
    %%%     figure(g+1)
    figr = figure('Position',[scrsz(4) scrsz(4)/2 560 scrsz(4)/1.4]);
    axes1 = axes('Parent',figr,'Position',[0.13 0.45 0.775 0.48],'FontSize',14);
    box(axes1,'on'); hold(axes1,'all');
    %
    plot(t,cluster_g','color',colors{g})
    hold on
    % polynomial regressors
    plot(t,solution.polynomials(:,:,g),'k--','linewidth',1)
    
    hold on,
    %
    plot(t,solution.Ex_g(:,g),'color',colors_cluster_means{g},'linewidth',3)
    title(['Cluster ',int2str(g)])
    ylabel('y');
    xlim([0 m-1])
    %
    %
    axes2 = axes('Parent',figr,'Position',[0.13 0.08 0.775 0.3],'FontSize',14);
    box(axes2,'on'); hold(axes2,'all');
    plot(t,solution.param.pi_jgk(1:m,:,g),'linewidth',2);
    ylabel('Logistic proportions')%['\pi_{jk}( w^g) , g = ',int2str(g)])
    xlabel('Time')
    set(gca,'ytick',[0:0.2:1])
    
     xlim([0 m-1])
end
 