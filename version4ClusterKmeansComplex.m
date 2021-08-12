clear
close all
cd /home/dylananderson/projects/duckGeomorph/

%addpath(genpath('/media/dylananderson/Elements/SERDP'))

load('ceofsSouthernTransect.mat')

numberCenters = 12;
sqrNC = sqrt(numberCenters);

[x1,y1] = pol2cart(phiRadian(:,1),Rt(:,1));
[x2,y2] = pol2cart(phiRadian(:,2),Rt(:,2));
[x3,y3] = pol2cart(phiRadian(:,3),Rt(:,3));
[x4,y4] = pol2cart(phiRadian(:,4),Rt(:,4));

data(1:length(Rt),1) = x1*(percentV(1)/percentV(1));
data(:,2) = y1*(percentV(1)/percentV(1));
data(:,3) = x2*(percentV(2)/percentV(1)); %Rt(:,2);%*(percentV(2)/percentV(1));
data(:,4) = y2*(percentV(2)/percentV(1)); %phiDegrees(:,2);
data(:,5) = x3*(percentV(3)/percentV(1)); %Rt(:,3)*(percentV(3)/percentV(1));
data(:,6) = x3*(percentV(3)/percentV(1)); %phiDegrees(:,3);
data(:,7) = x4*(percentV(4)/percentV(1)); %Rt(:,3)*(percentV(3)/percentV(1));
data(:,8) = x4*(percentV(4)/percentV(1));

% type = 'kmeans';
% escalar = [1,2,3,4,5,6];%,5];%,7,9,11,13];
% direccional = [];%[2,4];%,6];%,8,10,12,14];
% initmethod = 'mda';
% clusters = makeClustering_initMDA(data, type, numberCenters, escalar, direccional, initmethod);
% 
% 
% 
% clust = zeros(size(data,1),35);
% c = 1;
% for i=1:35
%     clusters = makeClustering_initMDA(data, type, i, escalar, direccional, initmethod);
% 
%     clust(:,c) = clusters.PatternsGroup;
%     c = c+1;
% end
% 
% 
% evaS = evalclusters(data,clust,'Silhouette')
%[IDX,C,SUMD,K] = kmeans_opt(data,35)

clusterType = 'linkage'
evaS = evalclusters(data,clusterType,'Silhouette','klist',[1:35])

figure
plot(evaS.CriterionValues)
title('Silhouette')
% 
% evaCH = evalclusters(data,clust,'CalinskiHarabasz')
evaCH = evalclusters(data,clusterType,'CalinskiHarabasz','klist',[1:35])

figure
plot(evaCH.CriterionValues)
title('Calinski-Harabasz')
% 
% evaDB = evalclusters(data,clust,'DaviesBouldin')
evaDB = evalclusters(data,clusterType,'DaviesBouldin','klist',[1:35])
figure
plot(evaDB.CriterionValues)
title('DaviesBouldin')

evaG = evalclusters(data,clusterType,'Gap','klist',[1:35])
figure
plot(evaG.CriterionValues)
title('Gap')

bmu = clusters.PatternsGroup;

[theta1,rho1] = cart2pol(clusters.Centers(:,1),clusters.Centers(:,2));
centers = theta1*180/pi;%clusters.Centers(:,2);

negcenter = find(centers < 0);
centers(negcenter) = centers(negcenter)+360;
[val, order] = sort(centers,'ascend');
colmap = jet(numberCenters);

% What do the clusters look like with respect to the first offshore
% migrating CEOF?
figure
polarscatter(phiDegrees(:,1)*pi/180,Rt(:,1))
hold on
[theta, rho] = cart2pol(clusters.Centers(:,2),clusters.Centers(:,1));
polarscatter(theta,rho,'r','filled')



%[val, order] = sort(clusters.Centers(:,2),'descend');
%colmap = jet(numberCenters);
figure
for pp = 1:numberCenters
    index = order(pp);%find(order==pp);    
    subplot(15,1,pp)

    plot(xinterp,mean(alllines(clusters.Group{index},:)), '-','color',colmap(pp,:))
    hold on
    text(300,-2,['',num2str(index),''])
end

% for qq = 1:numberCenters
%     figure
%     index = clusters.Group{qq};
%     th = phiRadian(index,1);
%     r = Rt(index,1);
%     polarscatter(th,r,'k')
%     hold on
%     polarscatter(clusters.Centers(qq,2)*pi/180,clusters.Centers(qq,1),'r','filled')
%     
% end

% % What do the clusters look like comparing the leading three magnitues?
% figure
% scatter(clusters.Centers(:,1)/(percentV(1)/percentV(1)),clusters.Centers(:,3)/(percentV(2)/percentV(1)),'o')

figure

for qq = 1:numberCenters
    subplot(4,4,qq)
    profiles = clusters.Group{qq};
    for xx = 1:length(profiles)
        plot(xinterp,alllines(profiles(xx),:),'--','color',[0.7 0.7 0.7])
        hold on
    end
    
    
    plot(xinterp,nanmean(alllines(profiles,:)),'k-','linewidth',1.5)
    hold on
    plot(xinterp,nanmean(alllines(profiles,:)) + nanstd(alllines(profiles,:)),'k--','linewidth',0.5)
    plot(xinterp,nanmean(alllines(profiles,:)) - nanstd(alllines(profiles,:)),'k--','linewidth',0.5)
    plot(xinterp,nanmean(alllines(profiles,:)) + 2*nanstd(alllines(profiles,:)),'r--','linewidth',0.5) %,'color',[0.5 0.5 0.5])
    plot(xinterp,nanmean(alllines(profiles,:)) - 2*nanstd(alllines(profiles,:)),'r--','linewidth',0.5) %,'color',[0.5 0.5 0.5])


end


% 
% 
% temp = clusters.Group{7};
% 
% p=find(diff(temp)==1);
% q=[p;p+1];
% temp(q)  % this gives all the pairs of consecutive numbers
% 
% s = find(diff(temp)>1);
% figure
% 
% for pp = 1:length(s);    
%     subplot(4,3,pp)
% 
%     if pp == 1
%         profiles = [1:s(pp)];
%     else
%         profiles = [s(pp-1)+1:s(pp)];
%     end
%     for xx = 1:length(profiles)
%         plot(xinterp,alllines(temp(profiles(xx)),:))
%         hold on
%     end
%     
% end



% What do the clusters look like with respect to the first offshore
% migrating CEOF?
figure
polarscatter(phiDegrees(:,1)*pi/180,Rt(:,1),10,'k','filled')
hold on
for qq = 1:numberCenters
    temp = order(qq);
    [theta,rho] = cart2pol(clusters.Centers(temp,1),clusters.Centers(temp,2));
    polarscatter(theta,rho,40,colmap(qq,:),'filled')
    
end

figure
%polarscatter(data(:,2)*pi/180,data(:,1),10,'k','filled')
%hold on
for qq = 1:numberCenters
    temp = order(qq);
    profiles = clusters.Group{temp};
    %polarscatter(clusters.Centers(temp,2)*pi/180,clusters.Centers(temp,1),40,colmap(qq,:),'filled')
    polarscatter(phiDegrees(profiles,1)*pi/180,Rt(profiles,1),20,colmap(qq,:),'filled')
    hold on
end
title('Clusters in CEOF R(t) and \theta(t)')


figure
for pp = 1:numberCenters
    index = order(pp);
    plot(xinterp,mean(alllines(clusters.Group{index},:)), '-','color',colmap(pp,:))
    hold on
end
xlim([0 500])
ylabel('Elevation (m)')
xlabel('Cross-shore (m)')
title('Cluster Centroids')
% for qq = 1:numberCenters
%     figure
%     index = clusters.Group{qq};
%     th = phiRadian(index,1);
%     r = Rt(index,1);
%     polarscatter(th,r,'k')
%     hold on
%     polarscatter(clusters.Centers(qq,2)*pi/180,clusters.Centers(qq,1),'r','filled')
%     
% end

% % What do the clusters look like comparing the leading three magnitues?
% figure
% scatter(clusters.Centers(:,1)/(percentV(1)/percentV(1)),clusters.Centers(:,3)/(percentV(2)/percentV(1)),'o')

figure

for qq = 1:numberCenters
    temp = order(qq);
        
    subplot(4,4,temp)

    profiles = clusters.Group{qq};
%     for xx = 1:length(profiles)
%         plot(xinterp,alllines(profiles(xx),:),'--','color',[0.7 0.7 0.7])
%         hold on
%     end
    
    
    plot(xinterp,nanmean(alllines(profiles,:)),'k-','linewidth',1.5)
    hold on
    plot(xinterp,nanmean(alllines(profiles,:)) + nanstd(alllines(profiles,:)),'k--','linewidth',0.5)
    plot(xinterp,nanmean(alllines(profiles,:)) - nanstd(alllines(profiles,:)),'k--','linewidth',0.5)

    if temp == 2
        title('Variability about mean profile')
    end

end



figure
for qq = 1:numberCenters
    index = order(qq);

    profiles = clusters.Group{qq};
    plot(ones(length(profiles),1)*index,time(profiles),'o','color',colmap(index,:))
    hold on
end


figure
for qq = 1:numberCenters
    subplot(4,4,qq)
    profiles = clusters.Group{qq};
    
    polarscatter(phiDegrees(profiles,2)*pi/180,Rt(profiles,2),'b')

end

figure
for qq = 1:numberCenters
    subplot(4,4,qq)
    profiles = clusters.Group{qq};
    
    polarscatter(phiDegrees(profiles,1)*pi/180,Rt(profiles,1),'r')

end

figure
for qq = 1:numberCenters
    subplot(4,4,qq)
    profiles = clusters.Group{qq};
    
    polarscatter(phiDegrees(profiles,3)*pi/180,Rt(profiles,3),'g')

end


figure
for qq = 1:numberCenters
    temp = order(qq);
    [theta,rho] = cart2pol(clusters.Centers(temp,1),clusters.Centers(temp,2));
    if theta < 0
        %plot(rho,(theta+2*pi)*180/pi,'.','color',colmap(qq,:),'markersize',30)
        plot((theta+2*pi)*180/pi,rho,'.','color',colmap(qq,:),'markersize',30)

    else
        %plot(rho,(theta)*180/pi,'.','color',colmap(qq,:),'markersize',30)
        plot((theta)*180/pi,rho,'.','color',colmap(qq,:),'markersize',30)
        
    end
    
    hold on
end
%  N = 1; % Required number of consecutive numbers following a first one
%  x = diff(temp)==1;
%  f = find([false,x]~=[x,false]);
%  g = find(f(2:2:end)-f(1:2:end-1)>=N,1,'first');
%  first_t = temp(f(2*g-1)); % First t followed by >=N consecutive numbers


% for qq = 1:numberCenters
%     figure
%     index = clusters.Group{qq};
%     th = phiRadian(index,1);
%     r = Rt(index,1);
%     polarscatter(th,r,'k')
%     hold on
%     polarscatter(clusters.Centers(qq,2)*pi/180,clusters.Centers(qq,1),'r','filled')
%     
% end






