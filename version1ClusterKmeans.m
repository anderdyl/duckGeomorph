clear
close all
cd /home/dylananderson/projects/duckGeomorph/

addpath(genpath('/media/dylananderson/Elements/SERDP'))

load('ceofsSouthernTransect.mat')

numberCenters = 15;
data(1:length(Rt),1) = Rt(:,1)*(percentV(1)/percentV(1));
data(:,2) = phiDegrees(:,1);
data(:,3) = Rt(:,2)*(percentV(2)/percentV(1));
data(:,4) = phiDegrees(:,2);
data(:,5) = Rt(:,3)*(percentV(3)/percentV(1));
data(:,6) = phiDegrees(:,3);
%data(:,7) = Rt(:,4)*(percentV(4)/percentV(1));
%data(:,8) = phiDegrees(:,4);
% data(:,9) = Rt(:,5);
% data(:,10) = phiDegrees(:,5)/(percentV(1));
% data(:,11) = Rt(:,6);
% data(:,12) = phiDegrees(:,6)/(percentV(1));
% data(:,13) = Rt(:,7);
% data(:,14) = phiDegrees(:,7)/(percentV(1));

type = 'kmeans';
escalar = [1,3,5];%,7,9,11,13];
direccional = [2,4,6];%,8,10,12,14];
initmethod = 'mda';
clusters = makeClustering_initMDA(data, type, numberCenters, escalar, direccional, initmethod);

% What do the clusters look like with respect to the first offshore
% migrating CEOF?
figure
polarscatter(data(:,2)*pi/180,data(:,1))
hold on
polarscatter(clusters.Centers(:,2)*pi/180,clusters.Centers(:,1),'r','filled')

[val, order] = sort(clusters.Centers(:,2),'descend');
colmap = jet(numberCenters);
figure
for pp = 1:numberCenters
    index = order(pp);
    plot(xinterp,mean(alllines(clusters.Group{index},:)), '-','color',colmap(pp,:))
    hold on
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

% What do the clusters look like comparing the leading three magnitues?
figure
scatter(clusters.Centers(:,1)/(percentV(1)/percentV(1)),clusters.Centers(:,3)/(percentV(2)/percentV(1)),'o')


temp = clusters.Group{7};

p=find(diff(temp)==1);
q=[p;p+1];
temp(q)  % this gives all the pairs of consecutive numbers

s = find(diff(temp)>1);
figure
for pp = 1:length(s);
    subplot(length(s),1,pp)
    if pp == 1
        profiles = [1:s(pp)];
    else
        profiles = [s(pp-1)+1:s(pp)];
    end
    for xx = 1:length(profiles)
        plot(xinterp,alllines(temp(profiles(xx)),:))
        hold on
    end
    
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






