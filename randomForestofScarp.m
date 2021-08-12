
clear




load('trialProfiles.mat')



%%

load('scarpPoints.mat')
load('hypoPoints.mat')
remove = find(isnan(scarpTopDist));%[1:44,46:63,65:121,123:161,163:209,211:309,311:359,361:370];

profBotDist(remove) = [];
profBotElev(remove) = [];
scarpBotDist(remove) = [];
scarpTopDist(remove) = [];
scarpTopElev(remove) = [];
scarpBotElev(remove) = [];
hypos(remove,:) = [];

nonScarpInds = find(isnan(profBotDist));
scarpInds = find(isfinite(profBotDist));


design1(:,1) = hypos(:,9);
design1(:,2) = hypos(:,10);  
design1(:,3) = hypos(:,13);  
design1(:,4) = hypos(:,11);  
design1(:,5) = hypos(:,12);  
design1(:,6) = hypos(:,1); 
design1(:,7) = hypos(:,2);  
design1(:,8) = hypos(:,3);  
design1(:,9) = hypos(:,4);  
design1(:,10) = hypos(:,5);  
design1(:,11) = hypos(:,6); 
design1(:,12) = hypos(:,7); 
design1(:,13) = hypos(:,8); 
design1(:,14) = hypos(:,14);  
design1(:,15) = hypos(:,15); 
design1(:,16) = hypos(:,16); 
design1(:,17) = hypos(:,17); 



scarpDesign = design1(scarpInds,:);
nonScarpDesign = design1(nonScarpInds,:);

nonScarp_Y = zeros(length(nonScarpDesign),1);
scarp_Y = ones(length(scarpInds),1);

Ycal = [nonScarp_Y;scarp_Y];
cal = [nonScarpDesign;scarpDesign];


scarpmdl = fitcsvm(cal,Ycal,'standardize',true,'KernelFunction','gaussian','KernelScale','auto');



nonScarplength = [160, 168, 176, 192, 200, 208, 216, 224, 240, 260, 272];%,120,120,120,120,120,120,120,120];
scarplength = [200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200];%,270,300,330,360,420,480,540,600];
%overlength = [132,132,132,132,132,132,130,132,132,132,132,132,132,132,132,132,132];
%runuplength = [33, 66, 88, 110,132,156,176,198,231,264,297,330,363,396,462,528,594];
%overlength = [100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100];
%runuplength = [25, 50, 62, 75, 88,100,112,125,150,175,200,250,300,350,400,450,500,550]; 
%overlength = [150,150,150,150,150,150,150,150,150,150,150,150,150,150,150];%,150,150,150]
%runuplength = [50,75, 95,125,150,175,200,225,250,300,350,400,450,525,600];%

for ggg = 1:30



for ppp = 1:length(nonScarplength)

    sizenonScarp = length(nonScarpDesign);
    [nonScarp_reorder] = randperm(sizenonScarp);
    temp_nonScarp_design = nonScarpDesign(nonScarp_reorder,:);
    temp_nonScarp_Y = nonScarp_Y(nonScarp_reorder);
    
    sizeScarp = length(scarpDesign);
    [scarp_reorder] = randperm(sizeScarp);
    temp_scarp_design = scarpDesign(scarp_reorder,:);
    temp_scarp_Y = scarp_Y(scarp_reorder);
    

    cal = [temp_nonScarp_design(1:nonScarplength(ppp),:);temp_scarp_design(1:scarplength(ppp),:)];
    Ycal = [temp_nonScarp_Y(1:nonScarplength(ppp));temp_scarp_Y(1:scarplength(ppp))];
    
    val = [temp_nonScarp_design(nonScarplength(ppp)+1:end,:);temp_scarp_design(scarplength(ppp)+1:end,:)];
    Yval = [temp_nonScarp_Y(nonScarplength(ppp)+1:end);temp_scarp_Y(scarplength(ppp)+1:end)];

mdl = generic_random_forests(cal,Ycal,100,'classification');

%mdl = fitcsvm(cal,Ycal,'standardize',true,'KernelFunction','gaussian','KernelScale','auto');

%mdl = fitglm(cal(:,1:21),logical(Ycal));
%close all
Ypred = predict(mdl,val);
Ypred = cellfun(@str2double,Ypred);
%Ypred = [Ypred{:}];
cm = zeros(2,2);
for qqq = 1:length(Yval)
    
%    p = str2num(Ypred{qqq});
     p = Ypred(qqq);
    if p == Yval(qqq)
        if p == 1
            cm(1,1) = cm(1,1)+1;
        else
            cm(2,2) = cm(2,2)+1;
        end
    else
        if p == 1
            cm(1,2) = cm(1,2) +1;
        else
            cm(2,1) = cm(2,1) +1;
        end
    end

end

cm(1,3) = cm(1,1)/(cm(1,1)+cm(1,2));
cm(3,1) = cm(1,1)/(cm(1,1)+cm(2,1));
cm(3,2) = cm(2,2)/(cm(2,2)+cm(1,2));
cm(2,3) = cm(2,2)/(cm(2,2)+cm(2,1));

accuracy(ggg,ppp) = (cm(1,1)+cm(2,2))/length(Yval);
sensitivity(ggg,ppp) = cm(3,1);
precision(ggg,ppp) = cm(1,3);
specificity(ggg,ppp) = cm(3,2);
npv(ggg,ppp) = cm(2,3);

fnr(ggg,ppp) = cm(1,2)/(cm(1,1)+cm(1,2));
fpr(ggg,ppp) = cm(2,1)/(cm(2,1)+cm(2,2));

misclass(ggg,ppp) = (cm(1,2)+cm(2,1))/length(Yval);


Fscore(ggg,ppp) = 2.*(cm(1,3)*cm(3,1))/(cm(1,3)+cm(3,1));

end

percent = scarplength./nonScarplength;



end

figure
%col = rgb('orange');
errorbar(percent,mean(sensitivity),mean(sensitivity)-min(sensitivity),max(sensitivity)-mean(sensitivity),'-o','color','r','markerfacecolor','r')
hold on
%col = rgb('darkgoldenrod');
errorbar(percent,mean(specificity),mean(specificity)-min(specificity),max(specificity)-mean(specificity),'-o','color','b','markerfacecolor','b')
% col = rgb('navy');
% errorbar(percent,mean(fpr),mean(fpr)-min(fpr),max(fpr)-mean(fpr),'-o','markerfacecolor',col)
% %col = rgb('purple');
% %errorbar(percent,mean(fnr),mean(fnr)-min(fnr),max(fnr)-mean(fnr),'-o','markerfacecolor',col)
% col = rgb('purple');
% errorbar(percent,mean(misclass),mean(misclass)-min(misclass),max(misclass)-mean(misclass),'-o','markerfacecolor',col)
% 


%errorbar(percent,mean(precision),mean(precision)-min(precision),max(precision)-mean(precision),'-o','markerfacecolor','b')
%hold on
%errorbar(percent,mean(npv),mean(npv)-min(npv),max(npv)-mean(npv),'-o','markerfacecolor',col)

errorbar(percent,mean(accuracy),mean(accuracy)-min(accuracy),max(accuracy)-mean(accuracy),'-ok','markerfacecolor','k')
%legend('Sensitivity','Specificity','False Positive Rate','Misclassification Rate','Accuracy')
legend('Sensitivity','Specificity','Accuracy')

%legend('Precision','Sensitivity','Specificity','Negative Predictive Value')
ylabel('Prediction skills')
xlabel('Scarp Cases / Non-Scarp Cases')
title('Effect of # of cases chosen for calibration (Random Forest)')




%%

clear
load('trialsAB.mat')
load('scarpPoints.mat')
load('hypoPoints.mat')

remove = find(isnan(a));
a(remove) = [];
b(remove) = [];
xLoc(remove) = [];
zLoc(remove) = [];
hypos(remove,:) = [];






%%



%{


%%

casenumbers = [1:204];
%load('/home/server/pi/homes/danderso/COSMOS_lite/LHS_1.mat')
load('profileTrainingLibrary2Updated3.mat')

designParams = predictors;%(casenumbers,:);

load('profileBreaks2.mat')
breakBotDist = [middleBreakBotDist{:}]; %cellfun(@str2double,middleBreakBotDist);
breakTopDist = [middleBreakTopDist{:}]; %cellfun(@str2double,middleBreakTopDist);
breakBotElev = [middleBreakBotElev{:}]; %cellfun(@str2double,middleBreakBotElev);
breakTopElev = [middleBreakTopElev{:}]; %cellfun(@str2double,middleBreakTopElev);
breakSlope = [middleBreakSlope{:}]; %cellfun(@str2double,middleBreakSlope);

keepers = [1:44,46:63,65:121,123:161,163:209,211:309,311:359,361:370];

midBreakDist = breakTopDist(keepers);
midBreakHeight = breakTopElev(keepers)-breakBotElev(keepers);
%midBreakHeight = middleBreakHeight(keepers);
midBreakDist = midBreakDist - 80;

design1(:,1) = designParams(:,9);
design1(:,2) = designParams(:,10);  
design1(:,3) = designParams(:,13);  
design1(:,4) = designParams(:,11);  
design1(:,5) = designParams(:,12);  
design1(:,6) = designParams(:,1); 
design1(:,7) = designParams(:,2);  
design1(:,8) = designParams(:,3);  
design1(:,9) = designParams(:,4);  
design1(:,10) = designParams(:,5);  
design1(:,11) = designParams(:,6); 
design1(:,12) = designParams(:,7); 
design1(:,13) = designParams(:,8); 
design1(:,14) = designParams(:,14);  
design1(:,15) = designParams(:,15); 
design1(:,16) = designParams(:,16); 
design1(:,17) = designParams(:,17); 


tbl = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
    design1(:,17));
tbl.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2'};
%gprMdl = fitrgp(tbl,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)



nonScarpInds = find(isnan(midBreakDist));
scarpInds = find(~isnan(midBreakDist));

% midBreakDist(distanceNans) = [];
% midBreakHeight(distanceNans) = [];
% trials = [1:245];
% 
% trials(distanceNans) = [];
% designParams(distanceNans,:) = [];


scarpDesign = design1(scarpInds,:);
nonScarpDesign = design1(nonScarpInds,:);

nonScarp_Y = zeros(length(nonScarpDesign),1);
scarp_Y = ones(length(scarpInds),1);

Ycal = [nonScarp_Y;scarp_Y];
cal = [nonScarpDesign;scarpDesign];


scarpmdl = fitcsvm(cal,Ycal,'standardize',true,'KernelFunction','gaussian','KernelScale','auto');




%%


nonScarplength = [35,35,35,35,35,35,35,35,35,35,35];%,120,120,120,120,120,120,120,120];
scarplength = [13, 21, 28, 35, 42, 50, 58, 66, 88,110,121];%,270,300,330,360,420,480,540,600];
%overlength = [132,132,132,132,132,132,130,132,132,132,132,132,132,132,132,132,132];
%runuplength = [33, 66, 88, 110,132,156,176,198,231,264,297,330,363,396,462,528,594];
%overlength = [100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100];
%runuplength = [25, 50, 62, 75, 88,100,112,125,150,175,200,250,300,350,400,450,500,550]; 
%overlength = [150,150,150,150,150,150,150,150,150,150,150,150,150,150,150];%,150,150,150]
%runuplength = [50,75, 95,125,150,175,200,225,250,300,350,400,450,525,600];%

for ggg = 1:30



for ppp = 1:length(nonScarplength)

    sizenonScarp = length(nonScarpDesign);
    [nonScarp_reorder] = randperm(sizenonScarp);
    temp_nonScarp_design = nonScarpDesign(nonScarp_reorder,:);
    temp_nonScarp_Y = nonScarp_Y(nonScarp_reorder);
    
    sizeScarp = length(scarpDesign);
    [scarp_reorder] = randperm(sizeScarp);
    temp_scarp_design = scarpDesign(scarp_reorder,:);
    temp_scarp_Y = scarp_Y(scarp_reorder);
    

    cal = [temp_nonScarp_design(1:nonScarplength(ppp),:);temp_scarp_design(1:scarplength(ppp),:)];
    Ycal = [temp_nonScarp_Y(1:nonScarplength(ppp));temp_scarp_Y(1:scarplength(ppp))];
    
    val = [temp_nonScarp_design(nonScarplength(ppp)+1:end,:);temp_scarp_design(scarplength(ppp)+1:end,:)];
    Yval = [temp_nonScarp_Y(nonScarplength(ppp)+1:end);temp_scarp_Y(scarplength(ppp)+1:end)];

mdl = generic_random_forests(cal,Ycal,100,'classification');

%mdl = fitcsvm(cal,Ycal,'standardize',true,'KernelFunction','gaussian','KernelScale','auto');

%mdl = fitglm(cal(:,1:21),logical(Ycal));
%close all
Ypred = predict(mdl,val);
Ypred = cellfun(@str2double,Ypred);
%Ypred = [Ypred{:}];
cm = zeros(2,2);
for qqq = 1:length(Yval)
    
%    p = str2num(Ypred{qqq});
     p = Ypred(qqq);
    if p == Yval(qqq)
        if p == 1
            cm(1,1) = cm(1,1)+1;
        else
            cm(2,2) = cm(2,2)+1;
        end
    else
        if p == 1
            cm(1,2) = cm(1,2) +1;
        else
            cm(2,1) = cm(2,1) +1;
        end
    end

end

cm(1,3) = cm(1,1)/(cm(1,1)+cm(1,2));
cm(3,1) = cm(1,1)/(cm(1,1)+cm(2,1));
cm(3,2) = cm(2,2)/(cm(2,2)+cm(1,2));
cm(2,3) = cm(2,2)/(cm(2,2)+cm(2,1));

accuracy(ggg,ppp) = (cm(1,1)+cm(2,2))/length(Yval);
sensitivity(ggg,ppp) = cm(3,1);
precision(ggg,ppp) = cm(1,3);
specificity(ggg,ppp) = cm(3,2);
npv(ggg,ppp) = cm(2,3);

fnr(ggg,ppp) = cm(1,2)/(cm(1,1)+cm(1,2));
fpr(ggg,ppp) = cm(2,1)/(cm(2,1)+cm(2,2));

misclass(ggg,ppp) = (cm(1,2)+cm(2,1))/length(Yval);


Fscore(ggg,ppp) = 2.*(cm(1,3)*cm(3,1))/(cm(1,3)+cm(3,1));

end

percent = scarplength./nonScarplength;



end

figure
%col = rgb('orange');
errorbar(percent,mean(sensitivity),mean(sensitivity)-min(sensitivity),max(sensitivity)-mean(sensitivity),'-o','color','r','markerfacecolor','r')
hold on
%col = rgb('darkgoldenrod');
errorbar(percent,mean(specificity),mean(specificity)-min(specificity),max(specificity)-mean(specificity),'-o','color','b','markerfacecolor','b')
% col = rgb('navy');
% errorbar(percent,mean(fpr),mean(fpr)-min(fpr),max(fpr)-mean(fpr),'-o','markerfacecolor',col)
% %col = rgb('purple');
% %errorbar(percent,mean(fnr),mean(fnr)-min(fnr),max(fnr)-mean(fnr),'-o','markerfacecolor',col)
% col = rgb('purple');
% errorbar(percent,mean(misclass),mean(misclass)-min(misclass),max(misclass)-mean(misclass),'-o','markerfacecolor',col)
% 


%errorbar(percent,mean(precision),mean(precision)-min(precision),max(precision)-mean(precision),'-o','markerfacecolor','b')
%hold on
%errorbar(percent,mean(npv),mean(npv)-min(npv),max(npv)-mean(npv),'-o','markerfacecolor',col)

errorbar(percent,mean(accuracy),mean(accuracy)-min(accuracy),max(accuracy)-mean(accuracy),'-ok','markerfacecolor','k')
%legend('Sensitivity','Specificity','False Positive Rate','Misclassification Rate','Accuracy')
legend('Sensitivity','Specificity','Accuracy')

%legend('Precision','Sensitivity','Specificity','Negative Predictive Value')
ylabel('Prediction skills')
xlabel('Scarp Cases / Non-Scarp Cases')
title('Effect of # of cases chosen for calibration (Random Forest)')




%%
% 
% maxHs=max(designParams(:,9));  minHs=min(designParams(:,9));
% maxTp=max(designParams(:,10));  minTp=min(designParams(:,10));
% maxDir=max(designParams(:,13));  minDir=min(designParams(:,13));
% maxNTR=max(designParams(:,11));  minNTR=min(designParams(:,11));
% maxDur=max(designParams(:,12));  minDur=min(designParams(:,12));
% maxEOF1=max(designParams(:,1));  minEOF1=min(designParams(:,1));
% maxEOF2=max(designParams(:,2));  minEOF2=min(designParams(:,2));
% maxEOF3=max(designParams(:,3));  minEOF3=min(designParams(:,3));
% maxEOF4=max(designParams(:,4));  minEOF4=min(designParams(:,4));
% maxEOF5=max(designParams(:,5));  minEOF5=min(designParams(:,5));
% maxEOF6=max(designParams(:,6));  minEOF6=min(designParams(:,6));
% maxEOF7=max(designParams(:,7));  minEOF7=min(designParams(:,7));
% maxEOFd=max(designParams(:,8));  minEOFd=min(designParams(:,8));
% maxM2=max(designParams(:,14));  minM2=min(designParams(:,14));
% maxN2=max(designParams(:,15));  minN2=min(designParams(:,15));
% maxK1=max(designParams(:,16));  minK1=min(designParams(:,16));
% maxS2=max(designParams(:,17));  minS2=min(designParams(:,17));
% 
% maxDist = max(designParams(:,18)); minDist = min(designParams(:,18));
% maxHeight = max(designParams(:,19)); minHeight = min(designParams(:,19));
% 
% subset_n(:,1)=(designParams(:,9)-minHs)./(maxHs-minHs);
% subset_n(:,2)=(designParams(:,10)-minTp)./(maxTp-minTp);
% subset_n(:,3)=designParams(:,13)*pi/180;
% subset_n(:,4)=(designParams(:,11)-minNTR)./(maxNTR-minNTR);
% subset_n(:,5)=(designParams(:,12)-minDur)./(maxDur-minDur);
% subset_n(:,6)=(designParams(:,1)-minEOF1)./(maxEOF1-minEOF1);
% subset_n(:,7)=(designParams(:,2)-minEOF2)./(maxEOF1-minEOF2);
% subset_n(:,8)=(designParams(:,3)-minEOF3)./(maxEOF3-minEOF3);
% subset_n(:,9)=(designParams(:,4)-minEOF4)./(maxEOF1-minEOF4);
% subset_n(:,10)=(designParams(:,5)-minEOF5)./(maxEOF5-minEOF5);
% subset_n(:,11)=(designParams(:,6)-minEOF6)./(maxEOF6-minEOF6);
% subset_n(:,12)=(designParams(:,7)-minEOF7)./(maxEOF5-minEOF7);
% subset_n(:,13)=(designParams(:,8)-minEOFd)./(maxEOFd-minEOFd);
% subset_n(:,14)=(designParams(:,14)-minN2)./(maxN2-minN2);
% subset_n(:,15)=(designParams(:,15)-minM2)./(maxM2-minM2);
% subset_n(:,16)=(designParams(:,16)-minK1)./(maxK1-minK1);
% subset_n(:,17)=(designParams(:,17)-minS2)./(maxS2-minS2);
% subset_n(:,18)=(designParams(:,18)-minDist)./(maxDist-minDist);
% subset_n(:,19)=(designParams(:,19)-minHeight)./(maxHeight-maxHeight);
% 
% 


%}

