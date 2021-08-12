
clear




casenumbers = [1:309];
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
midBreakDistBot = breakBotDist(keepers);
midBreakHeight = breakTopElev(keepers)-breakBotElev(keepers);
midBreakElev = breakTopElev(keepers);
midBreakElevBot = breakBotElev(keepers);
%midBreakHeight = middleBreakHeight(keepers);
midBreakDist = midBreakDist - 80;
midBreakDistBot = midBreakDistBot - 80;
distanceNans = find(isnan(midBreakDist));
% midBreakDist(distanceNans) = zeroCrossing;
% zerosNans = find(isnan(midBreakHeight));
% midBreakHeight(zerosNans) = 0;
% 

midBreakDist(distanceNans) = [];
midBreakDistBot(distanceNans) = [];
midBreakElev(distanceNans) = [];
midBreakElevBot(distanceNans) = [];
midBreakHeight(distanceNans) = [];
trials = [1:363];
trials(distanceNans) = [];
designParams(distanceNans,:) = [];


diffBreak = midBreakDistBot-midBreakDist;
tooBig = find(diffBreak > 10);

midBreakDist(tooBig) = [];
midBreakDistBot(tooBig) = [];
midBreakElev(tooBig) = [];
midBreakElevBot(tooBig) = [];
midBreakHeight(tooBig) = [];
trials(tooBig) = [];
designParams(tooBig,:) = [];

designParams(:,18) = midBreakDist;
designParams(:,19) = midBreakHeight;


lengthOfScarps = 363-length(distanceNans)-length(tooBig);


% WHATS IN WHICH COLUMN
% 1 - EOF1
% 2 - EOF2
% 3 - EOF3
% 4 - EOF4
% 5 - EOF5
% 6 - EOF6
% 7 - EOF7
% 8 - EOF differences N-S
% 9 - Hs
% 10 - Tp
% 11 - NTR
% 12 - duration
% 13 - direction
% 14 - M2
% 15 - N2
% 16 - K1
% 17 - S2
%Normalization of the selected cases with MaxDiss
maxHs=max(designParams(:,9));  minHs=min(designParams(:,9));
maxTp=max(designParams(:,10));  minTp=min(designParams(:,10));
maxDir=max(designParams(:,13));  minDir=min(designParams(:,13));
maxNTR=max(designParams(:,11));  minNTR=min(designParams(:,11));
maxDur=max(designParams(:,12));  minDur=min(designParams(:,12));
maxEOF1=max(designParams(:,1));  minEOF1=min(designParams(:,1));
maxEOF2=max(designParams(:,2));  minEOF2=min(designParams(:,2));
maxEOF3=max(designParams(:,3));  minEOF3=min(designParams(:,3));
maxEOF4=max(designParams(:,4));  minEOF4=min(designParams(:,4));
maxEOF5=max(designParams(:,5));  minEOF5=min(designParams(:,5));
maxEOF6=max(designParams(:,6));  minEOF6=min(designParams(:,6));
maxEOF7=max(designParams(:,7));  minEOF7=min(designParams(:,7));
maxEOFd=max(designParams(:,8));  minEOFd=min(designParams(:,8));
maxM2=max(designParams(:,14));  minM2=min(designParams(:,14));
maxN2=max(designParams(:,15));  minN2=min(designParams(:,15));
maxK1=max(designParams(:,16));  minK1=min(designParams(:,16));
maxS2=max(designParams(:,17));  minS2=min(designParams(:,17));

maxDist = max(designParams(:,18)); minDist = min(designParams(:,18));
maxHeight = max(designParams(:,19)); minHeight = min(designParams(:,19));

subset_n(:,1)=(designParams(:,9)-minHs)./(maxHs-minHs);
subset_n(:,2)=(designParams(:,10)-minTp)./(maxTp-minTp);
subset_n(:,3)=designParams(:,13)*pi/180;
subset_n(:,4)=(designParams(:,11)-minNTR)./(maxNTR-minNTR);
subset_n(:,5)=(designParams(:,12)-minDur)./(maxDur-minDur);
subset_n(:,6)=(designParams(:,1)-minEOF1)./(maxEOF1-minEOF1);
subset_n(:,7)=(designParams(:,2)-minEOF2)./(maxEOF1-minEOF2);
subset_n(:,8)=(designParams(:,3)-minEOF3)./(maxEOF3-minEOF3);
subset_n(:,9)=(designParams(:,4)-minEOF4)./(maxEOF1-minEOF4);
subset_n(:,10)=(designParams(:,5)-minEOF5)./(maxEOF5-minEOF5);
subset_n(:,11)=(designParams(:,6)-minEOF6)./(maxEOF6-minEOF6);
subset_n(:,12)=(designParams(:,7)-minEOF7)./(maxEOF5-minEOF7);
subset_n(:,13)=(designParams(:,8)-minEOFd)./(maxEOFd-minEOFd);
subset_n(:,14)=(designParams(:,14)-minN2)./(maxN2-minN2);
subset_n(:,15)=(designParams(:,15)-minM2)./(maxM2-minM2);
subset_n(:,16)=(designParams(:,16)-minK1)./(maxK1-minK1);
subset_n(:,17)=(designParams(:,17)-minS2)./(maxS2-minS2);
subset_n(:,18)=(designParams(:,18)-minDist)./(maxDist-minDist);
subset_n(:,19)=(designParams(:,19)-minHeight)./(maxHeight-maxHeight);






all = subset_n(:,1:17);

calibration = 1:lengthOfScarps-50;
validation = lengthOfScarps-50:lengthOfScarps;

target90 = predictands(casenumbers,1);%-predictors(casenumbers,1);
design1 = all(calibration,:);
design1(calibration,18) = target90(calibration);
tbl = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
    design1(:,17),design1(:,18));
tbl.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
gprMdl = fitrgp(tbl,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)

target902 = predictands(casenumbers,2);%-predictors(casenumbers,2);
design1 = all(calibration,:);
design1(calibration,18) = target902(calibration);
tbl2 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
    design1(:,17),design1(:,18));
tbl2.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
gprMdl2 = fitrgp(tbl2,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)

target903 = predictands(casenumbers,3);%-predictors(casenumbers,3);
design1 = all(calibration,:);
design1(calibration,18) = target903(calibration);
tbl3 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
    design1(:,17),design1(:,18));
tbl3.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
gprMdl3 = fitrgp(tbl3,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)

target904 = predictands(casenumbers,4);%-predictors(casenumbers,4);
design1 = all(calibration,:);
design1(calibration,18) = target904(calibration);
tbl4 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
    design1(:,17),design1(:,18));
tbl4.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
gprMdl4 = fitrgp(tbl4,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)

target905 = predictands(casenumbers,5);%-predictors(casenumbers,5);
design1 = all(calibration,:);
design1(calibration,18) = target905(calibration);
tbl5 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
    design1(:,17),design1(:,18));
tbl5.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
gprMdl5 = fitrgp(tbl5,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)




target906 = predictands(casenumbers,6);%-predictors(casenumbers,6);
design1 = all(calibration,:);
design1(calibration,18) = target906(calibration);
tbl6 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
    design1(:,17),design1(:,18));
tbl6.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
gprMdl6 = fitrgp(tbl6,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)

target907 = predictands(casenumbers,7);%-predictors(casenumbers,7);
design1 = all(calibration,:);
design1(calibration,18) = target907(calibration);
tbl7 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
    design1(:,17),design1(:,18));
tbl7.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
gprMdl7 = fitrgp(tbl7,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)

target90d = predictands(casenumbers,8);
design1 = all(calibration,:);
design1(calibration,18) = target90d(calibration);
tbld = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
    design1(:,17),design1(:,18));
tbld.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
gprMdld = fitrgp(tbld,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)



target90dist = midBreakDist(casenumbers);
design1 = all(calibration,:);
design1(calibration,18) = target90dist(calibration);
tbldist = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
    design1(:,17),design1(:,18));
tbldist.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
gprMdldist = fitrgp(tbldist,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)



target90height = midBreakHeight(casenumbers);
design1 = all(calibration,:);
design1(calibration,18) = target90height(calibration);
tblheight = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
    design1(:,17),design1(:,18));
tblheight.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
gprMdlheight = fitrgp(tblheight,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)



%%

%clear
%load('latest_gpr_model_wl.mat','tbl')
%gprMdl = fitrgp(tbl,'Hsig','KernelFunction','ardsquaredexponential','OptimizeHyperparameters','auto')
%load(fullfile(matlabroot,'examples','stats','gprdata2.mat'))
%Xtrain = [tbl.NTR,tbl.SLP,tbl.U,tbl.V,tbl.HsN,tbl.TpN,tbl.DirN,tbl.HsS,tbl.TpS,tbl.DirS,tbl.Hsea,tbl.Tsea,tbl.Dsea,tbl.M2,tbl.S2,tbl.N2,tbl.K2,tbl.K1,tbl.O1,tbl.P1,tbl.Q1];
%ytrain = tbl.Hsig;
%sigma0 = std(ytrain);
%sigmaF0 = sigma0;
%d = size(Xtrain,2);
%sigmaM0 = 10*ones(d,1);

%gprMdl = fitrgp(Xtrain,ytrain,'Basis','constant','FitMethod','exact','PredictMethod','exact','KernelFunction',...
%'ardsquaredexponential','KernelParameters',[sigmaM0:sigmaF0],'Sigma',sigma0,'Standardize',1);


%gprMdl = fitrgp(Xtrain,y,'KernelFunction','squaredexponential','OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%    struct('AcquisitionFunctionName','expected-improvement-plus'));

%%



ypred = predict(gprMdl,all(validation,:));
ypred2 = predict(gprMdl2,all(validation,:));
ypred3 = predict(gprMdl3,all(validation,:));
ypred4 = predict(gprMdl4,all(validation,:));
ypred5 = predict(gprMdl5,all(validation,:));
ypred6 = predict(gprMdl6,all(validation,:));
ypred7 = predict(gprMdl7,all(validation,:));
ypredd = predict(gprMdld,all(validation,:));
ypreddist = predict(gprMdldist,all(validation,:));
ypredheight = predict(gprMdlheight,all(validation,:));

% %%
figure();
plot(ypred,target90(validation),'o')
hold on
plot(ypred2,target902(validation),'o')
plot(ypred3,target903(validation),'o')
plot(ypred4,target904(validation),'o')
plot(ypred5,target905(validation),'o')
plot(ypred6,target906(validation),'o')
plot(ypred7,target907(validation),'o')


%plot(ypred,'b')
xlabel('true')
ylabel('predicted')
%legend({'data','predictions'},'Location','Best')
%axis([0 4300 0 30])
%hold off

final(1:length(ypred),1) = ypred;% + predictors(validation,1);
final(1:length(ypred),2) = ypred2;% + predictors(validation,2);
final(1:length(ypred),3) = ypred3;% + predictors(validation,3);
final(1:length(ypred),4) = ypred4;% + predictors(validation,4);
final(1:length(ypred),5) = ypred5;% + predictors(validation,5);
final(1:length(ypred),6) = ypred6;% + predictors(validation,6);
final(1:length(ypred),7) = ypred7;% + predictors(validation,7);
final(1:length(ypred),8) = ypredd;
final(1:length(ypred),9) = ypreddist;
final(1:length(ypred),10) = ypredheight;

valTrials = trials(lengthOfScarps-50:end);

save('predictions2.mat','valTrials','ypred','ypred2','ypred3','ypred4','ypred5','ypred6','ypred7','ypredd','final','ypreddist','ypredheight')
figure
plot(ypredd,target90d(validation),'o')
hold on
plot(ypreddist,target90dist(validation),'o')
plot(ypredheight,target90height(validation),'o')

%%

%sizedata = length(design);
%[reorder] = randperm(sizedata);
%design = design(reorder,:);
designParams = designParams(1:300,:);
midBreakDist = midBreakDist(1:300);
midBreakDistBot = midBreakDistBot(1:300);
midBreakElev = midBreakElev(1:300);
midBreakElevBot = midBreakElevBot(1:300);

midBreakHeight = midBreakHeight(1:300);
sizedata = length(designParams);
[reorder] = randperm(sizedata);

sets = [60]; %,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96];
totals = [300]; %,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480];

figure
set(gcf,'color','w');

for hh = 1:length(sets)
    size = sets(hh);
    total = totals(hh);
    design = all(reorder,:);
    targets = predictands(reorder,:);
    starts = predictors(reorder,:);
    distTargets = midBreakDist(reorder);
    distTargetsBot = midBreakDistBot(reorder);
    distHeights = midBreakHeight(reorder);
    
    distTargetsElev = midBreakElev(reorder);
    distTargetsElevBot = midBreakElevBot(reorder);
    newElev = distTargetsElev(1:total);
    newElevBot = distTargetsElevBot(1:total);
    
    newall = design(1:total,:);
    newtargets = targets(1:total,:);
    newdist = distTargets(1:total);
    newdistBot = distTargetsBot(1:total);
    
    newheight = distHeights(1:total);
    
    c = 1;
    for qq = 1:5

        design1 = newall;
        design1(c:c+size-1,:) = [];
        target1 = newtargets(:,1);%-starts(:,1);
        target1(c:c+size-1,:) = [];
        target2 = newtargets(:,2);%-starts(:,2);
        target2(c:c+size-1,:) = [];
        target3 = newtargets(:,3);%-starts(:,3);
        target3(c:c+size-1,:) = [];
        target4 = newtargets(:,4);%-starts(:,4);
        target4(c:c+size-1,:) = [];
        target5 = newtargets(:,5);%-starts(:,5);
        target5(c:c+size-1,:) = [];
        target6 = newtargets(:,6);%-starts(:,6);
        target6(c:c+size-1,:) = [];
        target7 = newtargets(:,7);%-starts(:,7);
        target7(c:c+size-1,:) = [];
        targetd = newtargets(:,8);
        targetd(c:c+size-1,:) = [];
        targetdist = newdist;
        targetdist(c:c+size-1) = [];
        targetdistBot = newdistBot;
        targetdistBot(c:c+size-1) = [];
        targetheight = newheight;
        targetheight(c:c+size-1) = [];
        
        targetElev = newElev;
        targetElev(c:c+size-1) = [];
        targetElevBot = newElevBot;
        targetElevBot(c:c+size-1) = [];        
        
        
        disp(['removing ',num2str(c),' to ',num2str(c+size-1),' from starting dataset'])
        data1 = newall(c:c+size-1,:);
        val1 = newtargets(c:c+size-1,1);%-starts(c:c+size-1,1);
        val2 = newtargets(c:c+size-1,2);%-starts(c:c+size-1,2);
        val3 = newtargets(c:c+size-1,3);%-starts(c:c+size-1,3);
        val4 = newtargets(c:c+size-1,4);%-starts(c:c+size-1,4);
        val5 = newtargets(c:c+size-1,5);%-starts(c:c+size-1,5);
        val6 = newtargets(c:c+size-1,6);%-starts(c:c+size-1,6);
        val7 = newtargets(c:c+size-1,7);%-starts(c:c+size-1,7);
        vald = newtargets(c:c+size-1,8);
        valdist = newdist(c:c+size-1);
        valdistBot = newdistBot(c:c+size-1);
        valElev = newElev(c:c+size-1);
        valElevBot = newElevBot(c:c+size-1);
        valheight = newheight(c:c+size-1);
        disp(['holding ',num2str(c),' to ',num2str(c+size-1),' aside for validation'])
        

%         disp(['working on EOF 1'])
%         
%         design1(:,18) = target1;
%         tbl1 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl1.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl1 = fitrgp(tbl1,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results1(1:size,qq) = predict(gprMdl1,data1);
%         Validation1(1:size,qq) = val1;
%         eval(['gprMdl',num2str(qq),' = gprMdl1;'])
%         subplot(2,10,1)
%         eval(['p',num2str(qq),' = scatter(val1,Results1(:,qq),''filled'')'])
%         hold on
%         xlabel('Xbeach PCs (m)')
%         ylabel('GPR PCs (m)')
        
%         
%         disp(['working on EOF 2'])
%         
%         design1(:,18) = target2;
%         tbl2 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl2.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl2 = fitrgp(tbl2,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results2(1:size,qq) = predict(gprMdl2,data1);
%         Validation2(1:size,qq) = val2;
%         eval(['gprMdl',num2str(qq),' = gprMdl1;'])
%         subplot(2,10,2)
%         eval(['p',num2str(qq),' = scatter(val2,Results2(:,qq),''filled'')'])
%         hold on
%         xlabel('Xbeach PCs (m)')
%         ylabel('GPR PCs (m)')        
        
        
%         disp(['working on EOF 3'])
%         
%         design1(:,18) = target3;
%         tbl3 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl3.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl3 = fitrgp(tbl3,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results3(1:size,qq) = predict(gprMdl3,data1);
%         Validation3(1:size,qq) = val3;
%         eval(['gprMdl',num2str(qq),' = gprMdl3;'])
%         subplot(2,10,3)
%         eval(['p',num2str(qq),' = scatter(val3,Results3(:,qq),''filled'')'])
%         hold on
%         xlabel('Xbeach PCs (m)')
%         ylabel('GPR PCs (m)')
        
        
        
%         disp(['working on EOF 4'])
%         
%         design1(:,18) = target4;
%         tbl4 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));        tbl4.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl4 = fitrgp(tbl4,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results4(1:size,qq) = predict(gprMdl4,data1);
%         Validation4(1:size,qq) = val4;
%         eval(['gprMdl',num2str(qq),' = gprMdl4;'])
%         subplot(2,10,4)
%         eval(['p',num2str(qq),' = scatter(val4,Results4(:,qq),''filled'')'])
%         hold on
%         xlabel('Xbeach PCs (m)')
%         ylabel('GPR PCs (m)')        
        
        
%         disp(['working on EOF 5'])
%         
%         design1(:,18) = target5;
%         tbl5 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl5.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl5 = fitrgp(tbl5,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results5(1:size,qq) = predict(gprMdl5,data1);
%         Validation5(1:size,qq) = val5;
%         eval(['gprMdl',num2str(qq),' = gprMdl5;'])
%         subplot(2,10,5)
%         eval(['p',num2str(qq),' = scatter(val5,Results5(:,qq),''filled'')'])
%         hold on
%         xlabel('Xbeach PCs (m)')
%         ylabel('GPR PCs (m)')        
        
        
%         disp(['working on EOF 6'])
%         
%         design1(:,18) = target6;
%         tbl6 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl6.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl6 = fitrgp(tbl6,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results6(1:size,qq) = predict(gprMdl6,data1);
%         Validation6(1:size,qq) = val6;
%         eval(['gprMdl',num2str(qq),' = gprMdl6;'])
%         subplot(2,10,6)
%         eval(['p',num2str(qq),' = scatter(val6,Results6(:,qq),''filled'')'])
%         hold on
%         xlabel('Xbeach PCs (m)')
%         ylabel('GPR PCs (m)')  
        
        
        
%         disp(['working on EOF 7'])
%         
%         design1(:,18) = target7;
%         tbl7 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl7.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl7 = fitrgp(tbl7,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results7(1:size,qq) = predict(gprMdl7,data1);
%         Validation7(1:size,qq) = val7;
%         eval(['gprMdl',num2str(qq),' = gprMdl7;'])
%         subplot(2,10,7)
%         eval(['p',num2str(qq),' = scatter(val7,Results7(:,qq),''filled'')'])
%         hold on
%         xlabel('Xbeach PCs (m)')
%         ylabel('GPR PCs (m)')          
        
%         disp(['working on EOF d'])
%         
%         design1(:,18) = targetd;
%         tbld = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbld.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdld = fitrgp(tbld,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Resultsd(1:size,qq) = predict(gprMdld,data1);
%         Validationd(1:size,qq) = vald;
%         eval(['gprMdl',num2str(qq),' = gprMdld;'])
%         subplot(2,10,8)
%         eval(['p',num2str(qq),' = scatter(vald,Resultsd(:,qq),''filled'')'])
%         hold on
%         xlabel('Xbeach PCs (m)')
%         ylabel('GPR PCs (m)')  
        
        
        disp('working on scarp #1 distance')
        
        design1(:,18) = targetdist;
        tbldist = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
            design1(:,17),design1(:,18));
        tbldist.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
        gprMdldist = fitrgp(tbldist,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsdist(1:size,qq) = predict(gprMdldist,data1);
        Validationdist(1:size,qq) = valdist;
        eval(['gprMdl',num2str(qq),' = gprMdldist;'])
        subplot(2,5,1)
        eval(['p',num2str(qq),' = scatter(valdist,Resultsdist(:,qq),''filled'')'])
        hold on
        xlabel('Xbeach (m)')
        ylabel('GPR (m)')          
        
        
        disp('working on scarp #2 distance')
        
        design1(:,18) = targetdistBot;
        tbldist2 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
            design1(:,17),design1(:,18));
        tbldist2.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
        gprMdldist2 = fitrgp(tbldist2,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsdist2(1:size,qq) = predict(gprMdldist2,data1);
        Validationdist2(1:size,qq) = valdistBot;
        eval(['gprMdl',num2str(qq),' = gprMdldist2;'])
        subplot(2,5,2)
        eval(['p',num2str(qq),' = scatter(valdistBot,Resultsdist2(:,qq),''filled'')'])
        hold on
        xlabel('Xbeach (m)')
        ylabel('GPR (m)')  
        
        
        
        disp('working on scarp #1 elev')
        
        design1(:,18) = targetElev;
        tbldist = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
            design1(:,17),design1(:,18));
        tbldist.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
        gprMdlelev = fitrgp(tbldist,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        ResultsElev(1:size,qq) = predict(gprMdlelev,data1);
        ValidationElev(1:size,qq) = valElev;
        eval(['gprMdl',num2str(qq),' = gprMdldist;'])
        subplot(2,5,3)
        eval(['p',num2str(qq),' = scatter(valElev,ResultsElev(:,qq),''filled'')'])
        hold on
        xlabel('Xbeach (m)')
        ylabel('GPR (m)')          
        
        
        disp('working on scarp #2 elev')
        
        design1(:,18) = targetElevBot;
        tbldist2 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
            design1(:,17),design1(:,18));
        tbldist2.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
        gprMdlelev2 = fitrgp(tbldist2,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        ResultsElev2(1:size,qq) = predict(gprMdlelev2,data1);
        ValidationElev2(1:size,qq) = valdistBot;
        eval(['gprMdl',num2str(qq),' = gprMdlelev2;'])
        subplot(2,5,4)
        eval(['p',num2str(qq),' = scatter(valElevBot,ResultsElev2(:,qq),''filled'')'])
        hold on
        xlabel('Xbeach (m)')
        ylabel('GPR (m)')  
        
        disp('working on scarp #1 height')
        
        design1(:,18) = targetheight;
        tblheight = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
            design1(:,17),design1(:,18));
        tblheight.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
        gprMdlheight= fitrgp(tblheight,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsheight(1:size,qq) = predict(gprMdlheight,data1);
        Validationheight(1:size,qq) = valheight;
        eval(['gprMdl',num2str(qq),' = gprMdlheight;'])
        subplot(2,5,5)
        eval(['p',num2str(qq),' = scatter(valheight,Resultsheight(:,qq),''filled'')'])
        hold on
        xlabel('Xbeach (m)')
        ylabel('GPR (m)')  
        
        %clear gprMdl1 gprMdl2 gprMdl3 gprMdl4 gprMdl5
        c = c +size;
        
    end
end



% subplot(2,10,1)
% [r m b] = regression(Results1(:)',Validation1(:)');
% title(['EOF 1 R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
% %plot([0.5 4.5],[0.5 4.5],'w--')
% p6 = plot([-11 11],[-11 11],'k--');
% lg = legend([p1 p2 p3 p4 p5 p6],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','1-to-1','location','northwest');
% % ylim([-4 3])
% % xlim([-4 3])
% subplot(2,10,11)
% errors = Results1-Validation1;
% binrng = -1:.1:1;
% counts1 = histc(errors(:,1),binrng);
% counts2 = histc(errors(:,2),binrng)+counts1;
% counts3 = histc(errors(:,3),binrng)+counts2;
% counts4 = histc(errors(:,4),binrng)+counts3;
% counts5 = histc(errors(:,5),binrng)+counts4;
% bar(binrng,counts5)
% hold on
% bar(binrng,counts4)
% bar(binrng,counts3)
% bar(binrng,counts2)
% bar(binrng,counts1)
% xlabel('PC error')
% title('Errors for all folds','color','k')

% 
% 
% subplot(2,10,2)
% [r m b] = regression(Results2(:)',Validation2(:)');
% title(['EOF 2 R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
% %plot([0.5 4.5],[0.5 4.5],'w--')
% plot([-10 10],[-10 10],'k--')
% %lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
% % ylim([-12 3])
% % xlim([-12 3])
% subplot(2,10,12)
% errors = Results2-Validation2;
% binrng = -1:0.1:1;
% counts1 = histc(errors(:,1),binrng);
% counts2 = histc(errors(:,2),binrng)+counts1;
% counts3 = histc(errors(:,3),binrng)+counts2;
% counts4 = histc(errors(:,4),binrng)+counts3;
% counts5 = histc(errors(:,5),binrng)+counts4;
% bar(binrng,counts5)
% hold on
% bar(binrng,counts4)
% bar(binrng,counts3)
% bar(binrng,counts2)
% bar(binrng,counts1)
% xlabel('PC error')
% title('Errors for all folds','color','k')


% subplot(2,10,3)
% [r m b] = regression(Results3(:)',Validation3(:)');
% title(['EOF 3 R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
% %plot([0.5 4.5],[0.5 4.5],'w--')
% plot([-6 6],[-6 6],'k--')
% %lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
% %ylim([-6 6])
% %xlim([-6 6])
% subplot(2,10,13)
% errors = Results3-Validation3;
% binrng = -1:0.1:1;
% counts1 = histc(errors(:,1),binrng);
% counts2 = histc(errors(:,2),binrng)+counts1;
% counts3 = histc(errors(:,3),binrng)+counts2;
% counts4 = histc(errors(:,4),binrng)+counts3;
% counts5 = histc(errors(:,5),binrng)+counts4;
% bar(binrng,counts5)
% hold on
% bar(binrng,counts4)
% bar(binrng,counts3)
% bar(binrng,counts2)
% bar(binrng,counts1)
% xlabel('PC error')
% title('Errors for all folds','color','k')


% 
% subplot(2,10,4)
% [r m b] = regression(Results4(:)',Validation4(:)');
% title(['EOF 4 R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
% %plot([0.5 4.5],[0.5 4.5],'w--')
% plot([-5 5],[-5 5],'k--')
% %lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
% %ylim([-5 5])
% %xlim([-5 5])
% subplot(2,10,14)
% errors = Results4-Validation4;
% binrng = -1:0.1:1;
% counts1 = histc(errors(:,1),binrng);
% counts2 = histc(errors(:,2),binrng)+counts1;
% counts3 = histc(errors(:,3),binrng)+counts2;
% counts4 = histc(errors(:,4),binrng)+counts3;
% counts5 = histc(errors(:,5),binrng)+counts4;
% bar(binrng,counts5)
% hold on
% bar(binrng,counts4)
% bar(binrng,counts3)
% bar(binrng,counts2)
% bar(binrng,counts1)
% xlabel('PC error')
% title('Errors for all folds','color','k')
% 

% 
% subplot(2,10,5)
% [r m b] = regression(Results5(:)',Validation5(:)');
% title(['EOF 5 R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
% %plot([0.5 4.5],[0.5 4.5],'w--')
% plot([-6 6],[-6 6],'k--')
% %lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
% %ylim([-6 6])
% %xlim([-6 6])
% subplot(2,10,15)
% errors = Results5-Validation5;
% binrng = -1:0.1:1;
% counts1 = histc(errors(:,1),binrng);
% counts2 = histc(errors(:,2),binrng)+counts1;
% counts3 = histc(errors(:,3),binrng)+counts2;
% counts4 = histc(errors(:,4),binrng)+counts3;
% counts5 = histc(errors(:,5),binrng)+counts4;
% bar(binrng,counts5)
% hold on
% bar(binrng,counts4)
% bar(binrng,counts3)
% bar(binrng,counts2)
% bar(binrng,counts1)
% xlabel('PC error')
% title('Errors for all folds','color','k')
% 


% 
% 
% 
% subplot(2,10,6)
% [r m b] = regression(Results6(:)',Validation6(:)');
% title(['EOF 6 R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
% %plot([0.5 4.5],[0.5 4.5],'w--')
% plot([-6 6],[-6 6],'k--')
% %lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
% %ylim([-6 6])
% %xlim([-6 6])
% subplot(2,10,16)
% errors = Results6-Validation6;
% binrng = -1:0.1:1;
% counts1 = histc(errors(:,1),binrng);
% counts2 = histc(errors(:,2),binrng)+counts1;
% counts3 = histc(errors(:,3),binrng)+counts2;
% counts4 = histc(errors(:,4),binrng)+counts3;
% counts5 = histc(errors(:,5),binrng)+counts4;
% bar(binrng,counts5)
% hold on
% bar(binrng,counts4)
% bar(binrng,counts3)
% bar(binrng,counts2)
% bar(binrng,counts1)
% xlabel('PC error')
% title('Errors for all folds','color','k')


% 
% 
% subplot(2,10,7)
% [r m b] = regression(Results7(:)',Validation7(:)');
% title(['EOF 7 R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
% %plot([0.5 4.5],[0.5 4.5],'w--')
% plot([-6 6],[-6 6],'k--')
% %lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
% %ylim([-6 6])
% %xlim([-6 6])
% subplot(2,10,17)
% errors = Results7-Validation7;
% binrng = -1:0.1:1;
% counts1 = histc(errors(:,1),binrng);
% counts2 = histc(errors(:,2),binrng)+counts1;
% counts3 = histc(errors(:,3),binrng)+counts2;
% counts4 = histc(errors(:,4),binrng)+counts3;
% counts5 = histc(errors(:,5),binrng)+counts4;
% bar(binrng,counts5)
% hold on
% bar(binrng,counts4)
% bar(binrng,counts3)
% bar(binrng,counts2)
% bar(binrng,counts1)
% xlabel('PC error')
% title('Errors for all folds','color','k')
% 
% 

% 
% subplot(2,10,8)
% [r m b] = regression(Resultsd(:)',Validationd(:)');
% title(['EOF Diff R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
% %plot([0.5 4.5],[0.5 4.5],'w--')
% plot([-6 6],[-6 6],'k--')
% %lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
% %ylim([-6 6])
% %xlim([-6 6])
% subplot(2,10,18)
% errors = Resultsd-Validationd;
% binrng = -1:0.1:1;
% counts1 = histc(errors(:,1),binrng);
% counts2 = histc(errors(:,2),binrng)+counts1;
% counts3 = histc(errors(:,3),binrng)+counts2;
% counts4 = histc(errors(:,4),binrng)+counts3;
% counts5 = histc(errors(:,5),binrng)+counts4;
% bar(binrng,counts5)
% hold on
% bar(binrng,counts4)
% bar(binrng,counts3)
% bar(binrng,counts2)
% bar(binrng,counts1)
% xlabel('PC error')
% title('Errors for all folds','color','k')
% 

subplot(2,5,1)
[r m b] = regression(Resultsdist(:)',Validationdist(:)');
title(['EOF Diff R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
plot([20 160],[20 160],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
%ylim([-6 6])
%xlim([-6 6])
subplot(2,5,6)
errors = Resultsdist-Validationdist;
binrng = -10:1:10;
counts1 = histc(errors(:,1),binrng);
counts2 = histc(errors(:,2),binrng)+counts1;
counts3 = histc(errors(:,3),binrng)+counts2;
counts4 = histc(errors(:,4),binrng)+counts3;
counts5 = histc(errors(:,5),binrng)+counts4;
bar(binrng,counts5)
hold on
bar(binrng,counts4)
bar(binrng,counts3)
bar(binrng,counts2)
bar(binrng,counts1)
xlabel('Top distance error')
title('Errors for all folds','color','k')


subplot(2,5,7)
[r m b] = regression(Resultsdist2(:)',Validationdist2(:)');
title(['EOF Diff R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
plot([-3 3],[-3 3],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
%ylim([-6 6])
%xlim([-6 6])
subplot(2,5,7)
errors = Resultsdist2-Validationdist2;
binrng = -1:0.1:1;
counts1 = histc(errors(:,1),binrng);
counts2 = histc(errors(:,2),binrng)+counts1;
counts3 = histc(errors(:,3),binrng)+counts2;
counts4 = histc(errors(:,4),binrng)+counts3;
counts5 = histc(errors(:,5),binrng)+counts4;
bar(binrng,counts5)
hold on
bar(binrng,counts4)
bar(binrng,counts3)
bar(binrng,counts2)
bar(binrng,counts1)
xlabel('Bottom distance error')
title('Errors for all folds','color','k')

subplot(2,5,3)
[r m b] = regression(ResultsElev(:)',ValidationElev(:)');
title(['EOF Diff R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
plot([-3 3],[-3 3],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
%ylim([-6 6])
%xlim([-6 6])
subplot(2,5,8)
errors = ResultsElev-ValidationElev;
binrng = -1:0.1:1;
counts1 = histc(errors(:,1),binrng);
counts2 = histc(errors(:,2),binrng)+counts1;
counts3 = histc(errors(:,3),binrng)+counts2;
counts4 = histc(errors(:,4),binrng)+counts3;
counts5 = histc(errors(:,5),binrng)+counts4;
bar(binrng,counts5)
hold on
bar(binrng,counts4)
bar(binrng,counts3)
bar(binrng,counts2)
bar(binrng,counts1)
xlabel('top elev error')
title('Errors for all folds','color','k')


subplot(2,5,4)
[r m b] = regression(ResultsElev2(:)',ValidationElev2(:)');
title(['EOF Diff R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
plot([-3 3],[-3 3],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
%ylim([-6 6])
%xlim([-6 6])
subplot(2,5,9)
errors = ResultsElev2-ValidationElev2;
binrng = -1:0.1:1;
counts1 = histc(errors(:,1),binrng);
counts2 = histc(errors(:,2),binrng)+counts1;
counts3 = histc(errors(:,3),binrng)+counts2;
counts4 = histc(errors(:,4),binrng)+counts3;
counts5 = histc(errors(:,5),binrng)+counts4;
bar(binrng,counts5)
hold on
bar(binrng,counts4)
bar(binrng,counts3)
bar(binrng,counts2)
bar(binrng,counts1)
xlabel('bottom elev error')
title('Errors for all folds','color','k')




subplot(2,5,5)
[r m b] = regression(Resultsheight(:)',Validationheight(:)');
title(['EOF Diff R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
plot([-3 3],[-3 3],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
%ylim([-6 6])
%xlim([-6 6])
subplot(2,5,10)
errors = Resultsheight-Validationheight;
binrng = -1:0.1:1;
counts1 = histc(errors(:,1),binrng);
counts2 = histc(errors(:,2),binrng)+counts1;
counts3 = histc(errors(:,3),binrng)+counts2;
counts4 = histc(errors(:,4),binrng)+counts3;
counts5 = histc(errors(:,5),binrng)+counts4;
bar(binrng,counts5)
hold on
bar(binrng,counts4)
bar(binrng,counts3)
bar(binrng,counts2)
bar(binrng,counts1)
xlabel('height error')
title('Errors for all folds','color','k')



% set(gcf,'Position',[100, 100, 1300, 500])
% 
% set(gcf,'color','k')
% subplot(121)
% %xlim([0.5 4.5])
% %ylim([0.5 4.5])
% xlim([-2 2.75])
% ylim([-2 2.75])
% set(gca,'xcolor','w')
% set(gca,'ycolor','w')
% set(lg,'color','k')
% set(lg,'textcolor','w')
% set(gca,'color','k')
% 
% subplot(122)
% set(gca,'xcolor','w')
% set(gca,'ycolor','w')
% set(gca,'color','k')
%%

%sizedata = length(design);
%[reorder] = randperm(sizedata);
%design = design(reorder,:);
sizedata = length(designParams);
%[reorder] = randperm(sizedata);

%reorder = [1:95];
sets = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60];%, 45, 49];
totals = [50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300];%, 225, 245];
%sets = [10, 15, 20, 25, 30]; %,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96];
%totals = [50, 75, 100, 125, 150]; %,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480];
    
%figure

for hh = 1:length(sets)
    
    size = sets(hh); % how many to be removed from the library

    total = totals(hh); % how big is the library in this iteration
    reorder = randperm(total); % randomize the order of the library
    
    % apply the reorder to all design and targets
    design = all(reorder,:);
    targets = predictands(reorder,:);
    starts = predictors(reorder,:);

    % make a new version of the reordered design and targers
    newall = design(1:total,:);
    newtargets = targets(1:total,:);
    newstarts = starts(1:total,:);
    
    distTargets = midBreakDist(reorder);
    distTargetsBot = midBreakDistBot(reorder);
    distHeights = midBreakHeight(reorder);
    
    distTargetsElev = midBreakElev(reorder);
    distTargetsElevBot = midBreakElevBot(reorder);
    newElev = distTargetsElev(1:total);
    newElevBot = distTargetsElevBot(1:total);
    
    newdist = distTargets(1:total);
    newdistBot = distTargetsBot(1:total);
    newheight = distHeights(1:total);
    
    c = 1;
    for qq = 1:5

        design1 = newall;
        design1(c:c+size-1,:) = [];
        target1 = newtargets(:,1);%-newstarts(:,1);
        target1(c:c+size-1,:) = [];
        target2 = newtargets(:,2);%-newstarts(:,2);
        target2(c:c+size-1,:) = [];
        target3 = newtargets(:,3);%-newstarts(:,3);
        target3(c:c+size-1,:) = [];
        target4 = newtargets(:,4);%-newstarts(:,4);
        target4(c:c+size-1,:) = [];
        target5 = newtargets(:,5);%-newstarts(:,5);
        target5(c:c+size-1,:) = [];
        target6 = newtargets(:,6);%-newstarts(:,6);
        target6(c:c+size-1,:) = [];
        target7 = newtargets(:,7);%-newstarts(:,7);
        target7(c:c+size-1,:) = [];
        targetd = newtargets(:,8);
        targetd(c:c+size-1,:) = [];        
        
        disp(['removing ',num2str(c),' to ',num2str(c+size-1),' from starting dataset'])
        data1 = newall(c:c+size-1,:);
        val1 = newtargets(c:c+size-1,1);%-newstarts(c:c+size-1,1);
        val2 = newtargets(c:c+size-1,2);%-newstarts(c:c+size-1,2);
        val3 = newtargets(c:c+size-1,3);%-newstarts(c:c+size-1,3);
        val4 = newtargets(c:c+size-1,4);%-newstarts(c:c+size-1,4);
        val5 = newtargets(c:c+size-1,5);%-newstarts(c:c+size-1,5);
        val6 = newtargets(c:c+size-1,6);%-newstarts(c:c+size-1,6);
        val7 = newtargets(c:c+size-1,7);%-newstarts(c:c+size-1,7);
        vald = newtargets(c:c+size-1,8);
        
        targetdist = newdist;
        targetdist(c:c+size-1) = [];
        targetdistBot = newdistBot;
        targetdistBot(c:c+size-1) = [];        
                
        targetElev = newElev;
        targetElev(c:c+size-1) = [];
        
        targetElevBot = newElevBot;
        targetElevBot(c:c+size-1) = [];  
        
        targetheight = newheight;
        targetheight(c:c+size-1) = [];
        valdist = newdist(c:c+size-1);
        valdistBot = newdistBot(c:c+size-1);
        valElev = newElev(c:c+size-1);
        valElevBot = newElevBot(c:c+size-1);
        valheight = newheight(c:c+size-1);
        disp(['holding ',num2str(c),' to ',num2str(c+size-1),' aside for validation'])
%         
%         disp(['working on EOF 1'])
%         
%         design1(:,18) = target1;
%         tbl1 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl1.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl1 = fitrgp(tbl1,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results1(1:size,qq) = predict(gprMdl1,data1);
%         Validation1(1:size,qq) = val1;
%         eval(['gprMdl',num2str(qq),' = gprMdl1;'])
% %         subplot(2,5,1)
% %         eval(['p',num2str(qq),' = scatter(val1,Results1(:,qq),''filled'')'])
% %         hold on
% %         xlabel('Xbeach PCs (m)')
% %         ylabel('GPR PCs (m)')
%         
%         
%         disp(['working on EOF 2'])
%         
%         design1(:,18) = target2;
%         tbl2 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl2.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl2 = fitrgp(tbl2,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results2(1:size,qq) = predict(gprMdl2,data1);
%         Validation2(1:size,qq) = val2;
%         eval(['gprMdl',num2str(qq),' = gprMdl1;'])
% %         subplot(2,5,2)
% %         eval(['p',num2str(qq),' = scatter(val2,Results2(:,qq),''filled'')'])
% %         hold on
% %         xlabel('Xbeach PCs (m)')
% %         ylabel('GPR PCs (m)')        
        
%         
%         disp(['working on EOF 3'])
%         
%         design1(:,18) = target3;
%         tbl3 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl3.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl3 = fitrgp(tbl3,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results3(1:size,qq) = predict(gprMdl3,data1);
%         Validation3(1:size,qq) = val3;
%         eval(['gprMdl',num2str(qq),' = gprMdl3;'])
% %         subplot(2,5,3)
% %         eval(['p',num2str(qq),' = scatter(val3,Results3(:,qq),''filled'')'])
% %         hold on
% %         xlabel('Xbeach PCs (m)')
% %         ylabel('GPR PCs (m)')
        
        
%         
%         disp(['working on EOF 4'])
%         
%         design1(:,18) = target4;
%         tbl4 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl4.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl4 = fitrgp(tbl4,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results4(1:size,qq) = predict(gprMdl4,data1);
%         Validation4(1:size,qq) = val4;
%         eval(['gprMdl',num2str(qq),' = gprMdl4;'])
% %         subplot(2,5,4)
% %         eval(['p',num2str(qq),' = scatter(val4,Results4(:,qq),''filled'')'])
% %         hold on
% %         xlabel('Xbeach PCs (m)')
% %         ylabel('GPR PCs (m)')        
        
        
%         disp(['working on EOF 5'])
%         
%         design1(:,18) = target5;
%         tbl5 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl5.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl5 = fitrgp(tbl5,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results5(1:size,qq) = predict(gprMdl5,data1);
%         Validation5(1:size,qq) = val5;
%         eval(['gprMdl',num2str(qq),' = gprMdl5;'])
% %         subplot(2,5,5)
% %         eval(['p',num2str(qq),' = scatter(val5,Results5(:,qq),''filled'')'])
% %         hold on
% %         xlabel('Xbeach PCs (m)')
% %         ylabel('GPR PCs (m)')        
        

%         disp(['working on EOF 6'])
%         
%         design1(:,18) = target6;
%         tbl6 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl6.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl6 = fitrgp(tbl6,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results6(1:size,qq) = predict(gprMdl6,data1);
%         Validation6(1:size,qq) = val6;
%         eval(['gprMdl',num2str(qq),' = gprMdl6;'])
% %         subplot(2,5,5)
% %         eval(['p',num2str(qq),' = scatter(val5,Results5(:,qq),''filled'')'])
% %         hold on
% %         xlabel('Xbeach PCs (m)')
% %         ylabel('GPR PCs (m)')    


%         disp(['working on EOF 7'])
%         
%         design1(:,18) = target7;
%         tbl7 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbl7.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdl7 = fitrgp(tbl7,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Results7(1:size,qq) = predict(gprMdl7,data1);
%         Validation7(1:size,qq) = val7;
%         eval(['gprMdl',num2str(qq),' = gprMdl7;'])
% %         subplot(2,5,5)
% %         eval(['p',num2str(qq),' = scatter(val5,Results5(:,qq),''filled'')'])
% %         hold on
% %         xlabel('Xbeach PCs (m)')
% %         ylabel('GPR PCs (m)')   
% 
%         disp(['working on EOF 8'])
%         
%         design1(:,18) = targetd;
%         tbld = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
%             design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
%             design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
%             design1(:,17),design1(:,18));
%         tbld.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
%         gprMdld = fitrgp(tbld,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
%         Resultsd(1:size,qq) = predict(gprMdld,data1);
%         Validationd(1:size,qq) = vald;
%         eval(['gprMdl',num2str(qq),' = gprMdld;'])
% %         subplot(2,5,5)
% %         eval(['p',num2str(qq),' = scatter(val5,Results5(:,qq),''filled'')'])
% %         hold on
% %         xlabel('Xbeach PCs (m)')
% %         ylabel('GPR PCs (m)')   

        disp('working on scarp #1 distance')
        
        design1(:,18) = targetdist;
        tbldist = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
            design1(:,17),design1(:,18));
        tbldist.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
        gprMdldist = fitrgp(tbldist,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsdist(1:size,qq) = predict(gprMdldist,data1);
        Validationdist(1:size,qq) = valdist;
        eval(['gprMdl',num2str(qq),' = gprMdldist;'])
   
        disp('working on scarp #2 distance')
        
        
        design1(:,18) = targetdistBot;
        tbldist2 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
            design1(:,17),design1(:,18));
        tbldist2.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
        gprMdldist2 = fitrgp(tbldist2,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsdist2(1:size,qq) = predict(gprMdldist2,data1);
        Validationdist2(1:size,qq) = valdistBot;
        eval(['gprMdl',num2str(qq),' = gprMdldist2;'])

        
         
        disp('working on scarp #1 elev')
        
        design1(:,18) = targetElev;
        tbldist = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
            design1(:,17),design1(:,18));
        tbldist.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
        gprMdlelev = fitrgp(tbldist,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        ResultsElev(1:size,qq) = predict(gprMdlelev,data1);
        ValidationElev(1:size,qq) = valElev;
        eval(['gprMdl',num2str(qq),' = gprMdlelev;'])
      
        
        
        disp('working on scarp #2 elev')
        
        design1(:,18) = targetElevBot;
        tbldist2 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
            design1(:,17),design1(:,18));
        tbldist2.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
        gprMdlelev2 = fitrgp(tbldist2,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        ResultsElev2(1:size,qq) = predict(gprMdlelev2,data1);
        ValidationElev2(1:size,qq) = valdistBot;
        eval(['gprMdl',num2str(qq),' = gprMdlelev2;'])

        
        design1(:,18) = targetheight;
        tblheight = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),...
            design1(:,17),design1(:,18));
        tblheight.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOFd','M2','N2','K1','S2','target'};
        gprMdlheight= fitrgp(tblheight,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsheight(1:size,qq) = predict(gprMdlheight,data1);
        Validationheight(1:size,qq) = valheight;
        eval(['gprMdl',num2str(qq),' = gprMdlheight;'])


        clear gprMdl1 gprMdl2 gprMdl3 gprMdl4 gprMdl5 gprMdl6 gprMdl7 gprMdld gprMdl8
        c = c +size;
        
    end
%     
%     errors1 = Results1-Validation1;
%     squarederrors1 = (errors1).^2;
%     mse1 = mean(squarederrors1);
%     RMSE1(1:5,hh) = sqrt(mse1);
%     
%     errors2 = Results2-Validation2;
%     squarederrors2 = (errors2).^2;
%     mse2 = mean(squarederrors2);
%     RMSE2(1:5,hh) = sqrt(mse2);
%     
%     errors3 = Results3-Validation3;
%     squarederrors3 = (errors3).^2;
%     mse3 = mean(squarederrors3);
%     RMSE3(1:5,hh) = sqrt(mse3);
%     
%     errors4 = Results4-Validation4;
%     squarederrors4 = (errors4).^2;
%     mse4 = mean(squarederrors4);
%     RMSE4(1:5,hh) = sqrt(mse4);
%     
%     errors5 = Results5-Validation5;
%     squarederrors5 = (errors5).^2;
%     mse5 = mean(squarederrors5);
%     RMSE5(1:5,hh) = sqrt(mse5);
%     
%     errors6 = Results6-Validation6;
%     squarederrors6 = (errors6).^2;
%     mse6 = mean(squarederrors6);
%     RMSE6(1:5,hh) = sqrt(mse6);
%     
%     errors7 = Results7-Validation7;
%     squarederrors7 = (errors7).^2;
%     mse7 = mean(squarederrors7);
%     RMSE7(1:5,hh) = sqrt(mse7);
%     
%     errorsd = Resultsd-Validationd;
%     squarederrorsd = (errorsd).^2;
%     msed = mean(squarederrorsd);
%     RMSEd(1:5,hh) = sqrt(msed);
%     
    errorsdist = Resultsdist-Validationdist;
    squarederrorsdist = (errorsdist).^2;
    msedist = mean(squarederrorsdist);
    RMSEdist(1:5,hh) = sqrt(msedist);
    
    errorsdist2 = Resultsdist2-Validationdist2;
    squarederrorsdist2 = (errorsdist2).^2;
    msedist2 = mean(squarederrorsdist2);
    RMSEdist2(1:5,hh) = sqrt(msedist2);    
    
    errorselev = ResultsElev-ValidationElev;
    squarederrorselev = (errorselev).^2;
    mseelev = mean(squarederrorselev);
    RMSEelev(1:5,hh) = sqrt(mseelev);
    
    errorselev2 = ResultsElev2-ValidationElev2;
    squarederrorselev2 = (errorselev2).^2;
    mseelev2 = mean(squarederrorselev2);
    RMSEelev2(1:5,hh) = sqrt(mseelev2);    
    
    errorsheight = Resultsheight-Validationheight;
    squarederrorsheight = (errorsheight).^2;
    mseheight = mean(squarederrorsheight);
    RMSEheight(1:5,hh) = sqrt(mseheight);
end


figure

set(gcf,'color','w');



subplot(151)
boxplot(RMSEdist,(totals-sets))
hold on
lines = findobj(gcf,'type','line','Tag','Median');
set(lines,'Color','r')
ylabel('RMSE (m)')
xlabel('# used to calibrate')
title('Distance 1')
% ylim([0.25 1.3])

subplot(152)
boxplot(RMSEdist2,(totals-sets))
hold on
lines = findobj(gcf,'type','line','Tag','Median');
set(lines,'Color','r')
ylabel('RMSE (m)')
xlabel('# used to calibrate')
title('Distance 2')

subplot(153)
boxplot(RMSEelev,(totals-sets))
hold on
lines = findobj(gcf,'type','line','Tag','Median');
set(lines,'Color','r')
ylabel('RMSE (m)')
xlabel('# used to calibrate')
title('Elev 1')

subplot(154)
boxplot(RMSEelev2,(totals-sets))
hold on
lines = findobj(gcf,'type','line','Tag','Median');
set(lines,'Color','r')
ylabel('RMSE (m)')
xlabel('# used to calibrate')
title('Elev 2')

subplot(155)
boxplot(RMSEheight,(totals-sets))
hold on
lines = findobj(gcf,'type','line','Tag','Median');
set(lines,'Color','r')
ylabel('RMSE (m)')
xlabel('# used to calibrate')
title('Height')
% ylim([0.25 1.3])

% 
% subplot(181)
% boxplot(RMSE1,(totals-sets))
% hold on
% lines = findobj(gcf,'type','line','Tag','Median');
% set(lines,'Color','r')
% ylabel('RMSE (m)')
% xlabel('# used to calibrate')
% title('EOF1')
% ylim([0.25 1.3])
% subplot(182)
% boxplot(RMSE2,(totals-sets))
% hold on
% lines = findobj(gcf,'type','line','Tag','Median');
% set(lines,'Color','r')
% ylabel('RMSE (m)')
% xlabel('# used to calibrate')
% title('EOF2')
% ylim([0.25 1.3])
% subplot(183)
% boxplot(RMSE3,(totals-sets))
% hold on
% lines = findobj(gcf,'type','line','Tag','Median');
% set(lines,'Color','r')
% ylabel('RMSE (m)')
% xlabel('# used to calibrate')
% title('EOF3')
% ylim([0.25 1.3])
% subplot(184)
% boxplot(RMSE4,(totals-sets))
% hold on
% lines = findobj(gcf,'type','line','Tag','Median');
% set(lines,'Color','r')
% ylabel('RMSE (m)')
% xlabel('# used to calibrate')
% title('EOF4')
% ylim([0.25 1.3])
% subplot(185)
% boxplot(RMSE5,(totals-sets))
% hold on
% lines = findobj(gcf,'type','line','Tag','Median');
% set(lines,'Color','r')
% ylabel('RMSE (m)')
% xlabel('# used to calibrate')
% title('EOF5')
% ylim([0.25 1.3])
% subplot(186)
% boxplot(RMSE6,(totals-sets))
% hold on
% lines = findobj(gcf,'type','line','Tag','Median');
% set(lines,'Color','r')
% ylabel('RMSE (m)')
% xlabel('# used to calibrate')
% title('EOF6')
% ylim([0.25 1.3])
% subplot(187)
% boxplot(RMSE7,(totals-sets))
% hold on
% lines = findobj(gcf,'type','line','Tag','Median');
% set(lines,'Color','r')
% ylabel('RMSE (m)')
% xlabel('# used to calibrate')
% title('EOF7')
% ylim([0.25 1.3])
% 
% subplot(188)
% boxplot(RMSEd,(totals-sets))
% hold on
% lines = findobj(gcf,'type','line','Tag','Median');
% set(lines,'Color','r')
% ylabel('RMSE (m)')
% xlabel('# used to calibrate')
% title('EOFDiff')
% ylim([0.25 1.3])