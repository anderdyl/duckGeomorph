
clear






casenumbers = [1:150];
%load('/home/server/pi/homes/danderso/COSMOS_lite/LHS_1.mat')
load('profileTrainingTemp2.mat')

designParams = predictors(casenumbers,:);

%Normalization of the selected cases with MaxDiss
maxHs=max(designParams(:,1));  minHs=min(designParams(:,1));
maxTp=max(designParams(:,2));  minTp=min(designParams(:,2));
maxDir=max(designParams(:,3));  minDir=min(designParams(:,3));
maxNTR=max(designParams(:,4));  minNTR=min(designParams(:,4));
maxDur=max(designParams(:,5));  minDur=min(designParams(:,5));
maxEOF1=max(designParams(:,6));  minEOF1=min(designParams(:,6));
maxEOF2=max(designParams(:,7));  minEOF2=min(designParams(:,7));
maxEOF3=max(designParams(:,8));  minEOF3=min(designParams(:,8));
maxEOF4=max(designParams(:,9));  minEOF4=min(designParams(:,9));
maxEOF5=max(designParams(:,10));  minEOF5=min(designParams(:,10));

subset_n(:,1)=(designParams(:,1)-minHs)./(maxHs-minHs);
subset_n(:,2)=(designParams(:,2)-minTp)./(maxTp-minTp);
subset_n(:,3)=designParams(:,3)*pi/180;
subset_n(:,4)=(designParams(:,4)-minNTR)./(maxNTR-minNTR);
subset_n(:,5)=(designParams(:,5)-minDur)./(maxDur-minDur);
subset_n(:,6)=(designParams(:,6)-minEOF1)./(maxEOF1-minEOF1);
subset_n(:,7)=(designParams(:,7)-minEOF2)./(maxEOF1-minEOF2);
subset_n(:,8)=(designParams(:,8)-minEOF3)./(maxEOF3-minEOF3);
subset_n(:,9)=(designParams(:,9)-minEOF4)./(maxEOF1-minEOF4);
subset_n(:,10)=(designParams(:,10)-minEOF5)./(maxEOF5-minEOF5);

all = subset_n;

calibration = 1:120;
validation = 121:150;

target90 = predictands(casenumbers,1)-predictors(casenumbers,6);
design1 = all(calibration,:);
design1(calibration,11) = target90(calibration);
tbl = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11));
tbl.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','target'};
gprMdl = fitrgp(tbl,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)

target902 = predictands(casenumbers,2)-predictors(casenumbers,7);
design1 = all(calibration,:);
design1(calibration,11) = target902(calibration);
tbl2 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11));
tbl2.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','target'};
gprMdl2 = fitrgp(tbl2,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)

target903 = predictands(casenumbers,3)-predictors(casenumbers,8);
design1 = all(calibration,:);
design1(calibration,11) = target903(calibration);
tbl3 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11));
tbl3.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','target'};
gprMdl3 = fitrgp(tbl3,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)

target904 = predictands(casenumbers,4)-predictors(casenumbers,9);
design1 = all(calibration,:);
design1(calibration,11) = target904(calibration);
tbl4 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11));
tbl4.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','target'};
gprMdl4 = fitrgp(tbl4,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)

target905 = predictands(casenumbers,5)-predictors(casenumbers,10);
design1 = all(calibration,:);
design1(calibration,11) = target905(calibration);
tbl5 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11));
tbl5.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','target'};
gprMdl5 = fitrgp(tbl5,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1)



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

%%
figure();
plot(ypred,target90(validation),'o')
hold on
plot(ypred2,target902(validation),'o')
plot(ypred3,target903(validation),'o')
plot(ypred4,target904(validation),'o')
plot(ypred5,target905(validation),'o')

%plot(ypred,'b')
xlabel('true')
ylabel('predicted')
%legend({'data','predictions'},'Location','Best')
%axis([0 4300 0 30])
%hold off

final(1:30,1) = ypred + predictors(validation,6);
final(1:30,2) = ypred2 + predictors(validation,7);
final(1:30,3) = ypred3 + predictors(validation,8);
final(1:30,4) = ypred4 + predictors(validation,9);
final(1:30,5) = ypred5 + predictors(validation,10);

%%

%sizedata = length(design);
%[reorder] = randperm(sizedata);
%design = design(reorder,:);
sizedata = length(designParams);
[reorder] = randperm(sizedata);

sets = [30]; %,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96];
totals = [150]; %,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480];

figure
set(gcf,'color','w');

for hh = 1:length(sets)
    size = sets(hh);
    total = totals(hh);
    design = all(reorder,:);
    targets = predictands(reorder,:);
    starts = predictors(reorder,:);
    
    newall = design(1:total,:);
    newtargets = targets(1:total,:);
    
    c = 1;
    for qq = 1:5

        design1 = newall;
        design1(c:c+size-1,:) = [];
        target1 = newtargets(:,1)-starts(:,6);
        target1(c:c+size-1,:) = [];
        target2 = newtargets(:,2)-starts(:,7);
        target2(c:c+size-1,:) = [];
        target3 = newtargets(:,3)-starts(:,8);
        target3(c:c+size-1,:) = [];
        target4 = newtargets(:,4)-starts(:,9);
        target4(c:c+size-1,:) = [];
        target5 = newtargets(:,5)-starts(:,10);
        target5(c:c+size-1,:) = [];
        
        disp(['removing ',num2str(c),' to ',num2str(c+size-1),' from starting dataset'])
        data1 = newall(c:c+size-1,:);
        val1 = newtargets(c:c+size-1,1)-starts(c:c+size-1,6);
        val2 = newtargets(c:c+size-1,2)-starts(c:c+size-1,7);
        val3 = newtargets(c:c+size-1,3)-starts(c:c+size-1,8);
        val4 = newtargets(c:c+size-1,4)-starts(c:c+size-1,9);
        val5 = newtargets(c:c+size-1,5)-starts(c:c+size-1,10);
        
        disp(['holding ',num2str(c),' to ',num2str(c+size-1),' aside for validation'])
        
        disp(['working on EOF 1'])
        
        design1(:,11) = target1;
        tbl1 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11));
        tbl1.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','target'};
        gprMdl1 = fitrgp(tbl1,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Results1(1:size,qq) = predict(gprMdl1,data1);
        Validation1(1:size,qq) = val1;
        eval(['gprMdl',num2str(qq),' = gprMdl1;'])
        subplot(2,5,1)
        eval(['p',num2str(qq),' = scatter(val1,Results1(:,qq),''filled'')'])
        hold on
        xlabel('Xbeach PCs (m)')
        ylabel('GPR PCs (m)')
        
        
        disp(['working on EOF 2'])
        
        design1(:,11) = target2;
        tbl2 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11));
        tbl2.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','target'};
        gprMdl2 = fitrgp(tbl2,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Results2(1:size,qq) = predict(gprMdl2,data1);
        Validation2(1:size,qq) = val2;
        eval(['gprMdl',num2str(qq),' = gprMdl1;'])
        subplot(2,5,2)
        eval(['p',num2str(qq),' = scatter(val2,Results2(:,qq),''filled'')'])
        hold on
        xlabel('Xbeach PCs (m)')
        ylabel('GPR PCs (m)')        
        
        
        disp(['working on EOF 3'])
        
        design1(:,11) = target3;
        tbl3 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11));
        tbl3.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','target'};
        gprMdl3 = fitrgp(tbl3,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Results3(1:size,qq) = predict(gprMdl3,data1);
        Validation3(1:size,qq) = val3;
        eval(['gprMdl',num2str(qq),' = gprMdl3;'])
        subplot(2,5,3)
        eval(['p',num2str(qq),' = scatter(val3,Results3(:,qq),''filled'')'])
        hold on
        xlabel('Xbeach PCs (m)')
        ylabel('GPR PCs (m)')
        
        
        
        disp(['working on EOF 4'])
        
        design1(:,11) = target4;
        tbl4 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11));
        tbl4.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','target'};
        gprMdl4 = fitrgp(tbl4,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Results4(1:size,qq) = predict(gprMdl4,data1);
        Validation4(1:size,qq) = val4;
        eval(['gprMdl',num2str(qq),' = gprMdl4;'])
        subplot(2,5,4)
        eval(['p',num2str(qq),' = scatter(val4,Results4(:,qq),''filled'')'])
        hold on
        xlabel('Xbeach PCs (m)')
        ylabel('GPR PCs (m)')        
        
        
        disp(['working on EOF 5'])
        
        design1(:,11) = target5;
        tbl5 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11));
        tbl5.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','EOF1','EOF2','EOF3','EOF4','EOF5','target'};
        gprMdl5 = fitrgp(tbl5,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Results5(1:size,qq) = predict(gprMdl5,data1);
        Validation5(1:size,qq) = val5;
        eval(['gprMdl',num2str(qq),' = gprMdl5;'])
        subplot(2,5,5)
        eval(['p',num2str(qq),' = scatter(val5,Results5(:,qq),''filled'')'])
        hold on
        xlabel('Xbeach PCs (m)')
        ylabel('GPR PCs (m)')        
        
        %clear gprMdl1 gprMdl2 gprMdl3 gprMdl4 gprMdl5
        c = c +size;
        
    end
end



subplot(2,5,1)
[r m b] = regression(Results1(:)',Validation1(:)');
title(['EOF 1 R^{2} = ',num2str(round(r*1000)/1000),' with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
p6 = plot([-11 11],[-11 11],'k--');
lg = legend([p1 p2 p3 p4 p5 p6],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','1-to-1','location','northwest');
ylim([-4 3])
xlim([-4 3])
subplot(2,5,6)
errors = Results1-Validation1;
binrng = -1:.1:1;
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
xlabel('PC error')
title('Errors for all folds','color','k')



subplot(2,5,2)
[r m b] = regression(Results2(:)',Validation2(:)');
title(['EOF 2 R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
plot([-10 10],[-10 10],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
ylim([-12 3])
xlim([-12 3])
subplot(2,5,7)
errors = Results2-Validation2;
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
xlabel('PC error')
title('Errors for all folds','color','k')


subplot(2,5,3)
[r m b] = regression(Results3(:)',Validation3(:)');
title(['EOF 3 R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
plot([-6 6],[-6 6],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
%ylim([-6 6])
%xlim([-6 6])
subplot(2,5,8)
errors = Results3-Validation3;
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
xlabel('PC error')
title('Errors for all folds','color','k')



subplot(2,5,4)
[r m b] = regression(Results4(:)',Validation4(:)');
title(['EOF 4 R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
plot([-5 5],[-5 5],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
%ylim([-5 5])
%xlim([-5 5])
subplot(2,5,9)
errors = Results4-Validation4;
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
xlabel('PC error')
title('Errors for all folds','color','k')



subplot(2,5,5)
[r m b] = regression(Results5(:)',Validation5(:)');
title(['EOF 5 R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
plot([-6 6],[-6 6],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
%ylim([-6 6])
%xlim([-6 6])
subplot(2,5,10)
errors = Results5-Validation5;
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
xlabel('PC error')
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

%{
clear
load(fullfile(matlabroot,'examples','stats','gprdata2.mat'))

gprMdl1 = fitrgp(x,y,'KernelFunction','squaredexponential');

sigma0 = 0.2;
kparams0 = [3.5, 6.2];
gprMdl2 = fitrgp(x,y,'KernelFunction','squaredexponential',...
     'KernelParameters',kparams0,'Sigma',sigma0);

ypred1 = resubPredict(gprMdl1);
ypred2 = resubPredict(gprMdl2);

figure();
plot(x,y,'r.');
hold on
plot(x,ypred1,'b');
plot(x,ypred2,'g');
xlabel('x');
ylabel('y');
legend({'data','default kernel parameters',...
'kparams0 = [3.5,6.2], sigma0 = 0.2'},...
'Location','Best');
title('Impact of initial kernel parameter values');
hold off

%}



