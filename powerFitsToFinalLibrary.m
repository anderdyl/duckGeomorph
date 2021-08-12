clear
load('trialsABFinal1250.mat')
load('scarpPointsFinal1250.mat')
load('hypoFinal500.mat')

trialNum = [1:1250];
hypos = hypo(trialNum,1:19);


% Bad fits

% 62 is a big berm removal
xLoc(63) = 279-217.5;
%a(63) = -0.7688733;
%b(63) = 0.84081691;

% 112 is a dune top erosion point...

% 134 eating away at the dune face
xLoc(135) = 235-164;
a(135) = -0.1447021;
b(135) = 0.64687113;

% 314 eating away at the dune face
xLoc(315) = 238-165;
a(315) = -0.16603765;
b(315) = 0.6616657;

% 319 looks like a normal scarp scenario...
xLoc(320) = 241-191.5;
a(320) = -0.77136149;
b(320) = 0.3802362;

%337 is a tiny scarp with more lower profile erosion
xLoc(338) = 224-210;
a(338) = -0.13217419;
b(338) = 0.69946215;

% 468 looks normal, maybe a bigger berm?
xLoc(469) = 252.5-211.5;
a(469) = -0.31718905;
b(469) = 0.54026136;

% 500 big summer berm on wide beach, inverted B
xLoc(501) = 271.5-218;
a(501) = -0.07538966;
b(501) = 0.78995607;

% 502 eats away to dune top
xLoc(503) = 180.5-134;
a(503) = -1.71586886;
b(503) = 0.27091895;

% 537 eats away at the dune face
xLoc(538) = 190.5-138.5;
a(538) = -2.12625674;
b(538) = 0.26282211;

% 549 dune overtopping
xLoc(550) = 183-134.5;
a(550) = -1.10235711;
b(550) = 0.38508868;

% 643 majority of erosion below 0 contour
xLoc(644) = 200-180.5;
a(644) = -0.08560981;
b(644) = 0.7948455;

% 648 is small erosion high above 0 contour
xLoc(649) = 211.5-168.5;
a(649) = -0.24661496;
b(649) = 0.53947659;

%658 is berm erosion
xLoc(659) = 247.5-196.5;
a(659) = -0.12493934;
b(659) = 0.7423418;

%662 is wide beach with scarp
xLoc(663) = 248.5-197.5;
a(663) = -0.33626497;
b(663) = 0.50353753;

%673 very flat removing a berm
xLoc(674) = 279.5-214.5;
a(674) = -0.11027043;
b(674) = 0.71955629;

% 693 small erosion high on beach above 0 contour
xLoc(694) = 207-167.5;
a(694) = -0.16779858;
b(694) = 0.6612156;

% 765 very flat slope removing a berm
xLoc(766) = 276-211.5;
a(766) = -0.11883329;
b(766) = 0.73470203;

% 802 very flat berm removal
xLoc(803) = 234.5-184;
a(803) = -0.10909526;
b(803) = 0.72666294;

% eroded into the dune
xLoc(812) = 216-160;
a(812) = -0.48472307;
b(812) = 0.39798088;

% 845 small erosion high on beach above 0 contour
xLoc(846) = 253.5-194.5;
a(846) = -0.18684069;
b(846) = 0.64304391;

% 958 removed the dune
%xLoc(959) = 
a(959) = -0.65819555;
b(959) = 0.50844227;

%1009 removes a lot of sediment under 0 contour
xLoc(1010) = 213 - 169;
a(1010) = -0.14395834;
b(1010) = 0.70480727;

% 1080 very flat slope
xLoc(1081) = 241.5-193;
a(1081) = -0.6628543;
b(1081) = 0.84010997;

% 1116 very flat slope
xLoc(1117) = 227-175;
a(1117) = -0.05429171;
b(1117) = 0.83676292;

% 1161 very flat slope berm removal
xLoc(1162) = 244-185;
a(1162) = -0.12285635;
b(1162) = 0.65809217;

% 1163 dune overtopping/erosion

% % remove = find(isnan(zeroCrossing));
% % 
% %oddProfiles = [3,9,49,55,90,133,167,181,304,319,388,408,411,416,419,459,481];
% %oddProfiles = [9, 32, 63, 90, 104, 117, 131, 150, 172, 182, 237, 432];
% remove = [134,314,319,337,500,693,845,1009,1080,1116]+1;
remove = [112,134,314,549,648,693,811,845,958,1080,1116,1161,1163]+1;
% oddProfiles = [112,502,537,549,958]+1;
% oddProfiles2 = [468,658,662,1161]+1;
% %remove = oddProfiles;
% remove = [remove,oddProfiles];
a(remove) = [];
b(remove) = [];
xLoc(remove) = [];
hypos(remove,:) = [];
trialNum(remove) = [];
% 

figure
subplot(1,3,1)
hist(xLoc,20)
subplot(1,3,2)
hist(a,30)
subplot(1,3,3)
hist(b,30)
xlim([0,1.7])

% weirdProfiles = find(b > 1);
% 
% remove2 = weirdProfiles;
% a(remove2) = [];
% b(remove2) = [];
% xLoc(remove2) = [];
% %zLoc(remove2) = [];
% hypos(remove2,:) = [];
% trialNum(remove2) = [];



designParams = hypos(1:length(a),:);


M2amp = 0.49;
N2amp = 0.114;
K1amp = 0.087;
S2amp = 0.088;

M2speed = 28.984104;
N2speed = 28.43973;
K1speed = 15.041069;
S2speed = 30.0;



M2comp = M2amp*cosd(designParams(:,14));
N2comp = N2amp*cosd(designParams(:,15));
K1comp = K1amp*cosd(designParams(:,16));
S2comp = S2amp*cosd(designParams(:,17));

totalWL = sum([M2comp';N2comp';K1comp';S2comp']);

designParams(:,20) = totalWL';

maxMSL = max(designParams(:,10)); minMSL = min(designParams(:,10));

maxHs=max(designParams(:,11));  minHs=min(designParams(:,11));
maxTp=max(designParams(:,12));  minTp=min(designParams(:,12));
maxDir=max(designParams(:,15)*pi/180);  minDir=min(designParams(:,15)*pi/180);

maxNTR=max(designParams(:,13));  minNTR=min(designParams(:,13));
maxDur=max(designParams(:,14));  minDur=min(designParams(:,14));

maxEOF1=max(designParams(:,1));  minEOF1=min(designParams(:,1));
maxEOF2=max(designParams(:,2));  minEOF2=min(designParams(:,2));
maxEOF3=max(designParams(:,3));  minEOF3=min(designParams(:,3));
maxEOF4=max(designParams(:,4));  minEOF4=min(designParams(:,4));
maxEOF5=max(designParams(:,5));  minEOF5=min(designParams(:,5));
maxEOF6=max(designParams(:,6));  minEOF6=min(designParams(:,6));
maxEOF7=max(designParams(:,7));  minEOF7=min(designParams(:,7));
maxEOF8=max(designParams(:,8));  minEOF8=min(designParams(:,8));
maxEOF9=max(designParams(:,9));  minEOF9=min(designParams(:,9));

maxM2=max(designParams(:,16)*pi/180);  minM2=min(designParams(:,16)*pi/180);
maxN2=max(designParams(:,17)*pi/180);  minN2=min(designParams(:,17)*pi/180);
maxK1=max(designParams(:,18)*pi/180);  minK1=min(designParams(:,18)*pi/180);
maxS2=max(designParams(:,19)*pi/180);  minS2=min(designParams(:,19)*pi/180);

maxTide=max(designParams(:,20));  minTide=min(designParams(:,20));

%maxDist = max(designParams(:,18)); minDist = min(designParams(:,18));
%maxHeight = max(designParams(:,(19)); minHeight = min(designParams(:,19));

subset_n(:,1)=(designParams(:,11)-minHs)./(maxHs-minHs);
subset_n(:,2)=(designParams(:,12)-minTp)./(maxTp-minTp);
subset_n(:,3)=((designParams(:,15)*pi/180)-minDir)./(maxDir-minDir);
subset_n(:,4)=(designParams(:,13)-minNTR)./(maxNTR-minNTR);
subset_n(:,5)=(designParams(:,14)-minDur)./(maxDur-minDur);
subset_n(:,6)=(designParams(:,10)-minMSL)./(maxMSL-minMSL);

subset_n(:,7)=(designParams(:,1)-minEOF1)./(maxEOF1-minEOF1);
subset_n(:,8)=(designParams(:,2)-minEOF2)./(maxEOF2-minEOF2);
subset_n(:,9)=(designParams(:,3)-minEOF3)./(maxEOF3-minEOF3);
subset_n(:,10)=(designParams(:,4)-minEOF4)./(maxEOF4-minEOF4);
subset_n(:,11)=(designParams(:,5)-minEOF5)./(maxEOF5-minEOF5);
subset_n(:,12)=(designParams(:,6)-minEOF6)./(maxEOF6-minEOF6);
subset_n(:,13)=(designParams(:,7)-minEOF7)./(maxEOF7-minEOF7);
subset_n(:,14)=(designParams(:,8)-minEOF8)./(maxEOF8-minEOF8);
subset_n(:,15)=(designParams(:,9)-minEOF9)./(maxEOF9-minEOF9);

subset_n(:,16)=(designParams(:,20)-minTide)./(maxTide-minTide);

% subset_n(:,14)=((designParams(:,14)*pi/180)-minM2)./(maxM2-minM2);
% subset_n(:,15)=((designParams(:,15)*pi/180)-minN2)./(maxN2-minN2);
% subset_n(:,16)=((designParams(:,16)*pi/180)-minK1)./(maxK1-minK1);
% subset_n(:,17)=((designParams(:,17)*pi/180)-minS2)./(maxS2-minS2);

%subset_n(:,18)=(designParams(:,18)-minDist)./(maxDist-minDist);
%subset_n(:,19)=(designParams(:,19)-minHeight)./(maxHeight-maxHeight);
%%

calWith = 1000;

disp(['working on distance to Erosion'])

design1 = subset_n(1:calWith,:);
trainNums = trialNum(1:calWith);
valNums = trialNum(calWith+1:end);

design1(:,17) = xLoc(1:calWith);
tbldist = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),design1(:,17));%,design1(:,16),...
    %design1(:,17),design1(:,18));
tbldist.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','MSL','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','tide','target'};
%gprMdldist = fitrgp(tbldist,'target','KernelFunction','ardsquaredexponential','FitMethod','sr','PredictMethod','fic','Standardize',1);
gprMdldist = fitrgp(tbldist,'target','BasisFunction','linear','KernelFunction','ardMatern52','FitMethod','sr','PredictMethod','fic','Standardize',1);

disp(['working on a parameter A'])
design1(:,17) = a(1:calWith);
tbl1 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),design1(:,17));%,design1(:,16),...
    %design1(:,17),design1(:,18));
tbl1.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','MSL','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','tide','target'};
gprMdla = fitrgp(tbl1,'target','BasisFunction','pureQuadratic','KernelFunction','ardMatern52','FitMethod','sr','PredictMethod','fic','Standardize',1);
     
disp(['working on a parameter B'])
design1(:,17) = b(1:calWith);
tbl1 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
    design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
    design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),design1(:,17));%,design1(:,16),...
    %design1(:,17),design1(:,18));
tbl1.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','MSL','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','tide','target'};
gprMdlb = fitrgp(tbl1,'target','BasisFunction','linear','KernelFunction','ardMatern52','FitMethod','sr','PredictMethod','fic','Standardize',1);



figure
subplot(131)
plot(xLoc(1:calWith),resubPredict(gprMdldist),'.')
title(['Loss on Dist Training model: ',num2str(resubLoss(gprMdldist)),''])
subplot(132)
plot(a(1:calWith),resubPredict(gprMdla),'.')
title(['Loss on A Training model: ',num2str(resubLoss(gprMdla)),''])
subplot(133)
plot(b(1:calWith),resubPredict(gprMdlb),'.')
title(['Loss on B Training model: ',num2str(resubLoss(gprMdlb)),''])



figure
d=16;
subplot(311)
sigmaMdist = gprMdldist.KernelInformation.KernelParameters(1:end-1,1);
plot((1:d)',log(sigmaMdist),'ro-');
%xlabel('Length scale number');
ylabel('Log of length scale');
xticks([1:16])
xlim([0 17])
xticklabels({'H_s','T_p','Dir','NTR','Dur','MSL','EOF_1','EOF_2','EOF_3','EOF_4','EOF_5','EOF_6','EOF_7','EOF_8','EOF_9','tide'})
subplot(312)
sigmaMa = gprMdla.KernelInformation.KernelParameters(1:end-1,1);
plot((1:d)',log(sigmaMa),'ro-');
%xlabel('Length scale number');
ylabel('Log of length scale');
xticks([1:16])
xlim([0 17])
xticklabels({'H_s','T_p','Dir','NTR','Dur','MSL','EOF_1','EOF_2','EOF_3','EOF_4','EOF_5','EOF_6','EOF_7','EOF_8','EOF_9','tide'})
subplot(313)
sigmaMb = gprMdlb.KernelInformation.KernelParameters(1:end-1,1);
plot((1:d)',log(sigmaMb),'ro-');
xlabel('Length scale number');
ylabel('Log of length scale');
xticks([1:16])
xlim([0 17])
xticklabels({'H_s','T_p','Dir','NTR','Dur','MSL','EOF_1','EOF_2','EOF_3','EOF_4','EOF_5','EOF_6','EOF_7','EOF_8','EOF_9','tide'})
set(gcf,'color','w')



ypreddist = predict(gprMdldist,subset_n(calWith+1:end,:));
ypreda = predict(gprMdla,subset_n(calWith+1:end,:));
ypredb = predict(gprMdlb,subset_n(calWith+1:end,:));


final(1:length(ypreddist),1) = ypreddist;% + predictors(validation,1);
final(1:length(ypreddist),2) = ypreda;% + predictors(validation,2);
final(1:length(ypreddist),3) = ypredb;% + predictors(validation,3);

valTrials = valNums;

%save('predictions3.mat','valTrials','ypreddist','ypreda','ypredb')


%%






%sizedata = length(design);
%[reorder] = randperm(sizedata);
%design = design(reorder,:);
sizeDataSet = 1235;
all = subset_n;
designParams = designParams(1:sizeDataSet,:);
a = a(1:sizeDataSet);
b = b(1:sizeDataSet);
xLoc = xLoc(1:sizeDataSet);
%midBreakElevBot = midBreakElevBot(1:sizeDataSet);
%midBreakHeight = midBreakHeight(1:300);

sizedata = length(designParams);
[reorder] = randperm(sizedata);

sets = [247]; %,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96];
totals = [1235]; %,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480];

figure
set(gcf,'color','w');

for hh = 1:length(sets)
    size = sets(hh);
    total = totals(hh);
    design = all(reorder,:);
%     targets = predictands(reorder,:);
%    starts = predictors(reorder,:);
    
    distTargets = xLoc(reorder);
    aTargets = a(reorder);
    bTargets = b(reorder);

   
    newall = design(1:total,:);
    newdist = distTargets(1:total);
    newa = aTargets(1:total);
    newb = bTargets(1:total);
    
    c = 1;
    for qq = 1:5

        design1 = newall;
        design1(c:c+size-1,:) = [];

        targetdist = newdist;
        targetdist(c:c+size-1) = [];
        targeta = newa;
        targeta(c:c+size-1) = [];
        targetb = newb;
        targetb(c:c+size-1) = [];
        
        
        disp(['removing ',num2str(c),' to ',num2str(c+size-1),' from starting dataset'])
        data1 = newall(c:c+size-1,:);

        valdist = newdist(c:c+size-1);
        vala = newa(c:c+size-1);
        valb = newb(c:c+size-1);

        disp(['holding ',num2str(c),' to ',num2str(c+size-1),' aside for validation'])
        

        disp(['working on distance to Erosion'])
        
        design1(:,17) = targetdist;
        tbldist = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),design1(:,17));%,design1(:,16),...
            %design1(:,17),design1(:,18));
        tbldist.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','MSL','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','tides','target'};
        gprMdldist = fitrgp(tbldist,'target','BasisFunction','linear','KernelFunction','ardMatern52','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsdist(1:size,qq) = predict(gprMdldist,data1);
        Validationdist(1:size,qq) = valdist;
        eval(['gprMdl',num2str(qq),' = gprMdldist;'])
        subplot(2,3,1)
        eval(['p',num2str(qq),' = scatter(valdist,Resultsdist(:,qq),''filled'');'])
        hold on
        xlabel('Xbeach (m)')
        ylabel('GPR (m)')          
        
        disp(['working on a parameter A'])
        design1(:,17) = targeta;
        tbl1 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),design1(:,17));%,design1(:,16),...
            %design1(:,17),design1(:,18));
        tbl1.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','MSL','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','tides','target'};
        gprMdla = fitrgp(tbl1,'target','BasisFunction','pureQuadratic','KernelFunction','ardMatern52','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsa(1:size,qq) = predict(gprMdla,data1);
        Validationa(1:size,qq) = vala;
        eval(['gprMdl',num2str(qq),' = gprMdl1;'])
        
        subplot(2,3,2)
        eval(['p',num2str(qq),' = scatter(vala,Resultsa(:,qq),''filled'');'])
        hold on
        xlabel('Xbeach A')
        ylabel('GPR A')
        
        disp(['working on a parameter B'])
        design1(:,17) = targetb;
        tbl1 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),design1(:,17));%,design1(:,16),...
            %design1(:,17),design1(:,18));
        tbl1.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','MSL','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','tides','target'};
        gprMdlb = fitrgp(tbl1,'target','BasisFunction','linear','KernelFunction','ardMatern52','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsb(1:size,qq) = predict(gprMdlb,data1);
        Validationb(1:size,qq) = valb;
        eval(['gprMdl',num2str(qq),' = gprMdl1;'])
        
        subplot(2,3,3)
        eval(['p',num2str(qq),' = scatter(valb,Resultsb(:,qq),''filled'');'])
        hold on
        xlabel('Xbeach B')
        ylabel('GPR B')
        
        
        
        c = c +size;
        
    end
end




subplot(2,3,1)
[r m2 b2] = regression(Resultsdist(:)',Validationdist(:)');
title(['xLoc R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
plot([100 260],[100 260],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
%ylim([-6 6])
%xlim([-6 6])
subplot(2,3,4)
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
xlabel('distance error')
title('Errors for all folds','color','k')


subplot(2,3,2)
[r m2 b2] = regression(Resultsa(:)',Validationa(:)');
title(['A R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
%plot([0.5 4.5],[0.5 4.5],'w--')
plot([-2.5 0.5],[-2.5 0.5],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
%ylim([-6 6])
%xlim([-6 6])
subplot(2,3,5)
errors = Resultsa-Validationa;
binrng = -0.75:.05:0.75;
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
xlabel('A error')
title('Errors for all folds','color','k')


subplot(2,3,3)
[r m2 b2] = regression(Resultsb(:)',Validationb(:)');
title(['B R^{2} = ',num2str(round(r*1000)/1000),''])% with K-fold of ',num2str(size),' with-held from ',num2str(total),''],'color','k')
plot([0.0 1.5],[0.0 1.5],'k--')
xlim([0.0,1.5])
ylim([0.0,1.5])
%plot([20 160],[20 160],'k--')
%lg = legend([p1 p2 p3 p4 p5],'1^{st} fold','2^{nd} fold','3^{rd} fold','4^{th} fold','5^{th} fold','location','northwest');
%ylim([-6 6])
%xlim([-6 6])
subplot(2,3,6)
errors = Resultsb-Validationb;
binrng = -0.75:.05:0.75;
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
xlabel('B error')
title('Errors for all folds','color','k')


    errorsdist = Resultsdist-Validationdist;
    squarederrorsdist = (errorsdist).^2;
    msedist = mean(squarederrorsdist);
    RMSEdist = sqrt(msedist);
    
errorsa = Resultsa-Validationa;
squarederrorsa = (errorsa).^2;
msea = mean(squarederrorsa);
RMSEa = sqrt(msea);

errorsb = Resultsb-Validationb;
squarederrorsb = (errorsb).^2;
mseb = mean(squarederrorsb);
RMSEb = sqrt(mseb);
%%


% sets = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,199]; %,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96];
% totals = [50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,995]; %,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480];

%sets = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250]; 
%totals = [50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250]; 
%sets = [20,40,60,80,100,120,140,160,180,200,220,240]; 
%totals = [100,200,300,400,500,600,700,800,900,1000,1100,1200]; 

sets = [25,50,75,100,125,150,175,200,225,234];
totals = [125,250,375,500,625,750,875,1000,1125,1234];
%targeting[50,100,200,
figure
set(gcf,'color','w');

for hh = 1:length(sets)
    size = sets(hh);
    total = totals(hh);
    design = all(reorder,:);
%     targets = predictands(reorder,:);
%    starts = predictors(reorder,:);
    
    distTargets = xLoc(reorder);
    aTargets = a(reorder);
    bTargets = b(reorder);

   
    newall = design(1:total,:);
    newdist = distTargets(1:total);
    newa = aTargets(1:total);
    newb = bTargets(1:total);
    
    c = 1;
    for qq = 1:5

        design1 = newall;
        design1(c:c+size-1,:) = [];

        targetdist = newdist;
        targetdist(c:c+size-1) = [];
        targeta = newa;
        targeta(c:c+size-1) = [];
        targetb = newb;
        targetb(c:c+size-1) = [];
        
        
        disp(['removing ',num2str(c),' to ',num2str(c+size-1),' from starting dataset'])
        data1 = newall(c:c+size-1,:);

        valdist = newdist(c:c+size-1);
        vala = newa(c:c+size-1);
        valb = newb(c:c+size-1);

        disp(['holding ',num2str(c),' to ',num2str(c+size-1),' aside for validation'])
        

        disp(['working on distance to Erosion'])
        
        design1(:,17) = targetdist;
        tbldist = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),design1(:,17));%,design1(:,16),...
            %design1(:,17),design1(:,18));
        tbldist.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','MSL','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','tides','target'};
        gprMdldist = fitrgp(tbldist,'target','BasisFunction','linear','KernelFunction','ardMatern52','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsdist(1:size,qq) = predict(gprMdldist,data1);
        Validationdist(1:size,qq) = valdist;
        eval(['gprMdl',num2str(qq),' = gprMdldist;'])

        disp(['working on a parameter A'])
        design1(:,17) = targeta;
        tbl1 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),design1(:,17));%,design1(:,16),...
            %design1(:,17),design1(:,18));
        tbl1.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','MSL','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','tides','target'};
        gprMdla = fitrgp(tbl1,'target','BasisFunction','pureQuadratic','KernelFunction','ardMatern52','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsa(1:size,qq) = predict(gprMdla,data1);
        Validationa(1:size,qq) = vala;
        eval(['gprMdl',num2str(qq),' = gprMdl1;'])
        

        
        disp(['working on a parameter B'])
        design1(:,17) = targetb;
        tbl1 = table(design1(:,1),design1(:,2),design1(:,3),design1(:,4),design1(:,5),design1(:,6),...
            design1(:,7),design1(:,8),design1(:,9),design1(:,10),design1(:,11),...
            design1(:,12),design1(:,13),design1(:,14),design1(:,15),design1(:,16),design1(:,17));%,design1(:,16),...
            %design1(:,17),design1(:,18));
        tbl1.Properties.VariableNames = {'Hs','Tp','Dir','NTR','Dur','MSL','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','tides','target'};
        gprMdlb = fitrgp(tbl1,'target','BasisFunction','linear','KernelFunction','ardMatern52','FitMethod','sr','PredictMethod','fic','Standardize',1);
        Resultsb(1:size,qq) = predict(gprMdlb,data1);
        Validationb(1:size,qq) = valb;
        eval(['gprMdl',num2str(qq),' = gprMdl1;'])
        

        
        
        c = c +size;
        
    end
    
    errorsdist = Resultsdist-Validationdist;
    squarederrorsdist = (errorsdist).^2;
    msedist = mean(squarederrorsdist);
    RMSEdist(1:5,hh) = sqrt(msedist);
    
    errorsa = Resultsa-Validationa;
    squarederrorsa = (errorsa).^2;
    msea = mean(squarederrorsa);
    RMSEa(1:5,hh) = sqrt(msea);    
    
    errorsb = Resultsb-Validationb;
    squarederrorsb = (errorsb).^2;
    mseb = mean(squarederrorsb);
    RMSEb(1:5,hh) = sqrt(mseb);
end






subplot(131)
boxplot(RMSEdist,(totals-sets))
hold on
lines = findobj(gcf,'type','line','Tag','Median');
set(lines,'Color','r')
ylabel('RMSE (m)')
xlabel('# used to calibrate')
title('Distance')
% ylim([0.25 1.3])

subplot(132)
boxplot(RMSEa,(totals-sets))
hold on
lines = findobj(gcf,'type','line','Tag','Median');
set(lines,'Color','r')
ylabel('RMSE')
xlabel('# used to calibrate')
title('A')

subplot(133)
boxplot(RMSEb,(totals-sets))
hold on
lines = findobj(gcf,'type','line','Tag','Median');
set(lines,'Color','r')
ylabel('RMSE')
xlabel('# used to calibrate')
title('B')