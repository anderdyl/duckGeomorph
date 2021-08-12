


clear

load('latestNourishmentFits.mat','gprMdldist','gprMdla','gprMdlb','designParams')


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
maxTide=max(designParams(:,20));  minTide=min(designParams(:,20));


load('init2019PCs.mat')

%'Hs','Tp','Dir','NTR','Dur','MSL','EOF1','EOF2','EOF3','EOF4','EOF5','EOF6','EOF7','EOF8','EOF9','tide',

%Storm1 = [2.885, 10.417, 41.8, 0.6, 72, 0];
Storm1 = [2.46, 6.84, -15, 0.6, 80, 0];

Storm1N(1) =(Storm1(1)-minHs)./(maxHs-minHs);
Storm1N(2)=(Storm1(2)-minTp)./(maxTp-minTp);
Storm1N(3)=((Storm1(3)*pi/180)-minDir)./(maxDir-minDir);
Storm1N(4)=(Storm1(4)-minNTR)./(maxNTR-minNTR);
Storm1N(5)=(Storm1(5)-minDur)./(maxDur-minDur);
Storm1N(6)=(Storm1(6)-minMSL)./(maxMSL-minMSL);
    
subset_n(:,1)=(initPCs(:,1)-minEOF1)./(maxEOF1-minEOF1);
subset_n(:,2)=(initPCs(:,2)-minEOF2)./(maxEOF2-minEOF2);
subset_n(:,3)=(initPCs(:,3)-minEOF3)./(maxEOF3-minEOF3);
subset_n(:,4)=(initPCs(:,4)-minEOF4)./(maxEOF4-minEOF4);
subset_n(:,5)=(initPCs(:,5)-minEOF5)./(maxEOF5-minEOF5);
subset_n(:,6)=(initPCs(:,6)-minEOF6)./(maxEOF6-minEOF6);
subset_n(:,7)=(initPCs(:,7)-minEOF7)./(maxEOF7-minEOF7);
subset_n(:,8)=(initPCs(:,8)-minEOF8)./(maxEOF8-minEOF8);
subset_n(:,9)=(initPCs(:,9)-minEOF9)./(maxEOF9-minEOF9);


Tide1 = [-0.36]; %[

for hhh = 1:1600

    

    dist1(hhh) = predict(gprMdldist,[Storm1N,subset_n(hhh,:),(Tide1-minTide)./(maxTide-minTide)]);
    a1(hhh) = predict(gprMdla,[Storm1N,subset_n(hhh,:),(Tide1-minTide)./(maxTide-minTide)]);
    b1(hhh) = predict(gprMdlb,[Storm1N,subset_n(hhh,:),(Tide1-minTide)./(maxTide-minTide)]);

end

save('storm1Fits.mat','dist1','a1','b1')




load('storm1PCs.mat')

%Storm2 = [6.289, 14.035, 82.6, 1.02, 72, 0];
Storm2 = [6.31, 10.0, 4.6, 1.02, 66, -0.05];

Tide2 = [0.25];
Storm2N(1) =(Storm2(1)-minHs)./(maxHs-minHs);
Storm2N(2)=(Storm2(2)-minTp)./(maxTp-minTp);
Storm2N(3)=((Storm2(3)*pi/180)-minDir)./(maxDir-minDir);
Storm2N(4)=(Storm2(4)-minNTR)./(maxNTR-minNTR);
Storm2N(5)=(Storm2(5)-minDur)./(maxDur-minDur);
Storm2N(6)=(Storm2(6)-minMSL)./(maxMSL-minMSL);
subset_n2(:,1)=(storm1PCs(:,1)-minEOF1)./(maxEOF1-minEOF1);
subset_n2(:,2)=(storm1PCs(:,2)-minEOF2)./(maxEOF2-minEOF2);
subset_n2(:,3)=(storm1PCs(:,3)-minEOF3)./(maxEOF3-minEOF3);
subset_n2(:,4)=(storm1PCs(:,4)-minEOF4)./(maxEOF4-minEOF4);
subset_n2(:,5)=(storm1PCs(:,5)-minEOF5)./(maxEOF5-minEOF5);
subset_n2(:,6)=(storm1PCs(:,6)-minEOF6)./(maxEOF6-minEOF6);
subset_n2(:,7)=(storm1PCs(:,7)-minEOF7)./(maxEOF7-minEOF7);
subset_n2(:,8)=(storm1PCs(:,8)-minEOF8)./(maxEOF8-minEOF8);
subset_n2(:,9)=(storm1PCs(:,9)-minEOF9)./(maxEOF9-minEOF9);


for hhh = 1:1600

    dist2(hhh) = predict(gprMdldist,[Storm2N,subset_n2(hhh,:),(Tide2-minTide)./(maxTide-minTide)]);
    a2(hhh) = predict(gprMdla,[Storm2N,subset_n2(hhh,:),(Tide2-minTide)./(maxTide-minTide)]);
    b2(hhh) = predict(gprMdlb,[Storm2N,subset_n2(hhh,:),(Tide2-minTide)./(maxTide-minTide)]);

end

save('storm2Fits.mat','dist2','a2','b2')






load('storm2PCs.mat')
%Storm3 = [4.073, 16.26, 62.1, 0.507, 108, 0];
Storm3 = [3.34, 9.17, -5.4, 0.507, 108-12, 0];

Tide3 = [-0.46];
Storm3N(1) =(Storm3(1)-minHs)./(maxHs-minHs);
Storm3N(2)=(Storm3(2)-minTp)./(maxTp-minTp);
Storm3N(3)=((Storm3(3)*pi/180)-minDir)./(maxDir-minDir);
Storm3N(4)=(Storm3(4)-minNTR)./(maxNTR-minNTR);
Storm3N(5)=(Storm3(5)-minDur)./(maxDur-minDur);
Storm3N(6)=(Storm3(6)-minMSL)./(maxMSL-minMSL);
subset_n3(:,1)=(storm2PCs(:,1)-minEOF1)./(maxEOF1-minEOF1);
subset_n3(:,2)=(storm2PCs(:,2)-minEOF2)./(maxEOF2-minEOF2);
subset_n3(:,3)=(storm2PCs(:,3)-minEOF3)./(maxEOF3-minEOF3);
subset_n3(:,4)=(storm2PCs(:,4)-minEOF4)./(maxEOF4-minEOF4);
subset_n3(:,5)=(storm2PCs(:,5)-minEOF5)./(maxEOF5-minEOF5);
subset_n3(:,6)=(storm2PCs(:,6)-minEOF6)./(maxEOF6-minEOF6);
subset_n3(:,7)=(storm2PCs(:,7)-minEOF7)./(maxEOF7-minEOF7);
subset_n3(:,8)=(storm2PCs(:,8)-minEOF8)./(maxEOF8-minEOF8);
subset_n3(:,9)=(storm2PCs(:,9)-minEOF9)./(maxEOF9-minEOF9);


for hhh = 1:1600

    dist3(hhh) = predict(gprMdldist,[Storm3N,subset_n3(hhh,:),(Tide3-minTide)./(maxTide-minTide)]);
    a3(hhh) = predict(gprMdla,[Storm3N,subset_n3(hhh,:),(Tide3-minTide)./(maxTide-minTide)]);
    b3(hhh) = predict(gprMdlb,[Storm3N,subset_n3(hhh,:),(Tide3-minTide)./(maxTide-minTide)]);

end

save('storm3Fits.mat','dist3','a3','b3')








load('storm3PCs.mat')
%Storm4 = [4.88, 15.209, 53.98, 0.77, 132, 0];
Storm4 = [4.02, 13.7, -19, 0.77, 132-12, 0];


Tide4 = [0.33];
Storm4N(1) =(Storm3(1)-minHs)./(maxHs-minHs);
Storm4N(2)=(Storm3(2)-minTp)./(maxTp-minTp);
Storm4N(3)=((Storm3(3)*pi/180)-minDir)./(maxDir-minDir);
Storm4N(4)=(Storm3(4)-minNTR)./(maxNTR-minNTR);
Storm4N(5)=(Storm3(5)-minDur)./(maxDur-minDur);
Storm4N(6)=(Storm3(6)-minMSL)./(maxMSL-minMSL);
subset_n4(:,1)=(storm2PCs(:,1)-minEOF1)./(maxEOF1-minEOF1);
subset_n4(:,2)=(storm2PCs(:,2)-minEOF2)./(maxEOF2-minEOF2);
subset_n4(:,3)=(storm2PCs(:,3)-minEOF3)./(maxEOF3-minEOF3);
subset_n4(:,4)=(storm2PCs(:,4)-minEOF4)./(maxEOF4-minEOF4);
subset_n4(:,5)=(storm2PCs(:,5)-minEOF5)./(maxEOF5-minEOF5);
subset_n4(:,6)=(storm2PCs(:,6)-minEOF6)./(maxEOF6-minEOF6);
subset_n4(:,7)=(storm2PCs(:,7)-minEOF7)./(maxEOF7-minEOF7);
subset_n4(:,8)=(storm2PCs(:,8)-minEOF8)./(maxEOF8-minEOF8);
subset_n4(:,9)=(storm2PCs(:,9)-minEOF9)./(maxEOF9-minEOF9);


for hhh = 1:1600

    dist4(hhh) = predict(gprMdldist,[Storm4N,subset_n4(hhh,:),(Tide4-minTide)./(maxTide-minTide)]);
    a4(hhh) = predict(gprMdla,[Storm4N,subset_n4(hhh,:),(Tide4-minTide)./(maxTide-minTide)]);
    b4(hhh) = predict(gprMdlb,[Storm4N,subset_n4(hhh,:),(Tide4-minTide)./(maxTide-minTide)]);

end

save('storm4Fits.mat','dist4','a4','b4')



