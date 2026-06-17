clear
mskFile='threeCellMsk_Stack.tif';
numImgs=numel(imfinfo(mskFile));
minCounts=1;
for u=numImgs:-1:1
    allMSK(:,:,u)=imread(mskFile,u);
end
numCells=max(allMSK,[],'all');
lowTRange=2:4;
highTRange=[1,5:numImgs];

folderName='/Users/schmick/Work/Microscopy.Data/Bijeesh_leicaSP8/Leica_EGFR_PTB_2604/defaultIRF_12N/';
matFileNames=dir(strcat(folderName,'0*mat'));

dcVal=zeros(size(allMSK));
rData=zeros(size(allMSK));
for u=1:numImgs
    load(strcat(folderName,matFileNames(u).name),'matrices');
    dcVal(:,:,u)=matrices.dcVals;
    rData(:,:,u)=matrices.rData;
    currMsk=allMSK(:,:,u);
    currMsk(matrices.dcVals<minCounts | isnan(matrices.rData))=0;
    allMSK(:,:,u)=currMsk;
end

%%
figure(1)
 set(gcf,'Position',[1,1,1700,500]);
 clf
 axes('Position',[0.025,0.03,0.25,0.9])
 imagesc(sum(allMSK,3)); axis image
 title('Overlay of all cell-masks')
 xticklabels([])
 yticklabels([])
 annotation('textbox',[0.145 0.46 .1 .1],'String','cell #1','FitBoxToText','on');
 annotation('textbox',[0.18 0.675 .1 .1],'String','cell #2','FitBoxToText','on');
 annotation('textbox',[0.075 0.145 .1 .1],'String','cell #3','FitBoxToText','on');

% redundant, sample for single experiment plot
% figure(2);
% clf;
% rRes=30;
% currCell=1;
% u=9;
%     currMsk=allMSK(:,:,u); currRData=rData(:,:,u); currDC=dcVal(:,:,u);
%     cellRData=currRData(currMsk==currCell);
%     cellDC=currDC(currMsk==currCell);
%     xCoords=real(cellRData); yCoords=imag(cellRData);
%     xRange=[min(xCoords),max(xCoords)]; deltaX=xRange(2)-xRange(1);
%     yRange=[min(yCoords),max(yCoords)]; deltaY=yRange(2)-yRange(1);
%     maxDelta=max(deltaX,deltaY);
%     dRes=rRes/maxDelta;
%     xBins=ceil(deltaX*dRes); yBins=ceil(deltaY*dRes);
%     xEdges=linspace(xRange(1),xRange(1)+xBins/dRes,xBins+1);
%     xCntrs=xEdges(1:end-1)+(xEdges(2)-xEdges(1))/2;
%     yEdges=linspace(yRange(1),yRange(1)+yBins/dRes,yBins+1);
%     yCntrs=yEdges(1:end-1)+(yEdges(2)-yEdges(1))/2;
%     cData=histcounts2(xCoords,yCoords,xEdges,yEdges);
%     alphaParam=fit(xCoords,yCoords,"poly1","Weights",cellDC);
% 
%     %imagesc(xCntrs,yCntrs,cData'); axis image
%     %set(gca,'YDir','normal')
%     scatter(real(cellRData),imag(cellRData),0.1,'b');
%     hold on
%         contour(xCntrs,yCntrs,cData','r');
%         plot(xCntrs,alphaParam.p1*xCntrs+alphaParam.p2,'k','LineWidth',2)
%     hold off


%% linear regression per temp for all cells
rRes=16;
allCellMSK=allMSK>0;
cellMSK37=allCellMSK;
cellMSK37(:,:,lowTRange)=false;

xCoords=real(rData(cellMSK37)); yCoords=imag(rData(cellMSK37));
xRange=[min(xCoords),max(xCoords)]; deltaX=xRange(2)-xRange(1);
yRange=[min(yCoords),max(yCoords)]; deltaY=yRange(2)-yRange(1);
maxDelta=max(deltaX,deltaY);
dRes=rRes/maxDelta;
xBins=ceil(deltaX*dRes); yBins=ceil(deltaY*dRes);
xEdges=linspace(xRange(1),xRange(1)+xBins/dRes,xBins+1);
yEdges=linspace(yRange(1),yRange(1)+yBins/dRes,yBins+1);

xCntrs37=xEdges(1:end-1)+(xEdges(2)-xEdges(1))/2;
yCntrs37=yEdges(1:end-1)+(yEdges(2)-yEdges(1))/2;
xCoords37=xCoords; yCoords37=yCoords;
cData37=histcounts2(xCoords,yCoords,xEdges,yEdges);
alphaParam37=fit(xCoords,yCoords,"poly1","Weights",dcVal(cellMSK37));

cellMSK05=allCellMSK;
cellMSK05(:,:,highTRange)=false;

xCoords=real(rData(cellMSK05)); yCoords=imag(rData(cellMSK05));
xRange=[min(xCoords),max(xCoords)]; deltaX=xRange(2)-xRange(1);
yRange=[min(yCoords),max(yCoords)]; deltaY=yRange(2)-yRange(1);
maxDelta=max(deltaX,deltaY);
dRes=rRes/maxDelta;
xBins=ceil(deltaX*dRes); yBins=ceil(deltaY*dRes);
xEdges=linspace(xRange(1),xRange(1)+xBins/dRes,xBins+1);
yEdges=linspace(yRange(1),yRange(1)+yBins/dRes,yBins+1);

xCntrs05=xEdges(1:end-1)+(xEdges(2)-xEdges(1))/2;
yCntrs05=yEdges(1:end-1)+(yEdges(2)-yEdges(1))/2;
xCoords05=xCoords; yCoords05=yCoords;
cData05=histcounts2(xCoords,yCoords,xEdges,yEdges);
alphaParam05=fit(xCoords,yCoords,"poly1","Weights",dcVal(cellMSK05));

function interSects=linHalfCircIntersect(m,y0)
    p=(0.5-m*y0)/(m^2+1);
    q=y0^2/(m^2+1);
    r=p^2-q;
    if (r>=0)
        interSects=[p-sqrt(r),p+sqrt(r)];
        interSects(2,:)=m*interSects+y0;
    else
        interSects=nan(2,2);
    end
end
alphaIntersects05=linHalfCircIntersect(alphaParam05.p1,alphaParam05.p2);
alphaIntersects37=linHalfCircIntersect(alphaParam37.p1,alphaParam37.p2);

%%
figure(1);
axes('Position',[0.325,0.1,0.3,0.8])
plot(0:0.002:1,sqrt(0.25-(-0.5:0.002:0.5).^2),'g','LineWidth',2);
hold on
    scatter(xCoords37,yCoords37,0.1,'r.')
    scatter(xCoords05,yCoords05,0.1,'b.')
    contour(xCntrs37,yCntrs37,cData37',5,'LineWidth',2);
    contour(xCntrs05,yCntrs05,cData05',5,'LineWidth',2);
    plot(xCntrs37,alphaParam37.p1*xCntrs37+alphaParam37.p2,'m','LineWidth',1.5);
    plot(xCntrs05,alphaParam05.p1*xCntrs05+alphaParam05.p2,'c','LineWidth',1.5);
    scatter(alphaIntersects05(1,:),alphaIntersects05(2,:),30,'k','filled');
    scatter(alphaIntersects37(1,:),alphaIntersects37(2,:),30,'k','filled');

    plot(0:0.002:1,sqrt(0.25-(-0.5:0.002:0.5).^2),'g--','LineWidth',2);
    hold off
    xlim([0.5,0.8]);
    ylim([0.375,0.575]);
    
    comVals=zeros(numImgs,numCells);
    cellCol=linspace(0.2,0.8,numCells);
    
    function tau=r2Tau(inCmplx)
        w=2.5634e-08/2/pi*1e9;
        normIn=(inCmplx-0.5)./abs(inCmplx-0.5)*0.5+0.5;
        tau=w*sqrt(1./abs(normIn).^2-1); %modulation lifetime
        %pLT=w*tan(angle(rD));        phase lifetime
    end
    
    twoTau05=[ r2Tau(complex(alphaIntersects05(1,1),alphaIntersects05(2,1))), ...
               r2Tau(complex(alphaIntersects05(1,2),alphaIntersects05(2,2)))];
    twoTau37=[ r2Tau(complex(alphaIntersects37(1,1),alphaIntersects37(2,1))), ...
               r2Tau(complex(alphaIntersects37(1,2),alphaIntersects37(2,2)))];

%%
hold on
for v=1:numCells
  for u=1:numImgs
    currMsk=allMSK(:,:,u);
    currDC=dcVal(:,:,u);
    currRData=rData(:,:,u);
    cellRData=currRData(currMsk==v);
    cellDC=currDC(currMsk==v);
    comVals(u,v)=sum(cellRData.*cellDC)/sum(cellDC);
  end
    xCoords=real(comVals(:,v));
    yCoords=imag(comVals(:,v));
    scatter(xCoords,yCoords,50,cellCol(v)*[1,1,1],'filled');
end
  scatter(real(comVals),imag(comVals),50,'white','LineWidth',1);
title(sprintf('phasor plot of half-circle-intersects (37°C:%1.3g/%1.3g; 5°C:%1.3g/%1.3g) [ns]',twoTau37(1),twoTau37(2),twoTau05(1),twoTau05(2)))
legend('','37°C','5°C','','','37°C','5°C','','','','cell #1','cell #2','cell #3','');
hold off

normTau=r2Tau(comVals);
 beta=zeros(size(normTau));
alpha=zeros(size(normTau));

 beta(highTRange,:)=(normTau(highTRange,:)-twoTau37(1))/(twoTau37(2)-twoTau37(1));
alpha(highTRange,:)=twoTau37(1)./(twoTau37(1)+twoTau37(2)./beta(highTRange,:)-twoTau37(2));
 beta(lowTRange,:)=(normTau(lowTRange,:)-twoTau05(1))/(twoTau05(2)-twoTau05(1));
alpha(lowTRange,:)=twoTau05(1)./(twoTau05(1)+twoTau05(2)./beta(lowTRange,:)-twoTau05(2));

timeVals=[-15,-10,-5,0,5,10,15,30,60];

figure(1);
axes('Position',[0.65,0.1,0.3,0.8])
    %plot(normTau,'*-');
    %yline(twoTau05,'Color','cyan');
    %yline(twoTau37,'Color','magenta');
    pH=plot(timeVals,alpha,'*-');
    for u=1:numCells; set(pH(u),'Color',cellCol(u)*ones(3,1)); end
    title('alpha pEGFR for individual cells')
    xticks(timeVals)
    xticklabels({'pre','5°C','EGF','wash','5min','10','15min','30min','60'})
    legend('cell #1', 'cell #2', 'cell #3')

%saveas(gcf,'Fig_pEGFR_FLIM.png','png')
%print('Fig_pEGFR_FLIM.svg','-dsvg','-vector')

%%
perimThick=5;
newMSK=zeros(size(allMSK),class(allMSK));
comValsPer=zeros(numImgs,numCells);

for u=1:numImgs
  currDC=dcVal(:,:,u);
  currRData=rData(:,:,u);
  for v=1:numCells
    currMsk=zeros(size(allMSK,1),size(allMSK,2),class(allMSK));
    currMsk(allMSK(:,:,u)==v)=v;
    newMsk=imdilate(currMsk,strel('diamond',1));
    newMsk=currMsk-imerode(newMsk,strel('diamond',perimThick+1));
    newMSK(:,:,u)=newMSK(:,:,u)+newMsk;
    cellRData=currRData(newMsk>0);
    cellDC=currDC(newMsk>0);
    comValsPer(u,v)=sum(cellRData.*cellDC)/sum(cellDC);
  end
end

figure(2);
  set(gcf,'Position',[1,1,1700,500]);
  clf
  axes('Position',[0.025,0.03,0.25,0.9])
  imagesc(sum(newMSK,3)); axis image
  title('Overlay of all cell-masks')
  xticklabels([])
  yticklabels([])
  annotation('textbox',[0.145 0.46 .1 .1],'String','cell #1','FitBoxToText','on');
  annotation('textbox',[0.18 0.675 .1 .1],'String','cell #2','FitBoxToText','on');
  annotation('textbox',[0.075 0.145 .1 .1],'String','cell #3','FitBoxToText','on');

  axes('Position',[0.325,0.1,0.3,0.8])
  plot(0:0.002:1,sqrt(0.25-(-0.5:0.002:0.5).^2),'g','LineWidth',2);
  hold on
    contour(xCntrs37,yCntrs37,cData37',5,'LineWidth',2);
    contour(xCntrs05,yCntrs05,cData05',5,'LineWidth',2);

    for v=1:numCells
      scatter(real(comValsPer(:,v)),imag(comValsPer(:,v)),50,cellCol(v)*[1,1,1],'filled');
    end
    scatter(real(comValsPer),imag(comValsPer),50,'white','LineWidth',1);
    xlim([0.5,0.8]);
    ylim([0.375,0.575]);
  hold off

normTauPer=r2Tau(comValsPer);
 betaPer=zeros(size(normTauPer));
alphaPer=zeros(size(normTauPer));

 betaPer(highTRange,:)=(normTauPer(highTRange,:)-twoTau37(1))/(twoTau37(2)-twoTau37(1));
alphaPer(highTRange,:)=twoTau37(1)./(twoTau37(1)+twoTau37(2)./betaPer(highTRange,:)-twoTau37(2));
 betaPer(lowTRange,:)=(normTauPer(lowTRange,:)-twoTau05(1))/(twoTau05(2)-twoTau05(1));
alphaPer(lowTRange,:)=twoTau05(1)./(twoTau05(1)+twoTau05(2)./betaPer(lowTRange,:)-twoTau05(2));

axes('Position',[0.65,0.1,0.3,0.8])
    pH=plot(alphaPer,'*-');
    for u=1:numCells; set(pH(u),'Color',cellCol(u)*ones(3,1)); end
    hold on
        pH=plot(alpha,'--');
        for u=1:numCells; set(pH(u),'Color',cellCol(u)*ones(3,1)); end
    hold off
    title('alpha pEGFR in periphery of individual cells')
    xticks(1:numImgs)
    xticklabels({'pre','5°C','EGF','wash','5min','10','15','30','60min'})
    legend('cell #1', 'cell #2', 'cell #3')

%saveas(gcf,'Fig_pEGFR_periphery_FLIM.png','png')

%% normalized rData -> two-tau-FRETFLIM
% This algorithm re-allocates all per-pixel rData complex phasor values
% into angular bins to reduce noise and impact of outliers for a more
% precise estimation of tau1 / tau2 (no-FRET, 100% FRET). Arbitrarily, 
% 256 angular bins were chosen.
photonsPerAngleBin=zeros(256,1);
allAngleBin=zeros(256,1);

for v=1:numCells
  for u=highTRange
    currMsk=allMSK(:,:,u);
    currDC=dcVal(:,:,u);
    currRData=rData(:,:,u);
    cellRData=currRData(currMsk==v);
    cellDC=currDC(currMsk==v);

    cellAngleBin=round((angle(cellRData-0.5)/2/pi+0.5)*255)+1;
    photonsPerAngleBin = photonsPerAngleBin + accumarray(cellAngleBin,cellDC,[256,1]);
    allAngleBin = allAngleBin + accumarray(cellAngleBin,cellRData.*cellDC,[256,1]);
  end
end
allAngleBin=allAngleBin./photonsPerAngleBin;
plotData=photonsPerAngleBin>0;
linReg=[photonsPerAngleBin(plotData).*ones(size(allAngleBin(plotData))), ...
        photonsPerAngleBin(plotData).*real(allAngleBin(plotData))]\...
       (photonsPerAngleBin(plotData).*imag(allAngleBin(plotData)));
alphaIntersects37_norm=linHalfCircIntersect(linReg(2),linReg(1));
twoTau37_norm=[ r2Tau(complex(alphaIntersects37_norm(1,1),alphaIntersects37_norm(2,1))), ...
                r2Tau(complex(alphaIntersects37_norm(1,2),alphaIntersects37_norm(2,2))) ];

tauData=r2Tau(rData);

 allBeta=zeros(size(tauData));
allAlpha=zeros(size(tauData));
 allBeta(:,:,highTRange)=(tauData(:,:,highTRange)-twoTau37_norm(1))/(twoTau37_norm(2)-twoTau37_norm(1));
allAlpha(:,:,highTRange)=twoTau37_norm(1)./(twoTau37_norm(1)+twoTau37_norm(2)./allBeta(:,:,highTRange)-twoTau37_norm(2));


 beta_norm=zeros(size(normTau));
alpha_norm=zeros(size(normTau));
 betaPer_norm=zeros(size(normTau));
alphaPer_norm=zeros(size(normTau));

 beta_norm(highTRange,:)=(normTau(highTRange,:)-twoTau37_norm(1))/(twoTau37_norm(2)-twoTau37_norm(1));
alpha_norm(highTRange,:)=twoTau37_norm(1)./(twoTau37_norm(1)+twoTau37_norm(2)./beta_norm(highTRange,:)-twoTau37_norm(2));

 betaPer_norm(highTRange,:)=(normTauPer(highTRange,:)-twoTau37_norm(1))/(twoTau37_norm(2)-twoTau37_norm(1));
alphaPer_norm(highTRange,:)=twoTau37_norm(1)./(twoTau37_norm(1)+twoTau37_norm(2)./betaPer_norm(highTRange,:)-twoTau37_norm(2));

currFig=[1,5,7,9];
for u=numel(currFig):-1:1
    currIm=dcVal(:,:,currFig(u))';
    currIm=currIm/max(currIm,[],"all");
    trnspIm(:,(u-1)*512+(1:512))=currIm;
    alphaIm(:,(u-1)*512+(1:512))=(double(allMSK(:,:,currFig(u))>0).*allAlpha(:,:,currFig(u)))';
end
colorLim=[0.25,0.6];
cntrstFactor=2.0;
cm=colormap("gray");
cm(:,1)=cm(end:-1:1,1);
rgbIm=(alphaIm-colorLim(1))/(colorLim(2)-colorLim(1));
rgbIm(rgbIm<0)=0;
rgbIm(rgbIm>1)=1;
rgbIm=cntrstFactor*ind2rgb(uint8(255*rgbIm+1),cm).*repmat(trnspIm,1,1,3);

figure(3)
  set(gcf,'Position',[1,1,1700,400]);
  clf
  aH=axes('Position',[0.025,0.1,0.75,0.9]);
  imshow(rgbIm);
  % imagesc(alphaIm,'AlphaData',2*trnspIm)
  % axis image
  clim(colorLim)
  % %cm=1-colormap("gray");
  % %cm(:,2)=cm(end:-1:1,2);
  % cm=colormap("gray");
  % cm(:,1)=cm(end:-1:1,1);
  colormap(cm);
  colorbar
  set(aH,'Color',[0,0,0])
  xticklabels([]);
  yticklabels([]);
  annotation('textbox',[0.09 0.625 .1 .1],'String','cell #1','FitBoxToText','on','Color','white','EdgeColor','white');
  annotation('textbox',[0.0375 0.55 .1 .1],'String','cell #2','FitBoxToText','on','Color','white','EdgeColor','white');
  annotation('textbox',[0.15 0.5 .1 .1],'String','cell #3','FitBoxToText','on','Color','white','EdgeColor','white');
  xline(512.5,'Color','white','LineWidth',3)
  xline(1024.5,'Color','white','LineWidth',3)
  xline(1536,'Color','white','LineWidth',3)
  annotation('textbox',[0.1,0.07,.1,.1],'String','pre','EdgeColor','none')
  annotation('textbox',[0.28,0.07,.1,.1],'String','5 min','EdgeColor','none')
  annotation('textbox',[0.46,0.07,.1,.1],'String','15 min','EdgeColor','none')
  annotation('textbox',[0.64,0.07,.1,.1],'String','60 min','EdgeColor','none')
  % axes('Position',[0.325,0.1,0.3,0.8])
  %   %contour(xCntrs37,yCntrs37,cData37',5,'LineWidth',2);
  %   imagesc(xCntrs37,yCntrs37,cData37'); axis image
  %   hold on
  %       contour(xCntrs37,yCntrs37,cData37',5,'LineWidth',2);
  %       plot(0:0.002:1,sqrt(0.25-(-0.5:0.002:0.5).^2),'g','LineWidth',2);
  %       %plot(xCntrs37,alphaParam37.p1*xCntrs37+alphaParam37.p2,'m','LineWidth',1.5);
  %       scatter(real(allAngleBin(plotData)),imag(allAngleBin(plotData)),photonsPerAngleBin(plotData)/1000,'filled')
  %       plot(xCntrs05,linReg(2)*xCntrs05+linReg(1),'k','LineWidth',1.5);
  %       scatter(alphaIntersects37_norm(1,:),alphaIntersects37_norm(2,:),30,'k','filled');
  %   hold off
  %   title(sprintf('phasor plot of half-circle-intersects (%1.3g/%1.3g) ns',twoTau37_norm(1),twoTau37_norm(2)))
  %   xlim([0.5,0.875]);
  %   ylim([0.3,0.6]);
  %   set(gca,'YDir','normal')

  timeVals=[0,5,10,15,30,60];
  axes('Position',[0.7825,0.17,0.165,0.76])
    %pH=plot(timeVals,alphaPer_norm(highTRange,:),'*-','LineWidth',2);
    %for u=1:numCells; set(pH(u),'Color',cellCol(u)*ones(3,1)); end
    %hold on
        pH=plot(timeVals,alpha_norm(highTRange,:),'*-','LineWidth',2);
        for u=1:numCells; set(pH(u),'Color',cellCol(u)*ones(3,1)); end
    %hold off
    %title('cumulative pEGFR-fraction of individual cells')
    xticks([0,5,10,15,30,60])
    xticklabels({'pre','5','10','15','30min','60min'})
    yticklabels([])
    ylim(colorLim)
    legend('cell #1', 'cell #2', 'cell #3','Location','northeast')
    ylabel('interacting fraction α','Rotation',270)
    set(gca,'YAxisLocation','right')
   
%saveas(gcf,strcat(folderName,'pEGFR_FLIM_green_mag.png'),'png')
saveas(gcf,strcat(folderName,'pEGFR_FLIM_cyan-red.png'),'png')
print(strcat(folderName,'pEGFR_FLIM_12N.svg'),'-dsvg','-vector')

%%
print(strcat(folderName,'pEGFR_FLIM_12N.eps'),'-depsc','-vector')
