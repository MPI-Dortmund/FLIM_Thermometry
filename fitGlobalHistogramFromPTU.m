% data import from SPC-file
fileName="SampleFile_01.ptu";
ptuData=parsePTUFile(fileName);

% unfiltered histogram of all photons
figure(1);
    semilogy(histcounts(ptuData.photons(:,5),'BinMethod','integers'));


% setting data ranges to determine background counts and for fitting
meanBGRange=(225:255);
fitRange=(3:200);
% for spatially resolved data: potentially limit analyzed photons to ROI 
% in x-y-direction (column 1+2) or discard partially scanned frame (col 3)
% cleanPhotons can be a meta-data flag of untrustworthy photons that is not
% utilized in SPC data format
xRange=[min(ptuData.photons(:,1)),max(ptuData.photons(:,1))];
yRange=[min(ptuData.photons(:,2)),max(ptuData.photons(:,2))];
fRange=[min(ptuData.photons(:,3)),max(ptuData.photons(:,3))];

filteredPhotons=( ptuData.cleanPhotons & ...
                  ptuData.photons(:,1)>=xRange(1) & ...
                  ptuData.photons(:,1)<=xRange(2) & ...
                  ptuData.photons(:,2)>=yRange(1) & ...
                  ptuData.photons(:,2)<=yRange(2) & ...
                  ptuData.photons(:,3)>=fRange(1) & ...
                  ptuData.photons(:,3)<=fRange(2) );

ArrTime=ptuData.photons(:,5)+1;
ArrTime=ArrTime(filteredPhotons);

binArrTime=histcounts(ArrTime,'BinMethod','integers');
binTau=(0:numel(binArrTime)-1)*ptuData.MTRes*1e9;

MTclock=ptuData.MTclock;
MTRes=ptuData.MTRes;
A=numel(ArrTime);

figure(2);
    plot(binTau,binArrTime,'k.');
    xlabel('arrival time [ns]')
    ylabel('# counts')
    %xlim([2.3 5])

% non-linear fit of histogram with Gaussian IRF (e.g. in Fig.1a)

bgCounts=mean(binArrTime(meanBGRange));
% 4-parameter generic model of a Gaussian IRF at p2 with width p3 convolved
% with a mono-exponential decay with decay time p1 and amplitude p4
% initial values for IRF t0 (p2) and width (p3) are extracted from
% histogram using 2ns as an initial mono-exponential decay
fitModel = @(p,xdata) bgCounts+p(4)*p(3)/p(1)*sqrt(pi/2)*exp(0.5*(p(3)/p(1))^2-(xdata-p(2))/p(1)).*erfc(sqrt(0.5)*(p(3)/p(1)-(xdata-p(2))/p(3)));

% initial guess for  tau
tau=2;
maxT=binTau(end);
t0=(max(binArrTime)+min(binArrTime))/2;
t0=find(binArrTime>t0,1);
t0=binTau(t0);
[~,sigma]=max(binArrTime);
if (isscalar(sigma))
    sigma=binTau(sigma)-t0;
else
    sigma=maxT/100;
end

lowerBound=[   0,   0,    0,0];
upperBound=[maxT,maxT, maxT,A];
initVals=  [ tau,  t0,sigma,A/2/exp(1)];

fitVals=lsqcurvefit(fitModel,initVals,binTau(fitRange),binArrTime(fitRange),lowerBound,upperBound);
yFitVals=fitModel(fitVals,binTau);
fprintf("Fitted tau=%gns\n",fitVals(1));

hold on
    plot(binTau,yFitVals,'g-','LineWidth',2);

% non-linear fit of histogram with asymmetric IRF (e.g. in Sup.Fig.5b)

A=A-bgCounts;
corrTCSPC = @(p,xdata) A*convTCSPCMonoDecay(p,xdata)+bgCounts;

lowerBound=[   0,    0,    0,  0, 0.05,    0,0];
upperBound=[maxT, maxT,    1,  1,    1, maxT,maxT];
initVals=  [  t0,sigma,sigma,0.5,sigma,sigma,3*sigma];

optimOptions = optimset('Display', 'off', 'TolFun', 1e-7, 'TolX', 1e-5,'MaxFunEvals',1e4,'MaxIter',1e4);
fitVals_2=lsqcurvefit(corrTCSPC,initVals,binTau,binArrTime,lowerBound,upperBound,optimOptions);
yFitVals_2=corrTCSPC(fitVals_2,binTau);
fprintf("Fitted tau=%gns\n",fitVals_2(end));

plot(binTau,yFitVals_2,'r-','LineWidth',2);
plot(binTau,bgCounts*ones(size(binTau)),'b--')
plot(binTau,0.1*A*tcspcIRF(fitVals_2(1:6),binTau),'r--')
plot(binTau,0.1*A*exp(-((binTau-fitVals(2))/fitVals(3)).^2),'g--')
hold off
ylim([1,0.1*A])
set(gca,'YScale','log');

dimDC=max(ptuData.photons(filteredPhotons,1:2))+1;
dcIm=accumarray(ptuData.photons(filteredPhotons,1:2)+1,...
                1,dimDC);
meanLTIm=accumarray(ptuData.photons(filteredPhotons,1:2)+1,...
                    binTau(ptuData.photons(filteredPhotons,5)+1),dimDC,@median)-t0;
figure(3)
  subplot(1,2,1)
    imagesc(dcIm);
    colorbar
    axis image
    colormap(gca,"gray")
    title('number of photons per pixel')
  subplot(1,2,2)
    imagesc(meanLTIm);
    colorbar
    axis image
    colormap(gca,"hot")
    title('mean photon-arrival time after IRF-pulse')