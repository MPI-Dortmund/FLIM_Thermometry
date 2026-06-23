% data import from SPC-file
fileName="SampleFile_02.spc";
spcData=parseSPCFile(fileName);

% unfiltered histogram of all photons
figure(1);
    semilogy(histcounts(spcData.photons(:,5),'BinMethod','integers'));


% setting data ranges to determine background counts and for fitting
meanBGRange=(2000:3500);
fitRange=(1:3500);
% for spatially resolved data: potentially limit analyzed photons to ROI 
% in x-y-direction (column 1+2) or discard partially scanned frame (col 3)
% cleanPhotons can be a meta-data flag of untrustworthy photons that is not
% utilized in SPC data format
xRange=[min(spcData.photons(:,1)),max(spcData.photons(:,1))];
yRange=[min(spcData.photons(:,2)),max(spcData.photons(:,2))];
fRange=[min(spcData.photons(:,3)),max(spcData.photons(:,3))];
% manually limit to first frame only
fRange=[0,0];

filteredPhotons=( spcData.cleanPhotons & ...
                  spcData.photons(:,1)>=xRange(1) & ...
                  spcData.photons(:,1)<=xRange(2) & ...
                  spcData.photons(:,2)>=yRange(1) & ...
                  spcData.photons(:,2)<=yRange(2) & ...
                  spcData.photons(:,3)>=fRange(1) & ...
                  spcData.photons(:,3)<=fRange(2) );

ArrTime=spcData.photons(:,5)+1;
ArrTime=ArrTime(filteredPhotons);

binArrTime=histcounts(ArrTime,'BinMethod','integers');
binTau=(0:numel(binArrTime)-1)*spcData.MTRes*1e9;

MTclock=spcData.MTclock;
MTRes=spcData.MTRes;
A=numel(ArrTime);

figure(2);
    plot(binTau,binArrTime,'k.');
    xlabel('arrival time [ns]')
    ylabel('# counts')
    xlim([2.3 5])

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

plot(binTau,yFitVals_2,'r-','LineWidth',2);
hold off
