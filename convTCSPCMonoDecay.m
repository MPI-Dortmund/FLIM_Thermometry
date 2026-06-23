% Taken from Walther, et. al "Precise measurement of protein interacting 
% fractions with fluorescence lifetime imaging microscopy",  Mol. BioSyst.
% 2011, 7, 322-336; DOI: 10.1039/c0mb00132e

function out=convTCSPCMonoDecay(p, t)
% calculates convolution of 6-param IRF model with asymmetric flanks with a
% mono-exponential decay with area normalized to 1 and no background counts
% input parameters: t     -- time, e.g. the detection bin
%                   tg    -- center of gaussian
%                   sigma -- width of gaussian
%                   dt    -- shift of asymmetric flanks
%                   rho   -- rel. amplitude of the asymmetric flanks
%                   th1   -- decay constant
%                   th2   -- rise constant
%                   tau   -- fluorescence decay
    sizeT=numel(t);
    numP=numel(p);
    if (numP<7)
        tau=0.75;
    else
        tau=p(7);
    end
    if (numP<6)
        th2=1;
    else
        th2=p(6);
    end
    if (numP<5)
        th1=0.5;
    else
        th1=p(5);
    end
    if (numP<4)
        rho=0.1;
    else
        rho=p(4);
    end
    if (numP<3)
        dt=0.5;
    else
        dt=p(3);
    end
    if (numP<2)
        sigma=0.2;
    else
        sigma=p(2);
    end
    if (numP<1)
        tg=2;
    else
        tg=p(1);
    end
    
    r2f=th2./th1; P=(1+r2f); r2f=(P./r2f).^(r2f); P=P.*r2f;
    Q=sqrt(2*pi)*sigma+rho.*th1.*r2f;
    r2f=th1.*th2./(th1+th2);
    Q=Q/(t(end)-t(1))*numel(t);
    t0=tg-t;
    sByTau=sigma./tau;

    %out=ones(size(t));
    out=tau./(th2.*(tau./r2f-1).*(tau./th1-1)).*exp((t0+dt)./th1);
    out=out+1./(tau./r2f-1).*exp((t0+dt)./r2f)-1./(tau./th1-1).*exp((t0+dt)./th1);
    msk=t<tg+dt;
    out(msk)=0;
    out=rho.*P.*out;
    out=out+sqrt(0.5*pi).*sByTau.*exp(0.5*sByTau.*sByTau)...
          .*erfc(sqrt(0.5)*(sByTau+t0./sigma))...
          .*exp(t0./tau);
    out=out./Q;
    out=out./sum(out(:));
end