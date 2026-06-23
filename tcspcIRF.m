% Taken from Walther, et. al "Precise measurement of protein interacting 
% fractions with fluorescence lifetime imaging microscopy",  Mol. BioSyst.
% 2011, 7, 322-336; DOI: 10.1039/c0mb00132e

function out=tcspcIRF(p, t)
% calculates the normalized IRF based on a gaussian at t+tg with asymmetric
% flanks by an exponential rise (left flank) + fall (right) at t+tg+dt
% input parameters: t     -- time, e.g. the detection bin
%                   tg    -- center of gaussian
%                   sigma -- width of gaussian
%                   dt    -- shift of asymmetric flanks
%                   rho   -- rel. amplitude of the asymmetric flanks
%                   th1   -- decay constant
%                   th2   -- rise constant
    sizeT=numel(t);
    numP=numel(p);
    if (numP<7)
        A=1;
    else
        A=p(7);
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
        dt=0;
    else
        dt=p(3);
    end
    if (numP<2)
        sigma=0.5;
    else
        sigma=p(2);
    end
    if (numP<1)
        tg=1;
    else
        tg=p(1);
    end
    
    r2f=th2./th1; P=(1+r2f); r2f=(P./r2f).^(r2f); P=P.*r2f;
    Q=A*sqrt(2*pi)*sigma+rho.*th1.*r2f;
    Q=Q/(t(end)-t(1))*numel(t);

    %out=ones(size(t));
    out=rho.*P.*exp((tg+dt-t)./th1);
    out=out.*(1-exp((tg+dt-t)./th2));
    msk=t<tg+dt;
    out(msk)=0;
    out=out+A*exp(-(t-tg).^2/2/sigma/sigma);
    out=out/Q;
end