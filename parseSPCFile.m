function out=parseSPCFile(fullFileString)

    s=dir(fullFileString);
    numBytes=s.bytes;
    numRecords=numBytes/4-1;
    
    fileID=fopen(fullfile(s.folder,s.name),'r');
    binData=fread(fileID,'*ubit16');
    fclose(fileID);
    clear s fileID fileName folderName
    
    out.numADCbins=4096;
    
    %first frame contains info about macro time clock in 0.1ns (byte 0-2),
    %should have the invalid bit set (byte 3, bit 7) and a "number of routing
    %bits" (byte 3, bit 3-6)
    out.MTclock=double(bitand(typecast(binData(1:2),'uint32'),0x00ffffff))*0.1e-9;
    
    if (bitget(binData(2),16)==0) disp('unreliable spc file...'); end
    %unused, should encode number of detectors, also, might contain leacy bug
    numRoutingBits=bitshift(bitshift(binData(2),1),-12);
    if (numRoutingBits==4) numRoutingBits=2; end
    
    isInvalid=bitget(binData(4:2:end),16)==1;
    isMarker=bitget(binData(4:2:end),13)==1;
    channel=uint8(bitshift(binData(3:2:end),-12));
    adcVal=(out.numADCbins-1)-uint16(bitand(binData(4:2:end),0x0fff));
    
    %for photon records, channel is the detector used (0 for FLIM, 2 for
    %confocal, 1+3 unused); for marker records: 1-pixel, 2-line, 4-frame
    
    out.photons=zeros(numRecords-numel(isInvalid(isInvalid)),5,'uint16');
    out.photons(:,4)=channel(~isInvalid);
    out.photons(:,5)=adcVal(~isInvalid);
    
    startFrame=find(channel==4 & isInvalid & isMarker);
    startLine=(channel==2 & isInvalid & isMarker);
    startPixel=(channel==1 & isInvalid & isMarker);
    out.dimZ=uint16(numel(startFrame));
    lastMarker=numel(isMarker);
    if (numel(startFrame)>1); lastMarker=startFrame(2); end
    out.dimY=uint16(floor(numel(isMarker(startLine(startFrame(1):lastMarker)))));
    out.dimX=uint16(floor(numel(isMarker(startPixel(startFrame(1):lastMarker)))/out.dimY));
    
    fprintf('parsing')
    curX=-1; curY=-1; curZ=-1; wrongMarker=0; numPhotons=0;
    for i=1:numRecords
        if (isInvalid(i))
            if (isMarker(i))
                switch channel(i)
                    case 4
                        curZ=curZ+1;
                        curY=-1;
                        curX=-1;
                        fprintf(".");                        
                    case 2
                        curY=curY+1;
                        curX=0;
                    case 1
                        curX=curX+1;
                    otherwise
                        wrongMarker=wrongMarker+1;
                end
            end
        else
            numPhotons=numPhotons+1;
            out.photons(numPhotons,1:3)=[curX,curY,curZ];
        end
    end
    fprintf("\n");
    clear curX curY curZ channel adcVal binData numBytes startFrame
    cleanPhotons=out.photons(:,1)>-1 & out.photons(:,2)>-1 & out.photons(:,1)<out.dimX & out.photons(:,2)<out.dimY;
    ch1Photons=cleanPhotons & out.photons(:,4)==0;
    ch2Photons=cleanPhotons & out.photons(:,4)==2;
    %check if data contained 1 or 2 channels 

    out.ignoreCh2=numel(ch1Photons(ch1Photons))==0;
    if (out.ignoreCh2) 
        ch1Photons=ch2Photons;
        clear out.ch2Photons;
        if (numel(ch1Photons(ch1Photons))==0)
            ch1Photons=cleanPhotons;
        end
    else
        if (numel(ch2Photons(ch2Photons))==0)
            out.ignoreCh2=true;
        end
    end
    if (out.ignoreCh2)
        out.cleanPhotons=ch1Photons;
    else
        out.cleanPhotons=[ch1Photons,ch2Photons];
    end

    fprintf('found %i photons in ch1 and %i photons in ch2 and there were %i wrong markers\n',...
        numel(ch1Photons(ch1Photons)),numel(ch2Photons(ch2Photons)),wrongMarker);

    out.syncRate=1/out.MTclock;
    out.MTRes=out.MTclock/double(out.numADCbins);
end