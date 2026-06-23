function out=parsePTUFile(fullFileString)
    % format types in header file
    tyEmpty8      = 0xFFFF0008;
    tyBool8       = 0x00000008;
    tyInt8        = 0x10000008;
    tyBitSet64    = 0x11000008;
    tyColor8      = 0x12000008;
    tyFloat8      = 0x20000008;
    tyTDateTime   = 0x21000008;
    tyFloat8Array = 0x2001FFFF;
    tyAnsiString  = 0x4001FFFF;
    tyWideString  = 0x4002FFFF;
    tyBinaryBlob  = 0xFFFFFFFF;

    s=dir(fullFileString);
    fileID=fopen(fullfile(s.folder,s.name),'r');
    fName=fullfile(s.folder,s.name);
    fName=[fName(1:end-4),'.meta'];
    metaID=fopen(fName,'w');
    header=fread(fileID,6,'char')';
    if (~all(header=='PQTTTR'))
        fprintf('Houston... %s\n',header);
    end
    fread(fileID,10,'char');
    while true    
        tagName=fread(fileID,32,'char');
        tagName=char(tagName(tagName~=0x00)');
        tagIdx=fread(fileID,1,'int32');
        tagTypeCode=fread(fileID,1,'int32');
        if (tagIdx>-1)
            tagName=[tagName,int2str(tagIdx)];
        end

        switch tagTypeCode
            case tyEmpty8
                tagData=fread(fileID,1,'int64');
                fprintf(metaID,'%s: %i\n',tagName,tagData);
            case tyInt8
                tagData=fread(fileID,1,'int64');
                fprintf(metaID,'%s: %i\n',tagName,tagData);
            case tyBitSet64
                tagData=fread(fileID,1,'int64');
                fprintf(metaID,'%s: %i\n',tagName,tagData);
            case tyColor8
                tagData=fread(fileID,1,'int64');
                fprintf(metaID,'%s: %i\n',tagName,tagData);
            case tyBool8
                tagData=logical(fread(fileID,1,'int64'));
                fprintf(metaID,'%s: %i\n',tagName,tagData);
            case tyFloat8
                tagData=fread(fileID,1,'float64');
                fprintf(metaID,'%s: %g (0x%s)\n',tagName,tagData,num2hex(tagData));
            case tyTDateTime
                tagData=fread(fileID,1,'float64');
                fprintf(metaID,'%s: %g (0x%s)\n',tagName,tagData,num2hex(tagData));
            case tyFloat8Array
                tagData=fread(fileID,1,'int64');
                numVals=tagData/8;
                tagData=fread(fileID,numVals,'float64');
                for u=1:numVals
                    fprintf(metaID,'%s: %g (0x%s)\n',tagName,tagData,num2hex(tagData));
                end
            case tyBinaryBlob
                tagData=fread(fileID,1,'int64');
                numVals=tagData/8;
                tagData=fread(fileID,numVals,'int64');
                for u=1:numVals
                    fprintf(metaID,'%s: 0x%s\n',tagName,num2hex(tagData));
                end
            case tyAnsiString
                tagData=fread(fileID,1,'int64');
                tagData=fread(fileID,tagData,'char');
                fprintf(metaID,'%s: %s\n',tagName,tagData);
            case tyWideString
                tagData=fread(fileID,1,'int64');
                tagData=fread(fileID,tagData,'char');
                fprintf(metaID,'%s: %s\n',tagName,tagData);
            otherwise
                tagData=fread(fileID,8,'uint8');
        end
        
        switch tagName
            case 'TTResult_SyncRate'
                out.syncRate=tagData;
            case 'TTResult_NumberOfRecords'
                numRecs=tagData;
            case 'MeasDesc_Resolution'
                out.MTRes=tagData;
            case 'MeasDesc_GlobalResolution'
                out.MTclock=tagData;
            case 'ImgHdr_PixX'
                out.dimX=tagData;
            case 'ImgHdr_PixY'
                out.dimY=tagData;
            case 'ImgHdr_Frame'
                markerFrameStart=tagData;
            case 'ImgHdr_LineStart'
                markerLineStart=tagData;
            case 'ImgHdr_LineStop'
                markerLineStop=tagData;
            case 'ImgHdr_BiDirect'
                biDirScan=tagData>0;
            case 'Header_End'
                break
            otherwise
%                disp(tagName);
        end
    end
    out.numADCbins=ceil(out.MTclock/out.MTRes);
    
    binData=fread(fileID,'*ubit16');
    numReadRecs=numel(binData)/2;
    if (numReadRecs~=numRecs)
        fprintf('found %i records, diff to expected: %i\n',numReadRecs,numReadRecs-numRecs);
    end
    
    fclose(fileID); fclose(metaID);
    clear s fileID fileName folderName fName
    

    nSync=typecast(binData(1:2:end),'uint16');
    adcVal=typecast(bitand(binData(2:2:end),0x0fff),'uint16');
    chan=uint8(bitshift(binData(2:2:end),-12));
    markerMsk=chan==15;

    adcVal(markerMsk)=bitand(adcVal(markerMsk),0x000f);
    photonMsk=chan>=1 & chan<=4;
    numPhotons=numel(chan(photonMsk));
    out.photons=zeros(numPhotons,5,'uint16');
    out.photons(:,5)=adcVal(photonMsk);
    out.photons(:,4)=chan(photonMsk);
% mark identical photons; use 3rd column temporarily to park nSync;
% currently mark both for deletion
    out.photons(:,3)=nSync(photonMsk);
    out.cleanPhotons=true(numPhotons,1);
    out.cleanPhotons(2:end  )=any(out.photons(1:end-1,3:5)~=out.photons(2:end,3:5),2);
    out.cleanPhotons(1:end-1)=any(out.photons(1:end-1,3:5)~=out.photons(2:end,3:5),2);
% set value of frame (3rd column)
    frameIdx=find(chan==15 & bitand(bitshift(adcVal,1-markerFrameStart),1));
    frame=zeros(numReadRecs,1);
    u=1;
    if ~isempty(frameIdx)
        if numel(frameIdx)>1
            for u=1:numel(frameIdx)-1
                frame(frameIdx(u):frameIdx(u+1)-1)=u;
            end
            u=u+1;
        end
        frame(frameIdx(u):end)=u;
    end    
    out.photons(:,3)=frame(photonMsk);
    out.dimZ=u;
    
    syncPeriod=1/double(out.syncRate);
    nSync=syncPeriod*double(nSync);
% only photons have a valid adc value that should be addet to the nSync
% time (jediFLIM bug: nSync in ns, but adc in s!)
    nSync(photonMsk)=nSync(photonMsk)+double(adcVal(photonMsk))*out.MTRes;
% add nSync-overflow of 16-bit counter to all subsequent events
    syncPeriod=syncPeriod*2^16;
    overflowIdx=find(chan==15 & adcVal==0);
    for u=1:numel(overflowIdx)-1
        nSync(overflowIdx(u):overflowIdx(u+1)-1)= ...
            nSync(overflowIdx(u):overflowIdx(u+1)-1) ...
           + u*syncPeriod; 
    end
    u=u+1;
    nSync(overflowIdx(u):end)=nSync(overflowIdx(u):end)+u*syncPeriod;
    out.TrueTimeFrames=nSync(frameIdx);
    
    lineStartIdx=find(chan==15 & bitand(bitshift(adcVal,1-markerLineStart),1));
    lineStopIdx =find(chan==15 & bitand(bitshift(adcVal,1-markerLineStop ),1));
    if (numel(lineStartIdx)>numel(lineStopIdx))
        out.cleanPhotons(lineStartIdx(end):end)=false;
        lineStartIdx=lineStartIdx(1:numel(lineStopIdx));
    end
    lineTimes=1./(nSync(lineStopIdx(:))-nSync(lineStartIdx(:))).*double(out.dimX);
    xPos=-1*ones(numReadRecs,1);
    yPos=-1*ones(numReadRecs,1);
    for u=1:numel(lineStartIdx)
        if (biDirScan) oddLine=mod(u,2);
        else           oddLine=1;
        end
        if (oddLine==1)
            xPos(lineStartIdx(u):lineStopIdx(u)) = ...
                (nSync(lineStartIdx(u):lineStopIdx(u))-nSync(lineStartIdx(u)))*lineTimes(u);
        else
            xPos(lineStartIdx(u):lineStopIdx(u)) = out.dimX-1 - ...
                (nSync(lineStartIdx(u):lineStopIdx(u))-nSync(lineStartIdx(u)))*lineTimes(u);
        end
        yPos(lineStartIdx(u):lineStopIdx(u)) = u-1-frame(lineStartIdx(u):lineStopIdx(u))*out.dimY;
    end
    out.photons(:,1)=uint16(xPos(photonMsk));
    out.photons(:,2)=uint16(yPos(photonMsk));
    out.cleanPhotons(out.photons(:,1)<0)=false;
    out.cleanPhotons(out.photons(:,1)>=out.dimX)=false;
    out.cleanPhotons(out.photons(:,2)<0)=false;
    out.cleanPhotons(out.photons(:,2)>=out.dimY)=false;
    out.cleanPhotons(out.photons(:,5)>=out.numADCbins)=false;
    
    numPhotonsPerChannel=zeros(4,1);
    chMsk=false(numPhotons,4);
    for u=1:4
        chMsk(:,u)=out.photons(:,4)==u & out.cleanPhotons;
        numPhotonsPerChannel(u)=numel(out.cleanPhotons(chMsk(:,u)));
    end
    if (numel(numPhotonsPerChannel(numPhotonsPerChannel>0))==1)
        out.ignoreCh2=true;
        out.photons(:,4)=1;
    else        
        out.ignoreCh2=false;
        chIdx=[1,2,3,4];
        out.cleanPhotons=chMsk(:,chIdx(numPhotonsPerChannel>0));
    end
end