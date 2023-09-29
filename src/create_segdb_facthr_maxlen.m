%--------------------------------
% function create_segdb_facthr_maxlen(strthr,facthr,seglen,fdrfile,facfile,segfile,matfile)
%--------------------------------
% This function takes an input binary stream network raster and partitions
%  the downstream stream network.  Using a user-provided flow accumulation
%  threshold ('facthr'), the network is first partitioned according to its
%  major confluences.  Each stream segment is subsequently partitioned into
%  subsegments of length <= 'seglen', using the maximum FAC jumps wihin a
%  segment to subdivide segments that are too long.
%
%  The subsegments are assigned a unique ID, which is written out to a BIL
%  file.  Also, a database table is created and stored in a MAT file, with
%  one record for each subsegment (see 'header' variable in the MAT file
%  for attribute labels).
%
% INPUT
%  strthr = flow accumulation (catchment) size used to define stream network
%  facthr = flow accumulation jump size used to define initial, natural breaks
%         = [] if no breaks are to be defined by FAC jumps
%  seglen = maximum segment length (in number of pixels)
%         = [] if no breaks are to be defined by segment length
%  fdrfile = path and filename for input flow direction BIL file ('uint8')
%  facfile = path and filename for input flow accumulation BIL file ('uint32')
%  segfile = path and filename for output segment ID BIL file (format will
%    be either 'uint8', 'uint16', or 'uint32', depending on number of segments)
%  matfile = path and filename where output 'seg_info' MATLAB table is to be stored
%--------------------------------
function create_segdb_facthr_maxlen(strthr,facthr,seglen,fdrf,facf,segf,matf)
%--------------------------------
%--------------------------------
%   3-by-3 numbering for analysis (determined by BIL ordering):
%           1		2 		3
%           4		x		5
%           6		7		8
%--------------------------------
%   3-by-3 numbering for flow direction (used in ArcHydroTools):
%           32		64 	128
%           16		x		1
%           8		4		2
%--------------------------------

%=============================
% Read FDR header & open FDR for reading
%-------------------------
fdrinfo = readbilheader(fdrf);
fdrfmt = fdrinfo.fmt;
fdrnumbyt = fdrinfo.nbyts;
if(fdrinfo.bytord)
   fdr = fopen(fdrf,'rb','b');
else
   fdr = fopen(fdrf,'rb');
end;
fdrfmt = fdrinfo.fmt;
r = fdrinfo.r;
c = fdrinfo.c;
pxsz = fdrinfo.pxszx;
ulx = fdrinfo.ulx;
uly = fdrinfo.uly;
%=============================

%=============================
% Read FAC header & open FAC for reading
%-------------------------
facinfo = readbilheader(facf);
facfmt = facinfo.fmt;
facnumbyt = facinfo.nbyts;
if(facinfo.bytord)
   fac = fopen(facf,'rb','b');
else
   fac = fopen(facf,'rb');
end;
%=============================

%=============================
% If break parameters are empty, set beyond UB value (so no breaks)
%=============================
if(isempty(facthr))
   facthr = r*c+1;
end;
if(isempty(seglen))
   seglen = r*c+1;
end;
%=============================

%=============================
% Step 0: count the number of stream pixels & read FAC value
%-----------------------------
if(~exist(matf,'file'))
   numpx = 0;
   strpx0 = zeros(round(r*c*.05),2); % pxid, fac value
   for j=1:r
      str0 = fread(fac,c,facfmt);
      f = find(str0>=strthr);
      len = length(f);
      if(len)
         strpx0(numpx+(1:len),:) = [(j-1)*c+f str0(f)];
         numpx = numpx+len;
      end;
   end;
   strpx0 = strpx0(1:numpx,:);
   save(matf,'strpx0');
else
   load(matf,'strpx0');
   numpx = size(strpx0,1);
end;
disp(sprintf('Stream network contains %d pixels',numpx));
%=============================

%=============================
% Step 1: identify headwater and confluence pixels
%-----------------------------
outflw = [32 64 128 16 1 8 4 2]'; % outflow
inflw = flipud(outflw); % inflow
idx0 = [-c+(-1:1)'; (-1:1)'; c+(-1:1)'];
wh = whos('-file',matf);
if(~ismember('strhd',{wh.name}))
   strhd = zeros(round(numpx*.01),1);
   strcf = strhd;
   hdct = 0; cfct = 0;
   for j = 1:numpx
      md = mod(strpx0(j,1),c);
      % check if pixel is not on image boundary
      if(~or(or(strpx0(j,1)<=c,strpx0(j,1)>(r-1)*c),ismember(md,[0 1])))
         buf0 = zeros(9,2);
         buf0(ismember(strpx0(j,1)+idx0,strpx0(:,1)),1) = 1;

         fseek(fdr,fdrnumbyt*(strpx0(j,1)-c-2),-1);
         buf0(1:3,2) = fread(fdr,3,fdrfmt);
         fseek(fdr,fdrnumbyt*(strpx0(j,1)-2),-1);
         buf0(4:6,2) = fread(fdr,3,fdrfmt);
         fseek(fdr,fdrnumbyt*(strpx0(j,1)+c-2),-1);
         buf0(7:9,2) = fread(fdr,3,fdrfmt);

         f = find(and(buf0([1:4,6:9],1),buf0([1:4,6:9],2)==inflw));
      elseif(strpx0(j,1)<=c) % on top boundary
         if(~ismember(strpx0(j,1),[1 c])) % not at a corner pixel
            buf0 = zeros(6,2);
            buf0(ismember(strpx0(j,1)+idx0(4:9),strpx0(:,1)),1) = 1;
            
            fseek(fdr,fdrnumbyt*(strpx0(j,1)-2),-1);
            buf0(1:3,2) = fread(fdr,3,fdrfmt);
            fseek(fdr,fdrnumbyt*(strpx0(j,1)+c-2),-1);
            buf0(4:6,2) = fread(fdr,3,fdrfmt);

            f = find(and(buf0([1,3:6],1),buf0([1,3:6],2)==inflw(4:8)));      
         elseif(strpx0==1) % UL pixel
            buf0 = zeros(4,2);
            buf0(ismember(strpx0(j,1)+idx0([5,6,8,9]),strpx0(:,1)),1) = 1;

            fseek(fdr,fdrnumbyt*(strpx0(j,1)-1),-1);
            buf0(1:2,2) = fread(fdr,2,fdrfmt);
            fseek(fdr,fdrnumbyt*(strpx0(j,1)+c-1),-1);
            buf0(3:4,2) = fread(fdr,2,fdrfmt);

            f = find(and(buf0(2:4,1),buf0(2:4,2)==inflw([5,7:8])));
         elseif(strpx0==c) % UR pixel
            buf0 = zeros(4,2);
            buf0(ismember(strpx0(j,1)+idx0([4,5,7,8]),strpx0(:,1)),1) = 1;

            fseek(fdr,fdrnumbyt*(strpx0(j,1)-2),-1);
            buf0(1:2,2) = fread(fdr,2,fdrfmt);
            fseek(fdr,fdrnumbyt*(strpx0(j,1)+c-2),-1);
            buf0(3:4,2) = fread(fdr,2,fdrfmt);

            f = find(and(buf0([1,3:4],1),buf0([1,3:4],2)==inflw([4,6:7])));
         end;
      elseif(strpx0(j,1)>(r-1)*c) % on bottom boundary
         if(~ismember(strpx0(j,1),[(r-1)*c+1 r*c])) % not at a corner pixel
            buf0 = zeros(6,2);
            buf0(ismember(strpx0(j,1)+idx0(1:6),strpx0(:,1)),1) = 1;

            fseek(fdr,fdrnumbyt*(strpx0(j,1)-c-2),-1);
            buf0(1:3,2) = fread(fdr,3,fdrfmt);
            fseek(fdr,fdrnumbyt*(strpx0(j,1)-2),-1);
            buf0(4:6,2) = fread(fdr,3,fdrfmt);

            f = find(and(buf0([1:4,6],1),buf0([1:4,6],2)==inflw(1:5)));
         elseif(strpx0==(r-1)*c+1) % LL pixel
            buf0 = zeros(4,2);
            buf0(ismember(strpx0(j,1)+idx0([2,3,5,6]),strpx0(:,1)),1) = 1;

            fseek(fdr,fdrnumbyt*(strpx0(j,1)-c-2),-1);
            buf0(1:2,2) = fread(fdr,2,fdrfmt);
            fseek(fdr,fdrnumbyt*(strpx0(j,1)-2),-1);
            buf0(3:4,2) = fread(fdr,2,fdrfmt);

            f = find(and(buf0([1:2,4],1),buf0([1:2,4],2)==inflw([2:3,5])));
         elseif(strpx0==r*c) % LR pixel
            buf0 = zeros(4,2);
            buf0(ismember(strpx0(j,1)+idx0([1,2,4,5]),strpx0(:,1)),1) = 1;

            fseek(fdr,fdrnumbyt*(strpx0(j,1)-c-2),-1);
            buf0(1:2,2) = fread(fdr,2,fdrfmt);
            fseek(fdr,fdrnumbyt*(strpx0(j,1)-2),-1);
            buf0(3:4,2) = fread(fdr,2,fdrfmt);

            f = find(and(buf0(1:3,1),buf0(1:3,2)==inflw([1:2,4])));
         end;
      elseif(and(md==1,~ismember(strpx0(j,1),[1,(r-1)*c+1]))) % in interior of left boundary
         buf0 = zeros(6,2);
         buf0(ismember(strpx0(j,1)+idx0([2,3,5,6,8,9]),strpx0(:,1)),1) = 1;

         fseek(fdr,fdrnumbyt*(strpx0(j,1)-c-1),-1);
         buf0(1:2,2) = fread(fdr,2,fdrfmt);
         fseek(fdr,fdrnumbyt*(strpx0(j,1)-1),-1);
         buf0(3:4,2) = fread(fdr,2,fdrfmt);
         fseek(fdr,fdrnumbyt*(strpx0(j,1)+c-1),-1);
         buf0(5:6,2) = fread(fdr,2,fdrfmt);

         f = find(and(buf0([1:2,4:6],1),buf0([1:2,4:6],2)==inflw([2:3,5,7:8])));
      elseif(and(~md,~ismember(strpx0(j,1),[c,r*c]))) % in interior of right boundary
         buf0 = zeros(6,2);
         buf0(ismember(strpx0(j,1)+idx0([1,2,4,5,7,8]),strpx0(:,1)),1) = 1;

         fseek(fdr,fdrnumbyt*(strpx0(j,1)-c-2),-1);
         buf0(1:2,2) = fread(fdr,2,fdrfmt);
         fseek(fdr,fdrnumbyt*(strpx0(j,1)-2),-1);
         buf0(3:4,2) = fread(fdr,2,fdrfmt);
         fseek(fdr,fdrnumbyt*(strpx0(j,1)+c-2),-1);
         buf0(5:6,2) = fread(fdr,2,fdrfmt);

         f = find(and(buf0([1:3,5:6],1),buf0([1:3,5:6],2)==inflw([1:2,4,6:7])));
      end;

      if(isempty(f))
         strhd(hdct+1) = strpx0(j,1);
         hdct = hdct+1;
      elseif(length(f)>1)
         strcf(cfct+1) = strpx0(j,1);
         cfct = cfct+1;
      end;
   end;
   strhd = strhd(1:hdct);
   strcf = strcf(1:cfct);
   save(matf,'strhd','strcf','-append');
else
   load(matf,'strhd','strcf');
   hdct = size(strhd,1);
   cfct = size(strcf,1);
end;
disp(sprintf('%d headwater pixels identified',hdct));
disp(sprintf('%d confluence pixels identified',cfct));
%=============================

%=============================
% Step 2: create initial 'seg_info0' table for all stream paths
%-----------------------------
% Columns of 'seg0':
% (1) start pixel
% (2) end pixel
% (3) FAC value @ start
% (4) FAC value @ end
% (5) segment length
% (6) network flag (0 if segment exits the image; 1 otherwise)

% Total number of segments will be hdct+cfct
sgct=hdct+cfct;

bdy0 = idx0([1:4,6:9]);
wh = whos('-file',matf);
if(~ismember('seg_info0',{wh.name}))
   strpx = zeros(numpx,3); % pxID, FAC, 'seg_info0' row ID
   pxct = 0;
   seg_info0 = [[strhd; strcf] zeros(sgct,5)];
   for j = 1:sgct
      v0 = seg_info0(j,1); % start pixel
      seg_info0(j,3) = strpx0(strpx0(:,1)==v0,2); % flow accumulation at segment start
      mv0 = mod(v0,c);
      fseek(fdr,fdrnumbyt*(v0-1),-1);
      dr0 = fread(fdr,1,fdrfmt); % flow direction at segment start
      fseek(fac,(v0-1)*facnumbyt,-1);
      strpx(pxct+1,:) = [v0 seg_info0(j,3) j];
      pxct = pxct+1;
      ct = 1; flg = 1;
      while(1>0)
         v1 = v0+bdy0(outflw==dr0);
         mv1 = mod(v1,c);
         % Stop if path leaves image extent
         if(or(or(v1<1,v1>r*c),all(ismember([0,1],[mv0,mv1]))))
            flg = 0;
            break;
         % Stop if path reaches a confluence
         elseif(ismember(v1,strcf))
            break;
         % Stop if path leaves study area extent, or continue downstream
         else
            fseek(fdr,fdrnumbyt*(v1-1),-1);
            dr0 = fread(fdr,1,fdrfmt);
            if(~ismember(dr0,outflw))
               flg =  0;
               break;
            else
               v0 = v1;
               mv0 = mv1;
               ct = ct+1;
               strpx(pxct+1,:) = [v0 strpx0(strpx0(:,1)==v0,2) j];
               pxct = pxct+1;
            end;
         end;
      end;
      seg_info0(j,2) = v0;
      seg_info0(j,4) = strpx0(strpx0(:,1)==v0,2);
      seg_info0(j,5) = ct; % segment length
      seg_info0(j,6) = flg; % 0 if segment flows out of image, 1 otherwise
   end;
   save(matf,'strpx','seg_info0','-append');
else
   load(matf,'strpx','seg_info0');
end;
%=============================

%=============================
% Step 3: subdivide segments based on 'facthr' and 'seglen'
%-----------------------------
% First determine 'natural' segmentation based of FAC jumps greater than or
%   equal to 'facthr'.
% Then bisect segments that are too long at the largest natural breaks
%   until all segments are no longer than 'seglen'.

% Columns of 'seg_info':
% (1) start pixel
% (2) end pixel
% (3) FAC value @ start
% (4) FAC value @ end
% (5) segment length
% (6) row in 'seg_info' for downstream segment (= 0 if none)
% (7) reach ID (= row of input variable 'pts')

wh = whos('-file',matf);
if(~ismember('seg_info',{wh.name}))
   seg_info = zeros(sgct*100,7);
   ct = 0;
   for j = 1:sgct
      p1 = strpx(strpx(:,3)==j,1:2); % p1: pxID, FAC value
      df = p1(2:end,2)-p1(1:end-1,2);

      f = find(df>=facthr);
      lf = length(f);
      if(lf)
         s0 = [p1(1,1) p1(f(1),1) p1(1,2) p1(f(1),2) f(1)];
         if(lf>1)
            s0 = [s0; [p1(f(1:end-1)+1,1) p1(f(2:end),1) p1(f(1:end-1)+1,2) p1(f(2:end),2) f(2:end)-f(1:end-1)]];
         end;
         s0 = [s0; [p1(f(end)+1,1) p1(end,1) p1(f(end)+1,2) p1(end,2) size(p1,1)-f(end)]];
      else
         s0 = [p1(1,1) p1(end,1) p1(1,2) p1(end,2) size(p1,1)];
      end;

      while(any(s0(:,5)>seglen))
         f = find(s0(:,5)>seglen); % identify subsegments that exceed max seg length
         if(~isempty(f))
            f0 = f(1);
%            n = size(s0,1);
            spth = p1(find(p1==s0(f0,1))+(0:s0(f0,5)-1),:);
            sdf = spth(2:end,2)-spth(1:end-1,2);
            f1 = find(sdf==max(sdf));
            lf1 = length(f1);
            % if 3 or more points, assume must have a straight segment through a flat area.
            %  In this case, just take the middle point & let the chips fall later.
            if(lf1>2)
               if(mod(lf1,2)) % odd number of points
                  f1 = median(f1);
               else % even
                  f1 = f1(lf1/2);
               end;
               lf1 = 1;
            end;
            s1 = [spth(1,1) spth(f1(1),1) spth(1,2) spth(f1(1),2) f1(1)];
            if(lf1>1)
               s1 = [s1; [spth(f1(1:end-1)+1,1) spth(f1(2:end),1) spth(f1(1:end-1)+1,2) spth(f1(2:end),2) f1(2:end)-f1(1:end-1)]];
            end;
            s1 = [s1; [spth(f1(end)+1,1) spth(end,1) spth(f1(end)+1,2) spth(end,2) size(spth,1)-f1(end)]];
            s0 = [s0(1:f0-1,:); s1; s0(f0+1:end,:)];
         end;
      end;
      s0 = [s0 [ct+(2:size(s0,1))'; 0] j*ones(size(s0,1),1)];
      sz0 = size(s0,1);
      seg_info(ct+(1:sz0),:) = s0;
      ct = ct+sz0;
   end;
   seg_info = seg_info(1:ct,:);
   f = find(~seg_info(:,6));
   lf = length(f);
   for j = 1:lf
      v0 = seg_info(f(j),2); % start pixel (end of a natural segment)
      fseek(fdr,fdrnumbyt*(v0-1),-1);
      dr0 = fread(fdr,1,fdrfmt); % flow direction at segment start
      v1 = v0+bdy0(outflw==dr0);
      f1 = find(seg_info(:,1)==v1);
      if(~isempty(f1))
         seg_info(f(j),6) = f1;
      end;
   end;
   
   % Add new column to 'strpx' for new segID
   strpx = [strpx zeros(numpx,1)]; 
   for j = 1:sgct
      f = find(seg_info(:,7)==j);
      seg0 = seg_info(f,:);
      sgct1 = size(seg0,1);
      fpx = find(strpx(:,3)==j);
      p1 = strpx(fpx,1:2); % p1: pxID, FAC value
      for k=1:sgct1
         f1=find(p1(:,1)==seg0(k,1));
         f2=find(p1(:,1)==seg0(k,2));
         p1(f1:f2,2)=f(k);
      end;
      strpx(fpx,4) = p1(:,2);
   end;
   
   header={'start pixel','end pixel','FAC at start','FAC at end','segment length',...
      'downstream segment (row index; 0 if none)','reach ID'};
   save(matf,'strpx','seg_info','header','-append');
else
   load(matf,'seg_info');
   sgct1 = size(seg_info,1);
end;
disp(sprintf('%d stream segments identified using threshold & length criteria',size(seg_info,1)));
fclose('all');
%=============================

%=============================
% Step 4: create BIL with segment IDs assigned to stream pixels
%-----------------------------
fhdr0 = fopen([facf(1:end-3),'hdr'],'rt');
txt0 = fgetl(fhdr0);
biltyp = 0; % ESRI
if(any(txt0==':'))
   biltyp = 1; % ERDAS
end;   

sgct = size(seg_info,1);
if(sgct<256)
   fmt = 'uint8';
   hdrfmt = 'U8';
   flg = 1;
elseif(sgct<65536)
   fmt = 'uint16';
   hdrfmt = 'U16';
   flg = 2;
else
   fmt = 'uint32';
   hdrfmt = 'U32';
   flg = 4;
end;

fhdr0 = fopen([facf(1:end-3),'hdr'],'rt');
fhdr = fopen([segf(1:end-3),'hdr'],'wt');

pr = ceil(strpx(:,1)/c);
pc = strpx(:,1)-(pr-1)*c;
urw = unique(pr);
zro = zeros(c,1);
fid = fopen(segf,'wb','b');
for j=1:r
   if(ismember(j,urw))
      f = find(pr==j);
      zro0 = zro;
      zro0(pc(f)) = strpx(f,4);
      fwrite(fid,zro0,fmt);
   else
      fwrite(fid,zro,fmt);
   end;
end;
fclose(fid);

%-----------------------------
% Write BIL header file (ERDAS or ESRI, depending on input)
%-----------------------------
fhdr0 = fopen([facf(1:end-3),'hdr'],'rt');
fhdr = fopen([segf(1:end-3),'hdr'],'wt');

if(biltyp) % ERDAS BIL
   while(~feof(fhdr0))
      txt0 = fgetl(fhdr0);
      if(any(txt0==':'))
         f = find(txt0==':');
         txt = txt0(f+1:end);
         if(isequal(txt0(1:f-1),'DATATYPE'))
            f1 = find(txt~=' ',1,'first');
            str = [txt0(1:f+f1-1),hdrfmt,10];
            fprintf(fhdr,str);
         elseif(isequal(txt0(1:f-1),'BYTE_ORDER'))
            f1 = find(txt0(f+1:end)~=' ',1,'first');
            if(flg==1)
               str = [txt0(1:f+f1-1),'NA',10];
            else
               str = [txt0(1:f+f1-1),'INTEL',10];
            end;
            fwrite(fhdr,str);
         else
            fprintf(fhdr,[txt0,10]);
         end;
      else
         fprintf(fhdr,[txt0,10]);
      end;
   end;
else % ESRI BIL
   fhdr0 = fopen([facf(1:end-3),'hdr'],'rt');
   while(~feof(fhdr0))
      txt0 = fgetl(fhdr0);
      f = find(txt0==' ',1,'first');
      txt = txt0(f+1:end);
      if(isequal(txt0(1:f-1),'NBITS'))
         nbits = str2num(txt(ismember(txt,'0123456789')));
      elseif(isequal(txt0(1:f-1),'BANDROWBYTES'))
         bandrowbytes = str2num(txt(ismember(txt,'0123456789')));
      elseif(isequal(txt0(1:f-1),'TOTALROWBYTES'))
         totalrowbytes = str2num(txt(ismember(txt,'0123456789')));
      end;
   end;
   bandrowbytes = num2str(bandrowbytes/(nbits/8)*flg);
   totalrowbytes = num2str(totalrowbytes/(nbits/8)*flg);
   
   fhdr0 = fopen([facf(1:end-3),'hdr'],'rt');
   while(~feof(fhdr0))
      txt0 = fgetl(fhdr0);
      f = find(txt0==' ',1,'first');
      txt = txt0(f+1:end);
      if(isequal(txt0(1:f-1),'BYTEORDER'))
         f1 = find(txt~=' ',1,'first');
         str = [txt0(1:f+f1-1),'M',10];
         fprintf(fhdr,str);
      elseif(isequal(txt0(1:f-1),'NBITS'))
         f1 = find(txt~=' ',1,'first');
         str = [txt0(1:f+f1-1),num2str(8*flg),10];
         fprintf(fhdr,str);
      elseif(isequal(txt0(1:f-1),'BANDROWBYTES'))
         f1 = find(txt~=' ',1,'first');
         str = [txt0(1:f+f1-1),bandrowbytes,10];
         fprintf(fhdr,str);
      elseif(isequal(txt0(1:f-1),'TOTALROWBYTES'))
         f1 = find(txt~=' ',1,'first');
         str = [txt0(1:f+f1-1),totalrowbytes,10];
         fprintf(fhdr,str);
      elseif(isequal(txt0(1:f-1),'PIXELTYPE'))
         f1 = find(txt~=' ',1,'first');
         str = [txt0(1:f+f1-1),'UNSIGNEDINT',10];
         fprintf(fhdr,str);
      elseif(isequal(txt0(1:f-1),'NODATA'))
         f1 = find(txt~=' ',1,'first');
         str = [txt0(1:f+f1-1),'0',10];
         fprintf(fhdr,str);
      else
         fprintf(fhdr,[txt0,10]);
      end;
   end;
end;
fclose('all');
if(exist([fdrf(1:end-3),'prj'],'file'))
   system(['cp ',fdrf(1:end-3),'prj ',segf(1:end-3),'prj']);
end;
%=============================
disp(sprintf('Segmented stream network file creation completed'));

