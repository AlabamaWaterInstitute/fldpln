%--------------------------------
% function create_segdb_facthr_maxlen_seedpts(pts,facthr,seglen,fdrfile,facfile,segfile,matfile)
%--------------------------------
% This function takes stream start points ('pts') and identifies the
%  downstream stream network.  Using a user-provided flow accumulation
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
%  pts = n-by-2 vector with <x,y> coordinates; OR
%        n-by-4 vector with last two columns containing stream endpoints; OR
%        path and filename for shapefile containing stream start points
%  facthr = flow accumulation jump size used to define initial, natural breaks
%         = [] if no breaks are to be defined by FAC jumps
%  seglen = maximum segment length (in number of pixels)
%         = [] if no breaks are to be defined by segment length
%  fdrfile = path and filename for input flow direction BIL file
%  facfile = path and filename for input flow accumulation BIL file
%     NOTE:  If facthr = [], seglen = [], and facfile = [], then no FAC
%            values are written to the output table
%  segfile = path and filename for output segment ID BIL file; format will
%    be either 'uint8', 'uint16', or 'uint32', depending on number of segments
%  matfile = path and filename where output 'seg_info' MATLAB table is to be stored
%
% NOTE: For 'pts', order the stream start points from main stem to tributary.
%       E.g., in a forked system, the first row should give the headwater 
%       coordinate for the branch that will eventually become the main stem;
%       the second row the next major branch; ...; and the last row for the
%       least significant tributary.
%  
%--------------------------------
function create_segdb_facthr_maxlen_seedpts(pts0,facthr,seglen,fdrf,facf,segf,matf)
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

if(ischar(pts0))
   shpf = pts0;
   pts0 = readshapefield(shpf,{'X','Y'});
end;
[npts,p2] = size(pts0);
if(and(p2~=2,p2~=4))
   error('Seed point matrix must have 2 or 4 columns (4 if path endpoints are specified)');
end;
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
r = fdrinfo.r;
c = fdrinfo.c;
pxsz = fdrinfo.pxszx;
ulx = fdrinfo.ulx;
uly = fdrinfo.uly;
%=============================

%=============================
% Read FAC header & open FAC for reading
%-------------------------
facflg = 1;
if(~isempty(facf))
    facinfo = readbilheader(facf);
    facfmt = facinfo.fmt;
    facnumbyt = facinfo.nbyts;
    if(facinfo.bytord)
       fac = fopen(facf,'rb','b');
    else
       fac = fopen(facf,'rb');
    end;
else
   if(and(isempty(facthr),isempty(seglen)))
      facflg = 0;
   else
      error('Both ''facthr'' and ''seglen'' must = [] if ''facfile'' = []');
   end;
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
% Determine precise (decimal) row and column values for each point
% Integer values correspond with pixel centroids [UL = (0,0)]
pts = zeros(size(pts0));
pts(:,1) = round((uly-pts0(:,2))/pxsz);
pts(:,2) = round((pts0(:,1)-ulx)/pxsz);

eflg = 0; % flag indicating whether or not endpoints are specified
if(any(or(or(pts(:,1)<0,pts(:,1)>r-1),or(pts(:,2)<0,pts(:,2)>c-1))))
   error('At least one start point file coordinate is outside the image dimensions');
else
   pts0 = pts(:,1)*c+pts(:,2)+1; % convert file coords to vector
end;

if(p2==4)
   if(any(or(or(pts(:,3)<0,pts(:,3)>r-1),or(pts(:,4)<0,pts(:,4)>c-1))))
      error('At least one end point file coordinate is outside the image dimensions');
   else
      pts0 = [pts0 pts(:,3)*c+pts(:,4)+1]; % convert file coords to vector
      eflg = 1;
   end;
end;
% pts0 = BIL pixel IDs
% pts = <row,col> IDs (UL = [0,0])
%=============================

%=============================
% Step 1: Identify & count the number of stream pixels
%-----------------------------
bdyflw = [32 64 128 16 1 8 4 2]'; % outflow
idx0 = [-c+(-1:1)'; [-1;1]; c+(-1:1)'];

if(~exist(matf,'file'))
   for j = 1:npts
      js = int2str(j);
      v0 = pts0(j,1);
      mv0 = mod(v0,c);
      fseek(fdr,fdrnumbyt*(v0-1),-1);
      dr0 = fread(fdr,1,fdrfmt);

      if(facflg)
         fseek(fac,facnumbyt*(v0-1),-1);
         f0 = fread(fac,1,facfmt);
      else
         f0 = 0;
      end;

      pth = zeros(100*max(r,c),2); % BIL pixel ID, FAC value

      ct = 1;
      bkflg = 0;
      while(1>0)
         pth(ct,:) = [v0 f0];
         didx = find(bdyflw==dr0);
         if(~isempty(didx))
            v0 = v0+idx0(didx);
            mv1 = mod(v0,c);
            % Stop if path leaves image extent
            if(or(bkflg,or(or(v0<1,v0>r*c),all(ismember([0,1],[mv0,mv1])))))
               break;
            end;
            if(eflg)
               if(v0==pts0(j,2)) % or if specified endpoint is reached
                  bkflg = 1;
               end;
            end;
            mv0 = mv1;
            ct = ct+1;
            fseek(fdr,fdrnumbyt*(v0-1),-1);
            dr0 = fread(fdr,1,fdrfmt);
            if(facflg)
               fseek(fac,facnumbyt*(v0-1),-1);
               f0 = fread(fac,1,facfmt);
            end;
         % Stop if path leaves study area extent
         else
            ct = ct-1;
            break;
         end;
      end;
      eval(['pth',js,' = pth(1:ct,:);']);

      if(j>1)
         save(matf,['pth',js],'-append');
      else
         save(matf,'pth1');
      end;
      disp(sprintf('%d of %d stream paths identified (%d pixels long)',j,npts,ct))
   end;
else
   load(matf);
end;
fclose(fdr);
if(facflg)
   fclose(fac);
end;
clear('pth');
%=============================

%=============================
% Step 2: Look for stream path overlaps and assign to highest priority path only.
%         Also, save confluence points as break points, in case they don't
%         exceed the FAC threshold.
%
%         If segment endpoints are specified, trim the downstream ends as well.
%------------------------
wh = whos('-file',matf);
if(~ismember('seg_info',{wh.name}))
   if(eflg)
      for j = 1:npts
         js = int2str(j);
         p1 = eval(['pth',js]);
         if(any(p1(:,1)==pts0(j,2)))
            f = find(p1(:,1)==pts0(j,2));
            p1 = p1(1:f,:);
         end;
         eval(['pth',js,' = p1;']);
      end;
   end;

   conpts = [];
   for j = 2:npts
      js = int2str(j);
      p1 = eval(['pth',js]);
      for k = 1:j-1
         ks = int2str(k);
         p2 = eval(['pth',ks]);
         if(any(ismember(p1(:,1),p2(:,1))))
            f = find(ismember(p1(:,1),p2(:,1)),1);
            conpts = [conpts; [p1(f,1) j]];
            p1 = p1(1:f-1,:);
         end;
      end;
      eval(['pth',js,' = p1;']);
   end;
   %=============================

   %=============================
   % Step 3: Create 'seg_info' tables for all stream paths
   %-----------------------------

   % First determine 'natural' segmentation based of FAC jumps greater than or
   %   equal to 'facthr'.
   % Then bisect segments that are too long at the largets natural breaks
   %   until all segments are no longer than 'seglen'.

   % Columns of 'seg_info':
   % (1) start pixel
   % (2) end pixel
   % (3) FAC value @ start
   % (4) FAC value @ end
   % (5) segment length
   % (6) row in 'seg_info' for downstream segment (= 0 if none)
   % (7) reach ID (= row of input variable 'pts')

   seg_info = [];
   for j = 1:npts
      js = int2str(j);
      p1 = eval(['pth',js]);
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
            n = size(s0,1);
            spth = p1(find(p1==s0(f0,1))+(0:s0(f0,5)-1),:);
            sdf = spth(2:end,2)-spth(1:end-1,2);
            f1 = find(sdf==max(sdf));
            lf1 = length(f1);
            s1 = [spth(1,1) spth(f1(1),1) spth(1,2) spth(f1(1),2) f1(1)];
            if(lf1>1)
               s1 = [s1; [spth(f1(1:end-1)+1,1) spth(f1(2:end),1) spth(f1(1:end-1)+1,2) spth(f1(2:end),2) f1(2:end)-f1(1:end-1)]];
            end;
            s1 = [s1; [spth(f1(end)+1,1) spth(end,1) spth(f1(end)+1,2) spth(end,2) size(spth,1)-f1(end)]];
   %         s0(f0,:)
   %         s1
            s0 = [s0(1:f0-1,:); s1; s0(f0+1:end,:)];
         end;
      end;
      s0 = [s0 [size(seg_info,1)+(2:size(s0,1))'; 0] j*ones(size(s0,1),1)];
      seg_info = [seg_info; s0];
   end;
   %f = find(seg_info(:,7)==1);
   %[size(pth1,1) sum(seg_info(f,5))]
   %f = find(seg_info(:,7)==2);
   %[size(pth2,1) sum(seg_info(f,5))]

   %-----------------------------
   % Check to see if all confluence points occur as natural breaks, and add
   %   new breaks as needed
   if(and(facflg,~isempty(conpts)))
      f = (~ismember(conpts(:,1),seg_info(:,1)));
      if(sum(f)>0)
         error('Need to add code for forcing a new segment break');
         conpts = conpts(f,:);
         len = size(conpts,1);
         for j = 1:len
            for k = 1:npts
               ks = int2str(k);
               p1 = eval(['pth',ks]);
               if(ismember(conpts(j,1),p1))
                  %***************************
                  % add code here
                  error('Congratulations!  You tripped an error!');
                  %***************************
               end;
            end;
         end;
      end;
      %-----------------------------

      f = find(ismember(conpts(:,1),seg_info(:,1)));
      lf = length(f);
      for j = 1:lf
         f0 = find(seg_info(:,7)==conpts(f(j),2));
         f1 = find(seg_info(:,1)==conpts(f(j),1));
         seg_info(f0(end),6) = f1;
      end;
   end;

   header = {'start pixel','end pixel','FAC at start','FAC at end','segment length',...
      'downstream segment (row index; 0 if none)','reach ID'};
   save(matf,'seg_info','header','-append');
else
   load(matf,'seg_info');
end;
%=============================

%=============================
% Step 4: create BIL with segment IDs assigned to stream pixels
%-----------------------------
fhdr0 = fopen([fdrf(1:end-3),'hdr'],'rt');
txt0 = fgetl(fhdr0);
biltyp = 0; % ESRI
if(any(txt0==':'))
   biltyp = 1; % ERDAS
end;   

% Assign all stream pixels to a segment
ct = 0;
pth = [];
for j = 1:npts
   js = int2str(j);
   p1 = eval(['pth',js]);
   f = find(seg_info(:,7)==j);
   seg0 = seg_info(f,:);
   sgct = size(seg0,1);
   for k = 1:sgct
      f1 = find(p1(:,1)==seg0(k,1));
      f2 = find(p1(:,1)==seg0(k,2));
      p1(f1:f2,2) = f(k);
   end;
   ct = ct+sgct;
   pth = [pth; p1];
end;

if(ct<256)
   fmt = 'uint8';
   fmthdr = 'U8';
   flg = 1;
elseif(ct<65536)
   fmt = 'uint16';
   fmthdr = 'U16';
   flg = 2;
else
   fmt = 'uint32';
   fmthdr = 'U32';
   flg = 4;
end;

pr = ceil(pth(:,1)/c);
pc = pth(:,1)-(pr-1)*c;
urw = unique(pr);
zro = zeros(c,1);
fid = fopen(segf,'wb','b');
for j = 1:r
   if(ismember(j,urw))
      f = find(pr==j);
      zro0 = zro;
      zro0(pc(f)) = pth(f,2);
      fwrite(fid,zro0,fmt);
   else
      fwrite(fid,zro,fmt);
   end;
end;
fclose(fid);

%-----------------------------
% Write BIL header file (ERDAS or ESRI, depending on input)
%-----------------------------
fhdr0 = fopen([fdrf(1:end-3),'hdr'],'rt');
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
   fhdr0 = fopen([fdrf(1:end-3),'hdr'],'rt');
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
   
   fhdr0 = fopen([fdrf(1:end-3),'hdr'],'rt');
   while(~feof(fhdr0))
      txt0 = fgetl(fhdr0);
      f = find(txt0==' ',1,'first');
      txt = txt0(f+1:end);
      if(isequal(txt0(1:f-1),'NBITS'))
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
