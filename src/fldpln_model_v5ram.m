%--------------------------------
% function fldpln_model_erdasbil(seg,inp,fldmn,fldmx,dh,mxht,fmt,ssflg,bg*)
%--------------------------------
% This function uses a filled DEM and a flow direction map (both ERDAS Imagine
%  or ESRI generic binary (BIL) rasters with header files) to map the floodplain.
%
% "Depth To Flood" (DTF) and "Flood Source Pixel" (FSP) information is
%  recorded for each pixel in the identified floodplain.  This info is
%  written to a MATLAB file that provides a floodplain database for the segment.
%
% INPUT
%  seg = segment ID number; if seg is a 2-vector, then this is assumed to be
%        a starting coordinate for a partial segment that is to be flooded;
%        NOTE: use Imagine file coordinate numbering (i.e., UL = (0,0))
%  inp = structure with input directory and file names
%     FIELDS:  bildr <directory where temporary FLDMAP files will be written>
%              segdr <directory where output FLDPLN MAT files will be saved>
%              dem <DEM BIL path & filename>
%              fil <FILLED DEM BIL path & filename>
%              fdr <Flow Direction BIL path & filename>
%              seg <Segment Information MAT path & filename>
%              segf <Segment ID BIL path & filename> (only need if start
%                    points are used rather than segment ID)
%  fldmn = minimum flood stage height (all floodplain heights lower than
%        this value get mapped to this value at the end)
%  fldmx = maximum flood stage height (all DTF values will be less than or
%        equal to this value)
%  dh = depth step for backfill-forward fill iterations
%  mxht = cap height = max(dem+flood height) to cease flooding;
%       enter 0 for no cap height
%  ssflg = 1 if the 'steady state' model is to be evaluated, where
%            spillover is performed exhastively during each iteration
%            (recommended for high resolution DEM, e.g., 10-m NED and 2-m LiDAR)
%        = 0 if only a single spillover flood step is to be used during
%            each iteration (recommended for low resolution DEM, e.g., 30-m
%            NED and 90-m NED)
%
% OPTIONAL INPUT(*)
%  bg = background value (default is 0); enter [] to skip 'bg' but define 'segf'
%
% ALGORITHM (applied to each pixel in stream network)
%  (0) floodplain is initialized to contain the stream segment (DTF = 0 at
%      all points)
%  (1) determine backfill flood reach for each floodplain interior boundary
%      pixel up to the current flood stage depth
%  (2) determine spillover flood reach for each floodplain interior boundary
%      pixel up to the current flood stage depth
%  (3) (ssflg = 1) repeat (1)-(2) until no new pixels are added
%  (4) update the flood stage depth
%  (5) go back to (1) until the maximum flood stage depth is attained
%
% OUTPUT properties for each identified floodplain pixel:
%  (i)  0<=DTF<=fldmx
%  (ii) FSP is a stream segment pixel indicating the stream pixel from
%       which floodwaters can originate that inundate the pixel at the
%       shallowest flood depth (which is the DTF value for the pixel)
%
% NOTE: Stream pixels in the trajectory of the most downstream segment pixel are
%       treated somewhat as background, to avoid redundant processing.
%--------------------------------
function y = fldpln_model_v5ram(seg0,inp,fildat,fdrdat,fldmn,fldmx,dh,mxht0,ssflg,varargin)
%--------------------------------
%   3-by-3 numbering for analysis (determined by BIL ordering):
%           1		2 		3
%           4		x		5
%           6		7		8
%--------------------------------
%   3-by-3 Jenson & Domingue (1988) numbering for flow direction:
%           64		128   1
%           32		x		2
%           16		8		4
%
%   What we use:
%           32		64 	128
%           16		x		1
%           8		4		2
%--------------------------------
y = 1;
bildr = inp.bildr;

%=============================
% Read FIL file specs from header
%-------------------------
ffil = inp.fil;
filinfo = readbilheader(ffil);
r = filinfo.r;
c = filinfo.c;
fmt = filinfo.fmt;
%=============================

bdy0 = [-1-c -c 1-c -1 1 c-1 c c+1]';
bdy0_row = bdy0';
inflw = [2 4 8 1 16 128 64 32]';
outflw = flipud(inflw);
ub = realmax('single');

%---------------------------
% Set background value (default is 0)
bg = 0;
if(nargin==10)
   bg = varargin{1};
end
%---------------------------

%----------------------------
% Make strings to use in file names
sh = num2str(fldmx);
f = find(sh=='.');
if(~isempty(f))
   sh(f) = 'p';
end
sdh = num2str(dh);
f = find(sdh=='.');
if(~isempty(f))
   sdh(f) = 'p';
end
sseg = int2str(seg0);
outf = ['h',sh,'_dh',sdh,'_seg',sseg];
%----------------------------
% Create 'segment_files' directory if one does not exist
segdr = inp.segdr;
exflg = 0;
if(~exist(segdr,'dir'))
   mkdir(segdr);
else
   % Check for segment completion
   if(exist([segdr,outf,'_tmp.mat'],'file'))
      load([segdr,outf,'_tmp.mat'],'fldht');
      if(fldht>=fldmx)
         if(exist([segdr,outf,'.mat'],'file'))
            return; % exit program if segment already processed through FLDMX
         else
            exflg = 1;
         end
      end
   end
end
%----------------------------
if(~exflg)

%=============================
fseg = inp.seg;
if(ischar(fseg))
   load(fseg,'seg_info'); %
else
    seg_info = fseg;
end
% Columns of 'seg_info' (row corresponds with segment ID):
% Columns of 'seg_info':
% (1) start pixel
% (2) end pixel
% (3) FAC value @ start
% (4) FAC value @ end
% (5) segment length
% (6) row in 'seg_info' for downstream segment (= 0 if none)
% (7) reach ID (= row of input variable 'pts')

% (5) segment length, (6) network flag (0 if segment exits the image; 1 otherwise),
% (7) immediately downstream segment (assigned in Step 4; = 0 if segment flows
% out of study area)

%----------------------------
% Determine segment pixels (create 'strpts', 'lstr')
bkflg = 0;
if(length(seg0)==1)
   ispt = 0;
   lstr = seg_info(seg0,5);
   strpts = zeros(lstr,1);
   v0 = seg_info(seg0,1); % start pixel
   d0 = double(fdrdat(v0)); % flow direction at segment start
   strpts(1) = v0;
   for k = 2:lstr-1
      v1 = strpts(k-1)+bdy0(outflw==d0);
      if(isempty(v1))
         lstr = k-1;
         strpts = strpts(1:lstr);
         bkflg = 1;
         break;
      end
      strpts(k) = v1;
      v0 = v1;
      d0 = double(fdrdat(v0));
   end
   if(and(lstr>1,~bkflg))
      v1 = v0+bdy0(outflw==d0);
      strpts(lstr) = v1;
   end
elseif(length(seg0)==2) % no need to adapt 'segf' reading to RAM
   ispt = 1;
   v0 = seg0(1)*c+seg0(2)+1;
   
   %------------------------
   % Open segment ID file for reading
   segf = inp.segf;
   filinfo = readbilheader(segf);
   segfmt = filinfo.fmt;
   filtyp = filinfo.bytord;
   numbyt = filinfo.nbyts;
   if(filtyp)
      seg = fopen(segf,'rb','b');
   else
      seg = fopen(segf,'rb');
   end
   %------------------------
   
   fseek(seg,(v0-1)*numbyt,-1);
   tmp = fread(seg,1,segfmt);
   if(~tmp)
      error('Input point is not in the stream network')
   end
   fclose(seg);
   seg0 = tmp;
   strpts = zeros(100*max(r,c),1);
   lstr = 0;
   while(~sum(v0==seg_info(:,2)))
      strpts(lstr+1) = v0;
      lstr = lstr+1;
      d0 = double(fdrdat(v0));
      v0 = strpts(lstr)+bdy0(outflw==d0);
   end
   strpts = strpts(1:lstr);
else
   error('''seg0'' must be either a 1-vector or 2-vector');
end
%----------------------------

%----------------------------
% Determine downstream trajectory pixels to be excluded (create 'excl')
if(~ispt)
   v0 = seg_info(seg0,2);
   mv0 = mod(v0,c);
   d0 = double(fdrdat(v0));

   excl = zeros(100*max(r,c),1);
   ct = 1;
   while(1>0)
      excl(ct) = v0;
      if(any(outflw==d0))
         v0 = v0+bdy0(outflw==d0);
         mv1 = mod(v0,c);
         % Stop if path leaves image extent
         if(or(or(v0<1,v0>r*c),and(sum(~[mv0,mv1])==1,sum(1==[mv0,mv1])==1)))
            break;
         end
         mv0 = mv1;
         ct = ct+1;
         d0 = double(fdrdat(v0));
         % Stop if path leaves study area extent
      else
         ct = ct-1;
         break;
      end
   end
   excl = excl(2:ct,:);
else % No exclude zone if a point is entered
%   v0 = strpts;
   excl = [];
end
%----------------------------

%----------------------------
% Set appropriate value for 'mxht0'
if(~mxht0)
   mxht0 = ub;
else
   mxht0 = 0.99999*mxht0;
end
%----------------------------
itnum = ceil(fldmx/dh);
sz = min(lstr,4000)*max(5000,max(r,c)); % hopefully this catches the largest of floodplains
%----------------------------
% In case of interruptus...
dr = dir(segdr);
if(~ismember([outf,'_tmp.mat'],{dr.name})) % unprocessed segment
   %------------------------------
   % Initialize variables from scratch
   %------------------------------
   % Initialize floodplain information matrices
   %  fldpln_info: reference stream pixel (FSP), floodplain pixel, flood height (DTF)
   fldpln_info = zeros(max(lstr,sz),3);
   %  fldpln_bdy: potential boundary pixel indicator (1 = yes, 0 = no)
   fldpln_bdy = ones(max(lstr,sz),1,'uint8');
   
   % Read FIL values to eliminate stream pixels in 'NoData' areas (such as
   %  lakes that have been masked to avoid processing)
   ct_tot = 0;
   for j = 1:lstr
      tmp = fildat(strpts(j));
      if(tmp~=bg)
         fldpln_info(ct_tot+1,1:2) = [strpts(j) strpts(j)];
         ct_tot = ct_tot+1;
      end
   end

   % bdy: reference stream pixel (FSP), floodplain pixel, flood height (DTF)
   bdy = fldpln_info(1:ct_tot,:); % initialized floodplain boundary to unmasked stream points
   bdct = ct_tot;
   %------------------------------
   % Initialize temporary fldmap vector
   flddat = zeros(r*c,1,'single');
   flddat(fldpln_info(1:ct_tot,1)) = fldmn;
   %------------------------------
   fldht = 0;
   g0 = 1;
else % partially processed segment
   %------------------------------
   % Initialize variables at the current iteration
   %------------------------------
   % bdct, bdy, ct_tot, fldht, fldpln_bdy, fldpln_info
   load([segdr,outf,'_tmp.mat']);
   g0 = fldht/dh+1;
   %------------------------------
   % Initialize temporary fldmap BIL file (used in program, deleted at the end)
   if(g0<=itnum)
      flddat = zeros(r*c,1,'single');
      flddat(fldpln_info(1:ct_tot,2)) = max(fldmn,fldpln_info(1:ct_tot,3));
   end
   %------------------------------
end
if(~exist('proctime','var'))
   proctime = [];
end
tic
%------------------------------
% Iterate on flood zone boundary to build successive flood zones
%------------------------------
for g = g0:itnum
   fldht = fldht+min(dh,fldmx-fldht); % current iteration max flood height
   %------------------------------
   % Create initial raw flood zone map for 'fldht', back-filling from the
   %  current flood plain boundary point set
   for j = 1:bdct
      %  col1 = reach pixel, col2 = flood height
      fldtmp = zeros(100000,2); % hopefully this catches the biggest catchments
      fj0 = bdy(j,2); % current boundary pixel
      stght = bdy(j,3); % current flood stage height
      
      fil0 = double(fildat(fj0)); % elevation at boundary pixel
      v0 = fil0+fldht-stght; % flood elevation at boundary pixel
      mxht = min(mxht0,v0); % temporary max height for backfilling
      %------------------------
      % Initialize pixel-specific reach

      % Determine if point is on an image boundary, and exclude 'out of bounds' neighbors
      r0 = ceil(fj0/c);
      c0 = fj0-(r0-1)*c;
      f0 = 1:8;
      bdy_rc = [sum(r0==[1 r]) sum(c0==[1 c])];
      if(sum(bdy_rc)==2)
         if(and(r0==1,c0==1)), f0 = [5,7,8];
         elseif(and(r0==1,c0==c)), f0 = [4,6,7];
         elseif(and(r0==r,c0==1)), f0 = [2,3,5];
         elseif(and(r0==r,c0==c)), f0 = [1,2,4];
         end
      elseif(bdy_rc(1)>0)
         if(r0==1), f0 = 4:8;
         elseif(r0==r), f0 = 1:5;
         end
      elseif(bdy_rc(2)>0)
         if(c0==1), f0 = [2,3,5,7,8];
         elseif(c0==c), f0 = [1,2,4,6,7];
         end
      end
      idx = fj0+bdy0(f0);
      lidx = length(idx);
      dat = zeros(lidx,3); % FDR, FIL, FLDMAP0

      % Retain neighbor points that are inflowing and in the new flood zone
      dat(:,1) = double(fdrdat(idx));
      dat(:,2) = double(fildat(idx));
      dat(:,3) = double(flddat(idx));
      f1 = and(dat(:,1)==inflw(f0),and(dat(:,2)<=mxht,dat(:,2)~=bg));
      idx = idx(f1);
      % Retain neighbor points that are not already in the old flood zone
      if(~isempty(idx))
         dat = dat(f1,:);
         f1 = ~dat(:,3);
         idx = idx(f1);
      end
      ct1 = 0;
      if(~isempty(idx))
         ct1 = sum(f1);
         fldtmp(1:ct1,:) = [idx stght+max(0,dat(f1,2)-fil0)];
      end
      %------------------------

      %------------------------
      % Generate boundary pixel-specific reach
      ct0 = 0;
      while(ct0<ct1)
         ct0 = ct0+1;
         fj = fldtmp(ct0,1);
         %------------------------
         % Initialize boundary pixel-specific reach

         % Determine if point is on a boundary, and exclude 'out of bounds' neighbors
         r0 = ceil(fj/c);
         c0 = fj-(r0-1)*c;
         f0 = 1:8;
         bdy_rc = [sum(r0==[1 r]) sum(c0==[1 c])];
         if(sum(bdy_rc)==2)
            if(and(r0==1,c0==1)), f0 = [5,7,8];
            elseif(and(r0==1,c0==c)), f0 = [4,6,7];
            elseif(and(r0==r,c0==1)), f0 = [2,3,5];
            elseif(and(r0==r,c0==c)), f0 = [1,2,4];
            end
         elseif(bdy_rc(1)>0)
            if(r0==1), f0 = 4:8;
            elseif(r0==r), f0 = 1:5;
            end
         elseif(bdy_rc(2)>0)
            if(c0==1), f0 = [2,3,5,7,8];
            elseif(c0==c), f0 = [1,2,4,6,7];
            end
         end
         idx = fj+bdy0(f0);
         lidx = length(idx);
         dat = zeros(lidx,3); % FDR, FIL, FLDMAP0
         
         % Retain neighbor points that are inflowing and in the new flood zone
         dat(:,1) = double(fdrdat(idx));
         dat(:,2) = double(fildat(idx));
         dat(:,3) = double(flddat(idx));
         f1 = and(dat(:,1)==inflw(f0),and(dat(:,2)<=mxht,dat(:,2)~=bg));
         idx = idx(f1);
         % Retain neighbor points that are not already in the old flood zone
         if(~isempty(idx))
            dat = dat(f1,:);
            f1 = ~dat(:,3);
            idx = idx(f1);
         end
         
         if(~isempty(idx))
            ct2 = sum(f1);
            fldtmp(ct1+(1:ct2),:) = [idx stght+max(0,dat(f1,2)-fil0)];
            ct1 = ct1+ct2;
         end
      end
      if(size(fldpln_info,1)<ct_tot+ct1)
         fldpln_info = [fldpln_info; zeros(max(ct1,10*max(r,c)),3)];
      end
      fldpln_info(ct_tot+(1:ct1),:) = [bdy(j,1)*ones(ct1,1) fldtmp(1:ct1,:)];
      % Update temporary flood vector
      flddat(fldpln_info(ct_tot+(1:ct1),2)) = max(fldmn,fldpln_info(ct_tot+(1:ct1),3));
      ct_tot = ct_tot+ct1;
   end

   %------------------------------
   % Identify current flood zone boundary pixels
   %------------------------------
   if(length(fldpln_bdy)<ct_tot+ct1)
      fldpln_bdy = [fldpln_bdy; ones(max(ct1,10*(ct_tot-length(fldpln_bdy))),1)];
   end
   f = find(fldpln_bdy(1:ct_tot));
   tmp = fldpln_info(f,2);
   tmp_r = ceil(tmp/c);
   tmp_c = tmp-(tmp_r-1)*c;
   f1 = find(~((tmp_r==1)+(tmp_r==r)+(tmp_c==1)+(tmp_c==c)));
   fbdr = find((tmp_r==1)+(tmp_r==r)+(tmp_c==1)+(tmp_c==c)>0);
   clear tmp_r tmp_c
   nfbdr = length(fbdr);
   tmp = tmp(f1);
   lf1 = length(f1);
   
   bdy_tot = zeros(lf1,8);
   for j = 1:lf1
      bdy_tot(j,:) = double(flddat(tmp(j)+bdy0)');
   end
   fldpln_bdy(f(f1(all(bdy_tot,2)))) = 0;
   clear bdy_tot
   for j = 1:nfbdr
      if(fldpln_bdy(f(fbdr(j))))
         fj = fldpln_info(f(fbdr(j)),2);
         % Determine if point is on a boundary, and exclude 'out of bounds' neighbors
         r0 = ceil(fj/c);
         c0 = fj-(r0-1)*c;
         f0 = 1:8;
         bdy_rc = [sum(r0==[1 r]) sum(c0==[1 c])];
         if(sum(bdy_rc)==2)
            if(and(r0==1,c0==1)), f0 = [5,7,8];
            elseif(and(r0==1,c0==c)), f0 = [4,6,7];
            elseif(and(r0==r,c0==1)), f0 = [2,3,5];
            elseif(and(r0==r,c0==c)), f0 = [1,2,4];
            end
         elseif(bdy_rc(1)>0)
            if(r0==1), f0 = 4:8;
            elseif(r0==r), f0 = 1:5;
            end
         elseif(bdy_rc(2)>0)
            if(c0==1), f0 = [2,3,5,7,8];
            elseif(c0==c), f0 = [1,2,4,6,7];
            end
         end
         idx = fj+bdy0(f0);
         if(all(ismember(idx,fldpln_info(1:ct_tot,2))))
            fldpln_bdy(f(fbdr(j))) = 0;
         end
      end
   end
   % col1 = reference stream pixel, col2 = boundary pixel, col3 = flood depth required to get to boundary pixel
   bdy = fldpln_info(fldpln_bdy(1:ct_tot)>0,:);
   bdct = size(bdy,1);

   %------------------------------
   % Identify points "spilled into" by boundary pixels
   %------------------------------
   spill_flg = 1;
   while(spill_flg)
      ct_tot2 = ct_tot;
      tmp = bdy(:,2);
      tmp_r = ceil(tmp/c);
      tmp_c = tmp-(tmp_r-1)*c;
      f1 = find(~((tmp_r==1)+(tmp_r==r)+(tmp_c==1)+(tmp_c==c)));
      lf1 = length(f1);
      fbdr = find((tmp_r==1)+(tmp_r==r)+(tmp_c==1)+(tmp_c==c)>0);
      clear tmp_r tmp_c
      nfbdr = length(fbdr);
      % Rearrange 'bdy' matrix so that entries on the image border are moved to the end
      bdy = bdy([f1;fbdr],:);
      tmp = bdy(:,2);
      % obdy_idx: off-boundary vector index location values (some will be
      %           mis-specified for image boundary pixels (rows lf1+1:bdct),
      %           but these will be ignored later)
      obdy_idx = tmp(:,ones(1,8))+bdy0_row(ones(bdct,1),:);
      clear tmp
      % obdy_fill: off-boundary DEM values; =ub if pixel is in the current
      %            flood zone, out of image extent, in the 'excl' set, or in
      %            an undefined area (bg FIL value)
      obdy_fill = ub*ones(bdct,8);
      % First look at 'bdy' pixels not on the image boundary
      for j = 1:lf1
         for k = 1:8
            if(~sum(obdy_idx(j,k)==excl)) % not in 'excl'
               if(~flddat(obdy_idx(j,k))) % not in current flood zone
                  if(fildat(obdy_idx(j,k))~=bg) % not in an undefined area
                     obdy_fill(j,k) = double(fildat(obdy_idx(j,k)));
                  end
               end
            end
         end
      end
      % Next look at 'bdy' pixels on the image boundary
      for j = 1:nfbdr
         fj = bdy(lf1+j,2);
         % Determine if point is on a boundary, and exclude 'out of bounds' neighbors
         r0 = ceil(fj/c);
         c0 = fj-(r0-1)*c;
         f0 = 1:8;
         bdy_rc = [sum(r0==[1 r]) sum(c0==[1 c])];
         if(sum(bdy_rc)==2)
            if(and(r0==1,c0==1)), f0 = [5,7,8];
            elseif(and(r0==1,c0==c)), f0 = [4,6,7];
            elseif(and(r0==r,c0==1)), f0 = [2,3,5];
            elseif(and(r0==r,c0==c)), f0 = [1,2,4];
            end
         elseif(bdy_rc(1)>0)
            if(r0==1), f0 = 4:8;
            elseif(r0==r), f0 = 1:5;
            end
         elseif(bdy_rc(2)>0)
            if(c0==1), f0 = [2,3,5,7,8];
            elseif(c0==c), f0 = [1,2,4,6,7];
            end
         end
         for k = f0 % ignore mis-specified 'out of bounds' entries of 'obdy_idx'
            if(~sum(obdy_idx(lf1+j,k)==excl)) % not in 'excl'
               if(~flddat(obdy_idx(lf1+j,k))) % not in current flood zone
                  if(fildat(obdy_idx(lf1+j,k))~=bg) % not in an undefined area
                     obdy_fill(lf1+j,k) = double(fildat(obdy_idx(lf1+j,k)));
                  end
               end
            end
         end
      end
      % Recode to 0 all 'obdy_idx' values that do not matter (i.e., all 'ub'
      %  entries in 'obdy_fill')
      obdy_idx(obdy_fill==ub) = 0;

      % Read in FIL values for boundary pixels
      bdy_fill = double(fildat(bdy(:,2)));

      % Determine if there are any candidate "spilled into" points for each
      %  boundary point; recode all non-candidates to 0 in 'obdy_idx'
      for j = 1:bdct
         f1 = find(obdy_fill(j,:)<=min(mxht0,bdy_fill(j)+fldht-bdy(j,3)));
         if(~isempty(f1))
            f2 = setxor(1:8,f1);
            obdy_idx(j,f2) = 0;
         else
            obdy_idx(j,:) = 0;
         end
      end

      % Determine unique, non-zero values in 'obdy_idx', and read FIL values
      unq = unique(obdy_idx(obdy_idx>0));
      obdct = length(unq);

      % obdy_fill:  obdy FIL (spilled into), bdy FIL (spilling pixel)
      obdy_fill = zeros(obdct,2);
      obdy_fill(:,1) = double(fildat(unq));
      
      % obdy: FSP, DTF for obdy "spilled into" pixel, "spilled into" pixel
      obdy = zeros(obdct,3);

      % Identify unique "spilled into" points, paired with the spilling-into
      %  boundary point with:
      % (1) [if DEM(obdy)<=DEM(bdy)] the shallowest required flood height to
      %     spill over, i.e., min{bdy(:,3)} (using highest elevation
      %     for secondary sorting, if needed--rather the flood come from
      %     upstream than downstream)
      % OR
      % (2) [if DEM(obdy)>DEM(bdy)] the deepest available flood depth for obdy,
      %     i.e., max{fldht-bdy(:,3)-max{0,DEM(obdy)-DEM(bdy)}} (using highest
      %     elevation for secondary sorting, if needed--rather the flood come
      %     from upstream than downstream)
      for j = 1:obdct
         [r1,c1] = find(obdy_idx==unq(j));
         lr1 = length(r1);
         if(lr1==1)
            obdy(j,:) = [bdy(r1,1) bdy(r1,3)+max(0,obdy_fill(j,1)-bdy_fill(r1)) unq(j)];
            obdy_fill(j,2) = bdy_fill(r1);
         else
            f1 = find(fldht-bdy(r1,3)-max(0,obdy_fill(j)-bdy_fill(r1))==max(fldht-bdy(r1,3)-max(0,obdy_fill(j)-bdy_fill(r1))));
            if(length(f1)==1)
               obdy(j,:) = [bdy(r1(f1),1) bdy(r1(f1),3)+max(0,obdy_fill(j,1)-bdy_fill(r1(f1))) unq(j)];
               obdy_fill(j,2) = bdy_fill(r1(f1));
            else
               f2 = find(bdy_fill(r1(f1))==max(bdy_fill(r1(f1))),1);
               obdy(j,:) = [bdy(r1(f1(f2(1))),1) bdy(r1(f1(f2(1))),3)+max(0,obdy_fill(j,1)-bdy_fill(r1(f1(f2(1))))) unq(j)];
               obdy_fill(j,2) = bdy_fill(r1(f1(f2(1))));
            end
         end
      end
      clear obdy_idx
      % Sort "spilled into" points by increasing flood height and secondarily by decreasing
      %  elevation of associated boundary pixel in order to process the earlier flooded,
      %  more upstream pixels first
      tmp = [obdy obdy_fill];
      tmp = sortrows(tmp,[2 -5]);
      obdy = tmp(:,1:3);
      obdy_fill = tmp(:,4); % retain "spilled into" FIL column
      clear tmp

      %------------------------------
      % Forward-flood the "spilled into" pixels
      %------------------------------
      for j = 1:obdct
         pt = obdy(j,3);
         fldzn = obdy(j,2);
         tmp1 = double(obdy_fill(j));
         strpts0 = zeros(max(r,c),2); % pxid, FIL
         ct = 0;
         f0 = double(fdrdat(pt));
         while(1>0)
            % Break out if trajectory is self-intersecting (i.e., an Escher loop)
            if(sum(pt==strpts0(1:ct,1)))
               break;
            end
            strpts0(ct+1,:) = [pt double(tmp1)];
            r0 = ceil(pt/c);
            c0 = pt-(r0-1)*c;
            f = find(outflw==f0);
            % Break out if point flows out of the image
            if(or(sum(r0==[1 r]),sum(c0==[1 c])))
               if(or(or(and(c0==1,sum(f==[1,4,6])),and(c0==c,sum(f==[3,5,8]))),...
                     or(and(r0==1,sum(f==1:3)),and(r0==r,sum(f==6:8)))))
                  break;
               end
            end
            pt = pt+bdy0(f);
            tmp1 = double(fildat(pt));
            if(sum(pt==excl))
               % Break out if point flows into 'excl'
               strpts0(ct+2,:) = [pt double(tmp1)];
               ct = ct+1;
               break;
            elseif(tmp1==bg)
               % Break out if point flows out of bounds
               break;
            else
               tmp2 = double(flddat(pt));
            end
            if(and(tmp2,fldzn>=tmp2))
               % Break out if point flows back into the current
               % floodplain, into an area with lower depth-to-flood
               break;
            else
               f0 = double(fdrdat(pt));
               ct = ct+1;
            end
         end
         ct = ct+1;
         strpts0 = strpts0(1:ct,:);
         % set cap height at "spilled into" point elevation plus current max flood depth
         d0 = min(mxht0,obdy_fill(j)+fldht-fldzn); % maximum floodwater elevation on this path
         fldstg0 = d0-obdy_fill(j); % maximum flood depth on this path

         %------------------------------
         %------------------------------
         % Backfill flood the spilled-into stream path
         fldpln0 = zeros(sz,2);
         % fldpln0: floodplain pixel, flood depth (unadjusted)
         fldpln0(1:ct,:) = [strpts0(:,1) zeros(ct,1)];
         ct_tot0 = ct;
         for k = 1:ct
            %  col1 = reach pixel, col2 = flood height
            fldtmp = zeros(100*max(r,c),2); % hopefully this catches the biggest catchments
            fk0 = strpts0(k,1); % current boundary pixel
            fil1 = strpts0(k,2); % reference elevation
            v0 = fil1+fldstg0; % max elevation
            %------------------------
            % Determine if point is on an image boundary, and exclude 'out of bounds' neighbors
            r0 = ceil(fk0/c);
            c0 = fk0-(r0-1)*c;
            f0 = 1:8;
            bdy_rc = [sum(r0==[1 r]) sum(c0==[1 c])];
            if(sum(bdy_rc)==2)
               if(and(r0==1,c0==1)), f0 = [5,7,8];
               elseif(and(r0==1,c0==c)), f0 = [4,6,7];
               elseif(and(r0==r,c0==1)), f0 = [2,3,5];
               elseif(and(r0==r,c0==c)), f0 = [1,2,4];
               end
            elseif(bdy_rc(1)>0)
               if(r0==1), f0 = 4:8;
               elseif(r0==r), f0 = 1:5;
               end
            elseif(bdy_rc(2)>0)
               if(c0==1), f0 = [2,3,5,7,8];
               elseif(c0==c), f0 = [1,2,4,6,7];
               end
            end
            idx = fk0+bdy0(f0);
            lidx = length(idx);
            dat = zeros(lidx,2); % FDR, FIL

            % Retain neighbor points that are inflowing and in the new flood zone
            dat(:,1) = double(fdrdat(idx));
            dat(:,2) = double(fildat(idx));
            f1 = and(dat(:,1)==inflw(f0),and(dat(:,2)<=v0,dat(:,2)~=bg));
            idx = idx(f1);

            % Exclude neighbor points that are in the current "spillover" segment
            if(~isempty(idx))
               dat = dat(f1,:);
               f1 = ~ismember(idx,strpts0(:,1));
               idx = idx(f1);
            end

            ct1 = 0;
            if(~isempty(idx))
               ct1 = sum(f1);
               fldtmp(1:ct1,:) = [idx max(0,dat(f1,2)-fil1)];
            end
            %------------------------

            %------------------------
            % Generate reach for specific pixel along "spilled into" path
            ct0 = 0;
            while(ct0<ct1)
               ct0 = ct0+1;
               fk = fldtmp(ct0,1);
               %------------------------
               % Initialize boundary pixel-specific reach

               % Determine if point is on a boundary, and exclude 'out of bounds' neighbors
               r0 = ceil(fk/c);
               c0 = fk-(r0-1)*c;
               f0 = 1:8;
               bdy_rc = [sum(r0==[1 r]) sum(c0==[1 c])];
               if(sum(bdy_rc)==2)
                  if(and(r0==1,c0==1)), f0 = [5,7,8];
                  elseif(and(r0==1,c0==c)), f0 = [4,6,7];
                  elseif(and(r0==r,c0==1)), f0 = [2,3,5];
                  elseif(and(r0==r,c0==c)), f0 = [1,2,4];
                  end
               elseif(bdy_rc(1)>0)
                  if(r0==1), f0 = 4:8;
                  elseif(r0==r), f0 = 1:5;
                  end
               elseif(bdy_rc(2)>0)
                  if(c0==1), f0 = [2,3,5,7,8];
                  elseif(c0==c), f0 = [1,2,4,6,7];
                  end
               end
               idx = fk+bdy0(f0);
               lidx = length(idx);
               dat = zeros(lidx,2); % FDR, FIL

               % Retain neighbor points that are inflowing and in the new flood zone
               dat(:,1) = double(fdrdat(idx));
               dat(:,2) = double(fildat(idx));
               f1 = and(dat(:,1)==inflw(f0),and(dat(:,2)<=v0,dat(:,2)~=bg));
               idx = idx(f1);

               % Exclude neighbor points that are in the current "spillover" segment
               if(~isempty(idx))
                  dat = dat(f1,:);
                  f1 = ~ismember(idx,strpts0(:,1));
                  idx = idx(f1);
               end
               if(~isempty(idx))
                  ct2 = sum(f1);
                  fldtmp(ct1+(1:ct2),:) = [idx max(0,dat(f1,2)-fil1)];
                  ct1 = ct1+ct2;
               end
            end
            fldpln0(ct_tot0+(1:ct1),:) = fldtmp(1:ct1,:);
            ct_tot0 = ct_tot0+ct1;
         end
         fldpln0 = fldpln0(1:ct_tot0,:);
         %------------------------------
         %------------------------------
         % Integrate new flood information into 'fldpln_info'
         fldpln0(:,2) = fldpln0(:,2)+fldzn;
         if(any(ismember(fldpln0(:,1),fldpln_info(1:ct_tot,2))))
            fldtmp = double(flddat(fldpln0(1:ct_tot0,1)));
            if(any(fldtmp))
               fldpln1 = fldpln0(fldtmp>0,:);
               fldpln0 = fldpln0(~fldtmp,:);
               ct_tot0 = size(fldpln0,1);
               f = find(ismember(fldpln_info(1:ct_tot,2),fldpln1(:,1)));
               lf = length(f);
               for k = 1:lf
                  f1 = (fldpln_info(f,2)==fldpln1(k,1));
                  if(fldpln_info(f(f1),3)>fldpln1(k,2))
                     fldpln_info(f(f1),:) = [obdy(j,1) fldpln1(k,:)];
                     flddat(fldpln1(k,1)) = max(fldmn,fldpln1(k,2));
                  end
               end
            end
         end
         if(size(fldpln_info,1)<ct_tot+ct_tot0)
            fldpln_info = [fldpln_info; zeros(max(ct_tot0,10*max(r,c)),3)];
         end
         fldpln_info(ct_tot+(1:ct_tot0),:) = [obdy(j,1)*ones(ct_tot0,1) fldpln0];
         flddat(fldpln_info(ct_tot+(1:ct_tot0),2)) = max(fldmn,fldpln_info(ct_tot+(1:ct_tot0),3));
         ct_tot = ct_tot+ct_tot0;
      end

      %------------------------------
      % Identify new flood zone boundary pixels
      %------------------------------
      if(length(fldpln_bdy)<ct_tot)
         fldpln_bdy = [fldpln_bdy;ones(10*(ct_tot-length(fldpln_bdy)),1)];
      end
      f = find(fldpln_bdy(1:ct_tot));
      tmp = fldpln_info(f,2);
      tmp_r = ceil(tmp/c);
      tmp_c = tmp-(tmp_r-1)*c;
      f1 = find(~((tmp_r==1)+(tmp_r==r)+(tmp_c==1)+(tmp_c==c)));
      fbdr = find((tmp_r==1)+(tmp_r==r)+(tmp_c==1)+(tmp_c==c)>0);
      clear tmp_r tmp_c
      nfbdr = length(fbdr);
      tmp = tmp(f1);
      lf1 = length(f1);
      bdy_tot = zeros(lf1,8);
      for j = 1:lf1
         bdy_tot(j,:) = double(flddat(tmp(j)+bdy0)');
      end
      fldpln_bdy(f(f1(all(bdy_tot,2)))) = 0;
      clear bdy_tot
      for j = 1:nfbdr
         if(fldpln_bdy(f(fbdr(j))))
            fj = fldpln_info(f(fbdr(j)),2);
            % Determine if point is on a boundary, and exclude 'out of bounds' neighbors
            r0 = ceil(fj/c);
            c0 = fj-(r0-1)*c;
            f0 = 1:8;
            bdy_rc = [sum(r0==[1 r]) sum(c0==[1 c])];
            if(sum(bdy_rc)==2)
               if(and(r0==1,c0==1)), f0 = [5,7,8];
               elseif(and(r0==1,c0==c)), f0 = [4,6,7];
               elseif(and(r0==r,c0==1)), f0 = [2,3,5];
               elseif(and(r0==r,c0==c)), f0 = [1,2,4];
               end
            elseif(bdy_rc(1)>0)
               if(r0==1), f0 = 4:8;
               elseif(r0==r), f0 = 1:5;
               end
            elseif(bdy_rc(2)>0)
               if(c0==1), f0 = [2,3,5,7,8];
               elseif(c0==c), f0 = [1,2,4,6,7];
               end
            end
            idx = fj+bdy0(f0);
            if(all(ismember(idx,fldpln_info(1:ct_tot,2))))
               fldpln_bdy(f(fbdr(j))) = 0;
            end
         end
      end
      % col1 = reference stream pixel, col2 = boundary pixel, col3 = flood depth required to get to boundary pixel
      bdy = fldpln_info(fldpln_bdy(1:ct_tot)>0,:);
      % Sort boundary points by increasing flood height and secondarily by decreasing
      %  elevation in order to process the earlier flooded, upstream pixels first
      lf = size(bdy,1);
      tmp = [bdy zeros(lf,1)];
      tmp(:,4) = double(fildat(bdy(:,2)));
      tmp = sortrows(tmp,[3 -4]);
      bdy = tmp(:,1:3); clear tmp
      bdct = size(bdy,1);
      
      if(ssflg)
         % check for additional possible spillover when calculating the
         %  'steady state' floodplain
         if(any(bdy(:,3)<fldht))
            if(ct_tot2<ct_tot)
               bdy = bdy(bdy(:,3)<fldht,:);
               bdct = size(bdy,1);
            else
               spill_flg = 0;
            end
         else
            spill_flg = 0;
         end
      else
         spill_flg = 0;
      end
   end

   t = toc;
   proctime = [proctime; [fldht t]];
   save([segdr,outf,'_tmp.mat'],'fldpln_info','fldpln_bdy','bdy','bdct','fldht','ct_tot','proctime','-v7.3');
   disp(sprintf('%d of %d iterations completed',g,itnum))
   disp(sprintf('SegID %d, iteration %d completed in %g seconds',seg0,g,t));
   tic
end
fclose('all');
end % exflg 'if' statement close
%return;

% Write final output file
if(~exist([segdr,outf,'.mat'],'file'))
   load([segdr,outf,'_tmp.mat'],'ct_tot','fldpln_info');
   fldpln_info = fldpln_info(1:ct_tot,:);
   fldpln_info(:,3) = max(fldmn,fldpln_info(:,3));
   
   % fldpln:
   %  col1 = reference stream pixel, col2 = floodplain pixel, col3 = flood height,
   %  col4 = flood height + fill depth
   fldpln = [fldpln_info zeros(ct_tot,1)];
   clear fldpln_info
   header = {'refernce stream pixel','floodplain pixel','flood height','flood height + fill depth'};
   save([segdr,outf,'.mat'],'fldpln','header','-v7.3');
   
   %------------------
   % open DEM
   demf = inp.dem;
   deminfo = readbilheader(demf);
   fmt = deminfo.fmt;
   numbyt = deminfo.nbyts;
   if(deminfo.bytord)
      dem = fopen(demf,'rb','b');
   else
      dem = fopen(demf,'rb');
   end
   %------------------
   
   for j = 1:ct_tot
      dat1 = double(fildat(fldpln(j,2)));
      fseek(dem,(fldpln(j,2)-1)*numbyt,-1);
      dat2 = fread(dem,1,fmt);
      fldpln(j,4) = fldpln(j,3)+dat1-max(0,dat2);
   end
   fclose('all');
   save([segdr,outf,'.mat'],'fldpln','header','-v7.3');
end
clear all
