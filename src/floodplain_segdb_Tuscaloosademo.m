%--------------------------------
% function floodplain_segdb_cks_walnut()
%--------------------------------
% This function creates a stream network and applies the FLDPLN model to
%  create an SLIE.  Modify parameters in the code to customize.
%--------------------------------
%--------------------------------
function floodplain_segdb_Tuscaloosademo()
%--------------------------------

%=============================
% Enter a home directory
%-----------------------------

dr1 = 'Tuscaloosa/';
dr0 = ['../',dr1];

%-----------------------------
% Initialize flood parameters & the I/O file names
%-----------------------------
fldmn = .01; % minimum flood depth (in vertical DEM units; typically meters or feet)
%-----------
% iteration step size (in vertical DEM units)
%dh = 0.01; 
%dh = 0.1;
%dh = 0.25;
dh = 0.5;
%dh = 1;
%-----------
fldmx = 3; % maximum flood depth
mxht = 0; % 'not-to-exceed' WSE value (not typically used; default 0 means no action)

demf = [dr0,'bil/Tuscaloosa_10m.bil'];
filf = [dr0,'bil/Tuscaloosa_Fil.bil'];
fdrf = [dr0,'bil/Tuscaloosa_Fdr.bil'];
facf = [dr0,'bil/Tuscaloosa_Fac.bil'];
strf0 = 'str.bil';
strf = [dr0,'bil/',strf0];
segf = [dr0,'bil/str_segid.bil'];
matf = [dr0,'mat/seg_info.mat'];
%shpf = [dr0,'vector/Tuscaloosa_start_pt.shp'];
%shpf_all = [dr0,'vector/str_segid_polyline_Dissolve.shp'];

fdrinfo = readbilheader(fdrf);
pxsz = 10
% pxsz = fdrinfo.pxszx;

%=============================
% Identify synthetic stream network
%-----------------------------
mi2px = 5280*.3048/pxsz; % horizontal units assumed to be meters (otherwise modify as needed) 
seglen = round(5*mi2px); % 5 miles maximum segment length
sqmi2px = (5280^2)*(.3048^2)/(pxsz^2); % horizontal units assumed to be meters (otherwise modify as needed)
facthr = round(25*sqmi2px); % 25 sq. mi minimum catchment size (min for segment break)
strthr = round(70*sqmi2px); % 70 sq. mi minimum catchment size (min for defining a stream)

if(~exist(strf,'file'))
   make_stream_map(facf,strthr,strf);
end

if(~exist([dr0,'mat'],'dir'))
   mkdir([dr0,'mat']);
end

if(~exist(matf,'file'))
create_segdb_facthr_maxlen(strthr,facthr,seglen,fdrf,facf,segf,matf);
% create_segdb_facthr_maxlen_seedpts(shpf,facthr,seglen,fdrf,facf,segf,matf);
end
%return;
% Columns of 'seg_info' (row corresponds with segment ID):
% 1) start pixel
% 2) end pixel
% 3) FAC value @ start
% 4) FAC value @ end
% 5) segment length
% 6) network flag (0 if segment exits the image; 1 otherwise)
% 7) immediately downstream segment ID (=0 if segment flows out of study area)
%whos('-file',matf);
%return;
%-----------------------------
% Define background value
[mn,mx,bg] = scanbil(filf);
%-----------------------------

load(matf,'seg_info');
num = size(seg_info,1);
seg_list = (1:num);
%seg_list0 = sortrows(readshapefield(shpf_all,{'GRID_CODE','StrID'}),2);
%seg_list = sort(readshapefield(shpf_all,{'grid_code'}));
num = length(seg_list);
%return;
%------------------
% Identify completed segments
if(1>0)
dr = dir([dr0,'segment_files']);
len = size(dr,1);
done_segs = zeros(len-2,1);
ct = 0;
for j = 3:len
   f2 = find(dr(j).name=='.');
   if(dr(j).name(f2-1)=='p')
      load([dr0,'segment_files/',dr(j).name],'fldht')
      if(fldht==fldmx)
         f1 = find(dr(j).name=='g',1,'last');
         f2 = find(dr(j).name=='_',1,'last');
         done_segs(ct+1) = str2double(dr(j).name(f1+1:f2-1));
         ct = ct+1;
      end
   end
end
done_segs = done_segs(1:ct);
seg_list = intersect(seg_list,setxor(seg_list,done_segs));
num = length(seg_list);
%return;
end
%------------------

%------------------
if(num>0)
   inp = struct;
   inp.bildr = [dr0,'bil/'];
   inp.segdr = [dr0,'segment_files/'];
   inp.dem = demf;
   inp.fil = filf;
   inp.fdr = fdrf;
   inp.seg = matf;
   inp.segf = segf;
   
   %-------------------------
   % Make sure FIL & FDR have same row-col dimensions
   f = inp.fil;
   filinfo = readbilheader(f);
   r = filinfo.r;
   c = filinfo.c;

   f = inp.fdr;
   filinfo = readbilheader(f);
   r0 = filinfo.r;
   c0 = filinfo.c;
   
   if(or(r0~=r,c0~=c))
      error('FIL & FDR do not have the same row-col dimensions');
      % Write this message to a file in the home directory
   end
   %-------------------------
   
   %=============================
   % Read FIL file specs & data
   %-------------------------
   f = inp.fil;
   filinfo = readbilheader(f);
   fmt = filinfo.fmt;
   filtyp = filinfo.bytord;

   if(filtyp)
      fid = fopen(f,'rb','b');
   else
      fid = fopen(f,'rb');
   end
   fildat = zeros(r*c,1,fmt);
   fildat(1:r*c) = fread(fid,inf,fmt);
   fclose(fid);
   %=============================
   
   %=============================
   % Read FDR data
   %-------------------------
   f = inp.fdr;
   filinfo = readbilheader(f);
   fmt = filinfo.fmt;
   filtyp = filinfo.bytord;

   if(filtyp)
      fid = fopen(f,'rb','b');
   else
      fid = fopen(f,'rb');
   end
   fdrdat = zeros(r*c,1,fmt);
   fdrdat(1:r*c) = fread(fid,inf,fmt);
   fclose(fid);
   %=============================
   
   
if(1>0)
   g = gcp('nocreate');
   if(~isempty(g))
      p = gcp;
   else
      p = parpool;
   end
   for j = 1:num
      pf(j) = parfeval(p,@fldpln_model_v5ram,1,seg_list(j),inp,fildat,fdrdat,fldmn,fldmx,dh,mxht,1,bg);
   end
%   wait(pf)
   pfresults = cell(num,1);
   for idx = 1:num
  %    fetchNext blocks until next results are available.
     [completedIdx,value] = fetchNext(pf);
     pfresults{completedIdx} = value;
     fprintf('Segment %d completed/n', completedIdx);
   end
%   delete(gcp('nocreate'))
end
return;
%   parfor j = 1:num
%      fldpln_model_v5ram(seg_list(j),inp,fildat,fdrdat,fldmn,fldmx,dh,mxht,1,bg);
%   end
%   return;
end

