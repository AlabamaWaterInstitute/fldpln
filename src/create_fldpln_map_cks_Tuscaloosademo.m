function create_fldpln_map_cks_Tuscaloosademo()

aoi = 'Tuscaloosa';
dr0 = ['/',aoi,'/'];

%=============================
% Initialize file names
%-----------------------------
fdrf = [dr0,'bil/Tuscaloosa_Fdr.bil'];
segdr = [dr0,'segment_files/'];
shpf = [dr0,'vector/str_segid_polyline_Dissolve.shp'];
%if(~exist(shpf,'file'))
%   shpf = [dr0,'vector/',aoi,'_str_segid_nogen_dissolve.shp'];
%end;
%=============================

%========================================================
%========================================================
% Read spatial parameters from FDR
%-------------------------
fdrinfo = readbilheader(fdrf);
r = fdrinfo.r; % rows
c = fdrinfo.c; % cols
%=============================

dr = dir(segdr);
nm0 = dr(3).name; % segment 0 = total merged SLIE for AOI
f = find(nm0=='_');
nm0 = nm0(1:f(2)-1);
h = eval(nm0(2:f(1)-1)); % max DTF for SLIE
hs = int2str(h);
dhs = nm0(f(1)+3:end);
dhs(dhs=='p') = '.';
dh = eval(dhs); % DTF step size used for FLDPLN
dhs(dhs=='.') = 'p';

%outf = [dr0,'bil/fldpln_map_h',hs,'_dh',dhs,'_',aoi,'_2p5segs9_10_11_12.bil'];
%outf = [dr0,'bil/fldpln_map_h',hs,'_dh',dhs,'_',aoi,'_2p5.bil'];
outf = [dr0,'bil/fldpln_map_h',hs,'_dh',dhs,'_',aoi,'.bil'];

fldmn = 0.01;

%load([dr0,'mat/seg_info.mat']); % seg_info
%seg_list = (1:size(seg_info,1))';

%seg_list = sort(readshapefield(shpf,{'GRID_CODE'}));
seg_list = sort(readshapefield(shpf,{'grid_code'}));
num = length(seg_list);

custom_floodplain_map(r,c,h,dh,fldmn,seg_list,h,segdr,outf,fdrf)
%custom_floodplain_map(r,c,h,dh,fldmn,seg_list,2.5,segdr,outf,fdrf)
return;
hts = h*ones(num,1);
%hts(ismember(seg_list,9:12)) = 2.5; % segment 9 has too much overland flooding
custom_floodplain_map(r,c,h,dh,fldmn,seg_list,hts,segdr,outf,fdrf)
%custom_floodplain_map(r,c,h,dh,fldmn,0,h,segdr,outf,fdrf)

%==================================
return;






load([dr0,'neuse_segments.mat'],'seg_list');
load([dr0,'problem_segments.mat'],'prob_segs','fill_less','fill_none');

if(1<0)
seg_ids=seg_list(find(~ismember(seg_list,fill_none)));
len=length(seg_ids);
seg_dpth=5*ones(len,1);
fill_less(:,2)=-fill_less(:,2);
segs0=[prob_segs; fill_less];
n=size(segs0,1);
for j=1:n
   if(segs0(j,1)~=99)
      f=find(seg_ids==segs0(j,1));
      seg_dpth(f)=seg_dpth(f)+segs0(j,2);
   end;
end;
%[seg_ids seg_dpth]
%size([seg_ids seg_dpth])
end;
seg_ids=99;
seg_dpth=5;

fldmn=0.01;
outf='fldpln_map_neuse_custom_1oct10';

%-----------------------

fmt='single'; % change if using integer data, to save memory
numbyt=4;

num=length(seg_ids);

%--------------------------
% Initialize FLD with the floodplain from the first segment
if(1<0)
zro=zeros(c,1);
seg=int2str(seg_ids(1));
hs=int2str(seg_dpth(1));
segf=['h',hs,'_dh1_seg',seg];
load([dr0,'segment_files_',hs,'m/',segf]);
% col2 = floodplain pixel, col3 = flood depth
fldpln=fldpln(:,2:3);
fldpln=fldpln(find(fldpln(:,2)<=seg_dpth(1)),:);
fld=fopen([dr0,'bil/',outf,'.bil'],'wb','b');
r0=floor(fldpln(:,1)/c)-(~mod(fldpln(:,1),c))+1;
c0=fldpln(:,1)-(r0-1)*c;
for k=1:r
   if(any(r0==k))
      f=find(r0==k);
      zro1=zro;
      zro1(c0(f))=max(fldmn,fldpln(f,2));
      fwrite(fld,zro1,fmt);
   else
      fwrite(fld,zro,fmt);
   end;
end;
fclose(fld);
if(num>1)
   disp(sprintf('1 of %d segments completed',num));
end;
end;
fld=fopen([dr0,'bil/',outf,'.bil'],'rb+','b');
%--------------------------

%--------------------------
% Overwrite FLD as needed for other segment floodplains
%for j=2:num
for j=1:num
   fseek(fld,0,-1); % rewind the floodplain file
   seg=int2str(seg_ids(j));
   hs=int2str(seg_dpth(j));
   segf=['h',hs,'_dh1_seg',seg];
   load([dr0,'segment_files_',hs,'m/',segf]);
   % col2 = floodplain pixel, col3 = flood depth
   fldpln=fldpln(:,2:3);
   fldpln=fldpln(find(fldpln(:,2)<=seg_dpth(j)),:);
   r0=floor(fldpln(:,1)/c)-(~mod(fldpln(:,1),c))+1;
   c0=fldpln(:,1)-(r0-1)*c;
   for k=1:r
      flg=0;
      if(~any(r0==k))
         fseek(fld,c*numbyt,0); % move the file pointer down a line
      else
         ln0=fread(fld,c,fmt); % read a line from the floodplain file
         f=find(r0==k);
         f1=find(ln0(c0(f)));
         if(~isempty(f1))
            if(any(ln0(c0(f(f1)))>max(fldmn,fldpln(f(f1),2))))
               flg=1;
               ln0(c0(f(f1)))=min(ln0(c0(f(f1))),max(fldmn,fldpln(f(f1),2)));
            end;
         end;
         if(length(f)~=length(f1))
            flg=1;
            f2=find(~ln0(c0(f)));
            ln0(c0(f(f2)))=max(fldmn,fldpln(f(f2),2));
         end;
         if(flg)
            fseek(fld,-c*numbyt,0);
            fwrite(fld,ln0,fmt);
         end;
      end;
   end;
   disp(sprintf('%d of %d segments completed',j,num));
end;
fclose(fld);


return;
%==================================

%dos(['rename ',dr0,'segment_files_',hs,'m segment_files']);
%dos(['rename ',dr0,'segment_files segment_files_',hs,'m']);
