%-----------------------
% function custom_floodplain_map(dr0,r,c,h,dh,fldmn,seg_ids,seg_dpth,segdr,outf,reff*)
%-----------------------
% This function merges segment-specific flood zone maps into a single
%  floodplain map and database file.
%
% INPUT
%  r = number of rows in full extent
%  c = number of columns in full extent
%  h = maximum flood depth (for file identification)
%  dh = flood depth step size (for file identification)
%  fldmn = minimum flood depth (to distinguish from background zeros)
%  seg_ids = n-vector of segment IDs with floodplains to be merged
%  seg_dpth = 1- or n-vector of segment-specific depths for floodplain
%     merging; n-vector entries correspond with segments listed in
%     'seg_ids'; 1-vector indicates a uniform max flood depth; enter the
%     empty set ([]) to use the maximum flood depth
%  segdr = path where FLDPLN segment files are stored
%  outf = path & filename to be used for output BIL file
%
% OPTIONAL INPUT*
%  reff = path & filename for BIL with same spatial info (for reference
%         during header file write)
%-----------------------
function custom_floodplain_map(r,c,h,dh,fldmn,seg_ids,seg_dpth,segdr,outf,varargin)
%-----------------------

fmt = 'single'; % change if using integer data, to save memory
numbyt = 4;

num=length(seg_ids);
if(~isempty(seg_dpth))
   if(~ismember(length(seg_dpth),[1,num]))
      error('Improperly specified ''seg_dpth'' variable');
   end
end
%----------------------------
% Make strings to use in file names
hs=num2str(h);
f=find(hs=='.');
if(~isempty(f))
   hs(f)='p';
end
dhs=num2str(dh);
f=find(dhs=='.');
if(~isempty(f))
   dhs(f)='p';
end
%----------------------------
if(length(seg_dpth)==1)
   seg_dpth=seg_dpth*ones(length(seg_ids),1);
end

%--------------------------
% Initialize FLD with the floodplain from the first segment
zro = zeros(c,1);
seg = int2str(seg_ids(1));
segf = ['h',hs,'_dh',dhs,'_seg',seg];
if(exist([segdr,segf,'.mat'],'file'))
   load([segdr,segf]);
   % col2 = floodplain pixel, col3 = flood depth
   fldpln = fldpln(:,2:3);
else
   segf = ['h',hs,'_dh',dhs,'_seg',seg,'_tmp'];
   load([segdr,segf],'fldpln_info','ct_tot');
   % col2 = floodplain pixel, col3 = flood depth
   fldpln = fldpln_info(1:ct_tot,2:3);
end
fldpln = fldpln(fldpln(:,2)<=seg_dpth(1),:);
fld = fopen(outf,'wb','b');
r0 = ceil(fldpln(:,1)/c);
c0 = fldpln(:,1)-(r0-1)*c;
for k=1:r
   if(any(r0==k))
      f = find(r0==k);
      zro1 = zro;
      zro1(c0(f)) = max(fldmn,fldpln(f,2));
      fwrite(fld,zro1,fmt);
   else
      fwrite(fld,zro,fmt);
   end
end
fclose(fld);
fld = fopen(outf,'rb+','b');
if(num>1)
   disp(sprintf('1 of %d segments completed',num));
end
%--------------------------

%--------------------------
% Overwrite FLD as needed for other segment floodplains
for j=2:num
   fseek(fld,0,-1); % rewind the floodplain file
   seg=int2str(seg_ids(j));
   segf=['h',hs,'_dh',dhs,'_seg',seg];
   if(exist([segdr,segf,'.mat'],'file'))
      load([segdr,segf]);
      % col2 = floodplain pixel, col3 = flood depth
      fldpln = fldpln(:,2:3);
   else
      segf = ['h',hs,'_dh',dhs,'_seg',seg,'_tmp'];
      load([segdr,segf],'fldpln_info','ct_tot');
      % col2 = floodplain pixel, col3 = flood depth
      fldpln = fldpln_info(1:ct_tot,2:3);
   end
   fldpln = fldpln(fldpln(:,2)<=seg_dpth(j),:);
   r0 = ceil(fldpln(:,1)/c);
   c0 = fldpln(:,1)-(r0-1)*c;
   for k = 1:r
      flg = 0;
      if(~any(r0==k))
         fseek(fld,c*numbyt,0); % move the file pointer down a line
      else
         ln0 = fread(fld,c,fmt); % read a line from the floodplain file
         f = find(r0==k);
         f1 = find(ln0(c0(f)));
         if(~isempty(f1))
            if(any(ln0(c0(f(f1)))>max(fldmn,fldpln(f(f1),2))))
               flg = 1;
               ln0(c0(f(f1))) = min(ln0(c0(f(f1))),max(fldmn,fldpln(f(f1),2)));
            end
         end
         if(length(f)~=length(f1))
            flg = 1;
            f2 = find(~ln0(c0(f)));
            ln0(c0(f(f2))) = max(fldmn,fldpln(f(f2),2));
         end
         if(flg)
            fseek(fld,-c*numbyt,0);
            fwrite(fld,ln0,fmt);
         end
      end
   end
   disp(sprintf('%d of %d segments completed',j,num));
end
fclose(fld);

%=============================
% Write header file (if reference file specified)
%-----------------------------
if(nargin==10)
   reff = varargin{1};
   fhdr0 = fopen([reff(1:end-3),'hdr'],'rt');
   fhdr = fopen([outf(1:end-3),'hdr'],'wt');

   txt0 = fgetl(fhdr0);
   while(~feof(fhdr0))
      if(any(txt0==':'))
         f = find(txt0==':');
         txt = txt0(f+1:end);
         if(isequal(txt0(1:f-1),'DATATYPE'))
            f1 = find(txt~=' ',1,'first');
            str = [txt0(1:f+f1-1),'F32',10];
            fprintf(fhdr,str);
         elseif(isequal(txt0(1:f-1),'BYTE_ORDER'))
            f1 = find(txt~=' ',1,'first');
            str = [txt0(1:f+f1-1),'INTEL',10];
            fprintf(fhdr,str);
         else
            fprintf(fhdr,[txt0,10]);
         end
      else
         fprintf(fhdr,[txt0,10]);
      end
      txt0 = fgetl(fhdr0);
   end
end
fclose('all');
%=============================
%-----------------------------
% Write BIL header file (ERDAS or ESRI, depending on input)
%  (Input reference file must be specified)
%-----------------------------
if(nargin==10)
   fhdr0 = fopen([reff(1:end-3),'hdr'],'rt');
   txt0 = fgetl(fhdr0);
   biltyp = 0; % ESRI
   if(any(txt0==':'))
      biltyp = 1; % ERDAS
   end
   
   fhdr0 = fopen([reff(1:end-3),'hdr'],'rt');
   fhdr = fopen([outf(1:end-3),'hdr'],'wt');
   if(biltyp) % ERDAS BIL
      fmts0 = {'F32','F64','U8','S8','U16','S16','U32','S32'};
      fmts = {'single','double','uint8','int8','uint16','int16','uint32','int32'};
      %-----------------------------
      % Need to test the ERDAS code
      %-----------------------------
      while(~feof(fhdr0))
         txt0 = fgetl(fhdr0);
         if(any(txt0==':'))
            f = find(txt0==':');
            txt = txt0(f+1:end);
            if(isequal(txt0(1:f-1),'DATATYPE'))
               f1 = find(txt~=' ',1,'first');
               str = [txt0(1:f+f1-1),fmts0{ismember(fmts,fmt)},10];
               fprintf(fhdr,str);
            elseif(isequal(txt0(1:f-1),'BYTE_ORDER'))
               f1 = find(txt0(f+1:end)~=' ',1,'first');
               str = [txt0(1:f+f1-1),'INTEL',10];
               fwrite(fhdr,str);
            else
               fprintf(fhdr,[txt0,10]);
            end
         else
            fprintf(fhdr,[txt0,10]);
         end
      end
   else % ESRI BIL
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
            str = [txt0(1:f+f1-1),num2str(8*numbyt),10];
            fprintf(fhdr,str);
         elseif(isequal(txt0(1:f-1),'BANDROWBYTES'))
            f1 = find(txt~=' ',1,'first');
            bandrowbytes = num2str(c*numbyt);
            str = [txt0(1:f+f1-1),bandrowbytes,10];
            fprintf(fhdr,str);
         elseif(isequal(txt0(1:f-1),'TOTALROWBYTES'))
            f1 = find(txt~=' ',1,'first');
            totalrowbytes = num2str(c*numbyt);
            str = [txt0(1:f+f1-1),totalrowbytes,10];
            fprintf(fhdr,str);
         elseif(isequal(txt0(1:f-1),'PIXELTYPE'))
            f1 = find(txt~=' ',1,'first');
            if(ismember(fmt,{'single','double'}))
               pxtyp = 'FLOAT';
            elseif(isequal(fmt(1:4),'uint'))
               pxtyp = 'UNSIGNEDINT';
            elseif(isequal(fmt(1:3),'int'))
               pxtyp = 'SIGNEDINT';
            else
               error('Unrecognized data type');
            end
            str = [txt0(1:f+f1-1),pxtyp,10];
            fprintf(fhdr,str);
         elseif(isequal(txt0(1:f-1),'NODATA'))
            f1 = find(txt~=' ',1,'first');
            str = [txt0(1:f+f1-1),'0',10];
            fprintf(fhdr,str);
         else
            fprintf(fhdr,[txt0,10]);
         end
      end
   end
   fclose('all');
   if(exist([reff(1:end-3),'prj'],'file'))
      system(['cp ',reff(1:end-3),'prj ',outf(1:end-3),'prj']);
   end
end
%=============================
