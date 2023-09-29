%------------------------
% function make_stream_map(facf,thr,strf)
%------------------------
% This function takes a FAC BIL file, and creates a stream network ERDAS BIL
% using the input threshold value
%
% INPUT
%  facf = path and filename for FAC ERDAS BIL file
%  thr = minimum number of pixels to use for stream pixel designation
%  strf = path and filename for output STR ERDAS BIL file ('uint8' format)
%------------------------
function make_stream_map(facf,thr,strf)
%------------------------

%=============================
% Read spatial parameters from FAC
%-------------------------
facinfo = readbilheader(facf);
r = facinfo.r; % rows
c = facinfo.c; % cols
bytord = facinfo.bytord; % byte order; if 1 then open file with 'b' option
fmt = facinfo.fmt; % FAC BIL file data format
%=============================

if(bytord)
   fac = fopen(facf,'rb','b');
else
   fac = fopen(facf,'rb');
end;
str=fopen(strf,'wb');

zro=zeros(c,1);
for j=1:r
   ln=fread(fac,c,fmt);
   if(any(ln>=thr))
      zro1=zro;
      zro1(ln>=thr)=1;
      fwrite(str,zro1,'uint8');
   else
      fwrite(str,zro,'uint8');
   end;
end;
fclose('all');

%=============================
% Write header file
%------------------------
fhdr0 = fopen([facf(1:end-3),'hdr'],'rt');
fhdr = fopen([strf(1:end-3),'hdr'],'wt');

txt0 = fgetl(fhdr0);
if(any(txt0==':')) % ERDAS BIL
   fhdr0 = fopen([facf(1:end-3),'hdr'],'rt');
   while(~feof(fhdr0))
      txt0 = fgetl(fhdr0);
      if(any(txt0==':'))
         f = find(txt0==':');
         txt = txt0(f+1:end);
         if(isequal(txt0(1:f-1),'DATATYPE'))
            f1 = find(txt~=' ',1,'first');
            str = [txt0(1:f+f1-1),'U8',10];
            fprintf(fhdr,str);
         elseif(isequal(txt0(1:f-1),'BYTE_ORDER'))
            f1 = find(txt~=' ',1,'first');
            str = [txt0(1:f+f1-1),'NA',10];
            fprintf(fhdr,str);
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
   bandrowbytes = num2str(bandrowbytes/(nbits/8));
   totalrowbytes = num2str(totalrowbytes/(nbits/8));
   
   fhdr0 = fopen([facf(1:end-3),'hdr'],'rt');
   while(~feof(fhdr0))
      txt0 = fgetl(fhdr0);
      f = find(txt0==' ',1,'first');
      txt = txt0(f+1:end);
      if(isequal(txt0(1:f-1),'NBITS'))
         f1 = find(txt~=' ',1,'first');
         str = [txt0(1:f+f1-1),'8',10];
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
%=============================
disp(sprintf('Binary stream network map creation completed'));
