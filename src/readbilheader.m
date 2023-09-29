%========================================================
% function filinfo = readbilheader(bilfil)
%========================================================
% This function reads header information from an Erdas BIL (direct read) or
%  an ESRI ArcMap 10.0+ BIL file and saves it in the structure 'filinfo'.
%
% INPUT
%  bilfil = path and file name for BIL file (must have an associated HDR
%           file containing header information to be read by this function)
%
% OUTPUT
%  filinfo = structure with the following fields:
%     r = number of rows
%     c = number of columns
%     bands = number of bands
%     ulx = upper left X coordinate
%     ulxstr = verbatim character string for 'ulx' value
%     uly = upper left Y coordinate
%     ulystr = verbatim character string for 'uly' value
%     pxszx = pixel width
%     pxszxstr = verbatim character string for 'pxszx' value
%     pxszy = pixel height
%     pxszystr = verbatim character string for 'pxszy' value
%     fmt = data format (MATLAB style)
%     bytord = 1 if Big Endian (open file with 'b' option), 0 if Little Endian
%     nbits = number of bits per value (Arc only; could be inferred from 'fmt')
%     nbyts = number of bytes per value
%     biltyp = 'ESRI' or 'ERDAS'
%     nodata = NODATA value (ESRI exclusive?)
%     nodatastr = verbatim character string for 'nodata' value
%     hdrascii = ASCII integer array for entire header file
%     hdrchar = character array for entire header file
%========================================================

%========================================================
function filinfo = readbilheader(bilfil)
%========================================================
fhdr = [bilfil(1:end-3),'hdr'];
fid = fopen(fhdr,'rt');

filinfo = struct;

txt0 = fgetl(fid);
fclose(fid);
fid = fopen(fhdr,'rt');
if(all(txt0~=':')) % ESRI BIL
   filinfo.biltyp = 'ESRI';
   bytords = {'I','M'};

   while(~feof(fid))
      txt0 = fgetl(fid);
      if(length(txt0)>4)
         f = find(txt0==' ',1,'first');
         txt = txt0(f+1:end);
         if(strcmpi(txt0(1:f-1),'BYTEORDER'))
            str = txt(txt~=' ');
            filinfo.bytord = find(ismember(bytords,str))-1;
         elseif(strcmpi(txt0(1:f-1),'NROWS'))
            str = txt(ismember(txt,'0123456789'));
            filinfo.r = str2double(str);
         elseif(strcmpi(txt0(1:f-1),'NCOLS'))
            str = txt(ismember(txt,'0123456789'));
            filinfo.c = str2double(str);
         elseif(strcmpi(txt0(1:f-1),'NBANDS'))
            str = txt(ismember(txt,'0123456789'));
            filinfo.bands = str2double(str);
         elseif(strcmpi(txt0(1:f-1),'NBITS'))
            str = txt(ismember(txt,'0123456789'));
            filinfo.nbits = str2double(str);
         elseif(strcmpi(txt0(1:f-1),'PIXELTYPE'))
            datfmt0 = txt(txt~=' ');
         elseif(strcmpi(txt0(1:f-1),'ULXMAP'))
            filinfo.ulxstr = txt(ismember(txt,'-.e0123456789'));
            filinfo.ulx = str2double(filinfo.ulxstr);
         elseif(strcmpi(txt0(1:f-1),'ULYMAP'))
            filinfo.ulystr = txt(ismember(txt,'-.e0123456789'));
            filinfo.uly = str2double(filinfo.ulystr);
         elseif(strcmpi(txt0(1:f-1),'XDIM'))
            filinfo.pxszxstr = txt(ismember(txt,'-.e0123456789'));
            filinfo.pxszx = str2double(filinfo.pxszxstr);
         elseif(strcmpi(txt0(1:f-1),'YDIM'))
            filinfo.pxszystr = txt(ismember(txt,'-.e0123456789'));
            filinfo.pxszy = str2double(filinfo.pxszystr);
         elseif(strcmpi(txt0(1:f-1),'NODATA'))
            filinfo.nodatastr = txt(ismember(txt,'-.e0123456789'));
            filinfo.nodata = str2double(filinfo.nodatastr);
         end;
      end;
   end;
   switch datfmt0
      case 'SIGNEDINT'
         filinfo.fmt = ['int',num2str(filinfo.nbits)];
      case 'UNSIGNEDINT'
         filinfo.fmt = ['uint',num2str(filinfo.nbits)];
      case 'FLOAT'
         if(filinfo.nbits==32)
            filinfo.fmt = 'single';
         elseif(filinfo.nbits==64)
            filinfo.fmt = 'double';
         end;
   end;
   filinfo.nbyts = filinfo.nbits/8;
else % Erdas BIL
   filinfo.biltyp = 'ERDAS';
   fmts0 = {'F32','F64','U8','S8','U16','S16','U32','S32'};
   fmts = {'single','double','uint8','int8','uint16','int16','uint32','int32'};
   bytords = {'MOTOROLA','INTEL'};
   % NOTE:  ERDAS has these backwards, or at least they say 'INTEL' but the
   %        file reads properly as big endien

   while(~feof(fid))
      txt0 = fgetl(fid);
      if(any(txt0==':'))
         f = find(txt0==':');
         txt = txt0(f+1:end);
         if(strcmpi(txt0(1:f-1),'BANDS'))
            str = txt(ismember(txt,'0123456789'));
            filinfo.bands = str2double(str);
         elseif(strcmpi(txt0(1:f-1),'ROWS'))
            str = txt(ismember(txt,'0123456789'));
            filinfo.r = str2double(str);
         elseif(strcmpi(txt0(1:f-1),'COLS'))
            str = txt(ismember(txt,'0123456789'));
            filinfo.c = str2double(str);
         elseif(strcmpi(txt0(1:f-1),'DATATYPE'))
            str = txt(txt~=' ');
            f = find(ismember(fmts0,str));
            filinfo.fmt = fmts{f};
            filinfo.nbyts = 1+ismember(f,5:6)+3*ismember(f,[1,7:8])+7*(f==2);
         elseif(strcmpi(txt0(1:f-1),'BYTE_ORDER'))
            str = txt(txt~=' ');
            filinfo.bytord = find(ismember(bytords,str))-1;
         elseif(strcmpi(txt0(1:f-1),'UL_X_COORDINATE'))
            filinfo.ulxstr = txt(ismember(txt,'-.e0123456789'));
            filinfo.ulx = str2double(filinfo.ulxstr);
         elseif(strcmpi(txt0(1:f-1),'UL_Y_COORDINATE'))
            filinfo.ulystr = txt(ismember(txt,'-.e0123456789'));
            filinfo.uly = str2double(filinfo.ulystr);
         elseif(strcmpi(txt0(1:f-1),'PIXEL_WIDTH'))
            filinfo.pxszxstr = txt(ismember(txt,'-.e0123456789'));
            filinfo.pxszx = str2double(filinfo.pxszxstr);
         elseif(strcmpi(txt0(1:f-1),'PIXEL_HEIGHT'))
            filinfo.pxszystr = txt(ismember(txt,'-.e0123456789'));
            filinfo.pxszy = str2double(filinfo.pxszystr);
         end;
      end;
   end;
end;
fclose(fid);

fid = fopen(fhdr,'rt');
filinfo.hdrascii = fread(fid,inf,'char');
fclose(fid);
filinfo.hdrchar = char(filinfo.hdrascii)';
