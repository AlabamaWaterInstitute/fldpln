%------------------------------
%function [mn,mx,bg] = scantif(fnam)
%------------------------------
% This function scans ERDAS & ESRI BIL files to extract the minimum,
% maximum, and NoData values.
%------------------------------
function [mn,mx,bg] = scanbil(fnam)
%------------------------------

filinfo = readbilheader(fnam);
rows = filinfo.r;
cols = filinfo.c;
datfmt = filinfo.fmt;
bytord = filinfo.bytord;

dv = floor(rows/100);
md = rows-dv*100;
mn = realmax;
mx = -realmax;
mn2 = mn;

if(bytord==1)
   fid = fopen(fnam,'rb','b');
else
   fid = fopen(fnam,'rb');
end;
%disp(sprintf('Scanning file... 0%%'))

%pct0 = 0;
for j=1:dv
   dat = fread(fid,cols*100,datfmt);
   mn = min(mn,min(dat));
   mx = max(mx,max(dat));
   if(any(dat>mn))
      mn2 = min(mn2,min(dat(dat>mn)));
   end;
%   if(100*j/dv>pct0+10)
%      disp(sprintf(['Scanning file... ',int2str(pct0+10),'%']));
%      pct0 = pct0+10;
%   end;
end;
if(md)
   dat = fread(fid,md,datfmt);
   mn = min(mn,min(dat));
   mx = max(mx,max(dat));
   if(any(dat>mn))
      mn2 = min(mn2,min(dat(dat>mn)));
   end;
end;
fclose(fid);

if(mn>0)
   bg = NaN;
else
   bg = mn; mn = mn2;
end;
