%---------------------
% function tabl = readshapefield(shpfil,flds,mflg (optional input))
%---------------------
% This function reads data from numeric fields in shapefiles.
%
% INPUT
%  shpfil = string with path & filename for shapefile
%  flds = cell array with field names to be read
%
% OPTIONAL INPUT
%  mflg = 1 if a *.mat file with the same name as 'shpfil' is to be saved
%         in (or read from if the same shapefile is accessed later) the
%         same loction as 'shpfil' (variables s & a are saved; type
%         "help shaperead" at the prompt to see what these are)
%
% OUTPUT
%  tabl = nrec-by-nfld table of values from the shapefile, where {nrec =
%         number of records in the shapefile} and {nfld = length of 'flds'
%         input}
%---------------------
function outp = readshapefield(shpfil,flds,varargin)

if(nargin==3)
   if(varargin{1}==1)
      if(~exist([shpfil(1:end-3),'mat'],'file'))
         [s,a] = shaperead(shpfil);
         save([shpfil(1:end-3),'mat'],'s','a','-v7.3');
      else
         load([shpfil(1:end-3),'mat']);
      end;
   else
      error('Optional input (3rd argument) must be either 1 or not used');
   end;
else
   [s,a] = shaperead(shpfil);
end;

num = length(a);
nfld = length(flds);
outp = zeros(num,nfld);
for k = 1:nfld
   if(or(isequal(flds{k},'X'),isequal(flds{k},'Y')))
      for j = 1:num
         outp(j,k) = eval(['s(j).',flds{k}]);
      end;
   else
      for j = 1:num
         outp(j,k) = eval(['a(j).',flds{k}]);
      end;
   end;
end;
