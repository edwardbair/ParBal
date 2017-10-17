function grib_struct=read_grib2(gribname,irec,varargin)
%
%grib_struct=read_grib2(gribname,irec,'HeaderFlag,[0|1],'DataFlag',[0|1]);
%
%read_grib for grib2 files (version 0.1)
%
%HeaderFlag and DataFlag are initiated, as is irec='inventory'
%
%Single and multiple parameter access (irec={'Parameter1",'Parameter2'})
%and irec=vector of record numbers (irec=1:10:100)
%
%Requires that wgrib2 version 0.1.5f or later is installed
%
%If named pipes are available (mkfifo) then they will be used.  This
%is a substantial speedup over writing to temp files
%

if ispc
  disp('This decoder only works on a unix-like system (Linux, BSD, Ubuntu, SUSE, etc)')
  return
end

% where is your local copy of wgrib2?
%[status,result]=system('which wgrib2');
%if status ~= 0
%  disp('wgrib2 is not in your path')
%  return
%end

%wgrib2=strtrim(result);

wgrib2='/usr1/tspindler/bin/wgrib2';

% version check (needs at least 0.1.5f to work)
base='v0.1.5f';
cmd=[wgrib2,' -version'];
[status,result]=system(cmd);
version = strtok(result);
if version(2) < base(2) | version(4) < base(4) | version(6) < base(6) 
  disp('wgrib2 version 0.1.5f or later is required for read_grib2')
  disp('you can download the latest version of wgrib2 from') 
  disp('http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/index.html')
  return
else  
  if length(version) == 7 && version(7) < base(7)
    disp('wgrib2 version 0.1.5f or later is required for read_grib2')
    disp('you can download the latest version of wgrib2 from') 
    disp('http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/index.html')
    return
  end
end

% pipes check (nice speedup if available)
[status,result]=system('which mkfifo');
if status==0
  UsePipes=1;
else
  UsePipes=0;
end
UsePipes=0;

if isempty(gribname)
   [fname,fpath]=uigetfile('*','Which GRiB File?');
   if fname==0,return,end
   gribname=[fpath fname];
else
   [fpath,fname,fext,fver]=fileparts(gribname);
   fname=[fname fext];
end

grib_struct=struct([]);
HeaderFlag=1;
DataFlag=1;

% Process propertyname/value pairs 
k=1;
while k<length(varargin),
  switch lower(varargin{k}),
    case 'headerflag',
      HeaderFlag=varargin{k+1};
      if ~(HeaderFlag==1 | HeaderFlag==0)
         error('Invalid HeaderFlag to READ_GRIB.')
      end
      varargin([k k+1])=[];
    case 'dataflag',
      DataFlag=varargin{k+1};
      if ~(DataFlag==1 | DataFlag==0)
         error('Invalid DataFlag to READ_GRIB.')
      end
      varargin([k k+1])=[];
    case 'screendiag',
      ScreenDiag=varargin{k+1};
      if ~(ScreenDiag==1 | ScreenDiag==0)
         error('Invalid ScreenDiag to READ_GRIB.')
      end
      varargin([k k+1])=[];
    case 'paramtable',
      ParamTable=varargin{k+1};
      %if ~(any(strcmp(ParamTable,{'NCEPOPER','ECMWF160','NCEPREAN','ECMWF128'})))
      %   error('Invalid Parameter Table to READ_GRIB.')
      %end
      varargin([k k+1])=[];
    otherwise
      k=k+2;
  end
end

% do the inventory here
if strncmp(lower(irec),'inv',3)
  [status,result]=system([wgrib2,' -v2 -t -var ',gribname]);
  inventory_str=sprintf('###################################################\n');
  inventory_str=[inventory_str sprintf('Inventory for GRiB file %s \n',fname)];
  inventory_str=[inventory_str sprintf('###################################################\n')];
  inventory_str=[inventory_str sprintf('\n')];
  gribrecords=strread(result,'%s','delimiter','\n');
%  tab=sprintf('\t');
%  inventory=[];
%  for i=1:length(gribrecords)
%    [t1,t2,t3,t4,t5,t6] = strread(char(gribrecords(i)),'%s%s%s%s%s%s','delimiter',':');
%    inventory=strvcat(inventory,[char(t1),tab,char(t2),tab,char(t3),tab,char(t4),tab,char(t5),tab,char(t6)]);
%  end
  helpwin([inventory_str,result]);
  return
end

if ~isnumeric(irec)
  
  nout = 0;
  for nparm = 1:length(irec)
    parmname = char(irec(nparm));
    cmd=[wgrib2,' -match '':',parmname,':'' ',gribname];
    [status,result]=system([cmd,' -v2 -t -nxny -var -grid -one_line ']);
    gribrecords=strread(result,'%s','delimiter','\n');
    
    for i = 1:length(gribrecords)
      nout = nout + 1;
      gribrecord = char(gribrecords(i));
      [nrec,vt,ninj,parm,grid_template,grid]=strread(gribrecord,'%n%*s%s%s%s%s%s%*[^\n]','delimiter',':');
      grib_struct(nout).record = nrec;
      [ni,nj] = strread(char(ninj),'%*c%n%*s%n%*c','delimiter',' ');
      vt = upper(char(strread(char(vt),'%s')));
      %this next is because wgrib2 returns dates as hhZddmmmyyyy
      grid_template = strread(char(grid_template),'%*s%n','delimiter','=');
      if grid_template == 0  %ww3 lat-lon grid decode
        [la1,la2,dila,lo1,lo2,dilo] = strread(char(grid),'%*s%*s%*s%*s%n%*s%n%*s%n%*s%n%*s%n%*s%n%*[^\n]','delimiter',' ');        
        if abs(1/15 - dila) < 1e-6
          dila = 1/15;  % roundoff error in grid definition 
        end
        grib_struct(nout).gds.La1 = la1;
        grib_struct(nout).gds.La2 = la2;
        grib_struct(nout).gds.Di = dila;
        grib_struct(nout).gds.Dj = dilo;
        grib_struct(nout).gds.Lo1 = lo1;
        grib_struct(nout).gds.Lo2 = lo2;        
      end        
      [parmname,description]=strtok(parm);
      grib_struct(nout).gds.Ni = ni;
      grib_struct(nout).gds.Nj = nj;
      grib_struct(nout).gds.template = grid_template;
      
      grib_struct(nout).stime = vt(4:end);
      grib_struct(nout).parameter = char(parmname);
      grib_struct(nout).description = strtrim(char(description));
      
      if DataFlag
        tempfile = tempname;
        cmd = [wgrib2,' -d ',num2str(grib_struct(nout).record),' ',gribname,' -no_header -bin ',tempfile];
        % a very nice speedup for linux and mac boxes
        if UsePipes
          system(['mkfifo ',tempfile]);
          [status,result] = system([cmd,' &']);
	else
          [status,result] = system(cmd);
        end
        fid=fopen(tempfile,'r');
        grib_struct(nout).fltarray = fread(fid,'single');
        fclose(fid);
        delete(tempfile);
      end
    end
  end
else
  for nout = 1:length(irec)
    cmd = [wgrib2,' -v2 -t -nxny -var -grid -one_line -d ',num2str(irec(nout)),' ',gribname];
    [status,result] = system(cmd);
    
    [nrec,vt,ninj,parm,grid_template,grid]=strread(result,'%n%*s%s%s%s%s%s%*[^\n]','delimiter',':');
    grib_struct(nout).record = nrec;
    [ni,nj] = strread(char(ninj),'%*c%n%*s%n%*c','delimiter',' ');
    vt = upper(char(strread(char(vt),'%s')));
    grid_template = strread(char(grid_template),'%*s%n','delimiter','=');
    if grid_template == 0  %ww3 lat-lon grid decode
      [la1,la2,dila,lo1,lo2,dilo] = strread(char(grid),'%*s%*s%*s%*s%n%*s%n%*s%n%*s%n%*s%n%*s%n%*[^\n]','delimiter',' ');        
      if abs(1/15 - dila) < 1e-6
        dila = 1/15;  % roundoff error in grid definition 
      end
      grib_struct(nout).gds.La1 = la1;
      grib_struct(nout).gds.La2 = la2;
      grib_struct(nout).gds.Di = dila;
      grib_struct(nout).gds.Dj = dilo;
      grib_struct(nout).gds.Lo1 = lo1;
      grib_struct(nout).gds.Lo2 = lo2;        
    end        
    [parmname,description]=strtok(parm);
    grib_struct(nout).gds.Ni = ni;
    grib_struct(nout).gds.Nj = nj;
    grib_struct(nout).gds.template = grid_template;
    
    grib_struct(nout).stime = vt(4:end);
    grib_struct(nout).parameter = char(parmname);
    grib_struct(nout).description = strtrim(char(description));
    
    tempfile = tempname;
    cmd = [wgrib2,' -d ',num2str(irec(nout)),' ',gribname,' -no_header -bin ',tempfile];
    % a very nice speedup for linux and mac boxes
    if UsePipes
      system(['mkfifo ',tempfile]);
      [status,result] = system([cmd,' &']);
    else
      [status,result] = system(cmd);
    end
    fid=fopen(tempfile,'r');
    grib_struct(nout).fltarray = fread(fid,'single');
    fclose(fid);
    delete(tempfile);
  end
end  

return
