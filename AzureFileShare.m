function [ drive,varargout ] = AzureFileShare( )
% drive = AzureFileShare( )
%returns drive letter for the Azure File Share

driveLetter = 'Z';
drive = [driveLetter ':\'];
command = ['net use ' driveLetter' ': '...
    '\\quinaya.file.core.windows.net\jeff /u:quinaya '...
    'pAH43zAdfJCNoxYCKvNdFQ51PY3WHKhsFrgSE6dgwXA/Uf7/ytW7EN30Jplm3gnyAzefKovC2Uoz7CTE6I+S/Q=='];
[status,cmdout] = system(command);
% make sure drive is there now
assert(exist(drive,'dir')==7,'drive %s should exist but doesn''t',drive);

if nargout>1
    varargout{1}=command;
end

end