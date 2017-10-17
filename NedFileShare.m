function [ drive,varargout ] = NedFileShare( )
% drive = NedFileShare( )
%returns drive letter for the Azure File Share

driveLetter = 'N';
drive = [driveLetter ':\'];
command = ['net use ' driveLetter' ': '...
    '\\tliyel.file.core.windows.net\tliyel /u:tliyel '...
    'JA7kwcFgHGLWurbz62ksXWnQJu9LdDxYEOqATkm4GwlSvPpCWu6Z0yMX/xvNq4TvaR188xCar6nf8OiFa/mwIA=='];
[status,cmdout] = system(command);
% make sure drive is there now
assert(exist(drive,'dir')==7,'drive %s should exist but doesn''t',drive);

if nargout>1
    varargout{1}=command;
end

end