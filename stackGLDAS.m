function GLDAS = stackGLDAS(gldas_filelist)
%function GLDAS = stackGLDAS(gldas)
% Load GLDAS files into a structure
%INPUT
% glads_filelist - structure of gldas filenames for a single day
%OUTPUT
%  GLDAS structure

%check for primary/secondary location
% if isfield(gldas_filelist,'vars_locs')
%     secondary_flag=true;
% else
%     secondary_flag=false;
% end
% Read files
for n=1:length(gldas_filelist.filenames)
    for t=1:length(gldas_filelist.var)
%         if secondary_flag && strcmp(gldas_filelist.var,'')
%             field_name='secondary_filenames';
%         else
            field_name='filenames';
%         end
        fname=gldas_filelist.(field_name){n};
        LDAS = read_LDAS(gldas_filelist.var{t},fname);
        GLDAS.(gldas_filelist.var{t})(:,:,n)=single(LDAS.matrix);
    end
end
% Datevals
GLDAS.datevalsUTC=gldas_filelist.datevalsUTC;
GLDAS.datevalsLocal=gldas_filelist.datevalsLocal;