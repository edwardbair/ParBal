function keep_only_melt(in)

i=str2double(in);
Adrive = AzureFileShare( );
recon_dir=[Adrive,'reconstruction\energy\h23v05'];
Ndrive = NedFileShare( );
new_dir=[Ndrive,'reconstruction\energy\h23v05'];
d=dir(fullfile(recon_dir,'*.mat'));
tic
fname=fullfile(recon_dir,d(i).name);
newfname=fullfile(new_dir,d(i).name);
try
    m=matfile(fname);
catch
    warning('mfile %s not readable');
%     delete(fname);
    return
end
if any(strcmp(fieldnames(m),'divisors'))
    disp('divisor found');
    [M,MATLABdates]=read_energy(fname,'M');
else
    disp('no divisor found');
    M=single(m.M)./(single(intmax('int32'))./10000);
    MATLABdates=m.MATLABdates;
end
M=int16(M);
save(newfname,'M','MATLABdates','-v7.3')
% delete(fname);
fprintf('done w %s\n',fname);
toc
