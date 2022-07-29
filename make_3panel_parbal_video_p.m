function make_3panel_parbal_video_p(infiles,vidname,rrb,varargin)
%parallel implementation of function to create reprojected parbal MODIS video
%input: infile - input h5 files, cell Nx1
%vidname : output vidname
%rrb : target rasterref w/ CRS
%optional: mask matching rrb, logical

mask0=[];
if nargin==4
    mask0=varargin{1};
end

vars={'melt','sweHybrid','swe'};

[x,y]=worldGrid(rrb);
[lat,lon]=projinv(rrb.ProjectedCRS,x,y);
dy=-2;
dx=2;
y1=ceil(lat(1,1));
y2=floor(lat(end,1));
x1=ceil(lon(1,1));
x2=floor(lon(1,end));
lat_l=y1:dy:y2;
lon_l=x1:dx:x2;

[x,y]=projfwd(rrb.ProjectedCRS,mean(lat_l)*ones(size(lon_l)),lon_l);
[clon,~]=worldToIntrinsic(rrb,x,y);

[x,y]=projfwd(rrb.ProjectedCRS,lat_l,mean(lon_l)*ones(size(lat_l)));
[~,rlat]=worldToIntrinsic(rrb,x,y);
cm=colormap(parula);
close;
cm(1,:)=[0.5 0.5 0.5];
%ltr={'(a)','(b)','(c)','(d)'};

spmd
    figure('Position',[1   1   1900 1200],'Color','w','Visible','off');
    set(gcf,'toolbar','none');
    tiledlayout(2,2,'TileSpacing','none','padding','tight');
    for j=1:length(vars)
        nexttile(j)
        imagesc;
        axis image;
        colormap(cm);
        set(gca,'YDir','reverse','Box','on','XTick',clon,'XTicklabel',[],...
            'YTick',rlat,'YTickLabel',[],'box','off');
        set(gca,'NextPlot','replaceChildren');
        if j==1 || j==3
            set(gca,'YTick',rlat,'YTickLabel',num2str(lat_l'));
        end
        if j>2
            set(gca,'XTick',clon,'XTickLabel',num2str(lon_l'))
        end
        c=colorbar('Location','south');
        c.Position(1)=c.Position(1)+0.25;
        c.Position(2)=c.Position(2)+0.05;
        c.Position(3)=c.Position(3)-0.20;
        if j==1
            %colorbar('Location','south');
            c.Label.String='daily melt, mm';
            caxis([0 100]);
        elseif j==2
            %colorbar('Location','south');
            c.Label.String='swe hybrid, mm';
            caxis([0 2000])
        elseif j==3
            %colorbar('Location','south');
            c.Label.String='recon swe, mm';
            caxis([0 2000]);
        end
    end
end
frames={};

for ii=1:size(infiles,1)
    fname=infiles{ii};
    for j=1:3
        if j==1
            [melt,hdr,matdates]=getMelt(fname,vars{j});
        elseif j==2
            sweHybrid=getMelt(fname,vars{j});
        elseif j==3
            swe=getMelt(fname,vars{j});
        end
    end

    parfor i=1:length(matdates)
        %need to initialize as empty temp var for parfor
        x=[];
        meltf=[];
        sweHybridf=[];
        swef=[];
       
        mask=[];

        melt_i=melt(:,:,i);
        sweHybrid_i=sweHybrid(:,:,i);
        swe_i=swe(:,:,i);

        matdates_i=matdates(i);

        for j=1:length(vars)
            if j==1
                meltf=rasterReprojection(melt_i,...
                    hdr.RasterReference,'InProj',hdr.ProjectionStructure,...
                    'rasterref',rrb);
                meltf(meltf<0)=0;
                if ~isempty(mask0)
                    mask=mask0 & ~isnan(meltf);
                else
                    mask=~isnan(meltf);
                end
                x=meltf;
            elseif j==2
                sweHybridf=rasterReprojection(sweHybrid_i,...
                    hdr.RasterReference,'InProj',hdr.ProjectionStructure,...
                    'rasterref',rrb);
                sweHybridf(sweHybridf<0)=0;
                sweHybridf(isnan(sweHybridf))=0;
                x=sweHybridf;
            elseif j==3
                swef=rasterReprojection(swe_i,...
                    hdr.RasterReference,'InProj',hdr.ProjectionStructure,...
                    'rasterref',rrb);
                swef(swef<0)=0;
                swef(isnan(swef))=0;
                x=swef;
            end
            nexttile(j)
            imagesc(x,'AlphaData',mask);

%             text(1,0.5,ltr{j},'FontSize',25,'VerticalAlignment','bottom',...
%                 'HorizontalAlignment','right','Units','normalized');
            if j==1
                text(0,0.25,datestr(matdates_i),'units','normalized',...
                    'FontSize',25);
            end

        end
        child=get(gcf,'children');
        child=child.Children;
        for y = 1:length(child)
            chi=child(y);
            set(chi, 'fontsize', 25);
        end
        frame=getframe(gcf);
        %parfor ensures frame order will be correct
        frames= [frames, frame];
        fprintf('done w/ %s\n', datestr(matdates_i));
    end
end

f=VideoWriter(vidname);
f.FrameRate=10;
f.Quality=90;
open(f);

for idx=1:numel(frames)
    writeVideo(f,frames{idx})
end

close(f);