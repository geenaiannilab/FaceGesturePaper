
figure(6)
x = round(rand(200,1)*1);
y = randn(200,1);
subplot(1,3,1)
beeswarm(x,y);
x = round(rand(400,1)*3);
y = randn(400,1);
subplot(1,3,2)
beeswarm(x,y);
x = round(rand(600,1)*7);
y = randn(600,1);
subplot(1,3,3)
beeswarm(x,y);
%% sort methods
x = round(rand(80,1)*2);
y = randn(80,1);
y(x==0)=y(x==0)+2;
y(x==2)=y(x==2)-2;
figure(1);
styles={'up','fan','down','rand','square','hex'};
for i=1:length(styles)
    subplot(2,3,i)
    beeswarm(x,y,'sort_style',styles{i},'dot_size',.5,'overlay_style','sd');
    title(styles{i})
    set(gca,'TitleFontSizeMultiplier',1.75)
    if i==1
        yl=ylim();
    else
        ylim(yl);
    end
end
%% dot size
x = round(rand(50,1)*2);
y = randn(50,1);
y(x==0)=y(x==0)+5;
figure(2)
dvec = [0.25 0.5 0.75 1]*2;
for i=1:length(dvec)
    subplot(1,length(dvec),i)
    beeswarm(x,y,'sort_style','hex','dot_size',dvec(i),'overlay_style','ci');
    title(num2str(dvec(i)))
end
%% overlay styles
x = round(rand(100,1)*2);
y = randn(100,1);
y(x==0)=y(x==0)+5;
figure(3)
styles={'box','sd','ci'};
for i=1:length(styles)
    subplot(1,length(styles),i)
    beeswarm(x,y,'sort_style','up','dot_size',1,'overlay_style',styles{i});
    title(styles{i})
end
%% corral styles
x = round(rand(500,1)*2);
y = randn(500,1);
y(x==0)=y(x==0)+2;
figure(4)
styles = {'none','gutter','omit','random'};
for i=1:length(styles)
    subplot(1,length(styles),i)
    beeswarm(x,y,'sort_style','up','corral_style',styles{i},'dot_size',.5,'overlay_style','ci');
    title(styles{i})
end
%% changing the axes will change results - we need to know the aspect ratio to detect "collisions"
figure(5)
subplot(1,3,1)
beeswarm(x,y,'sort_style','up','dot_size',.5,'overlay_style','sd');
ylim([-3 5])
subplot(1,3,2)
beeswarm(x,y,'sort_style','up','dot_size',.5,'overlay_style','sd');
axis tight
ylim([-3 5])
beeswarm(x,y,'sort_style','up','dot_size',.5,'overlay_style','sd','use_current_axes',true);
subplot(1,3,3)
ylim([-10 10])
xlim([-1 3])
beeswarm(x,y,'sort_style','up','dot_size',.5,'overlay_style','sd','use_current_axes',true);
ylim([-3 5])
%% color can be changed
figure(7)
subplot(1,3,1)
beeswarm(x,y,'MarkerFaceAlpha',1);
subplot(1,3,2)
beeswarm(x,y,'colormap','prism','MarkerFaceAlpha',1);
subplot(1,3,3)
cmap = [1 0 0; 0 1 0; 0 0 1];
xlim([-.5 2.5])
ylim([-3 4])
beeswarm(x,y,'colormap',cmap,'MarkerFaceAlpha',1,'MarkerEdgeColor','k');
