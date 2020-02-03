%% @authors: Mara Scussolini & Sara Garbarino 
% contact: garbarino@dima.unige.it

function  compartmental_gui(action,varargin)
clc
warning('off')

if nargin<1
    
    global analysis_gui
    
    if ispc==1, analysis_gui.slash='\'; else analysis_gui.slash='/'; end
    
    analysis_gui.start_path = pwd;
    analysis_gui.path = [];
    
    pathfile=[analysis_gui.start_path analysis_gui.slash 'CONFIG' analysis_gui.slash 'path.txt'];
    if exist(pathfile,'file')
        fo=fopen(pathfile,'r');
        gui_pyr.path=fgetl(fo);
        fclose(fo);
    else
        SetPath_Callback;
    end
    
    if exist(pathfile,'file')
        comp_g;
    end
    
else
    feval(action,varargin{:});
end

end

function analysis_gui = comp_g()

global analysis_gui

analysis_gui.OUTPUTfolder = [];
analysis_gui.DATAfolder = [];

analysis_gui.fig = figure('Visible','off','Position',[200,200,450,600],'MenuBar','none'); %,'Toolbar','figure');

analysis_gui.setdata = uicontrol('Style','pushbutton','String','Load PET Data',...
    'Position',[80,540,100,40],...
    'HorizontalAlignment','Right',...
    'Callback',{@setDATA_Callback});

analysis_gui.setIF = uicontrol('Style','pushbutton','String','Load BLOOD IF',...
    'Position',[80,490,100,40],...
    'HorizontalAlignment','Right',...
    'Callback',{@setIF_Callback});

analysis_gui.set3D = uicontrol('Style','pushbutton','String','View 3D',...
    'Position',[280,540,80,30],...
    'HorizontalAlignment','Right',...
    'Callback',{@set3D_Callback});

analysis_gui.setplotIF = uicontrol('Style','pushbutton','String','View IF',...
    'Position',[280,490,80,30],...
    'HorizontalAlignment','Right',...
    'Callback',{@setplotIF_Callback});

analysis_gui.PETslicetxt = uicontrol('Style','text','String','Select PET slice',...
    'Position',[80,400,100,40]);
analysis_gui.PETslice = uicontrol('Style','edit','String','',...
    'Position',[180,422,80,20]);

analysis_gui.timetxt = uicontrol('Style','text','String','Select time-frame:',...
    'Position',[80,370,100,40]);
analysis_gui.timetoview = uicontrol('Style','edit','String','',...
    'Position',[180,392,80,20]);

analysis_gui.setviewslice = uicontrol('Style','pushbutton','String','View',...
    'Position',[280,400,80,30],...
    'HorizontalAlignment','Right',...
    'Callback',{@setVIEWSLICE_Callback});

analysis_gui.modeltxt = uicontrol('Style','text','String','Select model',...
    'Position',[80,300,100,40]);
analysis_gui.setmodel = uicontrol('Style','popup','String','|Tumor|Kidney',...%'|Tumor|Kidney|Liver',...
    'Position',[175,320,90,20],...
    'HorizontalAlignment','Right'); %,...
% 'Callback',{@setMODEL_Callback});
analysis_gui.setviewmodel = uicontrol('Style','pushbutton','String','View',...
    'Position',[280,315,80,30], ...
    'HorizontalAlignment','Right',...
    'Callback',{@setVIEWMODEL_Callback});

analysis_gui.Vbtxt = uicontrol('Style','text','String','Blood volume fraction Vb=',...
    'Position',[80,270,100,40]);
analysis_gui.Vb = uicontrol('Style','edit','String','',...
    'Position',[180,290,80,20]);
% Vb tumor = 0.15;
% Vb Kidney = 0.2;

analysis_gui.setsegment = uicontrol('Style','pushbutton','String','Segment PET image',...
    'Position',[80,200,150,40],...
    'HorizontalAlignment','Right',...
    'Callback',{@setSEGMENT_Callback});

analysis_gui.setviewsegmentation = uicontrol('Style','pushbutton','String','View',...
    'Position',[280,200,80,40],...
    'HorizontalAlignment','Right',...
    'Callback',{@setVIEWSEGMENTATION_Callback});

analysis_gui.setREC = uicontrol('Style','pushbutton','String','Start RECONSTRUCTION',...
    'Position',[80,150,150,40],...
    'HorizontalAlignment','Right',...
    'Callback',{@setREC_Callback});

analysis_gui.setviewrec = uicontrol('Style','pushbutton','String','View',...
    'Position',[280,150,80,40],...
    'HorizontalAlignment','Right',...
    'Callback',{@setVIEWREC_Callback});

analysis_gui.setexit = uicontrol('Style','pushbutton','String','EXIT',...
    'Position',[280,30,80,60],...
    'HorizontalAlignment','Right',...
    'Callback',{@setEXIT_Callback});

% analysis_gui.setoutput = uicontrol('Style','pushbutton','String','set OUTPUT',...
%     'Position',[110,241.5,225,25],...
%     'HorizontalAlignment','Right',...
%     'Callback',{@setOUTPUT_Callback});

%-------------------------------------------------------------------------%
set([analysis_gui.fig,analysis_gui.setdata,analysis_gui.setIF,analysis_gui.set3D,analysis_gui.setplotIF,...
    analysis_gui.PETslicetxt,analysis_gui.PETslice,analysis_gui.timetxt,analysis_gui.timetoview,analysis_gui.setviewslice,...
    analysis_gui.modeltxt,analysis_gui.setmodel,analysis_gui.setviewmodel,analysis_gui.Vbtxt,analysis_gui.Vb,...
    analysis_gui.setsegment,analysis_gui.setviewsegmentation,analysis_gui.setREC,analysis_gui.setviewrec,...
    analysis_gui.setexit],...
    'Units','normalized');

set(analysis_gui.fig,'Name','Compartmental Analysis - PARAMETRIC IMAGING GUI');
set(analysis_gui.fig,'NumberTitle','off');
movegui(analysis_gui.fig,'center');
%-------------------------------------------------------------------------%
set(analysis_gui.fig,'Visible','on');

end

%% Callback

%-------------------------------------------------------------------------%
function SetPath_Callback(hObject, evendata, handles)

global analysis_gui;

if isempty(analysis_gui.path),
    str='No path selected. Click Yes to procede...';
else
    str=sprintf('Current path:\n %s\n\n Click Yes to change it...',analysis_gui.path);
end

choice = questdlg(str,'path setting:', ...
    'Yes','No','Yes');
switch choice
    case 'Yes'
        path=uigetdir(analysis_gui.start_path, 'Choose path');
        
        if path == 0
        else
            analysis_gui.path=path;
            config_dir=[analysis_gui.start_path analysis_gui.slash 'CONFIG'];
            if ~exist(config_dir,'dir'), mkdir(config_dir); end
            
            fo=fopen([config_dir,analysis_gui.slash 'path.txt'],'w');
            fprintf(fo,'%s',path);
            fclose(fo);
        end
    case 'No'
end

end

function setDATA_Callback(hObject, evendata, handles)

global analysis_gui

[namefile,location]=uigetfile('*.hdr','Load PET Data');

set([analysis_gui.setdata],'String',namefile,'Units','normalized');

if namefile ~= 0
    if ~iscell(namefile)
        if namefile==0
        else
            if isfield(analysis_gui,'namefile'), analysis_gui=rmfield(analysis_gui,'namefile'); end
            if isfield(analysis_gui,'location'), analysis_gui=rmfield(analysis_gui,'location');  end
        end
    end
    
    analysis_gui.pathDATA = strcat(location,namefile);
    
    analysis_gui.DATA4D = albira_readfiles(namefile,location);
    cd(analysis_gui.start_path)
    
    [analysis_gui.dim_x,analysis_gui.dim_y,analysis_gui.dim_z,analysis_gui.dim_t] = ...
        size(analysis_gui.DATA4D);
    
else  set([analysis_gui.setdata],'String','Load PET Data','Units','normalized');
end

end

function setIF_Callback(hObject, evendata, handles)

global analysis_gui

[namefile_if,location_if]=uigetfile('*.VOISTAT','Load BLOOD IF');

if namefile_if ~= 0
    
    set([analysis_gui.setIF],'String',namefile_if,'Units','normalized');
    
    if ~iscell(namefile_if)
        if namefile_if==0
        else
            if isfield(analysis_gui,'namefile_if'), analysis_gui=rmfield(analysis_gui,'namefile_if'); end
            if isfield(analysis_gui,'location_if'), analysis_gui=rmfield(analysis_gui,'location_if');  end
        end
    end
    
    analysis_gui.pathIF = strcat(location_if,namefile_if);
    
    delimiter = '\t';
    startRow = 4;
    % %s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%s%*s%*s%*s%[^\n\r]
    formatSpec = '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%[^\n\r]';
    fileID = fopen(analysis_gui.pathIF,'r');
    dataArray = textscan(fileID,formatSpec,'Delimiter',delimiter,'HeaderLines',startRow-1,'ReturnOnError',false);
    fclose(fileID);
    
    analysis_gui.time = dataArray{:, 4}./60;
    analysis_gui.IF = dataArray{:, 7};
    analysis_gui.IF(21) = (analysis_gui.IF(20)+analysis_gui.IF(22))*1/2; % #0 bad acquisition
    analysis_gui.Volume_IF = dataArray{:, 11}; analysis_gui.Volume_IF = analysis_gui.Volume_IF(1,:);
    analysis_gui.Ca = @(tt)(interp1([0;analysis_gui.time],[0;analysis_gui.IF],tt,'linear',0)).';

    
else  set([analysis_gui.setif],'String','Load BLOOD IF','Units','normalized');
end

end

function set3D_Callback(hObject, evendata, handles)

global analysis_gui

meanDATA = squeeze(mean(analysis_gui.DATA4D(:,:,:,end),4));

figure('numbertitle', 'off'); %('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'Color', 'white');
h_data = vol3d('cdata',meanDATA,'texture','3D');
view(3);
axis tight; axis off; daspect([1 1 1])
alphamap('rampup');
alphamap(2.* alphamap);

end

function setplotIF_Callback(hObject, evendata, handles)

global analysis_gui

figure('numbertitle', 'off'); %('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'Color', 'white'); %white bckgr
plot(analysis_gui.time,analysis_gui.IF,'r','LineWidth',2);
axis([min(analysis_gui.time) max(analysis_gui.time) min(analysis_gui.IF) max(analysis_gui.IF)+50])
axis square
xlabel('time [min]','FontSize',15); ylabel('concentration [kBq/mL]','FontSize',15);
legend({'BLOOD IF'},'FontSize',15);

end

function setVIEWMODEL_Callback(hObject, evendata, handles)

global analysis_gui

aux=get(analysis_gui.setmodel,'value');
pathfigure=[analysis_gui.start_path analysis_gui.slash 'model' analysis_gui.slash];

if aux == 2
    
    I = imread([pathfigure 'tumor2.jpg']);
    figure('numbertitle', 'off');
    imshow(I);
    
elseif aux == 3
    
    I = imread([pathfigure 'kidney2.jpg']);
    figure('numbertitle', 'off');
    imshow(I);
    
    % elseif aux == 4
    %
    %     I = imread([pathfigure 'liver.jpg']);
    %     figure('numbertitle', 'off');
    %     imshow(I);
    
    %
end

end

function setVIEWSLICE_Callback(hObject, evendata, handles)
addpath('./scripts')
global analysis_gui

slice = str2double(get(analysis_gui.PETslice,'String'));
tview = str2double(get(analysis_gui.timetoview,'String'));

C = analysis_gui.DATA4D(:,:,slice,:); C = squeeze(C);
C(:,:,21) = (C(:,:,20)+C(:,:,22)).*(1/2); % #0 bad acquisition

% Gaussian smoothing filter
filter_size=3; filter_sigma=1;
analysis_gui.DATA_filt = zeros(analysis_gui.dim_x,analysis_gui.dim_y,analysis_gui.dim_t);
for n=1:analysis_gui.dim_t
    analysis_gui.DATA_filt(:,:,n)=gaussian_filtering(C(:,:,n),filter_size,filter_sigma);
end

figure('numbertitle', 'off'); %('units','normalized','outerposition',[0 0 1 1]);
imagesc(analysis_gui.DATA_filt(:,:,tview));
colorbar; colormap(hot);
axis square; axis off;
title(['PET slice: ', num2str(slice),'    -    Time (min): ', num2str(analysis_gui.time(tview))],...
    'Fontsize',15)

end

function setSEGMENT_Callback(hObject, evendata, handles)
addpath('./scripts')
global analysis_gui

slice = str2double(get(analysis_gui.PETslice,'String'));
tview = str2double(get(analysis_gui.timetoview,'String'));

C = analysis_gui.DATA4D(:,:,slice,:); C = squeeze(C);
C(:,:,21) = (C(:,:,20)+C(:,:,22)).*(1/2); % #0 bad acquisition

% Gaussian smoothing filter
filter_size=3; filter_sigma=1;
analysis_gui.DATA_filt = zeros(analysis_gui.dim_x,analysis_gui.dim_y,analysis_gui.dim_t);
for n=1:analysis_gui.dim_t
    analysis_gui.DATA_filt(:,:,n)=gaussian_filtering(C(:,:,n),filter_size,filter_sigma);
end

aux=get(analysis_gui.setmodel,'value');

% PET image end-time
C_end = analysis_gui.DATA_filt(:,:,end);% mean(analysis_gui.DATA_filt,3);

% we fit a Gaussian mixture model with 2-peaks gaussian

[i_max,~] = find(C_end==max(max(C_end)));
mp_max = C_end(i_max,:);
f = fit([1:analysis_gui.dim_y]',mp_max','gauss2');
stationary_points = find(diff(sign(diff(f([1:analysis_gui.dim_y]'))))~=0); % now we are ensured to have 3 stationary points; t
cut = stationary_points(2); %the second is the cut
% j_max_kidneys = stationary_points(1); % the first is the point of max intensity for kidneys
% j_max_tumor = stationary_points(3); % the third is the point of max intensity for tumor

analysis_gui.col_limit = cut;

% Consider pixels with intensity upper the bound
bound = C_end(i_max,analysis_gui.col_limit);
index = (C_end>=bound);
if aux == 2
    index(:,1:cut)=0;
elseif aux == 3
    index(:,cut:end)=0;
end

% Segmented image
analysis_gui.DATA_segmented = zeros(size(analysis_gui.DATA_filt));

h = waitbar(0,'Image Segmentation. Please wait...');
steps = analysis_gui.dim_t;
c = 0;
for n=1:analysis_gui.dim_t
    analysis_gui.DATA_segmented(:,:,n) = double(index).*analysis_gui.DATA_filt(:,:,n);
    c = c + 1;
    waitbar(c / steps, h)%, sprintf('%d of %d', c, steps));
end
close(h)

end

function setVIEWSEGMENTATION_Callback(hObject, evendata, handles)

global analysis_gui

aux=get(analysis_gui.setmodel,'value');

tview = str2double(get(analysis_gui.timetoview,'String'));

figure('numbertitle', 'off'); %('units','normalized','outerposition',[0 0 1 1]);
imagesc(analysis_gui.DATA_segmented(:,:,tview));
colorbar; colormap(hot);
axis square; axis off;
if aux == 2
    title(['Segmented Tumor    -    Time (min):  ', num2str(analysis_gui.time(tview))],...
        'Fontsize',15)
elseif aux == 3
    title(['Segmented Kidney    -    Time (min):  ', num2str(analysis_gui.time(tview))],...
        'Fontsize',15)
end

end

function setREC_Callback(hObject, evendata, handles)
addpath('./scripts')

global analysis_gui

aux=get(analysis_gui.setmodel,'value');
% C_end = analysis_gui.DATA_filt(:,:,end);% mean(analysis_gui.DATA_filt,3);
C_end = analysis_gui.DATA_segmented;

if aux == 2
    analysis_gui.K1x = zeros(analysis_gui.dim_x,analysis_gui.dim_y); analysis_gui.K2x = zeros(analysis_gui.dim_x,analysis_gui.dim_y);
    analysis_gui.K3x = zeros(analysis_gui.dim_x,analysis_gui.dim_y); analysis_gui.K4x = zeros(analysis_gui.dim_x,analysis_gui.dim_y);
    
    relerr_ctr = zeros(analysis_gui.dim_x,analysis_gui.dim_y);
    nit_ctr = zeros(analysis_gui.dim_x,analysis_gui.dim_y);
        
    % Lower bound for radioactive activity: if the norm of Cdata (for a pixel)
    % is under this bound, we can consider no relevant radioactivity (in that pixel)
    activity = 100;
    h = waitbar(0,'Parametric reconstruction. Please wait...');
    steps = analysis_gui.dim_x * analysis_gui.dim_y;
    c = 0;
%     analysis_gui.index_i = [];
%     analysis_gui.index_j = [];
    for i = 1:analysis_gui.dim_x
        for j = 1:analysis_gui.dim_y
            Cdata = squeeze(C_end(i,j,:));
            c = c + 1;
            waitbar(c / steps, h)%, sprintf('%d of %d', c, steps));
            if norm(Cdata)>=activity
                [analysis_gui.K1x(i,j),analysis_gui.K2x(i,j),analysis_gui.K3x(i,j),analysis_gui.K4x(i,j),relerr_ctr(i,j),nit_ctr(i,j)] = ...
                    reconstruction_2C(Cdata,analysis_gui.Ca,analysis_gui.time,0,[0;0]);  
            end
        end
    end
    close(h)
    
elseif aux == 3
    analysis_gui.K1x_kid = zeros(analysis_gui.dim_x,analysis_gui.dim_y); analysis_gui.K2x_kid = zeros(analysis_gui.dim_x,analysis_gui.dim_y);
    analysis_gui.K3x_kid = zeros(analysis_gui.dim_x,analysis_gui.dim_y); analysis_gui.K4x_kid = zeros(analysis_gui.dim_x,analysis_gui.dim_y);
    analysis_gui.K5x_kid = zeros(analysis_gui.dim_x,analysis_gui.dim_y); analysis_gui.K6x_kid = zeros(analysis_gui.dim_x,analysis_gui.dim_y);
    analysis_gui.K7x_kid = zeros(analysis_gui.dim_x,analysis_gui.dim_y);
    
    relerr_ctr = zeros(analysis_gui.dim_x,analysis_gui.dim_y);
    nit_ctr = zeros(analysis_gui.dim_x,analysis_gui.dim_y);
        
    % Lower bound for radioactive activity: if the norm of Cdata (for a pixel)
    % is under this bound, we can consider no relevant radioactivity (in that pixel)
    activity = 100;
    h = waitbar(0,'Parametric reconstruction. Please wait...');
    steps = analysis_gui.dim_x * analysis_gui.dim_y;
    c = 0;
    analysis_gui.index_i = [];
    for i = 1:analysis_gui.dim_x
        for j = 1:analysis_gui.dim_y
            Cdata = squeeze(C_end(i,j,:));
            c = c + 1;
            waitbar(c / steps, h)%, sprintf('%d of %d', c, steps));
            if norm(Cdata)>=activity
                [analysis_gui.K1x_kid(i,j),analysis_gui.K2x_kid(i,j),analysis_gui.K3x_kid(i,j),analysis_gui.K4x_kid(i,j),...
                    analysis_gui.K5x_kid(i,j),analysis_gui.K6x_kid(i,j),analysis_gui.K7x_kid(i,j),relerr_ctr(i,j),nit_ctr(i,j)] = ...
                    reconstruction_3C_ctr(Cdata,analysis_gui.Ca,analysis_gui.time);
            end
        end
    end
    close(h)
end
end

function setVIEWREC_Callback(hObject, evendata, handles)

global analysis_gui

aux=get(analysis_gui.setmodel,'value');

if aux == 2
    
    [ii,jj]=find(analysis_gui.K1x>0);

    figure('units','normalized','outerposition',[0 0 1 1],'numbertitle', 'off');
    set(gcf,'Color','white');
    
    subplot(2,2,1)
    imagesc(analysis_gui.K1x(min(ii)-3:max(ii)+3,min(jj)-3:max(jj)+3),[min(min(analysis_gui.K1x(analysis_gui.K1x>0))) max(max(analysis_gui.K1x(analysis_gui.K1x>0)))]);
    axis image; axis off;
    colorbar; colormap(hot);
    title({'$k_{fb}$'},'FontSize',20,'Interpreter','Latex')
    
    subplot(2,2,2)
    imagesc(analysis_gui.K2x(min(ii)-3:max(ii)+3,min(jj)-3:max(jj)+3),[min(min(analysis_gui.K2x(analysis_gui.K2x>0))) max(max(analysis_gui.K2x(analysis_gui.K2x>0)))]);
    axis image; axis off;
    colorbar; colormap(hot);
    title({'$k_{bf}$'},'FontSize',20,'Interpreter','Latex')
    
    subplot(2,2,3)
    imagesc(analysis_gui.K3x(min(ii)-3:max(ii)+3,min(jj)-3:max(jj)+3),[min(min(analysis_gui.K3x(analysis_gui.K3x>0))) max(max(analysis_gui.K3x(analysis_gui.K3x>0)))]);
    axis image; axis off;
    colorbar; colormap(hot);
    title({'$k_{mf}$'},'FontSize',20,'Interpreter','Latex')
    
    subplot(2,2,4)
    imagesc(analysis_gui.K4x(min(ii)-3:max(ii)+3,min(jj)-3:max(jj)+3),[min(min(analysis_gui.K4x(analysis_gui.K4x>0))) max(max(analysis_gui.K4x(analysis_gui.K4x>0)))]);
    axis image; axis off;
    colorbar; colormap(hot);
    title({'$k_{fm}$'},'FontSize',20,'Interpreter','Latex')
    
elseif aux == 3
    
    [ii,jj]=find(analysis_gui.K1x_kid>0);

    figure('units','normalized','outerposition',[0 0 1 1],'numbertitle', 'off');
    set(gcf,'Color','white');
    
    subplot(2,4,1.5)
    imagesc(analysis_gui.K1x_kid(min(ii)-3:max(ii)+3,min(jj)-3:max(jj)+3),...
        [min(min(analysis_gui.K1x_kid(analysis_gui.K1x_kid>0))) max(max(analysis_gui.K1x_kid(analysis_gui.K1x_kid>0)))]);
    axis image; axis off;
    colorbar; colormap(hot);
    title({'$k_{fa}$'},'FontSize',20,'Interpreter','Latex')
    
    subplot(2,4,2.5)
    imagesc(analysis_gui.K2x_kid(min(ii)-3:max(ii)+3,min(jj)-3:max(jj)+3),...
        [min(min(analysis_gui.K2x_kid(analysis_gui.K2x_kid>0))) max(max(analysis_gui.K2x_kid(analysis_gui.K2x_kid>0)))]);
    axis image; axis off;
    colorbar; colormap(hot);
    title({'$k_{af}$'},'FontSize',20,'Interpreter','Latex')
    
    subplot(2,4,3.5)
    imagesc(analysis_gui.K1x_kid(min(ii)-3:max(ii)+3,min(jj)-3:max(jj)+3),...
        [min(min(analysis_gui.K1x_kid(analysis_gui.K1x_kid>0))) max(max(analysis_gui.K1x_kid(analysis_gui.K1x_kid>0)))]);
    axis image; axis off;
    colorbar; colormap(hot);
    title({'$k_{ma}$'},'FontSize',20,'Interpreter','Latex')
    
    subplot(2,4,5)
    imagesc(analysis_gui.K1x_kid(min(ii)-3:max(ii)+3,min(jj)-3:max(jj)+3),...
        [min(min(analysis_gui.K1x_kid(analysis_gui.K1x_kid>0))) max(max(analysis_gui.K1x_kid(analysis_gui.K1x_kid>0)))]);
    axis image; axis off;
    colorbar; colormap(hot);
    title({'$k_{mf}$'},'FontSize',20,'Interpreter','Latex')
    
    subplot(2,4,6)
    imagesc(analysis_gui.K5x_kid(min(ii)-3:max(ii)+3,min(jj)-3:max(jj)+3),...
        [min(min(analysis_gui.K5x_kid(analysis_gui.K5x_kid>0))) max(max(analysis_gui.K5x_kid(analysis_gui.K5x_kid>0)))]);
    axis image; axis off;
    colorbar; colormap(hot);
    title({'$k_{fm}$'},'FontSize',20,'Interpreter','Latex')
    
    subplot(2,4,7)
    imagesc(analysis_gui.K6x_kid(min(ii)-3:max(ii)+3,min(jj)-3:max(jj)+3),...
        [min(min(analysis_gui.K6x_kid(analysis_gui.K6x_kid>0))) max(max(analysis_gui.K6x_kid(analysis_gui.K6x_kid>0)))]);
    axis image; axis off;
    colorbar; colormap(hot);
    title({'$k_{tm}$'},'FontSize',20,'Interpreter','Latex')
    
    subplot(2,4,8)
    imagesc(analysis_gui.K7x_kid(min(ii)-3:max(ii)+3,min(jj)-3:max(jj)+3),...
        [min(min(analysis_gui.K7x_kid(analysis_gui.K7x_kid>0))) max(max(analysis_gui.K7x_kid(analysis_gui.K7x_kid>0)))]);
    axis image; axis off;
    colorbar; colormap(hot);
    title({'$k_{ut}$'},'FontSize',20,'Interpreter','Latex')
    
end


end

function setEXIT_Callback(hObject, evendata, handles)
clc
clear all
close all
% exit
end







