function [ sampleArea ] = SampleAreaCalcFunc( SampleFname, varargin )
%SAMPLEAREACALCFUNC is used to determine the size of samples testing in the
%SQUID / VSM. Place the sample on a piece of graph paper / engineering
%paper where each square is
%   SampleFname is the file name of the picture being used of the sample on
%   a piece of graph paper. The size of the cells on the graph paper needs
%   to be known before using the program
%   varargin: cell array of optional parameters
%               {1} - Median Number
%               {2} - Median Resolution
%               {3} - Black and White Threshold
%               {4} - Refinement

%Use this code to determine appropriate threshold value, validate by visual
%inspection yourself.

set(0,'DefaultFigureWindowStyle','Docked')
sampleFig = 1;

%% Set Options
if isempty(varargin)    
    median_num = 10;
    median_res = 5;
    threshold = .40;
    refinement = 5;
else 
    median_num = varargin{1};
    median_res = varargin{2};
    threshold = varargin{3};
    refinement = varargin{4};
end

%% Load Images
%Color images
S_uint8 = imread(SampleFname); % image of sample on grid
S_dbl = im2double(S_uint8);

[red, green, blue] = imageLayers(S_dbl);
[red_inv, green_inv, blue_inv] = inverseImages(red, green, blue);

%% Segment Grid Only Region for Rotaiton / Cell Size Calculation
%Segment out small grid region without sample in it
%User input - select two points to define a box (must include at least 4
%complete boxes from the image)

% Only change if debugging (hard code in region to study)

% rows = 132:132+1366;
% cols = 1511:1511+1522;
% S_crop = S_dbl(rows,cols,:);

[S_crop,~,~] = UserCrop(S_dbl, sampleFig,'Don''t enclose the sample');

%% Rotate image
%find the peak values of the grid spacing, and rotate to maximize the value
%average value
ave = mean(mean(mean(S_dbl)));

%Anonymous function that is minimized to find optimal angle
%minFunc is the actuall angle being called. Since we want to maximize the
%function, we take the negative, and feed it into a min finder
prominence = @(theta) -minFunc(S_crop,theta, ave);

%Bound domain minimum solver (sweep from -45 degree rotation to zero)
%this small domain works since a cartesian grid has 4-fold symmetry
theta = fminbnd(prominence,-45,0);

%Rotate Image (apply optimal angle
S_rot = imrotate(S_dbl,theta);
ave = mean(mean(mean(S_dbl)));
S_rot(S_rot == 0) = NaN;

S_crop_rot = imrotate(S_crop,theta);
ave = mean(mean(mean(S_crop)));
S_crop_rot(S_crop_rot == 0) = ave;

%% Determine Grid Spacing
factor = 3;
[pksX,  locsX, wX, pX, pksY,  locsY, wY, pY, spacingX, spacingY, spacing] = GridSpacing(S_crop_rot,ave,factor);

%% Calculate Pixel Area
spacings = [spacingX, spacingY, spacing];
[pixelArea,pixelX,pixelY] = UserInfo(spacings, sampleFig, 4);
pause(2)

%% Crop Images to Determine Sample Size
%User input - select two points to define a box (must include at least 4
%complete boxes from the image)
[S_crop,rows,cols] = UserCrop(S_rot, sampleFig,'Surround the sample');


%% Locate Sample Edges and Determine Area
pixels = SamplePixelCount(S_crop, sampleFig, median_num, median_res, threshold, refinement);

Area = pixels * pixelArea;
clc
disp(['Pixel Area = ', num2str(pixelArea),' mm^2'])
disp('')
disp(['Sample is ', num2str(pixels),' pixels'])
disp('')
disp(['Sample Area = ',num2str(Area),' mm^2'])

close all

end



%% Subfunctions

function [pks,  locs, w, p] = LineExtrema(img, scale)
flag = true; %use to make sure r,g and b have the same number of peaks

%image layers
[r,g,b] = imageLayers(img);
r = r(:)';
g = g(:)';
b = b(:)';
I = r+g+b;  %image intensity

while flag
    %Minimum Allowable Peak Prominence
    [rp, gp, bp] = Prominence(img, scale);
    Ip = 4*nanstd(I)/scale;
    
    %Find Peaks
    %     [rPks,  rLocs, rW, rP] = findpeaks(-r,'MinPeakProminence',rp,'Annotate','extents');
    %     [gPks,  gLocs, gW, gP] = findpeaks(-g,'MinPeakProminence',gp,'Annotate','extents');
    %     [bPks,  bLocs, bW, bP] = findpeaks(-b,'MinPeakProminence',bp,'Annotate','extents');
    [IPks,  ILocs, IW, IP] = findpeaks(-I,'MinPeakProminence',Ip,'Annotate','extents');
    if numel(IPks)<3
        scale = scale+1;
    else
        flag = false;
    end
    
    
    %     check1 = size(rPks)==size(gPks);
    %     check2 = size(rPks)==size(bPks);
    %     if all([check1,check2])
    %         flag = false;
    %     elseif scale > 6
    %         flag = false;
    %     else
    %         scale = scale + 1;
    %     end
end
% pks = abs([rPks(:);gPks(:);bPks(:)]);
% locs = [rLocs(:);gLocs(:);bLocs(:)];
% w = [rW(:);gW(:);bW(:)];
% p = [rP(:);gP(:);bP(:)];

pks = abs(IPks);
locs = ILocs;
w = IW;
p = IP;

end


function [rP,gP,bP] = Prominence(img,scale)
%Segment Image Into Layers
[r, g, b] = imageLayers(img);
r = r(:)';
g = g(:)';
b = b(:)';

%Minimum Allowable Peak Prominence
rP = (max(r) - min(r))./scale;
gP = (max(g) - min(g))./scale;
bP = (max(b) - min(b))./scale;
end


function [Area, CC, L, L2] = ImageProperties(BW)
% Feed in a black and white image (BW) and retrieve information for the
% varoius regions comprising the image
Area = regionprops(BW, ones(size(BW)),'area');    %Region Area
CC = bwconncomp(BW);    %Connected components list (pixel list for each region)
L = labelmatrix(CC);    %Label matrix
L2 = bwlabel(BW);       %Label matrix
end


function [red, green, blue] = imageLayers(Img)
red = Img(:,:,1);
green = Img(:,:,2);
blue = Img(:,:,3);
end


function [red_inv,green_inv,blue_inv] = inverseImages(red,green,blue)
%Inverse Images
Rmax = max(max(red));
Gmax = max(max(green));
Bmax = max(max(blue));
MAX = max([Rmax Gmax Bmax]);

red_inv = MAX-red;
green_inv = MAX-green;
blue_inv = MAX-blue;

end


function [BW_out,rows,cols] = UserCrop(BW, Fig, supString)
% Show the input image, then have the user crop the region of interest.
% After this is done, store only the cropped portion of the image, and
% redisplay it
cropFlag = true;

while cropFlag
    makeFigures
    uiwait(Fig)
end

    function makeFigures
        figure(Fig)
        clf
        subplot(1,2,1)
        iptsetpref('ImshowAxesVisible','on')
        iptsetpref('ImshowBorder','loose')
        imshow(BW, 'InitialMagnification', 'fit')
        title({'Select bottom left point, then top right point';supString})
        S_corners = ginput(2);
        title('Good Job!')
        % Get the x and y corner coordinates as integers
        sp(1) = min(floor(S_corners(1)), floor(S_corners(2))); %xmin
        sp(2) = min(floor(S_corners(3)), floor(S_corners(4))); %ymin
        sp(3) = max(ceil(S_corners(1)), ceil(S_corners(2)));   %xmax
        sp(4) = max(ceil(S_corners(3)), ceil(S_corners(4)));   %ymax
        w = sp(3)-sp(1);
        h = sp(4)-sp(2);
        
        h_rect = rectangle('Position',[sp(1),sp(2),w,h]);
        h_rect.EdgeColor = 'r';
        h_rect.LineWidth = 2;
        pause(.1)
        
        rows = sp(2):sp(4);
        cols = sp(1):sp(3);
        BW_out = BW(rows, cols,:);
        
        subplot(1,2,2)
        h_im = imshow(BW_out, 'InitialMagnification', 'fit');
        
        title('Did it work?')
        % Creat buttons
        h_ax = gca;
        ax_pos = h_ax.Position;
        %Button Group
        bGroup = uibuttongroup('Parent',gcf,...
            'Units','normalized',...
            'Position',[ax_pos(1) ax_pos(2)-.05 ax_pos(3) .1]);
        %Buttons
        bYes = uicontrol(bGroup,'Style','pushbutton',...
            'String','Yes, move on!',...
            'Units','normalized',...
            'Position',[.1 .1 .3 .8],...
            'Callback',@yes);
        bNo = uicontrol(bGroup,'Style','pushbutton',...
            'String','No, repeat...',...
            'Units','normalized',...
            'Position',[.6 .1 .3 .8],...
            'Callback',@no);
    end

    function yes(varargin)
        cropFlag = false;
        uiresume(gcbf)
    end

    function no(varargin)
        cropFlag = true;
        uiresume(gcbf)
    end

end


function [pixelArea,pixelX,pixelY] = UserInfo(spacings, Fig, subfig)
% Allow the user to enter the line spacing of the grid. After this, the
% area of each pixel may be calculated

infoFlag = true;

while infoFlag
    [editX, editY, txtArea] = makeFigures;
    uiwait(Fig)
end

    function [editX, editY, txtArea] = makeFigures
        figure(Fig)
        subplot(2,2,4)
        % Creat buttons
        h_ax = gca;
        h_ax.Visible = 'off';
        ax_pos = h_ax.Position;
        %Button Group
        bGroup = uibuttongroup('Parent',gcf,...
            'Units','normalized',...
            'Position',[.505 .04 .4 .45]);
        %Buttons
        bEnter = uicontrol(bGroup,'Style','pushbutton',...
            'String','Calculate Area!',...
            'Units','normalized',...
            'Position',[.1 .4 .8 .1],...
            'Callback',@Enter);
        
        %Text Box
        txt = uicontrol(bGroup,'Style','text',...
            'String','Enter Grid Dimensions',...
            'Units','normalized',...
            'FontSize',14,...
            'Position',[.2 .85 .6 .1]);
        txtX = uicontrol(bGroup,'Style','text',...
            'String','X(mm):',...
            'Units','normalized',...
            'FontSize',14,...
            'Position',[.2 .7 .3 .1]);
        txtY = uicontrol(bGroup,'Style','text',...
            'String','Y(mm):',...
            'Units','normalized',...
            'FontSize',14,...
            'Position',[.2 .6 .3 .1]);
        
        txtPxl = uicontrol(bGroup,'Style','text',...
            'String','Pixel Area (mm^2): ',...
            'Units','normalized',...
            'FontSize',14,...
            'Position',[.15 .15 .3 .1]);
        txtArea = uicontrol(bGroup,'Style','text',...
            'String','???',...
            'Units','normalized',...
            'FontSize',14,...
            'Position',[0.450 0.1500 0.400 0.1000]);
        
        %Edit Box
        editX = uicontrol(bGroup,'Style','edit',...
            'String','5.08',...
            'Units','normalized',...
            'FontSize',14,...
            'Position',[.5 .7 .3 .1]);
        editY = uicontrol(bGroup,'Style','edit',...
            'String','5.08',...
            'Units','normalized',...
            'FontSize',14,...
            'Position',[.5 .6 .3 .1]);
        
    end

    function Enter(varargin)
        infoFlag = false;
        mmX = str2num(editX.String);
        mmY = str2num(editY.String);
        PixelArea = spacings(1)*spacings(2);
        mmArea = mmX*mmY;
        
        pixelArea = mmArea / PixelArea;     %area of each pixel
        pixelX = mmX/spacings(1);           %X dim of each pixel
        pixelY = mmY/spacings(2);           %Y dim of each pixel
        txtArea.String = {[num2str(pixelArea,'%-5.3E')]};
        
        uiresume(gcbf)
    end

end


function [pksX,  locsX, wX, pX, pksY,  locsY, wY, pY, spacingX, spacingY, spacing] = GridSpacing(img,ave, factor)

imgDC = img;
% imgDC(imgDC == ave)= NaN;

%Calculate line average
lineAveX = nanmean(imgDC,1);
lineAveY = nanmean(imgDC,2);

%Determine location of local extrema
[pksX,  locsX, wX, pX] = LineExtrema(lineAveX, factor);
[pksY,  locsY, wY, pY] = LineExtrema(lineAveY, factor);

%Determine pixel spacing between extrema
spacingX = round(mean2(diff(locsX',1)));
spacingY = round(mean2(diff(locsY',1)));
spacing = round(mean([spacingX spacingY]));

end


function [prominence] = minFunc(img,theta,ave)

S_rot_func = @(theta) imrotate(img,theta);
S_rot = S_rot_func(theta);
S_rot(S_rot == 0) = ave;
% S_rot(S_rot == 0) = min(min(min(img)));

factor = 3;
[pksX,  locsX, wX, pX, pksY,  locsY, wY, pY, spacingX, spacingY, spacing] =  GridSpacing(S_rot,ave, factor);
% peak = mean([pksX pksY]);
numPts = numel([pX pY]);
prominence = mean([pX, pY])*numPts;
% width = mean([wX wY]);

PlotGrids(S_rot, ave, spacingX, spacingY, locsX, locsY, theta, 1)
stdX = std(diff(locsX));
stdY = std(diff(locsY));
check1 = stdX/spacingX;
check2 = stdY/spacingY;
if (check1>1) || (check2>1) || (abs(spacingX-spacingY)>100)
    prominence = prominence / max([check1,check2,1]);
end
% prominence
pause(.1)
end


function PlotGrids(img, ave, spacingX, spacingY, locsX, locsY, theta, sampleFig)
spacing = round(mean([spacingX spacingY]));

figure(sampleFig)
clf
hAx1 = subplot(2,2,1);
h1 = imagesc(img);
lims = [hAx1.XLim hAx1.YLim];
hAx1.XTick = [];
hAx1.YTick = [];
hAx1.Position = [0.1    0.5    0.4    0.45];

X = hAx1.XLim;
Y = hAx1.YLim;

hold on
for i = 1:numel(locsX)
    x = round(mean(locsX(:,1)))+(i-1)*spacing;
    plot([x,x], Y,'r-')
end
for i = 1:numel(locsY)
    y = round(mean(locsY(:,1)))+(i-1)*spacing;
    plot(X, [y,y],'r-')
end
title(['\theta = ',num2str(theta)])
%Image Sizes
[M,N,~] = size(img);
imgDC = img;
imgDC(imgDC == ave)= NaN;

%Calculate line average
lineAveX = nanmean(imgDC,1);
lineAveY = nanmean(imgDC,2);

%Plot 2
hAx2 = subplot(2,2,2);
h2 = plot(1:M,lineAveY(:,:,1),'r',1:M,lineAveY(:,:,2),'g',1:M,lineAveY(:,:,3),'b');
view(90,90)
hAx2.XLim = hAx1.YLim;
hAx2.XTick = [];
hAx2.YTick = [];
hAx2.Position = [0.505    0.5    0.4    0.45];
grid on

X = hAx2.XLim;
Y = hAx2.YLim;

hold on
for i = 1:size(locsY,2)
    x = round(mean(locsY(:,1)))+(i-1)*spacing;
    plot([x,x],Y,'r-')
end
axis tight

%Plot 3
hAx3 = subplot(2,2,3);
h3 = plot(1:N,lineAveX(:,:,1),'r',1:N,lineAveX(:,:,2),'g',1:N,lineAveX(:,:,3),'b');
hAx3.XLim = hAx1.XLim;
hAx3.XTick = [];
hAx3.YTick = [];
hAx3.Position = [0.1    0.04    0.4    0.45];
grid on

X = hAx3.XLim;
Y = hAx3.YLim;

hold on
for i = 1:size(locsX,2)
    x = round(mean(locsX(:,1)))+(i-1)*spacing;
    plot([x,x],Y,'r-')
end
axis tight

end


function [pixels] = SamplePixelCount(img, Fig, median_num, median_res, threshold, refinement)
% To start, the current figure should be displaying a 1x2 subplot, with the
% rotated image in the first plot, and the zoomed in image of the sample in
% the second. The user will input a threshold value which then converts the
% rgb image to a bw image. Once done, the edges of the binary image are
% detected and plotted. The user can update the threshold value if not
% satistfied.

xy_bound = [];
finish_enable = 'off';
updatePlot = true;

while updatePlot
    makeFigures
    uiwait(Fig)
end

    function makeFigures
        figure(Fig)
        hold off
        subplot(1,2,2)
        imshow(img, 'InitialMagnification', 'fit')
        hold on
        if logical(numel(xy_bound))
            plot(xy_bound(:,2),xy_bound(:,1),'r','LineWidth',2)
        end
        title({'Adjust image analysis parameters';'then press "Update"';'When satisfied, press "Finish"'})

        % Create buttons
        h_ax = gca;
        ax_pos = h_ax.Position;
        %Button Group
        bGroup = uibuttongroup('Parent',gcf,...
            'Units','normalized',...
            'Position',[ax_pos(1) ax_pos(2)-.05 ax_pos(3) .18]);
        %Buttons
        bUpdate = uicontrol(bGroup,'Style','pushbutton',...
            'String','Update',...
            'Units','normalized',...
            'Position',[.1 .1 .3 .2],...
            'Callback',@Update);
        bFinish = uicontrol(bGroup,'Style','pushbutton',...
            'String','Finish',...
            'Units','normalized',...
            'Enable',finish_enable,...
            'Position',[.6 .1 .3 .2],...
            'Callback',@Finish);
        
        %Text
        txtMedNum = uicontrol(bGroup,'Style','text',...
            'String','Filter #',...
            'Units','normalized',...
            'FontSize',12,...
            'Position',[.0 .4 .25 .2]);
        
        txtMedRes = uicontrol(bGroup,'Style','text',...
            'String','Resolution',...
            'Units','normalized',...
            'FontSize',12,...
            'Position',[.25 .4 .25 .2]);
        
        txtThresh = uicontrol(bGroup,'Style','text',...
            'String','BW Threshold',...
            'Units','normalized',...
            'FontSize',12,...
            'Position',[.50 .4 .25 .2]);
        
        txtRefine = uicontrol(bGroup,'Style','text',...
            'String','Refinement',...
            'Units','normalized',...
            'FontSize',12,...
            'Position',[.75 .4 .25 .2]);
        
        %Edit
        txtMedNum = uicontrol(bGroup,'Style','edit',...
            'String',num2str(median_num),...
            'Tag','Median Number',...
            'Units','normalized',...
            'FontSize',12,...
            'Position',[0 .7 .25 .2]);
        
        txtMedRes = uicontrol(bGroup,'Style','edit',...
            'String',num2str(median_res),...
            'Tag','Median Resolution',...
            'Units','normalized',...
            'FontSize',12,...
            'Position',[.25 .7 .25 .2]);
        
        txtThresh = uicontrol(bGroup,'Style','edit',...
            'String',num2str(threshold),...
            'Tag','Threshold',...
            'Units','normalized',...
            'FontSize',12,...
            'Position',[.5 .7 .25 .2]);
        
        txtRefine = uicontrol(bGroup,'Style','edit',...
            'String',num2str(refinement),...
            'Tag','Refinement',...
            'Units','normalized',...
            'FontSize',12,...
            'Position',[.75 .7 .25 .2]);
        
    end

    function Update(varargin)
        actionData = varargin{2};
        buttonGroup = actionData.Source.Parent;
        children = buttonGroup.Children;
        median_num = str2num(children(strcmp({children.Tag},'Median Number')).String);
        median_res = str2num(children(strcmp({children.Tag},'Median Resolution')).String);
        threshold = str2num(children(strcmp({children.Tag},'Threshold')).String);
        refinement = str2num(children(strcmp({children.Tag},'Refinement')).String);
        
        %use a median filter to sharpen / define the edges and define an initial mask
        median = rgb2gray(img);
        for i = 1:median_num
            median = medfilt2(median,[median_res median_res],'symmetric');
            %     disp(i)
            %     if any(i == [10:10:100])
            %         subplot(1,2,2)
            %         imshow(median, 'InitialMagnification', 'fit')
            %         colormap jet
            %         pause(.1)
            %     end
        end
        
        %Binary Image
        S_BW = im2bw(median,threshold);
        
        %Convex Hull
        CH = bwconvhull(~S_BW);
        
        %Make Mask (set pixels outside convex hull to zero)
        [M,N] = size(CH);
        S_mask = median;
        for i = 1:M
            for j = 1:N
                if CH(i,j)==0
                    S_mask(i,j,:)=0;
                end
            end
        end
        %Convert Mask to Binary
        S_mask_BW = im2bw(S_mask,.0001);
        %Get active contour (outline) of region inside of mask
        bw = activecontour(S_mask,S_mask_BW,refinement ,'edge');
        
        %Edge Detection
        S_boundaries = bwboundaries(bw,'noholes');
        xy_bound = S_boundaries{1};
        
        %Determine cross sectional area
        pixels = sum(double(bw(:)));
        
        %find / delete old line
        h = findobj(gca, 'color','r');
        delete(h)
        %plot
        plot(xy_bound(:,2),xy_bound(:,1),'r.','LineWidth',2)
        
        children(strcmp({children.String},'Finish')).Enable = 'on';
        finish_enable = 'on';
        uiresume(gcbf)
    end

    function Finish(varargin)
        updatePlot = false;
        uiresume(gcbf)
    end

end


