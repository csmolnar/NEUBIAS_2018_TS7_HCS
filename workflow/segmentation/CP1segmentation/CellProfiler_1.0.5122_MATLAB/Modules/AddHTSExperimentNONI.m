function handles = AddHTSExperiment(handles)
% Help for the AddHTSExperiment module:
% Category: Other
%
% SHORT DESCRIPTION:
% Create web pages from High-throughput screening images.
% For more information visit: 
% http://www.lmc.ethz.ch/People/PeterHorvath
%
% Author: Peter Horvath, Light Microscopy Centre, ETH Zurich
% *************************************************************************
% $Revision: 1000 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%pathnametextVAR01 = Please selet the folder contains main.html?
%defaultVAR01 = .
BaseFolder = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Enter the name of the experiment
%defaultVAR02 = Experiment1
ProjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});


%textVAR03 = Choose the well type you used?
%choiceVAR03 = 96 (12x8)
%choiceVAR03 = 384 (24x16)
WellType = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu custom

%textVAR04 = Splitting for visualization purposes. Choose 1x1 if you want to see your plate on one page 2x2 if you want to split it. 2x2 recommended for 384 well plates of for smaller displays. 
%choiceVAR04 = 1x1
%choiceVAR04 = 2x2
VisualizationType = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu custom

%pathnametextVAR05 = Selet the folder contains your experiment?
%defaultVAR05 = .
BaseExpFolder = char(handles.Settings.VariableValues{CurrentModuleNum,5});


%textVAR06 = Type the text that the first (red) channel of image has in common?
%defaultVAR06 = 
ChannelCommon1 = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Type the text that the second (green) channel of image has in common?
%defaultVAR07 = tritc
ChannelCommon2 = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = Type the text that the third (blue) channel of image has in common?
%defaultVAR08 = dapi
ChannelCommon3 = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR9 = Choose resize method (nearest is two orders of magnitude faster). 
%choiceVAR09 = nearest
%choiceVAR09 = bilinear
%choiceVAR09 = bicubic
ImageResizeMethod = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu custom

%%%VariableRevisionNumber = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY INPUT CHECK                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

if strcmp(BaseFolder, '.')
    disp('No folder had been given.');
end;

switch WellType
    case '96 (12x8)'
        WellTypeCode = 1;
    case '384 (24x16)'
        WellTypeCode = 2;
    otherwise
        disp('Unknown well type.')
end;

switch VisualizationType
    case '1x1'
        VisualizationTypeCode = 1;
    case '2x2'
        VisualizationTypeCode = 2;
    otherwise
        disp('Unknown visualization type.')
end;

warning off;

% check whether we are in batch mode
if ~isfield(handles.Current, 'BatchInfo'),
    
    % Add link to the main page
    add_html_link([BaseFolder '\main.html'], ProjectName, ['_data\' ProjectName '\plates.html']);

    % create a new directory
    mkdir([BaseFolder '\_data\' ProjectName]);

    % copy plate page
    copyfile([BaseFolder '\source\plate_template.html'], [BaseFolder '\_data\' ProjectName '\plates.html']);
    copyfile([BaseFolder '\source\platelist.css'], [BaseFolder '\_data\' ProjectName '\platelist.css']);

    % list the projectpath folrder, add plates to plate html
    dir_list = dir(BaseExpFolder);
    
    for i=3:size(dir_list, 1)

        % time left information
        StartTime = clock;           
        
        if (dir_list(i).isdir)
            % create folder
            mkdir([BaseFolder '\_data\' ProjectName '\' dir_list(i).name]);

            if VisualizationTypeCode == 1
                add_html_link([BaseFolder '\_data\' ProjectName '\plates.html'], dir_list(i).name, [dir_list(i).name '\thumbnails.html' ]);
                % copy thumbnail html and css
                copyfile([BaseFolder '\source\thumbnails_template.html'], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails.html']);
                copyfile([BaseFolder '\source\plate.css'], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\plate.css']);
                copyfile([BaseFolder '\source\back.gif'], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\back.gif']);

                add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails.html'], '<table border="0">');

            elseif VisualizationTypeCode == 2
                add_html_link([BaseFolder '\_data\' ProjectName '\plates.html'], dir_list(i).name, [dir_list(i).name '\thumbnails_1.html' ]);
                % copy thumbnail html and css
                for wellnum=1:4
                    copyfile([BaseFolder '\source\thumbnails_template.html'], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_' int2str(wellnum) '.html']);
                end;
                % cross links
                if WellTypeCode == 1
                    copyfile([BaseFolder '\source\well_6_4_w.gif '], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\well_6_4_w.gif']);
                    copyfile([BaseFolder '\source\well_6_4_r.gif '], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\well_6_4_r.gif']);
                    for wellnum =1:4
                        for imagenum=1:4
                            if imagenum == wellnum
                                add_html_image([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_' int2str(wellnum) '.html'], 'well_6_4_r.gif', ['thumbnails_' int2str(imagenum) '.html'], 2);
                            else
                                add_html_image([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_' int2str(wellnum) '.html'], 'well_6_4_w.gif', ['thumbnails_' int2str(imagenum) '.html'], 2);
                            end;
                            if imagenum == 2
                                add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_' int2str(wellnum) '.html'], '<br>');
                            end;
                        end;
                    end;
                elseif WellTypeCode == 2
                    copyfile([BaseFolder '\source\well_12_8_w.gif '], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\well_12_8_w.gif']);
                    copyfile([BaseFolder '\source\well_12_8_r.gif '], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\well_12_8_r.gif']);
                    for wellnum =1:4
                        for imagenum=1:4
                            if imagenum == wellnum
                                add_html_image([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_' int2str(wellnum) '.html'], 'well_12_8_r.gif', ['thumbnails_' int2str(imagenum) '.html'], 2);
                            else
                                add_html_image([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_' int2str(wellnum) '.html'], 'well_12_8_w.gif', ['thumbnails_' int2str(imagenum) '.html'], 2);
                            end;
                            if imagenum == 2
                                add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_' int2str(wellnum) '.html'], '<br>');
                            end;
                        end;
                    end;
                end;

                for wellnum=1:4
                    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_' int2str(wellnum) '.html'], '<table border="0">');
                end;
                copyfile([BaseFolder '\source\plate.css'], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\plate.css']);
                copyfile([BaseFolder '\source\back.gif'], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\back.gif']);
            end;

            if WellTypeCode == 1
                if VisualizationTypeCode == 1
                    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails.html'], '<tr> <td></td> <td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td></tr>');
                elseif VisualizationTypeCode == 2
                    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_1.html'], '<tr> <td></td> <td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td></tr>');
                    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_3.html'], '<tr> <td></td> <td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td></tr>');
                    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_2.html'], '<tr> <td></td> <td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td></tr>');
                    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_4.html'], '<tr> <td></td> <td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td></tr>');
                end;
            elseif WellTypeCode == 2
                if VisualizationTypeCode == 1
                    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails.html'], '<tr> <td></td> <td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td><td>13</td><td>14</td><td>15</td><td>16</td><td>17</td><td>18</td><td>19</td><td>20</td><td>21</td><td>22</td><td>23</td><td>24</td></tr>');
                elseif VisualizationTypeCode == 2
                    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_1.html'], '<tr> <td></td> <td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td></tr>');
                    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_3.html'], '<tr> <td></td> <td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td></tr>');
                    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_2.html'], '<tr> <td></td> <td>13</td><td>14</td><td>15</td><td>16</td><td>17</td><td>18</td><td>19</td><td>20</td><td>21</td><td>22</td><td>23</td><td>24</td></tr>');
                    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_4.html'], '<tr> <td></td> <td>13</td><td>14</td><td>15</td><td>16</td><td>17</td><td>18</td><td>19</td><td>20</td><td>21</td><td>22</td><td>23</td><td>24</td></tr>');
                end;
            end;

             handles = add_experiment(handles, WellTypeCode, VisualizationTypeCode, BaseFolder, ProjectName, dir_list(i).name, BaseExpFolder, ChannelCommon1, ChannelCommon2, ChannelCommon3, ImageResizeMethod);

            if VisualizationTypeCode == 1
                add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails.html'], '</table>');
            elseif VisualizationTypeCode == 2
                add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_1.html'], '</table>');
                add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_3.html'], '</table>');
                add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_2.html'], '</table>');
                add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails_4.html'], '</table>');
            end;

        end;
        SpentTime = etime(clock, StartTime);
         timertext = sprintf('Estimated time left: %f sec.', SpentTime * (size(dir_list, 1)-i));
         disp(timertext);
    end;
    % decrease the value of the sets by one, because it starts with 1
    % instead of 0.
else
    % batch mode is active, read lists of images and convert them
    
    % find AddHTSExperiment in the pipeline
    
    ModulePosition = find(strcmp(handles.Settings.ModuleNames ,{'AddHTSExperiment'}));
    
    % change path name to the one defined in Batch module
    RemoveLength = length(handles.Settings.VariableValues{ModulePosition, 5});
    path = handles.Pipeline.ImageFileNames{handles.Current.SetBeingAnalyzed}.Path;
    path = path(RemoveLength+1:length(path));
    path = [handles.Current.DefaultImageDirectory path];
    % change backslash to unix mode
    BSPos = strfind(path, '\'); path(BSPos) = '/';
       
    %change output name for Bach purpose
    RemoveLength = length(handles.Settings.VariableValues{ModulePosition, 1});
    outPrefix = handles.Pipeline.ImageFileNames{handles.Current.SetBeingAnalyzed}.OutPrefix;
    outPrefix = outPrefix(RemoveLength+1:length(outPrefix));
    outPrefix = [handles.Current.DefaultOutputDirectory outPrefix];
    % change backslash to unix mode
    BSPos = strfind(outPrefix, '\'); outPrefix(BSPos) = '/';

    imgprefix   = handles.Pipeline.ImageFileNames{handles.Current.SetBeingAnalyzed}.ImagePrefix;
    dir_list_r  = handles.Pipeline.ImageFileNames{handles.Current.SetBeingAnalyzed}.ImageListRed;
    dir_list_g  = handles.Pipeline.ImageFileNames{handles.Current.SetBeingAnalyzed}.ImageListGreen;
    dir_list_b  = handles.Pipeline.ImageFileNames{handles.Current.SetBeingAnalyzed}.ImageListBlue;

    create_thumbnail(path, outPrefix, dir_list_r, dir_list_g, dir_list_b, imgprefix, ImageResizeMethod);
end;

%% add_html_link
function add_html_link(fileName, text, linkto)

fid1=fopen(fileName, 'r');
lines = '';
comment = '<! --MATLAB-->';

while 1
    tline = fgets(fid1);
    if ~ischar(tline),   break,   end
    lines = [lines  tline];
end

fclose(fid1);
a = strfind(lines, comment);

if length(a) ~= 1
    disp('error');
end;
newtext = '';


if length(linkto) == 0
    newtext = ['<p>' text '</p>'];
else
    newtext = ['<li><a href = "' linkto '">' text '</a>'];
end;

lines = [lines(1:a-1) newtext lines(a : length(lines))];
fid1=fopen(fileName, 'w');
fprintf(fid1, '%s', lines);
fclose(fid1);

%% add_html_text
function add_html_text(fileName, text);

fid=fopen(fileName);

lines = '';

comment = '<! --MATLAB-->';

while 1

    tline = fgets(fid);
    if ~ischar(tline),   break,   end
    
    lines = [lines  tline];
    
end

fclose(fid);

a = strfind(lines, comment);

if length(a) ~= 1
   
    disp('error');
    
end;

newtext = text;
    
lines = [lines(1:a-1) newtext lines(a : length(lines))];

fid=fopen(fileName, 'w');

fprintf(fid, '%s', lines);

fclose(fid);

%% add html image
function add_html_image(fileName, imagelocation, linkto, newpage);

fid=fopen(fileName);

lines = '';

comment = '<! --MATLAB-->';

while 1

    tline = fgets(fid);
    if ~ischar(tline),   break,   end
    
    lines = [lines  tline];
    
end

fclose(fid);

a = strfind(lines, comment);

if length(a) ~= 1
   
    disp('error');
    
end;

newtext = '';


if length(linkto) == 0

    newtext = ['<img border="0" src="' imagelocation '">'];

else
    
    if newpage == 1
    
        newtext = ['<a href="' linkto '" target="_blank"><img border="0" src="' imagelocation '"></a>'];

    elseif newpage == 2

        newtext = ['<a href="' linkto '" target="_self"><img border="0" src="' imagelocation '"></a>'];
        
    elseif newpage == 3

        newtext = ['<a href="' linkto '" target="_parent"><img border="0" src="' imagelocation '"></a>'];
        
    end;
        
end;

lines = [lines(1:a-1) newtext lines(a : length(lines))];

fid=fopen(fileName, 'w');

fprintf(fid, '%s', lines);

fclose(fid);


function handles = add_experiment(handles, WellTypeCode, VisualizationTypeCode, BaseFolder, ProjectName, dir_list_name, BaseExpFolder, ChannelCommon1, ChannelCommon2, ChannelCommon3, ImageResizeMethod);

if WellTypeCode == 1
    max_c = 'H';
    max_r = 12;
elseif WellTypeCode == 2
    max_c = 'P';
    max_r = 24;
end;



for y='A':max_c

    % infectX extension
    %y = 'wA';
    
    if VisualizationTypeCode == 1

        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails.html'], ['<tr><td>' y '</td>']);

    end;

    for x = 1:max_r

        if VisualizationTypeCode == 2 %% 2x2
            if WellTypeCode == 1 %% 96 well
                if y < 'E'
                    if x == 1
                        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_1.html'], ['<tr><td>' y '</td>']);
                    elseif x == 7
                        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_2.html'], ['<tr><td>' y '</td>']);
                    end;

                else
                    if x == 1
                        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_3.html'], ['<tr><td>' y '</td>']);
                    elseif x == 7
                        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_4.html'], ['<tr><td>' y '</td>']);
                    end;
                end;

            elseif WellTypeCode == 2 %% 384 well
                if y < 'I'
                    if x == 1
                        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_1.html'], ['<tr><td>' y '</td>']);
                    elseif x == 13
                        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_2.html'], ['<tr><td>' y '</td>']);
                    end;

                else
                    if x == 1
                        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_3.html'], ['<tr><td>' y '</td>']);
                    elseif x == 13
                        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_4.html'], ['<tr><td>' y '</td>']);
                    end;
                end;

            end;

        end;



        if x<10
            smallname = sprintf('%s_%s0%d.jpg', dir_list_name, y, x);
            imgprefix = sprintf('%s_%s0%d', dir_list_name, y, x);
            htmlname = sprintf('%s_%s0%d.html', dir_list_name, y, x);

        else

            smallname = sprintf('%s_%s%d.jpg', dir_list_name, y, x);
            imgprefix = sprintf('%s_%s%d', dir_list_name, y, x);
            htmlname = sprintf('%s_%s%d.html', dir_list_name, y, x);

        end;

        % create the small images
        % first read all the images in the well and make a global
        % statistics
        % create prefix for list

        path = [BaseExpFolder '\' dir_list_name '\'];

        if (x<10)

            outPrefix = [BaseFolder '\_data\' ProjectName '\' dir_list_name '\' sprintf('%s_%s0%d', dir_list_name, y, x)];

        else

            outPrefix = [BaseFolder '\_data\' ProjectName '\' dir_list_name '\' sprintf('%s_%s%d', dir_list_name, y, x)];

        end;

        % for this we should find a better naming since dapi trict are not nice, maybe the user should define them
        if (x<10)

            imageprefix_r = sprintf('%s%s_*%s0%d*_%s*.tif', path, dir_list_name, y, x, ChannelCommon1);

            imageprefix_g = sprintf('%s%s_*%s0%d*_%s*.tif', path, dir_list_name, y, x, ChannelCommon2);

            imageprefix_b = sprintf('%s%s_*%s0%d*_%s*.tif', path, dir_list_name, y, x, ChannelCommon3);

        else

            imageprefix_r = sprintf('%s%s_*%s%d*_%s*.tif', path, dir_list_name, y, x, ChannelCommon1);

            imageprefix_g = sprintf('%s%s_*%s%d*_%s*.tif', path, dir_list_name, y, x, ChannelCommon2);

            imageprefix_b = sprintf('%s%s_*%s%d*_%s*.tif', path, dir_list_name, y, x, ChannelCommon3);

        end;

        dir_list_r = dir(imageprefix_r);

        dir_list_g = dir(imageprefix_g);

        dir_list_b = dir(imageprefix_b);

        exist_r = size(dir_list_r, 1);

        exist_g = size(dir_list_g, 1);

        exist_b = size(dir_list_b, 1);


        % create medium images to see
        copyfile([BaseFolder '\source\wellview_template.html'], [BaseFolder '\_data\' ProjectName '\' dir_list_name '\' htmlname]);
        
        % if CreateBatchFiles module is present than add images to a list
        % of tasks otherwise process them
        
        IsCreateBatchFilesPresent = 0;
        for i=1:handles.Current.NumberOfModules
            if strcmp(handles.Settings.ModuleNames{i}, 'CreateBatchFiles') == 1
                IsCreateBatchFilesPresent = 1;
            end;
        end;

        if IsCreateBatchFilesPresent 
            handles.Current.NumberOfImageSets = handles.Current.NumberOfImageSets + 1;
            handles.Pipeline.ImageFileNames{handles.Current.NumberOfImageSets}.Path             = path;
            handles.Pipeline.ImageFileNames{handles.Current.NumberOfImageSets}.OutPrefix        = outPrefix;            
            handles.Pipeline.ImageFileNames{handles.Current.NumberOfImageSets}.ImagePrefix      = imgprefix;
            handles.Pipeline.ImageFileNames{handles.Current.NumberOfImageSets}.ImageListRed     = dir_list_r;
            handles.Pipeline.ImageFileNames{handles.Current.NumberOfImageSets}.ImageListGreen   = dir_list_g;
            handles.Pipeline.ImageFileNames{handles.Current.NumberOfImageSets}.ImageListBlue    = dir_list_b;            
        else
            create_thumbnail(path, outPrefix, dir_list_r, dir_list_g, dir_list_b, imgprefix, ImageResizeMethod);
        end;

        if x<10
            basename = sprintf('%s_%s0%d', dir_list_name, y, x);
        else
            basename = sprintf('%s_%s%d', dir_list_name, y, x);
        end;

        % add thumbnail image
        if VisualizationTypeCode == 1

            add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails.html'], '<td>');
            add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails.html'], ['<div id="pic"><a class="thumbnail"  href="' htmlname '" target="_blank" ><img src="' smallname '" /><span><img  src=" ' basename 'montage.jpg"    /></span></a></div>']);
            add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails.html'], '</td>');

        elseif VisualizationTypeCode == 2
            if WellTypeCode == 1 %% 96 well
                if y < 'E'
                    if x < 7
                        pagenum=1;
                    else
                        pagenum=2;
                    end;
                else
                    if x < 7
                        pagenum=3;
                    else
                        pagenum=4;
                    end;
                end;
            else
                if y < 'I'
                    if x < 13
                        pagenum=1;
                    else
                        pagenum=2;
                    end;
                else
                    if x < 13
                        pagenum=3;
                    else
                        pagenum=4;
                    end;
                end;
            end;
            add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_' num2str(pagenum) '.html'], '<td>');
            add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_' num2str(pagenum) '.html'], ['<div id="pic"><a class="thumbnail"  href="' htmlname '" target="_blank" ><img src="' smallname '" /><span><img  src=" ' basename 'montage.jpg" /><br></span></a></div>']);
            add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_' num2str(pagenum) '.html'], '</td>');

        end;

        % end table row
        if VisualizationTypeCode == 1
            if x == max_r
                add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails.html'], '</tr>');
            end;
        elseif VisualizationTypeCode == 2
            pagenum = 0;
            if WellTypeCode == 1 %% 96 well
                if x == 6
                    if y < 'E'
                        pagenum=1;
                    else
                        pagenum=3;
                    end;
                elseif x == 12
                    if y > 'E'
                        pagenum=2;
                    else
                        pagenum=4;
                    end;
                end;
            elseif WellTypeCode == 2 %% 384 well
                if x == 12
                    if y < 'I'
                        pagenum=1;
                    else
                        pagenum=3;
                    end;
                elseif x == 24
                    if y > 'I'
                        pagenum=2;
                    else
                        pagenum=4;
                    end;
                end;

            end;
            if pagenum ~= 0
                add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails_' num2str(pagenum) '.html'],  '</tr>');
            end;
        end;
    end;
end;

function out = create_thumbnail(prefix, outPrefix, r_list, g_list, b_list, imgprefix, ImageResizeMethod);

% load images
% find the nonzero channels length
exist_r = size(r_list, 1); 
exist_g = size(g_list, 1); 
exist_b = size(b_list, 1);

list_length = max([exist_r exist_g exist_b]);

% read all the images 
for i=1:list_length    
    if exist_r        
        rch(:, :, i) = imread([prefix r_list(i).name])';        
    else
        rch = 0;        
    end;

    if exist_g        
        gch(:, :, i) = imread([prefix g_list(i).name])';        
    else        
        gch=0;        
    end;
    
    if exist_b        
        bch(:, :, i) = imread([prefix b_list(i).name])';        
    else
        bch=0;        
    end;
end;

% normalize the images
% min max values

if exist_r               
    rch = rch - min(rch(:));
    rmax = max(rch(:));
    if rmax ~= 0
        shiftBits = floor(log2(double(rmax))-7);
        rch = bitshift(rch, -shiftBits);
    end;
end;

if exist_g               
    gch = gch - min(gch(:));
    gmax = max(gch(:));
    if gmax ~= 0           
        shiftBits = floor(log2(double(gmax))-7);
        gch = bitshift(gch, -shiftBits);
    end;
end;

if exist_b
    bch = bch - min(bch(:));
    bmax = max(bch(:));
    if bmax ~= 0
        shiftBits = floor(log2(double(bmax))-7);
        bch = bitshift(bch, -shiftBits);
    end;
end;

% prepare montage images
xdim = max([size(rch,1) size(gch,1) size(bch,1)]);
ydim = max([size(rch,2) size(gch,2) size(bch,2)]);

if exist_r
    rch = reshape(rch, xdim, ydim, 1, list_length);    
    [h, r_mon] = montageClone(rch);    
    r_mon = imresize(r_mon, [350 NaN], ImageResizeMethod);
    [xdim_mon, ydim_mon]  = size(r_mon);
end;

if exist_g
    gch = reshape(gch, xdim, ydim, 1, list_length);    
    [h, g_mon] = montageClone(gch);    
    g_mon = imresize(g_mon, [350 NaN], ImageResizeMethod);    
    [xdim_mon, ydim_mon]  = size(g_mon);    
end;

if exist_b
    bch = reshape(bch, xdim, ydim, 1, list_length);    
    [h, b_mon] = montageClone(bch);    
    b_mon = imresize(b_mon, [350 NaN], ImageResizeMethod);    
    [xdim_mon, ydim_mon]  = size(b_mon);    
end;

if xdim_mon >1
    fullmontage = uint8(zeros(xdim_mon, ydim_mon, 3));
    if exist_r
        fullmontage(:,:,1) = bitshift(r_mon, 1);    
    end;
    if exist_g
        fullmontage(:,:,2) = bitshift(g_mon, 1);    
    end;
    if exist_b
        fullmontage(:,:,3) = bitshift(b_mon, 1);    
    end;
    % save it
    imwrite(fullmontage, [outPrefix 'montage.jpg']);
end;

% make the small image

rand_num = floor(rand * list_length) + 1;

vimg = uint8(zeros(ydim, xdim, 3));
if exist_r
    vimg(:,:,1) = bitshift(rch(:,:,rand_num)', 1);
end;
if exist_g
    vimg(:,:,2) = bitshift(gch(:,:,rand_num)', 1);
end;
if exist_b
    vimg(:,:,3) = bitshift(bch(:,:,rand_num)', 1);
end;
vimg = imresize(vimg, [100 NaN], ImageResizeMethod);   

imwrite(vimg, [outPrefix '.jpg']);

% make the view images
for i=1:list_length

    % create image
    outimg = uint8(zeros(ydim, xdim, 3));
    if exist_r
        outimg(:,:,1) = bitshift(rch(:,:,i)', 1); 
    end;
    if exist_g
        outimg(:,:,2) = bitshift(gch(:,:,i)', 1);
    end;
    if exist_b
        outimg(:,:,3) = bitshift(bch(:,:,i)', 1);
    end;
    outimg = imresize(outimg, [300 NaN], ImageResizeMethod);   

    name = sprintf('%s_0%d.jpg', outPrefix, i);
        
    imwrite(outimg, name);

    name = sprintf('%s_0%d.jpg', imgprefix, i);
    
    add_html_image([outPrefix '.html'], name , '');
        
end;

function [h, val] = montageClone(V)

[r c s]=size(V);
h=ceil(sqrt(s));
w=ceil(s/h);
val=V;
fill=0;
if h*w~=s
    val(r,c,h*w)=fill; % pads empty region out with some value
end
val=reshape(val,[r c h w]);
val=permute(val,[1 3 2 4]);
val=reshape(val,[r*h c*w]);
val = val';

