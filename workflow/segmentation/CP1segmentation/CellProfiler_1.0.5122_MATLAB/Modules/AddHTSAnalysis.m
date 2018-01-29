function handles = AddHTSAnalysis(handles)

% Help for the AddHTSAnalysis module:
% Category: Other
%
% SHORT DESCRIPTION:
% Add meta-data to web pages created with AddHTSExperiment module.
% For more information visit: 
% http://www.lmc.ethz.ch/People/PeterHorvath
%
% Author: Peter Horvath, Light Microscopy Centre, ETH Zurich
% *************************************************************************
% $Revision: 1000 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

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

%pathnametextVAR04 = Selet the folder contains your experiment?
%defaultVAR04 = .
BaseExpFolder = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Type the name of the folders contains the metadata?
%defaultVAR05 = .
MetadataFolder = char(handles.Settings.VariableValues{CurrentModuleNum,5});

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
        WellTypeCode = 1
    case '384 (24x16)'
        WellTypeCode = 2
    otherwise
        disp('Unknown well type.')
end;

% list the folder of the project in the viewer directory
BaseViewFolder = [BaseFolder '\_data\' ProjectName];
dir_list = dir(BaseViewFolder);

for i=3:size(dir_list, 1)
     
     if (dir_list(i).isdir)               
    
         % list html files in the folder
         CurrentBaseViewFolder = [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\*.html'];
         dir_html_list = dir(CurrentBaseViewFolder);
         
         for j=1:size(dir_html_list, 1)
             % list all files in the analysis folder having the same prefix
             % as the html file
             [pathstr, prefix, ext, versn] = fileparts(dir_html_list(j).name);
             OriginalFolderList = [BaseExpFolder '\' dir_list(i).name '\' MetadataFolder '\' prefix '*'] ;
             dir_analysis_list = dir(OriginalFolderList);
             for k=1:size(dir_analysis_list, 1)
                 add_html_link( [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\' dir_html_list(j).name] , dir_analysis_list(k).name ,[BaseExpFolder '\' dir_list(i).name '\' MetadataFolder '\' dir_analysis_list(k).name] );
                %add_html_image( [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\' dir_html_list(j).name] ,[BaseExpFolder '\' dir_list(i).name '\' MetadataFolder '\' dir_analysis_list(k).name] , '');                                 
             end;
             
         end;
         
         
%         % create folder
%         mkdir([BaseFolder '\_data\' ProjectName '\' dir_list(i).name]);
%         % add link to the platelist
%         add_html_link([BaseFolder '\_data\' ProjectName '\plates.html'], dir_list(i).name, [dir_list(i).name '\thumbnails.html' ]); 
%         % copy thumbnail html and css
%         copyfile([BaseFolder '\source\thumbnails_template.html'], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails.html']);        
%         copyfile([BaseFolder '\source\plate.css'], [BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\plate.css']);   
% 
%         add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails.html'], '<table border="0">');   
% 
%         if WellTypeCode == 1
%             add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails.html'], '<tr> <td></td> <td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td></tr>');
%         elseif WellTypeCode == 2
%             add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list(i).name '\thumbnails.html'], '<tr> <td></td> <td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td><td>13</td><td>14</td><td>15</td><td>16</td><td>17</td><td>18</td><td>19</td><td>20</td><td>21</td><td>22</td><td>23</td><td>24</td></tr>');
%         end;
        
%         switch MicroscopeType
%             case 'CellWorx'
%                 add_CellWorx(WellTypeCode, BaseFolder, ProjectName, dir_list(i).name, BaseExpFolder, ChannelCommon1, ChannelCommon2, ChannelCommon3)
%             case 'MD micro'
%                 add_MD_micro(WellTypeCode)
%             case 'BD pathway'
%                 add_BD_pathway(WellTypeCode)
%             otherwise
%                 disp('Unknown microscope.')                
%         end;
                        
     end;
     
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
function add_html_image(fileName, imagelocation, linkto);

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
    
    newtext = ['<a href="' linkto '" target="_blank"><img border="0" src="' imagelocation '"></a>'];
    
end;

lines = [lines(1:a-1) newtext lines(a : length(lines))];

fid=fopen(fileName, 'w');

fprintf(fid, '%s', lines);

fclose(fid);


function add_CellWorx(WellTypeCode, BaseFolder, ProjectName, dir_list_name, BaseExpFolder, ChannelCommon1, ChannelCommon2, ChannelCommon3)

if WellTypeCode == 1
    max_c = 'H';
    max_r = 12;
elseif WellTypeCode == 2
    max_c = 'P';
    max_r = 24;
end;

for y='A':max_c

    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails.html'], ['<tr><td>' y '</td>']);

    for x = 1:max_r

         if x<10
% %             dapiname = sprintf('%s_%s0%d_010_dapi.tif', dir_list_name, y, x);
% % 
% %             tritcname = sprintf('%s_%s0%d_010_tritc.tif', dir_list_name, y, x);
% % 
             smallname = sprintf('%s_%s0%d.jpg', dir_list_name, y, x);
% % 
             htmlname = sprintf('%s_%s0%d.html', dir_list_name, y, x);
% % 
         else
% %             dapiname = sprintf('%s_%s%d_010_dapi.tif', dir_list_name, y, x);
% % 
% %             tritcname = sprintf('%s_%s%d_010_tritc.tif', dir_list_name, y, x);
% % 
             smallname = sprintf('%s_%s%d.jpg', dir_list_name, y, x);
% % 
             htmlname = sprintf('%s_%s%d.html', dir_list_name, y, x);
% % 
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

            imageprefix_r = sprintf('%s%s_%s0%d*_%s.tif', path, dir_list_name, y, x, ChannelCommon1);

            imageprefix_g = sprintf('%s%s_%s0%d*_%s.tif', path, dir_list_name, y, x, ChannelCommon2);

            imageprefix_b = sprintf('%s%s_%s0%d*_%s.tif', path, dir_list_name, y, x, ChannelCommon3);

        else

            imageprefix_r = sprintf('%s%s_%s%d*_%s.tif', path, dir_list_name, y, x, ChannelCommon1);

            imageprefix_g = sprintf('%s%s_%s%d*_%s.tif', path, dir_list_name, y, x, ChannelCommon2);

            imageprefix_b = sprintf('%s%s_%s%d*_%s.tif', path, dir_list_name, y, x, ChannelCommon3);

        end;
        
        dir_list_r = dir(imageprefix_r);

        dir_list_g = dir(imageprefix_g);
        
        dir_list_b = dir(imageprefix_b);        
                
        exist_r = size(dir_list_r, 1);
        
        exist_g = size(dir_list_g, 1); 
        
        exist_b = size(dir_list_b, 1);

        
        % create medium images to see

        copyfile([BaseFolder '\source\wellview_template.html'], [BaseFolder '\_data\' ProjectName '\' dir_list_name '\' htmlname]);
                
        create_thumbnail(path, outPrefix, dir_list_r, dir_list_g, dir_list_b);
        
        % list
        
        % create small thumbnail
        %
        %                   small = create_thumbnail('', [projectpath '\' dir_list(i).name '\' tritcname], [projectpath '\' dir_list(i).name '\' dapiname], 0.2);
        %
        %                   imwrite(small, [databasePath '_data\' projectname '\' dir_list(i).name '\' smallname]);
        %

        % copy the wellview html
        
        %create the detailed view images
        %                 for counter = 0:19
        %
        %                     if x<10
        %
        %                         dapiname = sprintf('%s_%s0%d_0%d_dapi.tif', dir_list_name, y, x, counter);
        %
        %                         tritcname = sprintf('%s_%s0%d_0%d_tritc.tif', dir_list_name, y, x, counter);
        %
        %                         imagename = sprintf('%s_%s0%d_0%d.jpg', dir_list_name, y, x, counter);
        %
        %                     else
        %
        %                         dapiname = sprintf('%s_%s%d_0%d_dapi.tif', dir_list_name, y, x, counter);
        %
        %                         tritcname = sprintf('%s_%s%d_0%d_tritc.tif', dir_list_name, y, x, counter);
        %
        %                         imagename = sprintf('%s_%s%d_0%d.jpg', dir_list_name, y, x, counter);
        %
        %                     end;
        %
        %
        %
        % %                         img = create_thumbnail('', [projectpath '\' dir_list(i).name '\' tritcname], [projectpath '\' dir_list(i).name '\' dapiname], 0.4);
        % %
        % %                         imwrite(img, [databasePath '_data\' projectname '\' dir_list(i).name '\' imagename]);
        %
%%%                                 add_html_image([BaseFolder '\_data\' ProjectName '\' dir_list_name '\' htmlname], imagename , '');
        %
        %
        %                 end;
        if x<10
            basename = sprintf('%s_%s0%d', dir_list_name, y, x);
        else
            basename = sprintf('%s_%s%d', dir_list_name, y, x);
        end;
        %
        %                 list = dir([databasePath '_data\' projectname '\' dir_list(i).name '\' basename '_*.jpg']);
        %
        %                 for ccc = 1:size(list, 1)
        %
        %                     list(ccc).name = [databasePath '_data\' projectname '\' dir_list(i).name '\' list(ccc).name];
        %
        %                 end;
        %
        %                 fileNames = {list.name}';
        %
        %                 [h, mimg] = montage(fileNames);
        %
        %                 mimg = imresize(mimg, 0.4);
        %
        %                 imwrite(mimg, [databasePath '_data\' projectname '\' dir_list(i).name '\' basename 'montage.jpg']);

        % add thumbnail image
        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails.html'], '<td>');

        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails.html'], ['<div id="pic"><a class="thumbnail"  href="' htmlname '" target="_blank" ><img src="' smallname '" /><span><img  src=" ' basename 'montage.jpg"    /></span></a></div>']);

        add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails.html'], '</td>');


    end;
    add_html_text([BaseFolder '\_data\' ProjectName '\' dir_list_name '\thumbnails.html'], '</tr>');

end;

function out = create_thumbnail(prefix, outPrefix, r_list, g_list, b_list);

% load images
% find the nonzero channels length
exist_r = size(r_list, 1); 
exist_g = size(g_list, 1); 
exist_b = size(b_list, 1);

list_length = max([exist_r exist_g exist_b]);

% read all the images

for i=1:list_length
    
    if exist_r
        
        rch(:, :, i) = imread([prefix r_list(i).name]);
        
    else
        rch = 0;        
    end;

    if exist_g
        
        gch(:, :, i) = imread([prefix g_list(i).name]);
        
    else
        
        gch=0;
        
    end;

    if exist_b
        
        bch(:, :, i) = imread([prefix b_list(i).name]);
        
    else
        bch=0;        
    end;
end;

% normalize the images
% min max values

if exist_r               
    rch = double(rch);    
    rmax = max(max(max(rch)));
    rmin = min(min(min(rch)));
    rch = uint8(((rch - rmin)/(rmax - rmin)) * 255);
end;

if exist_g               
    gch = double(gch);    
    gmax = max(max(max(gch)));
    gmin = min(min(min(gch)));
    gch = uint8(((gch - gmin)/(gmax - gmin)) * 255);    
end;

if exist_b
    bch = double(bch);    
    bmax = max(max(max(bch)));
    bmin = min(min(min(bch)));
    bch = uint8(((bch - bmin)/(bmax - bmin)) * 255);    
end;

% prepare montage images

xdim = max([size(rch,1) size(gch,1) size(bch,1)]);
ydim = max([size(rch,2) size(gch,2) size(bch,2)]);

if exist_r
    rch = reshape(rch, xdim, ydim, 1, list_length);    
    [h, r_mon] = montage(rch);    
    r_mon = imresize(r_mon, [350 NaN]);
    [xdim_mon, ydim_mon]  = size(r_mon);
end;

if exist_g
    gch = reshape(gch, xdim, ydim, 1, list_length);    
    [h, g_mon] = montage(gch);    
    g_mon = imresize(g_mon, [350 NaN]);    
    [xdim_mon, ydim_mon]  = size(g_mon);    
end;

if exist_b
    bch = reshape(bch, xdim, ydim, 1, list_length);    
    [h, b_mon] = montage(bch);    
    b_mon = imresize(b_mon, [350 NaN]);    
    [xdim_mon, ydim_mon]  = size(b_mon);    
end;

if xdim_mon >1
    fullmontage = uint8(zeros(xdim_mon, ydim_mon, 3));
    if exist_r
        fullmontage(:,:,1) = r_mon*2;    
    end;
    if exist_g
        fullmontage(:,:,2) = g_mon*2;    
    end;
    if exist_b
        fullmontage(:,:,3) = b_mon*2;    
    end;
    % save it
    imwrite(fullmontage, [outPrefix 'montage.jpg']);
end;
    
% make the small image
rand_num = floor(rand * list_length) + 1;
vimg = uint8(zeros(xdim, ydim, 3));
if exist_r
    vimg(:,:,1) = rch(:,:,1,rand_num)*2;
end;
if exist_g
    vimg(:,:,2) = gch(:,:,1,rand_num)*2;
end;
if exist_b
    vimg(:,:,3) = bch(:,:,1,rand_num)*2;
end;
vimg = imresize(vimg, [100 NaN]);   

imwrite(vimg, [outPrefix '.jpg']);

% make the view images
for i=1:list_length

    % create image
    outimg = uint8(zeros(xdim, ydim, 3));
    if exist_r
        outimg(:,:,1) = rch(:,:,1,i)*2;
    end;
    if exist_g
        outimg(:,:,2) = gch(:,:,1,i)*2;
    end;
    if exist_b
        outimg(:,:,3) = bch(:,:,1,i)*2;
    end;
    outimg = imresize(outimg, [300 NaN]);   

    name = sprintf('%s_0%d.jpg', outPrefix, i);
    
    imwrite(outimg, name);

    [outPrefix '.html']
    
    add_html_image([outPrefix '.html'], name , '');
    
end;


% save it
% imwrite(fullmontage, [outPrefix 'montage.jpg']);



% for i=1:list_length
%     
%     out = zeros(xdim, ydim, 3);
%     
%     if exist_r
%         
%         
%         
%     end;
%     
% end;
    
    
    
% % rch = 0; gch=0; bch=0;
% % 
% % if length(r) > 0
% % 
% %     if (exist(r))
% %     
% %         rch = imread(r);
% %         rch = imadjust(rch,stretchlim(rch),[]);
% %     else
% %         
% %         r=zeros(20, 20);
% %         
% %     end;
% % 
% % end;
% % 
% % [rx, ry] = size(rch);
% % 
% % if length(g) > 0
% % 
% %     if (exist(g))
% %     
% %         gch = imread(g);
% %         gch = imadjust(gch,stretchlim(gch),[]);
% %         
% %     else
% %         
% %         g=zeros(20, 20);
% %         
% %     end;
% %         
% %     
% % end;
% % 
% % [gx, gy] = size(gch);
% % 
% % if length(g) > 0
% % 
% %     if exist(b)
% %     
% %         bch = imread(b);
% %         bch = imadjust(bch,stretchlim(bch),[]);
% %         
% %     else
% %         
% %         b = zeros(20, 20);
% %         
% %     end;
% %     
% % end;
% % 
% % [bx, by] = size(bch);
% % 
% % maxx = max([rx gx bx]);
% % maxy = max([ry gy by]);
% % 
% % out = zeros(maxx, maxy, 3);
% % 
% % 
% % out(1:rx, 1:ry, 1) = rch;
% % out(1:gx, 1:gy, 2) = gch;
% % out(1:bx, 1:by, 3) = bch;
% % 
% % out = imresize(out, scale);
% % 
% % minint = min(min(min(out)));
% % maxint = max(max(max(out)));
% % 
% % out = uint8(((out - minint)/(maxint - minint)) * 255);
