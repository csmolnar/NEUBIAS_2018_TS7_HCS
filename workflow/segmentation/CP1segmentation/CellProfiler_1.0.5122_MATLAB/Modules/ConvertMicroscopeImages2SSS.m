function handles = ConvertMicroscopeImages2SSS(handles)

% Help for the ConvertMicroscopeImages2SSS module:
% Category: File Processing
%
% SHORT DESCRIPTION:
% Converts images from screening microscope formats to standard screening
% structure
% *************************************************************************
%
%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by Peter Horvath.
% Copyright 2010.
%
% Please see the AUTHORS file for credits.
%
% Website: http://acc.ethz.ch
%
% $Revision: 5025 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Microscope type:
%choiceVAR01 = BD Pathway
%choiceVAR01 = MD Micro
MicType = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = Specify folder contains images. For BD Pathway define a folder contains Experiment.exp file
%defaultVAR02 = /
OriginalFolder = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Folder where you will save SSS
%defaultVAR03 = /
SSSFolder = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = What do you want to call the plate?
%defaultVAR04 = plate001
PlateName = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = 1st channel name:
%defaultVAR05 = w1
ChannelName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = 2nd channel name:
%defaultVAR06 = w2
ChannelName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = 3rd channel name:
%defaultVAR07 = w3
ChannelName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = 4th channel name:
%defaultVAR08 = w4
ChannelName{4} = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = Number of output folders:
%defaultVAR09 = 3
AnalFolderNumber = char(handles.Settings.VariableValues{CurrentModuleNum,9});


%%%VariableRevisionNumber = 2

drawnow;
% only in the first cycle
if handles.Current.SetBeingAnalyzed == 1

    % operating system specific slash
    if ispc == 1
        OSSlash = '\';
    else
        OSSlash = '/';
    end;

    % check folder names
    % check for OS slash
    if SSSFolder(length(SSSFolder)) ~= OSSlash
        SSSFolder = [SSSFolder OSSlash];
    end;
    if OriginalFolder(length(OriginalFolder)) ~= OSSlash
        OriginalFolder = [OriginalFolder OSSlash];
    end;

    %% create new folder        
    mkdir([SSSFolder PlateName]);
    
    % create analysis folders
    for i=1:str2num(AnalFolderNumber)
        mkdir([SSSFolder PlateName OSSlash 'anal' num2str(i)]);        
    end;
    
    
    % BD Pathway reader
    if strcmp(MicType, 'BD Pathway')
        fid=fopen([OriginalFolder 'Experiment.exp']);

        % Detect number of channels
        % Determine the name of the channels
        % Determine montage and pixel size
        while 1
            tline = fgetl(fid);
            if ~isempty(strfind(tline, 'Dyes='))
                numberOfChannels = str2num(tline(6:end));
            end;

            if ~isempty(strfind(tline, '[Dyes]'))
                for i=1:numberOfChannels
                    tline = fgetl(fid);
                    colorName{i} = tline(3:end);
                end;
            end;
                                    
            if ~isempty(strfind(tline, 'TilesX='))
                TilesX = str2num(tline(8:end));
            end;
            if ~isempty(strfind(tline, 'TilesY='))
                TilesY = str2num(tline(8:end));
            end;
            if ~isempty(strfind(tline, 'TilePixelsX='))
                TilePixelsX = str2num(tline(13:end));
            end;
            if ~isempty(strfind(tline, 'TilePixelsY='))
                TilePixelsY = str2num(tline(13:end));
            end;
            if ~ischar(tline),   break,   end
        end
                
        fclose(fid);

        % list wells
        wellList = dir([OriginalFolder 'Well*']);
        for i=1:length(wellList)
            
            wellName = wellList(i).name(6:end);
            for j=1:numberOfChannels
                % load image  
                counter=0;
                fileList = dir([OriginalFolder wellList(i).name OSSlash colorName{j} '*.tif']);
                
                in = imread([OriginalFolder wellList(i).name OSSlash fileList(1).name]);
                
                for ii=1:TilePixelsX:TilePixelsX*TilesX
                    for jj=1:TilePixelsY:TilePixelsY*TilesY
                        % sitename
                        if counter < 10
                            siteName = ['0' num2str(counter)];
                        else
                            siteName = [num2str(counter)];
                        end;                                       
                        out = in(jj:jj+TilePixelsY-1, ii:ii+TilePixelsX-1);
                        
                        imwrite(out,[SSSFolder PlateName OSSlash PlateName '_' wellName '_' siteName '_' ChannelName{j} '.tif']);
                        counter = counter+1;
                    end;
                end;
                
            end;            
            
        end;    
        
    end;
    
    % MD Micro reader
    if strcmp(MicType, 'MD Micro')

        fileList = dir([OriginalFolder '*.tif']);
        %create new name
        for i=1:length(fileList)
            isThumb = strfind(fileList(i).name, 'thumb');

            if isempty(isThumb)
                % delete the md database string
                newname = fileList(i).name;
                newname = newname(1:length(newname) - 40);
                % change stage position from sX or sXX to XX
                if newname(length(newname) - 4) == 's'
                    newname(length(newname) - 4) = '0';
                elseif newname(length(newname) - 5) == 's'
                    newname = [newname(1:(length(newname)-6)) newname((length(newname)-4) : length(newname))];
                else
                    disp('Unsupported conversion');
                    return;
                end;
                newname = [newname '.tif'];
                if ~strcmp(PlateName, '')
                    %find the name given by md and replace with the input name
                    underscorePos = strfind(newname, '_');
                    newname = [PlateName newname(underscorePos(length(underscorePos) - 2) : length(newname))];
                end;
                copyfile([OriginalFolder fileList(i).name], [SSSFolder PlateName OSSlash newname]);
            end;
        end;
    end;
    
    % change input/output folders
    handles.Current.DefaultOutputDirectory = [SSSFolder PlateName OSSlash];
    handles.Current.DefaultImageDirectory = [SSSFolder PlateName OSSlash];
    
end;

