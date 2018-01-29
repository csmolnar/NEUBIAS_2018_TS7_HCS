function handles = ExportToACC(handles)

% DATE: 
% May 03, 2016
%
% AUTHOR: 
% Peter Horvath
%
% Save data into ACC format.
% Category: File Processing
%
% DESCRIPTION:
% Save the measurements obtained with CellProfiler into text files that 
% can be used by Advanced Cell Classifier (ACC).
%
% WEBSITE: 
% http://www.cellclassifier.org/
%
% COPYRIGHT
% Advanced Cell Classifier (ACC) Toolbox. All rights reserved.
% Copyright © 2016 Peter Horvath
% Synthetic and System Biology Unit, Hungarian Academia of Sciences,
% Biological Research Center, Szeged, Hungary; Institute for Molecular
% Medicine Finland, University of Helsinki, Helsinki, Finland.
%
% *************************************************************************

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 5025 $


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);


%textVAR01 = Where do you want to save the measurements
%defaultVAR01 = anal2
OriginalFolder = char(handles.Settings.VariableValues{CurrentModuleNum,1});


%textVAR02 = Which objects do you want to export? Cell location will be calculted from the first object. Normally use here nuclei.
%infotypeVAR02 = objectgroup
%choiceVAR02 = /
%choiceVAR02 = Image
%choiceVAR02 = Experiment
Object{1} = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 =
%infotypeVAR03 = objectgroup
%choiceVAR03 = /
%choiceVAR03 = Image
%choiceVAR03 = Experiment
Object{2} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 =
%infotypeVAR04 = objectgroup
%choiceVAR04 = /
%choiceVAR04 = Image
%choiceVAR04 = Experiment
Object{3} = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%%%VariableRevisionNumber = 2

a = tic;

drawnow
tmp = {};
for n = 1:3
    if ~strcmp(Object{n}, '/')
        tmp{end+1} = Object{n};
    end
end
Object = tmp;

% last cycle
%if handles.Current.SetBeingAnalyzed == handles.Current.NumberOfImageSets
% in the new version we do it for every cycle
if 1
    % operating system specific slash
    if ispc == 1
        OSSlash = '\';
    else
        OSSlash = '/';
    end;

    % output folder
    outPutFolder = handles.Current.DefaultOutputDirectory;
    if outPutFolder(length(outPutFolder)) ~= OSSlash
        outPutFolder = [outPutFolder OSSlash];
    end;
    
    
    ObjectNames = unique(Object);
        
    % for each image sets
    % Note: position of the objects comes from the first measured object
    % type (usually nuclei)
    
    %for i=1:handles.Current.NumberOfImageSets
    % current cycle
    i = handles.Current.SetBeingAnalyzed;
        wdata = [];
        [~, baseOutputFileName, ~] = fileparts(char(handles.Measurements.Image.FileNames{i}(1)));
        outputFileName = [baseOutputFileName '.txt'];
        
        % get xy positions
        PosField = getfield(handles.Measurements, Object{1});
        wdata(:,1:2) = PosField.Location{i}(:,1:2); % nucleii location       
        header = {};
        header{1} = 'LocationX';
        header{2} = 'LocationY';        
        featureCounter = 2;
        % get measured fileds:        
        for j=1:length(ObjectNames)
            % get measurements
            measurementField = getfield(handles.Measurements, ObjectNames{j});
            measurementFieldNames = fieldnames(measurementField);
            measurementFieldNamesDescrPos = strfind(measurementFieldNames, 'Features');
            for k=1:length(measurementFieldNames)                
                if ~isempty(measurementFieldNamesDescrPos{k})
                    % avoid location
                    isLocationField = strfind(measurementFieldNames{k}, 'Location');  
                    isChildrenField = strfind(measurementFieldNames{k}, 'Children');  
                    isParentField = strfind(measurementFieldNames{k}, 'Parent');  
                    if isempty(isLocationField) & isempty(isChildrenField) & isempty(isParentField)
                        % put feature names                                                
                        featureNameList = getfield(measurementField, measurementFieldNames{k});
                        % put data
                        for l=1:length(featureNameList)
                            header = [header [ObjectNames{j} '.' measurementFieldNames{k} '.' featureNameList{l}]];                        
                        end;
                        featureCounter = featureCounter + length(featureNameList);
                        featureValues = getfield(measurementField, measurementFieldNames{k+1});
                        currentData = featureValues{i};
                        wdata = [wdata (currentData)];
                    end;
                end;
            end;                                                
        end;
        % save output
        save([outPutFolder OriginalFolder OSSlash outputFileName], 'wdata', '-ascii');
        % save feature description for further analisys
        if ~exist([outPutFolder OriginalFolder OSSlash 'featureNames.acc'],'file') %TODO wrong if we change the features to save!!!
            outFile = fopen([outPutFolder OriginalFolder OSSlash 'featureNames.acc'], 'w');
            for l=1:length(header);
                fprintf(outFile, '%s\n', header{l});
            end;
            fclose(outFile);
        end;
    %end;    
end;

toc(a);