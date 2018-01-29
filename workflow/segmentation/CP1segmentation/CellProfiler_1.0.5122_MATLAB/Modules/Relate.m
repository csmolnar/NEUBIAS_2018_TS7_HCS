function handles = Relate(handles)

% Help for the Relate module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Assigns relationships: All objects (e.g. speckles) within a parent object
% (e.g. nucleus) become its children.
% *************************************************************************
%
% Allows associating "children" objects with "parent" objects. This is
% useful for counting the number of children associated with each parent,
% and for calculating mean measurement values for all children that are
% associated with each parent. For every measurement that has been made of
% the children objects upstream in the pipeline, this module calculates the
% mean value of that measurement over all children and stores it as a
% measurement for the parent. For this reason, this module should be placed
% *after* all Measure modules that make measurements of the children
% objects.
%
% An object will be considered a child even if the edge is the only part
% touching a parent object. If an object is touching two parent objects,
% the objects parent will be the higher numbered parent.
%
% The minimum distances of each child to its parent are also calculated.
% These values are associated with the child objects. If an "Other" object
% is defined (e.g. Nuclei), then distances are calculated to this object
% too, as well as normalized distances.  Normalized distances for each
% child have a range [0 1] and are calculated as: 
% (distance to the Parent) / sum(distances to parent and Other object)

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
% $Revision: 5052 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What objects are the children objects (subobjects)?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
SubObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What are the parent objects?
%infotypeVAR02 = objectgroup
%inputtypeVAR02 = popupmenu
ParentName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What other object do you want to find distances to? (Must be one object per parent object, e.g. Nuclei)
%infotypeVAR03 = objectgroup
%choiceVAR03 = None
%inputtypeVAR03 = popupmenu
ParentName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%VariableRevisionNumber = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
SubObjectLabelMatrix = CPretrieveimage(handles,['Segmented', SubObjectName],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
ParentObjectLabelMatrix = CPretrieveimage(handles,['Segmented', ParentName{1}],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
if ~strcmp(ParentName{2},'None')
    
    %% Sanity checks
    if strcmp(SubObjectName,ParentName{1}) || strcmp(SubObjectName,ParentName{2})
        CPwarndlg('The Children and at least one of the Parent objects are the same.  Your results may be erroneous.','Relate module')
    end
    if strcmp(ParentName{1},ParentName{2})
        CPwarndlg('The Parent and Other Object are the same.  Your results may be erroneous.','Relate module')
    end

    StepParentObjectLabelMatrix = CPretrieveimage(handles,['Segmented', ParentName{2}],ModuleName,'MustBeGray','DontCheckScale');
   
    %% Sanity check
    if max(ParentObjectLabelMatrix(:)) ~= max(StepParentObjectLabelMatrix(:))
        CPwarndlg(['The number of Parent Objects (' num2str(max(ParentObjectLabelMatrix(:))) ...
            ') does not equal the number of Other objects (' num2str(max(StepParentObjectLabelMatrix(:))) ...
            ') in the Relate Module, Cycle#' num2str(handles.Current.SetBeingAnalyzed) ...
            '.  If the difference is large, this may cause the Relate Module output to be suspect.'],'Relate module')
    end
else
    ParentName = {ParentName{1}};
end

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

try
    [handles,NumberOfChildren,ParentsOfChildren] = CPrelateobjects(handles,SubObjectName,ParentName{1},...
        SubObjectLabelMatrix,ParentObjectLabelMatrix,ModuleName);
    handles.Measurements.(SubObjectName).SubObjectFlag=1;
    
    %% Save Distance 'Features'
    handles.Measurements.(SubObjectName).DistanceFeatures = ParentName;

    %% Initialize Distance
    handles.Measurements.(SubObjectName).Distance{handles.Current.SetBeingAnalyzed} = ...
        zeros(max([length(ParentsOfChildren) 1]),length(ParentName));
    
    %% Calcuate the smallest distance from each Child to their Parent
    %% If no parent exists, then Distance = NaN
        
    if isfield(handles.Measurements.(SubObjectName),'Location')
        iObj = 0;
        for thisParent = ParentName %% Will need to change if we add more StepParents
            iObj = iObj + 1;
            
            %% Calculate perimeters for all parents simultaneously
            DistTransAll = CPlabelperim(handles.Pipeline.(['Segmented' ParentName{iObj}]));
            
            for iParentsOfChildren = 1:max(ParentsOfChildren)
                %% Calculate distance transform to perimeter of Parent objects
                DistTrans = (bwdist(DistTransAll == iParentsOfChildren));
                
                %% Get location of each child object
                ChList = find(ParentsOfChildren == iParentsOfChildren);
                ChildrenLocations = handles.Measurements.(SubObjectName).Location{handles.Current.SetBeingAnalyzed}(ChList,:);

                roundedChLoc = round(ChildrenLocations);
                idx = sub2ind(size(DistTrans),roundedChLoc(:,2), roundedChLoc(:,1));
                Dist = DistTrans(idx);

                %% SAVE Distance to 'handles'
                handles.Measurements.(SubObjectName).Distance{handles.Current.SetBeingAnalyzed}(ChList,iObj) = Dist;
            end
        end
    else
        warning('There is no ''Location'' field with which to find subObj to Parent distances')
    end
    
    %% Calculate normalized distances
    %% All distances are relative to the *first* parent.
    if length(ParentName) > 1
        Dist = handles.Measurements.(SubObjectName).Distance{handles.Current.SetBeingAnalyzed};   
        %% Initialize NormDistance
        handles.Measurements.(SubObjectName).NormDistance{handles.Current.SetBeingAnalyzed} = zeros(size(Dist,1),1);
        NormDist = Dist(:,1) ./ sum(Dist,2);
        NormDist(isnan(NormDist)) = 0;  %% In case sum(Dist,2) == 0 for any reason (no parents/child, or child touching either parent
        
        %% Save Normalized Distances
        handles.Measurements.(SubObjectName).NormDistanceFeatures = {ParentName{1}}; %% outer curly brackets needed for correct length(MeasurementFeatures) calculation in inner loop below
        handles.Measurements.(SubObjectName).NormDistance{handles.Current.SetBeingAnalyzed} = NormDist;
    end

    %% Adds a 'Mean<SubObjectName>' field to the handles.Measurements structure
    %% which finds the mean measurements of all the subObjects that relate to each parent object
    MeasurementFieldnames = fieldnames(handles.Measurements.(SubObjectName))';
    NewObjectName=['Mean',SubObjectName];
    if isfield(handles.Measurements.(SubObjectName),'Parent')

        % Why is test line here? Isn't this always the case?  Or is it in case Relate is called twice?- Ray 2007-08-09
        if length(handles.Measurements.(SubObjectName).Parent) >= handles.Current.SetBeingAnalyzed
            Parents=handles.Measurements.(SubObjectName).Parent{handles.Current.SetBeingAnalyzed};
            for RemainingMeasurementFieldnames = MeasurementFieldnames
                if isempty(strfind(char(RemainingMeasurementFieldnames),'Features')) || ~isempty(strfind(char(RemainingMeasurementFieldnames),'Parent')) || ~isempty(strfind(char(RemainingMeasurementFieldnames),'Children'))
                    continue
                else
                    FieldName=char(RemainingMeasurementFieldnames);
                    MeasurementFeatures=handles.Measurements.(SubObjectName).(FieldName);
                    handles.Measurements.(NewObjectName).(FieldName)=MeasurementFeatures;
                    if length(handles.Measurements.(SubObjectName).(FieldName(1:end-8))) >= handles.Current.SetBeingAnalyzed;
                        Measurements=handles.Measurements.(SubObjectName).(FieldName(1:end-8)){handles.Current.SetBeingAnalyzed};
                        % The loop over 'j' below will never be entered if
                        % there are no parents in the image, leading to a bug
                        % where the data is truncated if the last few images
                        % don't contain parents.  This next statement handles
                        % that case by ensuring that at least something is
                        % written at the correction location in the
                        % Measurements structure.
                        handles.Measurements.(NewObjectName).(FieldName(1:end-8)){handles.Current.SetBeingAnalyzed} = [];
                        for i=1:length(MeasurementFeatures)
                            for j=1:max(max(ParentObjectLabelMatrix))
                                index=find(Parents==j);
                                if isempty(index)
                                    handles.Measurements.(NewObjectName).(FieldName(1:end-8)){handles.Current.SetBeingAnalyzed}(j,i)=0;
                                else
                                    handles.Measurements.(NewObjectName).(FieldName(1:end-8)){handles.Current.SetBeingAnalyzed}(j,i)=mean(Measurements(index,i));
                                end
                            end
                        end
                    else
                        CPwarndlg('The Relate module is attempting to take the mean of a measurement downstream.  Be advised that unless the Relate module is placed *after* all Measurement modules, some ''Mean'' measurements will not be calculated.','Relate Module warning','replace')
                    end
                end
            end
        end
    end
catch
    error('The Relate Module errored.  This may be fixed by ensuring that all Relate Modules occur *after* all Measure Modules, and that two Relate Modules are not relating the same objects.')
end

%%% Since the label matrix starts at zero, we must include this value in
%%% the list to produce a label matrix image with children re-labeled to
%%% their parents values. This does not get saved and is only for display.
if ~isempty(ParentsOfChildren)
    ParentsOfChildrenLM = [0;ParentsOfChildren];
    NewObjectParentLabelMatrix = ParentsOfChildrenLM(SubObjectLabelMatrix+1);
else
    NewObjectParentLabelMatrix = SubObjectLabelMatrix;
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(ParentObjectLabelMatrix,'TwoByTwo',ThisModuleFigureNumber);
    end
    subplot(2,2,1);
    ColoredParentLabelMatrixImage = CPlabel2rgb(handles,ParentObjectLabelMatrix);
    CPimagesc(ColoredParentLabelMatrixImage,handles);
    title(['Parent Objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,2,2);
    ColoredSubObjectLabelMatrixImage = CPlabel2rgb(handles,SubObjectLabelMatrix);
    CPimagesc(ColoredSubObjectLabelMatrixImage,handles);
    title('Original Sub Objects');
    subplot(2,2,3);
    ColoredNewObjectParentLabelMatrix = CPlabel2rgb(handles,NewObjectParentLabelMatrix);
    CPimagesc(ColoredNewObjectParentLabelMatrix,handles);
    title('New Sub Objects');
end