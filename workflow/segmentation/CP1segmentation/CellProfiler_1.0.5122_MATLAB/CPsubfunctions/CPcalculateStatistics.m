function handles = CPcalculateStatistics(handles,DataName,LogOrLinear,FigureName,ModuleName,LicenseStats)

FigureIncrement = 1;

%%% Get all fieldnames in Measurements
ObjectFields = fieldnames(handles.Measurements);
GroupingStrings = handles.Measurements.Image.(DataName);
%%% Need column vector
GroupingValues = str2num(char(GroupingStrings')); %#ok Ignore MLint

%%% Get the handle to the waitbar and update the text in the waitbar
waitbarhandle = CPwaitbar(0,'');
CPwaitbar(0,waitbarhandle,'CPcalculateStatistics Progress');

for i = 1:length(ObjectFields)
    ObjectName = char(ObjectFields(i));
    %%% Filter out Experiment and Image fields
    if ~strcmpi(ObjectName,'Experiment')

        try
            %%% Get all fieldnames in Measurements.(ObjectName)
            MeasureFields = fieldnames(handles.Measurements.(ObjectName));
        catch %%% Must have been text field and ObjectName is class 'cell'
            continue
        end

        for j = 1:length(MeasureFields)
            MeasureFeatureName = char(MeasureFields(j));
            if length(MeasureFeatureName) > 7
                if strcmpi(MeasureFeatureName(end-7:end),'Features')

                    %%% Not placed with above if statement since
                    %%% MeasureFeatureName may not be 8 characters long
                    if ~strcmpi(MeasureFeatureName(1:8),'Location')

                        if strcmpi(MeasureFeatureName,'ModuleErrorFeatures')
                            continue;
                        end

                        %%% Get Features
                        MeasureFeatures = handles.Measurements.(ObjectName).(MeasureFeatureName);

                        %%% Get Measure name
                        MeasureName = MeasureFeatureName(1:end-8);
                        %%% Check for measurements
                        if ~isfield(handles.Measurements.(ObjectName),MeasureName)
                            CPwarndlg(['There is a problem in the ' ModuleName ' module becaue it could not find the measurements you specified. CellProfiler will proceed but this module will be skipped.']);
                            close(waitbarhandle)
                            return;
                        end

                        Ymatrix = zeros(handles.Current.NumberOfImageSets,length(MeasureFeatures));
                        for k = 1:handles.Current.NumberOfImageSets
                            for l = 1:length(MeasureFeatures)
                                if isempty(handles.Measurements.(ObjectName).(MeasureName){k})
                                    Ymatrix(k,l) = 0;
                                else
                                    Ymatrix(k,l) = mean(handles.Measurements.(ObjectName).(MeasureName){k}(:,l));
                                end
                            end
                        end

                        GroupingValueRows = size(GroupingValues,1);
                        YmatrixRows = size(Ymatrix,1);
                        if GroupingValueRows ~= YmatrixRows
                            CPwarndlg('There was an error in the Calculate Statistics module involving the number of text elements loaded for it.  CellProfiler will proceed but this module will be skipped.');
                            close(waitbarhandle)
                            return;
                        else
                            [v, z, OrderedUniqueDoses, OrderedAverageValues] = CP_VZfactors(GroupingValues,Ymatrix);
                            if LicenseStats == 1
                                if ~strcmpi(FigureName,'Do not save')
                                    PartialFigureName = fullfile(handles.Current.DefaultOutputDirectory,FigureName);
                                else PartialFigureName = FigureName;
                                end
                                try
                                    [FigureIncrement, ec50stats] = CPec50(OrderedUniqueDoses',OrderedAverageValues,LogOrLinear,PartialFigureName,ModuleName,DataName,FigureIncrement);
                                catch
                                    ec50stats = zeros(size(OrderedAverageValues,2),4);
                                end
                                ec = ec50stats(:,3);
                                if strcmpi(LogOrLinear,'Yes')
                                    ec = exp(ec);
                                end
                            end
                        end

                        measurefield = [ObjectName,'Statistics'];
                        featuresfield = [ObjectName,'StatisticsFeatures'];
                        if isfield(handles.Measurements,'Experiment')
                            if isfield(handles.Measurements.Experiment,measurefield)
                                OldEnd = length(handles.Measurements.Experiment.(featuresfield));
                            else OldEnd = 0;
                            end
                        else OldEnd = 0;
                        end
                        for a = 1:length(z)
                            handles.Measurements.Experiment.(measurefield){1}(1,OldEnd+a) = z(a);
                            handles.Measurements.Experiment.(measurefield){1}(1,OldEnd+length(z)+a) = v(a);
                            handles.Measurements.Experiment.(featuresfield){OldEnd+a} = ['Zfactor_',MeasureName,'_',MeasureFeatures{a}];
                            handles.Measurements.Experiment.(featuresfield){OldEnd+length(z)+a} = ['Vfactor_',MeasureName,'_',MeasureFeatures{a}];
                            if LicenseStats == 1
                                handles.Measurements.Experiment.(measurefield){1}(1,OldEnd+2*length(z)+a) = ec(a);
                                handles.Measurements.Experiment.(featuresfield){OldEnd+2*length(z)+a} = ['EC50_',MeasureName,'_',MeasureFeatures{a}];
                            end
                        end
                    end
                end
            end
        end
    end
    %%% Update waitbar
    CPwaitbar(i./length(ObjectFields),waitbarhandle);
end
close(waitbarhandle)