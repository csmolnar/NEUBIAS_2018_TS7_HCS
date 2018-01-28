function succFinished = evalCP(StartImage, EndImage, handles)

BatchFilePrefix = 'Batch_';

tic;

handles.Current.BatchInfo.Start = StartImage;
handles.Current.BatchInfo.End = EndImage;
for BatchSetBeingAnalyzed = StartImage:EndImage,
    handles.Current.SetBeingAnalyzed = BatchSetBeingAnalyzed;
    for SlotNumber = 1:handles.Current.NumberOfModules,        
        ModuleNumberAsString = sprintf('%02d', SlotNumber);
        ModuleName = char(handles.Settings.ModuleNames(SlotNumber));
        handles.Current.CurrentModuleNumber = ModuleNumberAsString;
        try
            handles = feval(ModuleName,handles);
        catch e
            handles.BatchError = [ModuleName ' ' e];
            disp(['Batch Error: ' ModuleName ' ' e]);
            rethrow(e);
            succFinished = 0;
            quit;
        end
        disp(['  -> Running module ' ModuleNumberAsString]);
    end
    fprintf('Set %d done. Elapsed time %f', BatchSetBeingAnalyzed, toc);
end

succFinished = 1;
%cd(folderName);
%handles.Pipeline = []; eval(['save ',sprintf('%s%d_to_%d_OUT', BatchFilePrefix, StartImage, EndImage), ' handles;']);
