function failedJobs = evalBatch(inputFolderName, outputFolderName, pipeLineName, poolSize, batchSize, numAttempts)

disp('- Analysis: started -');

warning('off', 'all');

prepareBatchAnalysis(pipeLineName, inputFolderName, outputFolderName);

% delete(gcp('nocreate'));
% parpool(poolSize);

t = load([outputFolderName 'CPBatchInfo.mat']);

loadImagesIndex = find(~cellfun('isempty', strfind(t.handles.Settings.ModuleNames, 'LoadImages')));

firstImageName = t.handles.Settings.VariableValues{loadImagesIndex,3};

first = 1;

% needed to set in the code for parfor syntax
last = 1;
eval( ['last = length(t.handles.Pipeline.FileList' firstImageName ');'] );

% clear image names for the first batch based on the pipeline settings
for i=3:2:9
    if ~strcmp(t.handles.Settings.VariableValues{loadImagesIndex,i},'/')
        t.handles.Pipeline = rmfield(t.handles.Pipeline, t.handles.Settings.VariableValues{loadImagesIndex,i});
    end
end

loopEnd = ceil((last-first+1)/batchSize);

succFinishedJobs = zeros(loopEnd, 1);

for i=1:loopEnd
        
    currentI = (i-1) * batchSize+first;
    
    status = 0;
    
    attempCounter = 0;
    
    while status == 0 && attempCounter < numAttempts
    
%         try
            if currentI +batchSize-1 > last
                fprintf('%d->%d (attemp %d)\n', currentI, last, attempCounter);
                status = feval(@evalCP, currentI, last, t.handles);
                
            else
                fprintf('%d->%d (attemp %d)\n', currentI, currentI+batchSize-1, attempCounter);
                status = feval(@evalCP, currentI, currentI+batchSize-1, t.handles);
                
            end
            if status == 1
                succFinishedJobs(i) = 1;
            end
%         catch errID
%             disp(errID);
%         end
    
        attempCounter = attempCounter + 1;
        
    end
            
end

% save(fullfile(folderName, 'results.mat'), 't');

% checking for non-completed jobs and report them

failedJobs = succFinishedJobs == 0;

if ~isempty(failedJobs)
    disp('There are jobs which could not be done after multiple attempts!');
else
    disp('Every job successfully finished!');
end

% delete(gcp('nocreate'));

warning('on', 'all');

disp('- Analysis: finished -');