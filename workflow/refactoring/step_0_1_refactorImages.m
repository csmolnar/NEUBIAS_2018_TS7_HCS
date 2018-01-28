function step_0_1_refactorImages()

% rename all images to the nomenclature of ACC (SSS structure)
% SSS structure: PlateName_w%Row%Column_s%Num_*.ext (e.g. Plate001_wB01_s03.tif)
%
% Example format:
% Perkin Elmer Operetta: r01c01f01p01-ch1sk1fk1fl1.tiff ->
%           'r%02dc%02df%02dp%02d-ch%dsk%dfk%dfl%d.tiff'

global options;

if options.refactor.doRefactor
    for i=1:length(options.data.rawDataChannels)
        [rawdataDir, ~, ~] = options.data.rawDataDir;
        pathFolders = strsplit(rawdataDir, filesep);
        plateName = pathFolders{end};
        imageList = dir(fullfile(options.data.rawDataDir, options.data.rawDataChannels{i}));
        refactorFileId = fopen(fullfile(options.data.rawDataDir),...
                    sprintf('%s_%s.txt',plateName,strrep(options.data.rawDataChannels{i},'*','')), 'w');
        if ~isempty(imageList(1).name, plateName)
            newNamePrefix = '';
        else
            newNamePrefix = plateName;
        end
        for j = 1:length(imageList)
            oldName = imageList(j).name;
%             [~, oldBaseName, ext] = fileparts(oldName);
            newName = sprintf('%s_%s', newNamePrefix, oldName);
            fprintf(refactorFileId, '%s,%s\n', oldName, newName);
            copyfile(fullfile(options.data.rawDataDir,oldName),fullfile(options.data.refactoredDataDir,newName));
        end
        fclose(refactorFileId);
    end
    options.data.rawDataDir;
end
    
end