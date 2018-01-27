%% this file contains the settings for the HCS workflow

fprintf('Loading options... ');

global options;

options.data.rawDataDir = 'D:\Projects\NEUBIAS\data\artificial\Week1_22123_unfocus\';
options.data.rawDataChannels = {'*w1*.tif','*w2*.tif','*w4*.tif'};
options.data.refactoredDataDir = 'd:\Projects\NEUBIAS\data\refactored\NEUBIAS_example_data_1_siRNA_screen\';
options.data.correctedDataDir = 'd:\Projects\NEUBIAS\data\corrected\Week1_22123_unfocus\';
options.data.focusedDataDir = 'd:\Projects\NEUBIAS\data\focused\Week1_22123_unfocus\';
options.data.segmentedDataDir = [];
options.data.featuresDir = [];
options.data.statisticsDir = [];
options.data.allowedFileTypes = {'.bmp', '.gif', '.jpg', '.jpeg', '.tif', '.tiff', '.png',...
                                '.BMP', '.GIF', '.JPG', '.JPEG', '.TIF', '.TIFF', '.PNG'};

options.refactor.doRefactor = 0;
options.refactor.inputFormat = [];
options.refactor.outputFormat = '%s_w%c%02d_s%02d_p%02'; % SSS Structure
                     
options.ffc.name = 'cidre';
options.ffc.vizualize = 1;

options.focusDetection = '';

options.segmentation.name = 'cp_identify_objects';
options.segmentation.pipe = fullfile(options.HCSWorkflowPath,'segmentation','pipelines','first_pipe.cppipe');
options.segmentation.commands = {sprintf('python CellProfiler.py -c -r -i %s -o %s -p %s',...
                        options.data.focusedDataDir, options.data.segmentedDataDir, options.segmentation.pipe)};
options.features = {};

fprintf('DONE.\n');
