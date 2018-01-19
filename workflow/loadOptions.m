%% this file contains the settings for the HCS workflow

fprintf('Loading options... ');

global options;

options = struct();

options.data.rawDataDir = 'd:\Projects\NEUBIAS\data\artificial\NEUBIAS_example_data_1_siRNA_screen_illuminated\';
options.data.rawDataChannels = {'*DAPI.tif','*GFP.tif'};
options.data.correctedDataDir = 'd:\Projects\NEUBIAS\data\corrected\NEUBIAS_example_data_1_siRNA_screen_illuminated\';
options.data.focusedDataDir = [];
options.data.segmentedDataDir = [];
options.data.featuresDir = [];
options.data.statisticsDir = [];
options.data.allowedFileTypes = {'.bmp', '.gif', '.jpg', '.jpeg', '.tif', '.tiff', '.png',...
                                '.BMP', '.GIF', '.JPG', '.JPEG', '.TIF', '.TIFF', '.PNG'};

options.ffc = [];
options.focusDetection = '';
options.segmentation = [];
options.features = {};

fprintf('DONE.\n');
