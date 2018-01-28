function step_1_flatFieldCorrection()
% flat field correction

global options;

fprintf('\n====================\nFlat field correction started...\n');

if ~isfield(options,'ffc') || isempty(options.ffc)
    options.ffc.name = 'cidre';
end

plateList = dir( fullfile(options.data.rawDataDir) );
% clear not plate folders
plateList(~cat(1,plateList.isdir)) = [];

for i=length(plateList):-1:1
    if strcmp(plateList(i).name,'.') || strcmp(plateList(i).name,'..')
        plateList(i) = [];
    end
end

if strcmp(options.ffc.name, 'cidre')
    % run CIDRE method
    
    for ii=1:length(plateList)
        for i=1:length(options.data.rawDataChannels)
            dataChannelRegExp =  options.data.rawDataChannels{i};
            [~,baseName,ext] = fileparts(dataChannelRegExp);
            if ~isempty(ext)
                if any(cellfun(@(x) strcmp(x,ext), options.data.allowedFileTypes)) || strcmp(ext,'*')
                    baseName(baseName=='*') = '';
                else
                    
                end
            else
                baseName(baseName=='*') = '';
            end
            
            cidre(fullfile(options.data.rawDataDir,plateList(ii).name, options.data.rawDataChannels{i}), 'destination',fullfile(options.data.correctedDataDir, plateList(ii).name));
            
            movefile(fullfile(options.data.correctedDataDir,plateList(ii).name,'cidre_model.mat'),...
                fullfile(options.data.correctedDataDir,plateList(ii).name,['cidre_model_' baseName '.mat']));
            if options.ffc.vizualize
                load(fullfile(options.data.correctedDataDir,plateList(ii).name,['cidre_model_' baseName '.mat']),'model');
                figure;
                subplot(1,2,1,'replace'); imagesc(model.v); colorbar; title(sprintf('Gain (v) (%s)',options.data.rawDataChannels{i}));
                subplot(1,2,2,'replace'); imagesc(model.z); colorbar; title(sprintf('Additive noise (z) (%s)',options.data.rawDataChannels{i}));
                clear model;
            end
        end
        
    end
    % TODO option for existing correction model
    
elseif strcmp(options.ffc.name, 'mean')
    % run mean image computation based ffc
    error('HCS:error','Method ''mean'' is not implemented yet.' );
else
    error('HCS:error','Method %s is not implemented.' , options.ffc.name);
end

fprintf('\nFlat field correction DONE.\n');

end

% function isElementOf( str, strCellArray )
%
% end