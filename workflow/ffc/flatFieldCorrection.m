function flatFieldCorrection()
% flat field correction

global options;

fprintf('\n====================\nFlat field correction started...\n');

if ~isfield(options,'ffc') || isempty(options.ffc)
    options.ffc.name = 'cidre';
end

if strcmp(options.ffc.name, 'cidre')
    % run CIDRE method
    for i=1:length(options.data.rawDataChannels)    
        dataChannelRegExp =  options.data.rawDataChannels{i};
        [~,baseName,ext] = fileparts(dataChannelRegExp);
        if ~isempty(ext)
            if any(cellfun(@(x) strcmp(x,ext), options.data.allowedFileTypes)) || strcmp(ext,'*')
                nonStars = strfind(baseName, '*');
                baseName(nonStars) = [];
            else
                
            end
        else
            
        end
        
        cidre(fullfile(options.data.rawDataDir, options.data.rawDataChannels{i}), 'destination',options.data.correctedDataDir);
        
        movefile(fullfile(options.data.correctedDataDir,'cidre_model.mat'),...
                 fullfile(options.data.correctedDataDir,['cidre_model_' baseName '.mat']));
        if options.ffc.vizualize
            load(fullfile(options.data.correctedDataDir,['cidre_model_' baseName '.mat']),'model');
            figure;
            subplot(1,2,1,'replace'); imagesc(model.v); colorbar; title('Gain (v)');
            subplot(1,2,2,'replace'); imagesc(model.z); colorbar; title('Additive noise (z)');
            clear model;
        end
    end
    
    % TODO option for existing correction model
    
elseif strcmp(options.ffc.name, 'mean')
    % run mean image computation based ffc
end

fprintf('\nFlat field correction DONE.\n');

end

% function isElementOf( str, strCellArray )
%     
% end