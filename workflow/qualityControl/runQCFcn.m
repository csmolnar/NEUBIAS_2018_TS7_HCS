% run quality control on corrected images
function runQCFcn(proj,plateList,focus_errdir,blue_errdir,regexpr,report,numChannels,thres,adjust,blueThres,blueSize,thresAdj)

% report=true;    % create .csv file of focus scores

for i=1:numel(plateList)
    plate_name=plateList(i).name;
    if exist([proj filesep focus_errdir filesep plate_name],'dir') ||...
            exist([proj filesep blue_errdir filesep plate_name],'dir')
        % plate already processed
        fprintf('plate %s already processed\n',plate_name);
        continue;
    end
    tic;
    filtall_func_BBBC_parfor_new(proj,plate_name,focus_errdir,regexpr,report,numChannels,thres);
    t1=toc;
    fprintf('plate %s focus done in %f seconds...\n',plate_name,t1);
    hibakeres_func_BBBC_parfor_new(proj,plate_name,blue_errdir,regexpr,numChannels,adjust,blueThres,blueSize,thresAdj);   % function blue is needed!
end
end


function out=filtall_func_BBBC_parfor_new(proj,pdir,errdir,regexpr,report,numChannels,thres)
% proj: project folder, e.g. 'SZBK letöltött cuccok'
% pdir: plates dir in src , e.g. 'Week1_22123'--> always a \Corrected folder
% errdir: folder to collect error images (copies only)
% regexp: regular expression to list files, e.g. '*.tif'
% report: flag to indicate writing a report file of the focus scores
% thres: threshold for std(grad(orig))/std(grad(filt)), default: 1.69
if nargin==6
    thres=1.69;
end

src=[proj filesep];
%     platesDir = [src pdir filesep 'Corrected'];
platesDir = [src pdir];
cd(platesDir);
asterisk=strfind(regexpr,'*');
if isempty(asterisk)
    dir_list=dir(['*' regexpr '*.tif']);
else
    dir_list = dir(regexpr);
end
%     dir_list = dir(regexp);
completeErrdir=fullfile(src,errdir,pdir);
if ~exist(completeErrdir)
    completeErrdir=fullfile(errdir,pdir);
end

mkdir(completeErrdir);
if report
    file=fopen([completeErrdir '_focus.csv'],'w');
end
tic;
n=size(dir_list, 1);
for i=1:n
    img=imread(dir_list(i).name);
    imgd=im2double(img);
    % creating the mask
    f=fspecial('average',[12 12]);
    % filter image to get blurry version
    filt=imfilter(img,f);
    filt=im2double(filt);
    %         figure; imshow([imgd,filt]),title('orig, filt');
    % gradients of orig & blurry
    [g,~]=gradient(filt);
    [gorig,~]=gradient(imgd);
    
    % count standard deviations of grad(orig) & grad(blurry)
    stdor=std(gorig(:));
    stdfilt=std(g(:));
    ratio=stdor/stdfilt;
    
    % if sd-s differ beyond a point, the image was not blurry in the
    % beginning, otherwise it was blurry --> out of focus image
    % detected and can be selected into a separate directory
    if ratio < thres         % ~1.9-2.5 for normal, ~ 1.5 for out of focus
        disp(['found ' dir_list(i).name ' : focus']);
        %% useful part
        movefile(dir_list(i).name,[completeErrdir filesep ...
            dir_list(i).name]);
        otherChannelFiles=dir([dir_list(i).name(1:end-length(regexpr)) '*.tif']);
        if ~isempty(otherChannelFiles)
            if size(otherChannelFiles,1)>numChannels-1
                % didn't retrieve corresponding files, something went
                % wrong --> do nothing with these files
            else
                % move other channel files too
                for j=1:numel(otherChannelFiles)
                    movefile([platesDir filesep otherChannelFiles(j).name], ...
                        [completeErrdir filesep otherChannelFiles(j).name]);
                end
            end
        end
        %% end of useful part
    end
    if report
        toWrite(i)=ratio;
    end
end
if report
    for i=1:n
        fprintf(file,'%s,%f\n',dir_list(i).name,toWrite(i));
    end
    fclose(file);
end
cd(src);
end

function hibakeres_func_BBBC_parfor_new(proj,pdir,errdir,regexpr,numChannels,adjust,thres,blueSize,thresAdj)
% proj: project folder, e.g. 'SZBK letöltött cuccok'
% pdir: plates dir in src , e.g. 'Week1_22123'--> always a \Corrected folder
% errdir: folder to collect error images (copies only)

src=[proj filesep];
platesDir = [src pdir];
%     platesDir = [src pdir filesep 'Corrected'];
cd(platesDir);
asterisk=strfind(regexpr,'*');
if isempty(asterisk)
    dir_list=dir(['*' regexpr '*.tif']);
else
    dir_list = dir(regexpr);
end

completeErrdir=fullfile(src,errdir,pdir);
if ~exist(completeErrdir)
    completeErrdir=fullfile(errdir,pdir);
end
mkdir(completeErrdir);

tic;
n=size(dir_list, 1);
for i=1:n
    [blueArea,blueMajax]=blue_adjust([platesDir filesep dir_list(i).name],adjust);
    if blueArea>thres || blueMajax>blueSize
        disp(['found ' dir_list(i).name ' : stain']);
        movefile([platesDir filesep dir_list(i).name], ...
            [completeErrdir filesep dir_list(i).name]);
        
        otherChannelFiles=dir([dir_list(i).name(1:end-length(regexpr)) '*.tif']);
        if ~isempty(otherChannelFiles)
            if size(otherChannelFiles,1)>numChannels-1
                % didn't retrieve corresponding files, something went
                % wrong --> do nothing with these files
            else
                % move other channel files too
                for j=1:numel(otherChannelFiles)
                    movefile([platesDir filesep otherChannelFiles(j).name], ...
                        [completeErrdir filesep otherChannelFiles(j).name]);
                end
            end
        end
    end
end
t=toc;
fprintf('time for plate %s: ',pdir); disp(t);
cd(src);
end
