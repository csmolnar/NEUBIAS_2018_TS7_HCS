function handles = CorrectIllumination_Apply_FP_v7(handles)

% Help for the Filippo Piccinini (FP) Correct Illumination Apply module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Applies an illumination function, created by a Filippo Piccinini 
% Matlab external codex, to an image in order to correct for uneven 
% illumination (uneven shading).
% *************************************************************************
%
% This module corrects for uneven illumination of each image. An
% illumination function image that represents the variation in
% illumination across the field of view is made with a Filippo Piccinini 
% Matlab external codex,
% This module applies the illumination function to each image coming 
% through the pipeline to produce the corrected image.
%
% They are loaded two curves relative to the illumination field. 
% The first come from the foreground, acquired from a image with a dense 
% uniform virtual object that covers the entire field of view of the 
% camera. Here is expressed the effect of the uneven illumination like 
% in the objects present in the images that must to be corrected. 
% The second curve comes from a dense "virtual" background, 
% that covers the entire field of view of the camera. 
% Here is expressed the effect of the uneven illumination and the 
% problems due to the camera noise and the border effect. 
%
% For the implemented formula of the flat field correction, 
% normalized to the mean value of the light field from the foreground, 
% seeing the article:
% Jericevic Z, Wiese B, Bryan J, Smith LC. Validation of an imaging
% system: steps to evaluate and validate a microscope imaging
% system for quantitative studies. Methods Cell Biol 1989;
% 30:47–83.

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the image to be corrected?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the corrected image?
%defaultVAR02 = CorrBlue
%infotypeVAR02 = imagegroup indep
CorrectedImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What did you call the illumination function detected from the FOREGROUND (loaded as a .mat format image using Load Single Image)?
%infotypeVAR03 = imagegroup
FG_IF_ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = What did you call the illumination function detected from the BACKGROUND (loaded as a .mat format image using Load Single Image)?
%infotypeVAR04 = imagegroup
BG_IF_ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = Enter the intensity from the original matrix that should be set to the highest value in the rescaled image.
%choiceVAR05 = 0.0039
%choiceVAR05 = 0.0625
%choiceVAR05 = 1
ScaleFactor = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,5}))*2^16;
%inputtypeVAR05 = popupmenu

%textVAR06 = Enter the value of the mean of the foreground of the reference.
%defaultVAR06 = 1
Mean_Reference_FG = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,6}))/ScaleFactor;

%textVAR07 = Enter the value of the mean of the background of the reference.
%defaultVAR07 = 0
Mean_Reference_BG = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,7}))/ScaleFactor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable.
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');

%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable.
FG_IF_Matrix = CPretrieveimage(handles,FG_IF_ImageName,ModuleName,'MustBeGray','DontCheckScale',size(OrigImage));

%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable.
BG_IF_Matrix = CPretrieveimage(handles,BG_IF_ImageName,ModuleName,'MustBeGray','DontCheckScale',size(OrigImage));

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

BG_IF_Matrix = BG_IF_Matrix./ScaleFactor;
FG_IF_Matrix = FG_IF_Matrix./ScaleFactor;

ImageLessBG = double(OrigImage) - double(BG_IF_Matrix);
FGLessBG = double(FG_IF_Matrix) - double(BG_IF_Matrix);
CorrectedImage = ((double(ImageLessBG)./double(FGLessBG)).*double(Mean_Reference_FG-Mean_Reference_BG))+Mean_Reference_BG; %Peter Horvath Law
%%%%CorrectedImage(CorrectedImage < 0) = 0;

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
    end
    %%% A subplot of the figure window is set to display the original
    %%% image, some intermediate images, and the final corrected image.
    subplot(2,1,1);
    CPimagesc(OrigImage,handles);
    title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% The mean image does not absolutely have to be present in order to
    %%% carry out the calculations if the illumination image is provided,
    %%% so the following subplot is only shown if MeanImage exists in the
    %%% workspace.
    subplot(2,1,2);
    CPimagesc(CorrectedImage,handles);
    title('Illumination Corrected Image');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the corrected image to the
%%% handles structure so it can be used by subsequent modules.
handles.Pipeline.(CorrectedImageName) = CorrectedImage;