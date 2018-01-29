function handles = MedianFilter(handles)
  
% Help for the Median Filter module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
%
%
% *************************************************************************
%
%
% See also ColorToGray.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 5025 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the input image?
%infotypeVAR01 = imagegroup
InputImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu custom


%textVAR02 = What do you want to call the resulting image?
%defaultVAR02 = FilteredImage
%infotypeVAR02 = imagegroup indep
FilteredImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Enter the size of the filter
%defaultVAR03 = 5
FilterSize = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%%%VariableRevisionNumber = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Determines whether the user has specified an image to be loaded in
%%% blue.
 if ~strcmp(InputImageName, '')
     %%% Read (open) the images and assign them to variables.
     InputImage = CPretrieveimage(handles,InputImageName,ModuleName,'MustBeGray','CheckScale');
     InputImageExists = 1;
 else
     InputImageExists = 0;
 end
 drawnow

%%% Check to see if the filter size is bigger than 0

if (str2double(FilterSize) < 0.0) || isnan(str2double(FilterSize))
     if isempty(findobj('Tag',['Msgbox_' ModuleName ', ModuleNumber ' num2str(CurrentModuleNum) ': Filter size invalid']))
         CPwarndlg(['The filter size you have entered in the ' ModuleName ' module is invalid or less than 0. It is being set to the default value of 5.'],[ModuleName ', ModuleNumber ' num2str(CurrentModuleNum) ': Filter size invalid']);
     end
     FilterSize = '5';
end

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% If any of the images are binary/logical format, they must be
%%% converted to a double first before immultiply.
filteredImage = medfilt2(InputImage, [str2double(FilterSize) str2double(FilterSize)]);

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

 ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
 if any(findobj == ThisModuleFigureNumber)
     %%% Activates the appropriate figure window.
     
      CPfigure(handles,'Image',ThisModuleFigureNumber);
      if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
          CPresizefigure(filteredImage,'TwoByOne',ThisModuleFigureNumber);
      end
 
     %%% A subplot of the original image

      subplot(2,1,1); 
      CPimagesc(InputImage,handles);
      title(['Original Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);

     %%% A subplot of the filtered image.

      subplot(2,1,2); 
      CPimagesc(filteredImage,handles);
     title('Filtered Image');

 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the filtered image to the handles structure so it can be used by
%%% subsequent modules.
handles.Pipeline.(FilteredImageName) = filteredImage;