function handles = PauseCellProfiler(handles)
% Help for the PauseCP module:
% Category: Other
%
% SHORT DESCRIPTION:
% Pauses CellProfiler interactively.
% *************************************************************************

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
% $Revision: 4436 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = This module pauses CellProfiler until the dialog box is clicked.

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% The figure window display is unnecessary for this module, so the figure
%%% window is closed the first time through the module.
CPclosefigure(handles,CurrentModule)

drawnow

%% ButtonName = questdlg(Question, Title, Btn1, Btn2,..., DEFAULT);
ButtonName = CPquestdlg('Continue processing?','PauseCP','Continue','Cancel','Continue');
%% TODO - add Modify
% ButtonName = CPquestdlg('Continue or Modify previous
% module?','PauseCP','Continue','Modify','Cancel','Continue');

switch ButtonName
    case 'Continue'
        return
        %% TODO - add Modify
%     case 'Modify'
%         handles.Current.CurrentModuleNumber = num2str(str2num(handles.Current.CurrentModuleNumber) - 1);
%         set(cat(2,handles.VariableBox{:}),'enable','on','foregroundcolor','black'); %% re-enable variable boxes
        
    case 'Cancel'

        %%% This should cause a cancel so no further processing is done
        %%% on this machine.
        set(handles.timertexthandle,'string','Canceling after current module')
end
