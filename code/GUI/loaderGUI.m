function varargout = loaderGUI(varargin)
% LOADERGUI MATLAB code for loaderGUI.fig
%      LOADERGUI, by itself, creates a new LOADERGUI or raises the existing
%      singleton*.
%
%      H = LOADERGUI returns the handle to a new LOADERGUI or the handle to
%      the existing singleton*.
%
%      LOADERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOADERGUI.M with the given input arguments.
%
%      LOADERGUI('Property','Value',...) creates a new LOADERGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before loaderGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to loaderGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help loaderGUI

% Last Modified by GUIDE v2.5 12-Apr-2013 16:10:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @loaderGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @loaderGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before loaderGUI is made visible.
function loaderGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to loaderGUI (see VARARGIN)

% Choose default command line output for loaderGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%keyboard
% Update settings structure to default
setappdata(handles.settingsPanel, 'params', defaultParams);
% UIWAIT makes loaderGUI wait for user response (see UIRESUME)
% uiwait(handles.mainFigure);


% --- Outputs from this function are returned to the command line.
function varargout = loaderGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%keyboard
varargout{1} = handles.output;


% --- Executes on button press in LoadSpike.
function LoadSpike_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSpike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[matFile, matpath] = uigetfile('*.mat','Please choose the song Spike2 file','data');
birdID = regexp(matpath, '\w{1,2}\d{1,3}','match'); 
if ~isempty(birdID)
    birdID = birdID{:};
else 
    birdID = '?';
end
sessionID = regexp(matFile, ...
    sprintf('(?<=%s_)\\d{1,2}_\\d{1,2}_\\d{1,3}',birdID),'match');
if ~isempty(sessionID)
    sessionID = sessionID{:};
else 
    sessionID = '?';
end
set(handles.birdField,'String', sprintf('Bird: %s', birdID));
set(handles.sessionField,'String', sprintf('Session: %s',sessionID));

% load file
songStruct = load([matpath matFile]);
% trick to get the main struct into a standard name, if there's only one
% variable in the file
fld=fieldnames(songStruct);
songStruct=songStruct.(fld{1});
fs = 1/songStruct.interval;
songStruct.fs = fs;
songStruct.title = matFile;
% set variables in the guiData
setappdata(handles.mainFigure, 'songStruct', songStruct);


% --- Executes on button press in FindSounds.
function FindSounds_Callback(hObject, eventdata, handles)
% hObject    handle to FindSounds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hfig = figure;
sounds = stepSpectrogram(getappdata(handles.mainFigure, 'songStruct'), ...
                         getappdata(handles.settingsPanel, 'params'),...
                         'Nsplits',400,'plot',true);
close(hfig);
divs = guidata(handles.divisionsPanel); 
setappdata(handles.divisionsPanel, 'sounds', divs);  

% update workspace
updateWorkspace(hObject, eventdata, handles);


% --- Executes on button press in keyboardAccess.
function keyboardAccess_Callback(hObject, eventdata, handles)
% hObject    handle to keyboardAccess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard
% update workspace
updateWorkspace(hObject, eventdata, handles);



% --- Executes on button press in loadVariables.
function loadVariables_Callback(hObject, eventdata, handles)
% hObject    handle to loadVariables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[tmpFile tmpPath] = uigetfile('*.mat','Select matfile to load');
tmpStruct = load([tmpPath filesep tmpFile]);
varNames = fieldnames(tmpStruct);
for ii = 1:numel(varNames);
    setappdata(handles.divisionsPanel, varNames{ii}, tmpStruct.(varNames{ii}));
end

% update workspace
updateWorkspace(hObject, eventdata, handles);

function updateWorkspace(hObject, eventdata, handles)
    values = getappdata(handles.divisionsPanel);
    fNames = fieldnames(values);
    fNames(strcmp('lastValidTag',fNames)) = [];
    set(handles.divisionsPanel,'String', fNames);
    
