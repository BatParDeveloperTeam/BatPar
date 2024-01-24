function varargout = User_Inputs(varargin)
% USER_INPUTS MATLAB code for User_Inputs.fig
%      USER_INPUTS, by itself, creates a new USER_INPUTS or raises the existing
%      singleton*.
%
%      H = USER_INPUTS returns the handle to a new USER_INPUTS or the handle to
%      the existing singleton*.
%
%      USER_INPUTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in USER_INPUTS.M with the given input arguments.
%
%      USER_INPUTS('Property','Value',...) creates a new USER_INPUTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before User_Inputs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to User_Inputs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help User_Inputs

% Last Modified by GUIDE v2.5 02-Sep-2022 17:30:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @User_Inputs_OpeningFcn, ...
                   'gui_OutputFcn',  @User_Inputs_OutputFcn, ...
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


% --- Executes just before User_Inputs is made visible.
function User_Inputs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to User_Inputs (see VARARGIN)

% Choose default command line output for User_Inputs

% background image
handles.output = hObject;

ha=axes('units','normalized','pos',[0 0 1 1]);

uistack(ha,'down');

ii=imread('Background image1.jpg');

image(ii);

colormap gray(1)

set(ha,'handlevisibility','off','visible','off');

disp('2022, Electrochemical Science and Engineering Group, Imperial College London, UK.')
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes User_Inputs wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = User_Inputs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

exp_path = get(handles.edit35, 'string');
assignin('base','exp_path',exp_path);

record_column = get(handles.edit65, 'string');
assignin('base','record_column',record_column);

time_column = get(handles.edit38, 'string');
assignin('base','time_column',time_column);

current_column = get(handles.edit39, 'string');
assignin('base','current_column',current_column);

voltage_column = get(handles.edit40, 'string');
assignin('base','voltage_column',voltage_column);

temp_column = get(handles.edit73, 'string');
assignin('base','temp_column',temp_column);

row_start = get(handles.edit36, 'string');
assignin('base','row_start',row_start);

row_end = get(handles.edit55, 'string');
assignin('base','row_end',row_end);

Nom_capacity = get(handles.edit42, 'string');
assignin('base','Nom_capacity',Nom_capacity);

voltage_up_limit = get(handles.edit43, 'string');
assignin('base','voltage_up_limit',voltage_up_limit);

voltage_low_limit = get(handles.edit56, 'string');
assignin('base','voltage_low_limit',voltage_low_limit);

OCVsheet_name = get(handles.edit47, 'string');
assignin('base','OCVsheet_name',OCVsheet_name);

R0sheet_name = get(handles.edit74, 'string');
assignin('base','R0sheet_name',R0sheet_name);

RCpair_number = get(handles.edit51, 'string');
assignin('base','RCpair_number',RCpair_number);

SOC_window = get(handles.edit50, 'string');
assignin('base','SOC_window',SOC_window);

SOC_param_lowlimit = get(handles.edit62, 'string');
assignin('base','SOC_param_lowlimit',SOC_param_lowlimit);

SOC_param_uplimit = get(handles.edit66, 'string');
assignin('base','SOC_param_uplimit',SOC_param_uplimit);

tau1_limit = get(handles.edit58, 'string');
assignin('base','tau1_limit',tau1_limit);

tau2_limit = get(handles.edit59, 'string');
assignin('base','tau2_limit',tau2_limit);

tau3_limit = get(handles.edit60, 'string');
assignin('base','tau3_limit',tau3_limit);

Ri_limit = get(handles.edit61, 'string');
assignin('base','Ri_limit',Ri_limit);

Currents_4_Denp = get(handles.edit67, 'string');
assignin('base','Currents_4_Denp',Currents_4_Denp);

CurrentSensi_4_Denp = get(handles.edit68, 'string');
assignin('base','CurrentSensi_4_Denp',CurrentSensi_4_Denp);

Temps_4_Denp = get(handles.edit71, 'string');
assignin('base','Temps_4_Denp',Temps_4_Denp);

TempSensi_4_Denp = get(handles.edit72, 'string');
assignin('base','TempSensi_4_Denp',TempSensi_4_Denp);


%OPTfactor = get(handles.edit63, 'string');
%assignin('base','OPTfactor',OPTfactor);

if get(handles.radiobutton3,'Value')==1
    assignin('base','Current_belwozero_tf','yes');
else
    assignin('base','Current_belwozero_tf','no'); 
end

if get(handles.radiobutton11,'Value')==1
    assignin('base','MaxSOCis1','yes');
else
    assignin('base','MaxSOCis1','no'); 
end

if get(handles.radiobutton17,'Value')==1
    assignin('base','Mode_anode_tf','yes');
else
    assignin('base','Mode_anode_tf','no'); 
end

if get(handles.radiobutton13,'Value')==1
    assignin('base','Mode_charge_tf','yes');
else
    assignin('base','Mode_charge_tf','no'); 
end

if get(handles.radiobutton15,'Value')==1
    assignin('base','Mode_cathode_tf','yes');OCVtable
else
    assignin('base','Mode_cathode_tf','no'); 
end

if get(handles.radiobutton18,'Value')==1
    assignin('base','Pulse_head_tf','yes');
else
    assignin('base','Pulse_head_tf','no'); 
end

if get(handles.radiobutton26,'Value')==1
    assignin('base','fixedtau_realy_tf','yes');
else
    assignin('base','fixedtau_realy_tf','no'); 
end

if get(handles.radiobutton27,'Value')==1
    assignin('base','OCVpseudo_tf','yes');
else
    assignin('base','OCVpseudo_tf','no'); 
end

if get(handles.radiobutton28,'Value')==1
    assignin('base','OCVtrue_tf','yes');
else
    assignin('base','OCVtrue_tf','no'); 
end

if get(handles.radiobutton29,'Value')==1
    assignin('base','R0only_tf','yes');
else
    assignin('base','R0only_tf','no'); 
end

if get(handles.checkbox1,'Value')==1
    assignin('base','OCVtable_tf','yes');
else
    assignin('base','OCVtable_tf','no');
end 

if get(handles.checkbox12,'Value')==1
    assignin('base','R0table_tf','yes');
else
    assignin('base','R0table_tf','no');
end 

if get(handles.checkbox5,'Value')==1
    assignin('base','datasheet_cleanup_tf','yes');
else
    assignin('base','datasheet_cleanup_tf','no');
end

if get(handles.checkbox6,'Value')==1
    assignin('base','Ri_taui_relimit_tf','yes');
else
    assignin('base','Ri_taui_relimit_tf','no');
end 

if get(handles.checkbox9,'Value')==1
    assignin('base','CurrentDenp_tf','yes');
else
    assignin('base','CurrentDenp_tf','no');
end 

if get(handles.checkbox11,'Value')==1
    assignin('base','TempDenp_tf','yes');
else
    assignin('base','TempDenp_tf','no');
end 

run("SingleRun.m")
%if get(handles.checkbox7,'Value')==1
%    assignin('base','deeperOP_tf','yes');
%else
%    assignin('base','deeperOP_tf','no');
%end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
if get(handles.checkbox1,'Value')==1
    set(handles.edit47,'visible','on');
else
    set(handles.edit47,'visible','off'); 
end


function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double


% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double


% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit40_Callback(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit40 as text
%        str2double(get(hObject,'String')) returns contents of edit40 as a double


% --- Executes during object creation, after setting all properties.
function edit40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit42_Callback(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit42 as text
%        str2double(get(hObject,'String')) returns contents of edit42 as a double


% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit43 as text
%        str2double(get(hObject,'String')) returns contents of edit43 as a double


% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit47 as text
%        str2double(get(hObject,'String')) returns contents of edit47 as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit50 as text
%        str2double(get(hObject,'String')) returns contents of edit50 as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit51 as text
%        str2double(get(hObject,'String')) returns contents of edit51 as a double


% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit55_Callback(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit55 as text
%        str2double(get(hObject,'String')) returns contents of edit55 as a double


% --- Executes during object creation, after setting all properties.
function edit55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit56_Callback(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit56 as text
%        str2double(get(hObject,'String')) returns contents of edit56 as a double


% --- Executes during object creation, after setting all properties.
function edit56_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
if get(handles.checkbox5,'Value')==1
    set(handles.edit36,'visible','on');
    set(handles.edit55,'visible','on');
else
    set(handles.edit36,'visible','off'); 
    set(handles.edit55,'visible','off');
end



% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
if get(handles.checkbox6,'Value')==1
    set(handles.edit58,'visible','on');
    set(handles.edit59,'visible','on');
    set(handles.edit60,'visible','on');
    set(handles.edit61,'visible','on');
else
    set(handles.edit58,'visible','off'); 
    set(handles.edit59,'visible','off');
    set(handles.edit60,'visible','off'); 
    set(handles.edit61,'visible','off');
end


function edit58_Callback(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit58 as text
%        str2double(get(hObject,'String')) returns contents of edit58 as a double


% --- Executes during object creation, after setting all properties.
function edit58_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit59_Callback(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit59 as text
%        str2double(get(hObject,'String')) returns contents of edit59 as a double


% --- Executes during object creation, after setting all properties.
function edit59_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit60_Callback(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit60 as text
%        str2double(get(hObject,'String')) returns contents of edit60 as a double


% --- Executes during object creation, after setting all properties.
function edit60_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit61_Callback(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit61 as text
%        str2double(get(hObject,'String')) returns contents of edit61 as a double


% --- Executes during object creation, after setting all properties.
function edit61_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton11.
function radiobutton11_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton11



function edit62_Callback(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit62 as text
%        str2double(get(hObject,'String')) returns contents of edit62 as a double


% --- Executes during object creation, after setting all properties.
function edit62_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit63_Callback(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit63 as text
%        str2double(get(hObject,'String')) returns contents of edit63 as a double


% --- Executes during object creation, after setting all properties.
function edit63_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text58_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in radiobutton12.
function radiobutton12_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton12


% --- Executes on button press in radiobutton17.
function radiobutton17_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton17


% --- Executes on button press in radiobutton13.
function radiobutton13_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton13



function edit65_Callback(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit65 as text
%        str2double(get(hObject,'String')) returns contents of edit65 as a double


% --- Executes during object creation, after setting all properties.
function edit65_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit66_Callback(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit66 as text
%        str2double(get(hObject,'String')) returns contents of edit66 as a double


% --- Executes during object creation, after setting all properties.
function edit66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton15.
function radiobutton15_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton15


% --- Executes on button press in radiobutton18.
function radiobutton18_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton18


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exp_path = get(handles.edit35, 'string');
assignin('base','exp_path',exp_path);

record_column = get(handles.edit65, 'string');
assignin('base','record_column',record_column);

time_column = get(handles.edit38, 'string');
assignin('base','time_column',time_column);

current_column = get(handles.edit39, 'string');
assignin('base','current_column',current_column);

voltage_column = get(handles.edit40, 'string');
assignin('base','voltage_column',voltage_column);

temp_column = get(handles.edit73, 'string');
assignin('base','temp_column',temp_column);

row_start = get(handles.edit36, 'string');
assignin('base','row_start',row_start);

row_end = get(handles.edit55, 'string');
assignin('base','row_end',row_end);

Nom_capacity = get(handles.edit42, 'string');
assignin('base','Nom_capacity',Nom_capacity);

voltage_up_limit = get(handles.edit43, 'string');
assignin('base','voltage_up_limit',voltage_up_limit);

voltage_low_limit = get(handles.edit56, 'string');
assignin('base','voltage_low_limit',voltage_low_limit);

OCVsheet_name = get(handles.edit47, 'string');
assignin('base','OCVsheet_name',OCVsheet_name);

R0sheet_name = get(handles.edit74, 'string');
assignin('base','R0sheet_name',R0sheet_name);

RCpair_number = get(handles.edit51, 'string');
assignin('base','RCpair_number',RCpair_number);

SOC_window = get(handles.edit50, 'string');
assignin('base','SOC_window',SOC_window);

SOC_param_lowlimit = get(handles.edit62, 'string');
assignin('base','SOC_param_lowlimit',SOC_param_lowlimit);

SOC_param_uplimit = get(handles.edit66, 'string');
assignin('base','SOC_param_uplimit',SOC_param_uplimit);

tau1_limit = get(handles.edit58, 'string');
assignin('base','tau1_limit',tau1_limit);

tau2_limit = get(handles.edit59, 'string');
assignin('base','tau2_limit',tau2_limit);

tau3_limit = get(handles.edit60, 'string');
assignin('base','tau3_limit',tau3_limit);

Ri_limit = get(handles.edit61, 'string');
assignin('base','Ri_limit',Ri_limit);

Currents_4_Denp = get(handles.edit67, 'string');
assignin('base','Currents_4_Denp',Currents_4_Denp);

CurrentSensi_4_Denp = get(handles.edit68, 'string');
assignin('base','CurrentSensi_4_Denp',CurrentSensi_4_Denp);

Temps_4_Denp = get(handles.edit71, 'string');
assignin('base','Temps_4_Denp',Temps_4_Denp);

TempSensi_4_Denp = get(handles.edit72, 'string');
assignin('base','TempSensi_4_Denp',TempSensi_4_Denp);


%OPTfactor = get(handles.edit63, 'string');
%assignin('base','OPTfactor',OPTfactor);

if get(handles.radiobutton3,'Value')==1
    assignin('base','Current_belwozero_tf','yes');
else
    assignin('base','Current_belwozero_tf','no'); 
end

if get(handles.radiobutton11,'Value')==1
    assignin('base','MaxSOCis1','yes');
else
    assignin('base','MaxSOCis1','no'); 
end

if get(handles.radiobutton17,'Value')==1
    assignin('base','Mode_anode_tf','yes');
else
    assignin('base','Mode_anode_tf','no'); 
end

if get(handles.radiobutton13,'Value')==1
    assignin('base','Mode_charge_tf','yes');
else
    assignin('base','Mode_charge_tf','no'); 
end

if get(handles.radiobutton15,'Value')==1
    assignin('base','Mode_cathode_tf','yes');OCVtable
else
    assignin('base','Mode_cathode_tf','no'); 
end

if get(handles.radiobutton18,'Value')==1
    assignin('base','Pulse_head_tf','yes');
else
    assignin('base','Pulse_head_tf','no'); 
end

if get(handles.radiobutton26,'Value')==1
    assignin('base','fixedtau_realy_tf','yes');
else
    assignin('base','fixedtau_realy_tf','no'); 
end

if get(handles.radiobutton27,'Value')==1
    assignin('base','OCVpseudo_tf','yes');
else
    assignin('base','OCVpseudo_tf','no'); 
end

if get(handles.radiobutton28,'Value')==1
    assignin('base','OCVtrue_tf','yes');
else
    assignin('base','OCVtrue_tf','no'); 
end

if get(handles.radiobutton29,'Value')==1
    assignin('base','R0only_tf','yes');
else
    assignin('base','R0only_tf','no'); 
end

if get(handles.checkbox1,'Value')==1
    assignin('base','OCVtable_tf','yes');
else
    assignin('base','OCVtable_tf','no');
end 

if get(handles.checkbox12,'Value')==1
    assignin('base','R0table_tf','yes');
else
    assignin('base','R0table_tf','no');
end 


if get(handles.checkbox5,'Value')==1
    assignin('base','datasheet_cleanup_tf','yes');
else
    assignin('base','datasheet_cleanup_tf','no');
end

if get(handles.checkbox6,'Value')==1
    assignin('base','Ri_taui_relimit_tf','yes');
else
    assignin('base','Ri_taui_relimit_tf','no');
end 

if get(handles.checkbox9,'Value')==1
    assignin('base','CurrentDenp_tf','yes');
else
    assignin('base','CurrentDenp_tf','no');
end 

if get(handles.checkbox11,'Value')==1
    assignin('base','TempDenp_tf','yes');
else
    assignin('base','TempDenp_tf','no');
end 


run("ConstantRun.m")


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

run("VariableRun.m")


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9
if get(handles.checkbox9,'Value')==1
    set(handles.edit67,'visible','on');
    set(handles.edit68,'visible','on');
else
    set(handles.edit67,'visible','off'); 
    set(handles.edit68,'visible','off');
end


function edit67_Callback(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit67 as text
%        str2double(get(hObject,'String')) returns contents of edit67 as a double


% --- Executes during object creation, after setting all properties.
function edit67_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit68_Callback(hObject, eventdata, handles)
% hObject    handle to edit68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit68 as text
%        str2double(get(hObject,'String')) returns contents of edit68 as a double


% --- Executes during object creation, after setting all properties.
function edit68_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11
if get(handles.checkbox11,'Value')==1
    set(handles.edit71,'visible','on');
    set(handles.edit72,'visible','on');
else
    set(handles.edit71,'visible','off'); 
    set(handles.edit72,'visible','off');
end


function edit71_Callback(hObject, eventdata, handles)
% hObject    handle to edit71 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit71 as text
%        str2double(get(hObject,'String')) returns contents of edit71 as a double


% --- Executes during object creation, after setting all properties.
function edit71_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit71 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit72_Callback(hObject, eventdata, handles)
% hObject    handle to edit72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit72 as text
%        str2double(get(hObject,'String')) returns contents of edit72 as a double


% --- Executes during object creation, after setting all properties.
function edit72_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit73_Callback(hObject, eventdata, handles)
% hObject    handle to edit73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit73 as text
%        str2double(get(hObject,'String')) returns contents of edit73 as a double


% --- Executes during object creation, after setting all properties.
function edit73_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox12
if get(handles.checkbox12,'Value')==1
    set(handles.edit74,'visible','on');
else
    set(handles.edit74,'visible','off'); 
end


function edit74_Callback(hObject, eventdata, handles)
% hObject    handle to edit74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit74 as text
%        str2double(get(hObject,'String')) returns contents of edit74 as a double


% --- Executes during object creation, after setting all properties.
function edit74_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton29.
function radiobutton29_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton29


% --- Executes during object creation, after setting all properties.
function uibuttongroup14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in radiobutton30.
function radiobutton30_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton30
