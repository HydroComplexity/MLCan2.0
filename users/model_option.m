function varargout = model_option(varargin)
% MODEL_OPTION M-file for model_option.fig
%      MODEL_OPTION, by itself, creates a new MODEL_OPTION or raises the existing
%      singleton*.
%
%      H = MODEL_OPTION returns the handle to a new MODEL_OPTION or the handle to
%      the existing singleton*.
%
%      MODEL_OPTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL_OPTION.M with the given input arguments.
%
%      MODEL_OPTION('Property','Value',...) creates a new MODEL_OPTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before model_option_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to model_option_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help model_option

% Last Modified by GUIDE v2.5 22-Jan-2014 11:21:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @model_option_OpeningFcn, ...
                   'gui_OutputFcn',  @model_option_OutputFcn, ...
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


% --- Executes just before model_option is made visible.
function model_option_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to model_option (see VARARGIN)

% Choose default command line output for model_option
handles.output = hObject;

load './Temps/temporary.mat' 'working_name';
load ('./Temps/temp_variable.mat', 'Sim_species', 'Sim_species_con', 'vanGen', 'RHC', 'CO2_Ambient', 'CO2_Elev', 'CO2_Elev_con', ...
    'Temp_Elev', 'Temp_Elev_con', 'Soil_nutrient', 'Soil_heat', 'Turbulence', 'HR', ...
    'Soil_C_pool1', 'N_denit', 'N_Fix', 'N_Fert', 'N_Fert_DOY', 'N_Fert_amm', 'N_Fert_nit', 'N_Fert_urea', ...
    'num_species', 'N_Uptake_RB', 'N_Remo', 'opt_root_para', 'N_Adepo', 'N_Adepo_amm', 'N_Adepo_nit');

if Sim_species == 1
    set(handles.rad_txt_single_species,'Enable','off');
else
    set(handles.rad_txt_single_species,'Enable','on');
end

if num_species == 1
    set(handles.rad_multi_species_sim,'Value',0);    
    set(handles.rad_txt_single_species,'String',1);
    set(handles.rad_multi_species_sim,'Enable','off'); 
    set(handles.rad_txt_single_species,'Enable','off');
else
    set(handles.rad_multi_species_sim,'Value',Sim_species);    
    set(handles.rad_txt_single_species,'String',Sim_species_con);
    %set(handles.rad_txt_single_species,'Enable','on');
end

set(handles.rad_vanGenuchten,'Value',vanGen);
set(handles.rad_rhc_linear,'Value',RHC);
%set(handles.rad_rnm_implicit,'Value',rnm);
set(handles.rad_txt_CO2_ambient,'String',CO2_Ambient);
set(handles.rad_CO2_elev_on,'Value',CO2_Elev);
set(handles.rad_txt_CO2_ele,'String',CO2_Elev_con);
if CO2_Elev == 0
    set(handles.rad_CO2_elev_off,'Value',1);
end
set(handles.rad_ele_temp_on,'Value',Temp_Elev);
set(handles.rad_txt_temp_ele,'String',Temp_Elev_con);
if Temp_Elev == 0
    set(handles.rad_ele_temp_off,'Value',1);
end

set(handles.chk_soilnutrient,'Value',Soil_nutrient);
set(handles.chk_soilheat,'Value',Soil_heat);
set(handles.chk_Turbulence,'Value',Turbulence);
set(handles.chk_HR,'Value',HR);
%set(handles.chk_entropy,'Value',Entropy);

set(handles.rad_carbon_1pool,'Value',Soil_C_pool1);
set(handles.options_root_tab_para,'Data',opt_root_para);
set(handles.rad_nuptake_rootbiomass,'Value',N_Uptake_RB);
set(handles.chk_fixation,'Value',N_Fix);
set(handles.chk_denitrification,'Value',N_denit);
set(handles.chk_atdeposition,'Value',N_Adepo);
set(handles.txt_amm_dep,'String',N_Adepo_amm);
set(handles.txt_nit_dep,'String',N_Adepo_nit);
set(handles.chk_fertilizer,'Value',N_Fert);
set(handles.txt_fertilizer_DOY,'String',N_Fert_DOY);
set(handles.txt_amm_fertilizer_amount,'String',N_Fert_amm);
set(handles.txt_nit_fertilizer_amount,'String',N_Fert_nit);
set(handles.txt_urea_fertilizer_amount,'String',N_Fert_urea);
    

if CO2_Elev == 1
    set(handles.rad_txt_CO2_ele,'Enable','on');
    set(handles.rad_txt_CO2_ambient,'Enable','off');
else
    set(handles.rad_txt_CO2_ele,'Enable','off');
    set(handles.rad_txt_CO2_ambient,'Enable','on');
end
if Temp_Elev == 1
    set(handles.rad_txt_temp_ele,'Enable','on');
else
    set(handles.rad_txt_temp_ele,'Enable','off');
end

if Soil_nutrient == 1
    set(handles.push_soil_CN_process,'Enable','on');
    set(handles.chk_soilheat,'Value',1);
    set(handles.chk_soilheat,'Enable','off');
else
    set(handles.push_soil_CN_process,'Enable','off');
end
    
if N_Uptake_RB == 0
    set(handles.rad_nuptake_carbon,'Value',1);
    set(handles.chk_fixation,'Enable','on');
    set(handles.chk_remobilization,'Enable','on');
    set(handles.options_root_tab_para,'Enable','off');
   
else
    set(handles.chk_fixation,'Enable','off');
    set(handles.chk_remobilization,'Enable','off');
    set(handles.options_root_tab_para,'Enable','on');
end



% Set constrain for remobilization doy should be 1:365 or 1:366
% and simulation year should be more than 2 years
load './Temps/temp_variable.mat' ...
    'DOY_start' 'DOY_end'
if DOY_start == 1 && (DOY_end == 365 || DOY_end == 366)
   if N_Uptake_RB == 0
       set(handles.chk_remobilization,'Enable','on');
       set(handles.chk_remobilization,'Value',N_Remo);
   end
else
    set(handles.chk_remobilization,'Enable','off');
    N_Remo=0;
    set(handles.chk_remobilization,'Value',N_Remo);
end


if N_Adepo == 1
    set(handles.txt_amm_dep,'Enable','on');
    set(handles.txt_nit_dep,'Enable','on');
else
    set(handles.txt_amm_dep,'Enable','off');
    set(handles.txt_nit_dep,'Enable','off');
end

if N_Fert == 1
    set(handles.txt_fertilizer_DOY,'Enable','on');
    set(handles.txt_amm_fertilizer_amount,'Enable','on');
    set(handles.txt_nit_fertilizer_amount,'Enable','on');
    set(handles.txt_urea_fertilizer_amount,'Enable','on');
else
    set(handles.txt_fertilizer_DOY,'Enable','off');
    set(handles.txt_amm_fertilizer_amount,'Enable','off');
    set(handles.txt_nit_fertilizer_amount,'Enable','off');
    set(handles.txt_urea_fertilizer_amount,'Enable','off');
end


set(handles.uipanel_soil_CN_process,'Visible','off');
set(handles.uipanel_soil_CN_model_process,'Visible','off');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes model_option wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = model_option_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function axes1_CreateFcn(hObject, eventdata, handles)
optioniconsmall = imread('./users/icons/option_small.png');
image(optioniconsmall);
axis off
%imshow('./users/icons/option_small.png');


function mod_option_ok_Callback(hObject, eventdata, handles)
load './Temps/temporary.mat' 'working_name';

Sim_species     = get(handles.rad_multi_species_sim,'Value');
Sim_species_con = str2num(get(handles.rad_txt_single_species,'String'));
vanGen          = get(handles.rad_vanGenuchten,'Value');
RHC             = get(handles.rad_rhc_linear,'Value');
%rnm             = get(handles.rad_rnm_implicit,'Value');
CO2_Ambient     = get(handles.rad_txt_CO2_ambient,'String');
CO2_Elev        = get(handles.rad_CO2_elev_on,'Value');
CO2_Elev_con    = str2num(get(handles.rad_txt_CO2_ele,'String'));
Temp_Elev       = get(handles.rad_ele_temp_on,'Value');
Temp_Elev_con   = str2num(get(handles.rad_txt_temp_ele,'String'));

Soil_nutrient   = get(handles.chk_soilnutrient,'Value');
Soil_heat       = get(handles.chk_soilheat,'Value');
Turbulence      = get(handles.chk_Turbulence,'Value');
HR              = get(handles.chk_HR,'Value');
%Entropy         = get(handles.chk_entropy,'Value');

Soil_C_pool1    = get(handles.rad_carbon_1pool,'Value'); % 1 = 1 pool, 0 = 2 pool
N_Uptake_RB     = get(handles.rad_nuptake_rootbiomass,'Value'); % 1 = Juan, 0 = Dongkook
opt_root_para   = get(handles.options_root_tab_para,'Data');
N_Fix           = get(handles.chk_fixation,'Value');
N_Remo           = get(handles.chk_remobilization,'Value');
N_denit         = get(handles.chk_denitrification,'Value');
N_Adepo         = get(handles.chk_atdeposition,'Value');
N_Adepo_amm     = str2num(get(handles.txt_amm_dep,'String'));
N_Adepo_nit     = str2num(get(handles.txt_nit_dep,'String'));
N_Fert          = get(handles.chk_fertilizer,'Value');
N_Fert_DOY      = str2num(get(handles.txt_fertilizer_DOY,'String'));
N_Fert_amm      = str2num(get(handles.txt_amm_fertilizer_amount,'String'));
N_Fert_nit      = str2num(get(handles.txt_nit_fertilizer_amount,'String'));
N_Fert_urea     = str2num(get(handles.txt_urea_fertilizer_amount,'String'));
N_Fert          = get(handles.chk_fertilizer,'Value');
N_Fert_DOY      = str2num(get(handles.txt_fertilizer_DOY,'String'));
N_Fert_lbha     = str2num(get(handles.txt_urea_fertilizer_amount,'String'));

save './Temps/temp_variable.mat' 'Sim_species' 'Sim_species_con' 'vanGen' 'RHC' 'CO2_Ambient' 'CO2_Elev' 'CO2_Elev_con' ...
    'Temp_Elev' 'Temp_Elev_con' 'Soil_nutrient' 'Soil_heat' 'Turbulence' 'HR' ...
    'Soil_C_pool1' 'N_denit' 'N_Fix' 'N_Fert' 'N_Fert_DOY' 'N_Fert_amm' 'N_Fert_nit' 'N_Fert_urea' ...
    'N_Uptake_RB'   'N_Remo' 'opt_root_para' 'N_Adepo' 'N_Adepo_amm' 'N_Adepo_nit' -append;

load ('./Temps/temp_variable.mat',...
    'num_species');
integerTest=~mod(Sim_species_con,1);
if (Sim_species_con < 1 || Sim_species_con > 4 || integerTest ~= 1)
    msgbox ('Intended single species must be a integer between 1 and 4 ','MLCan Error','error');
    return
else
    if (Sim_species_con > num_species);
        msgbox ('Intended single species must be <= number of species','MLCan Error','error');
        return
    else
    end
end

load ('./Temps/temporary.mat',...
    'working_name'  )
copyfile('./Temps/temp_variable.mat',working_name,'f');

close


function figure1_CreateFcn(hObject, eventdata, handles)


function mod_option_cancel_Callback(hObject, eventdata, handles)
close


function rad_but_C3_CreateFcn(hObject, eventdata, handles)


function rad_txt_CO2_ambient_Callback(hObject, eventdata, handles)


function rad_txt_CO2_ambient_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function rad_txt_CO2_ele_Callback(hObject, eventdata, handles)


function rad_txt_CO2_ele_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function chk_soilheat_Callback(hObject, eventdata, handles)


function chk_Turbulence_Callback(hObject, eventdata, handles)


function chk_HR_Callback(hObject, eventdata, handles)


function chk_carbon_Callback(hObject, eventdata, handles)


function chk_nitrogen_Callback(hObject, eventdata, handles)


function uipanel8_SelectionChangeFcn(hObject, eventdata, handles)
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'rad_CO2_elev_off'
        set(handles.rad_txt_CO2_ele,'Enable','off');
        set(handles.rad_txt_CO2_ambient,'Enable','on');
    case 'rad_CO2_elev_on'
        set(handles.rad_txt_CO2_ele,'Enable','on');
        set(handles.rad_txt_CO2_ambient,'Enable','off');
    otherwise
end
%updates the handles structure
guidata(hObject, handles);


function rad_CO2_elev_off_CreateFcn(hObject, eventdata, handles)


function rad_CO2_elev_on_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in chk_soilnutrient.
function chk_soilnutrient_Callback(hObject, eventdata, handles)
% hObject    handle to chk_soilnutrient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_soilnutrient
handles = guidata(hObject);
Soil_model_on=get(hObject,'Value');
if Soil_model_on == 1
    set(handles.push_soil_CN_process,'Enable','on');
    set(handles.chk_soilheat,'Value',1);
    set(handles.chk_soilheat,'Enable','off');    
else
    set(handles.push_soil_CN_process,'Enable','off');
    set(handles.chk_soilheat,'Enable','on');
end
guidata(hObject, handles);

% --- Executes on button press in chk_entropy.
function chk_entropy_Callback(hObject, eventdata, handles)
% hObject    handle to chk_entropy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_entropy



function rad_txt_temp_ele_Callback(hObject, eventdata, handles)
% hObject    handle to rad_txt_temp_ele (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rad_txt_temp_ele as text
%        str2double(get(hObject,'String')) returns contents of rad_txt_temp_ele as a double


% --- Executes during object creation, after setting all properties.
function rad_txt_temp_ele_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rad_txt_temp_ele (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rad_txt_single_species_Callback(hObject, eventdata, handles)
% hObject    handle to rad_txt_single_species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rad_txt_single_species as text
%        str2double(get(hObject,'String')) returns contents of rad_txt_single_species as a double


% --- Executes during object creation, after setting all properties.
function rad_txt_single_species_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rad_txt_single_species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_soil_CN_process.
function push_soil_CN_process_Callback(hObject, eventdata, handles)
% hObject    handle to push_soil_CN_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanel_canopy_process,'Visible','off');
set(handles.uipanel_canopy_model_process,'Visible','off');
set(handles.uipanel_soil_CN_process,'Visible','on');
set(handles.uipanel_soil_CN_model_process,'Visible','on');

% --- Executes on button press in push_canopy_process.
function push_canopy_process_Callback(hObject, eventdata, handles)
% hObject    handle to push_canopy_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanel_canopy_process,'Visible','on');
set(handles.uipanel_canopy_model_process,'Visible','on');
set(handles.uipanel_soil_CN_process,'Visible','off');
set(handles.uipanel_soil_CN_model_process,'Visible','off');

% --- Executes on button press in chk_denitrification.
function chk_denitrification_Callback(hObject, eventdata, handles)
% hObject    handle to chk_denitrification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_denitrification


% --- Executes on button press in chk_fixation.
function chk_fixation_Callback(hObject, eventdata, handles)
% hObject    handle to chk_fixation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of chk_fixation


% --- Executes on button press in chk_fertilizer.
function chk_fertilizer_Callback(hObject, eventdata, handles)
% hObject    handle to chk_fertilizer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of chk_fertilizer
handles = guidata(hObject);
N_Fert_on=get(hObject,'Value');
if N_Fert_on == 1
    set(handles.txt_fertilizer_DOY,'Enable','on');
    set(handles.txt_amm_fertilizer_amount,'Enable','on');
    set(handles.txt_nit_fertilizer_amount,'Enable','on');
    set(handles.txt_urea_fertilizer_amount,'Enable','on');
else
    set(handles.txt_fertilizer_DOY,'Enable','off');
    set(handles.txt_amm_fertilizer_amount,'Enable','off');
    set(handles.txt_nit_fertilizer_amount,'Enable','off');
    set(handles.txt_urea_fertilizer_amount,'Enable','off');
end
guidata(hObject, handles);




function txt_fixation_Callback(hObject, eventdata, handles)
% hObject    handle to txt_fixation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_fixation as text
%        str2double(get(hObject,'String')) returns contents of txt_fixation as a double


% --- Executes during object creation, after setting all properties.
function txt_fixation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_fixation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_fertilizer_DOY_Callback(hObject, eventdata, handles)
% hObject    handle to txt_fertilizer_DOY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_fertilizer_DOY as text
%        str2double(get(hObject,'String')) returns contents of txt_fertilizer_DOY as a double


% --- Executes during object creation, after setting all properties.
function txt_fertilizer_DOY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_fertilizer_DOY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_urea_fertilizer_amount_Callback(hObject, eventdata, handles)
% hObject    handle to txt_urea_fertilizer_amount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_urea_fertilizer_amount as text
%        str2double(get(hObject,'String')) returns contents of txt_urea_fertilizer_amount as a double


% --- Executes during object creation, after setting all properties.
function txt_urea_fertilizer_amount_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_urea_fertilizer_amount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel12.
function uipanel12_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel12 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'rad_ele_temp_off'
        set(handles.rad_txt_temp_ele,'Enable','off');
    case 'rad_ele_temp_on'
        set(handles.rad_txt_temp_ele,'Enable','on');
    otherwise
end
%updates the handles structure
guidata(hObject, handles);


% --- Executes when selected object is changed in uipanel16.
function uipanel16_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel16 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'rad_multi_species_sim'
        set(handles.rad_txt_single_species,'Enable','off');
    case 'rad_single_species_sim'
        set(handles.rad_txt_single_species,'Enable','on');
    otherwise
end
%updates the handles structure
guidata(hObject, handles);


% --- Executes on key press with focus on chk_soilnutrient and none of its controls.
function chk_soilnutrient_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to chk_soilnutrient (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chk_litter.
function chk_litter_Callback(hObject, eventdata, handles)
% hObject    handle to chk_litter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_litter



function chk_txt_litterM_Callback(hObject, eventdata, handles)
% hObject    handle to chk_txt_litterM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chk_txt_litterM as text
%        str2double(get(hObject,'String')) returns contents of chk_txt_litterM as a double


% --- Executes during object creation, after setting all properties.
function chk_txt_litterM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chk_txt_litterM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chk_remobilization.
function chk_remobilization_Callback(hObject, eventdata, handles)
% hObject    handle to chk_remobilization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_remobilization

function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to txt_fixation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_fixation as text
%        str2double(get(hObject,'String')) returns contents of txt_fixation as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_fixation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_remobilization_Callback(hObject, eventdata, handles)
% hObject    handle to txt_remobilization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_remobilization as text
%        str2double(get(hObject,'String')) returns contents of txt_remobilization as a double


% --- Executes during object creation, after setting all properties.
function txt_remobilization_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_remobilization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel44.
function uipanel44_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel44 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'rad_nuptake_rootbiomass'
        set(handles.chk_fixation,'Enable','off');
        set(handles.chk_remobilization,'Enable','off');
        N_Remo=0;
        set(handles.chk_remobilization,'Value',N_Remo);
        set(handles.options_root_tab_para,'Enable','on');
    case 'rad_nuptake_carbon' 
        set(handles.chk_fixation,'Enable','on');
        
        % Set constrain for remobilization doy should be 1:365 or 1:366
        % and simulation year should be more than 2 years
        load './Temps/temp_variable.mat' ...
            'DOY_start' 'DOY_end' 'N_Remo'
        if DOY_start == 1 && (DOY_end == 365 || DOY_end == 366)
            set(handles.chk_remobilization,'Enable','on');
        else
            set(handles.chk_remobilization,'Enable','off');
            N_Remo=0;
            set(handles.chk_remobilization,'Value',N_Remo);
        end

        set(handles.options_root_tab_para,'Enable','off');
    otherwise
end
%updates the handles structure
guidata(hObject, handles);


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



function edit41_Callback(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit41 as text
%        str2double(get(hObject,'String')) returns contents of edit41 as a double


% --- Executes during object creation, after setting all properties.
function edit41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
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



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit44 as text
%        str2double(get(hObject,'String')) returns contents of edit44 as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit46_Callback(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit46 as text
%        str2double(get(hObject,'String')) returns contents of edit46 as a double


% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_nit_dep_Callback(hObject, eventdata, handles)
% hObject    handle to txt_nit_dep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_nit_dep as text
%        str2double(get(hObject,'String')) returns contents of txt_nit_dep as a double


% --- Executes during object creation, after setting all properties.
function txt_nit_dep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_nit_dep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_amm_dep_Callback(hObject, eventdata, handles)
% hObject    handle to txt_amm_dep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_amm_dep as text
%        str2double(get(hObject,'String')) returns contents of txt_amm_dep as a double


% --- Executes during object creation, after setting all properties.
function txt_amm_dep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_amm_dep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chk_atdeposition.
function chk_atdeposition_Callback(hObject, eventdata, handles)
% hObject    handle to chk_atdeposition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_atdeposition
handles = guidata(hObject);
N_Adepo_on=get(hObject,'Value');
if N_Adepo_on == 1
    set(handles.txt_amm_dep,'Enable','on');
    set(handles.txt_nit_dep,'Enable','on');
else
    set(handles.txt_amm_dep,'Enable','off');
    set(handles.txt_nit_dep,'Enable','off');
end
guidata(hObject, handles);



function txt_nit_fertilizer_amount_Callback(hObject, eventdata, handles)
% hObject    handle to txt_nit_fertilizer_amount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_nit_fertilizer_amount as text
%        str2double(get(hObject,'String')) returns contents of txt_nit_fertilizer_amount as a double


% --- Executes during object creation, after setting all properties.
function txt_nit_fertilizer_amount_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_nit_fertilizer_amount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_amm_fertilizer_amount_Callback(hObject, eventdata, handles)
% hObject    handle to txt_amm_fertilizer_amount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_amm_fertilizer_amount as text
%        str2double(get(hObject,'String')) returns contents of txt_amm_fertilizer_amount as a double


% --- Executes during object creation, after setting all properties.
function txt_amm_fertilizer_amount_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_amm_fertilizer_amount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
