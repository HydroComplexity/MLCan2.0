function varargout = model_parameters(varargin)
% MODEL_PARAMETERS M-file for model_parameters.fig
%      MODEL_PARAMETERS, by itself, creates a new MODEL_PARAMETERS or raises the existing
%      singleton*.
%
%      H = MODEL_PARAMETERS returns the handle to a new MODEL_PARAMETERS or the handle to
%      the existing singleton*.
%
%      MODEL_PARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL_PARAMETERS.M with the given input arguments.
%
%      MODEL_PARAMETERS('Property','Value',...) creates a new MODEL_PARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before model_parameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to model_parameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help model_parameters

% Last Modified by GUIDE v2.5 20-Jan-2014 16:12:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @model_parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @model_parameters_OutputFcn, ...
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


% --- Executes just before model_parameters is made visible.
function model_parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to model_parameters (see VARARGIN)

% Choose default command line output for model_parameters
handles.output = hObject;
load    './Temps/temporary.mat' 'working_name';
load (  './Temps/temp_variable.mat',...
        'num_species'                , 'Soil_nutrient'             , 'N_Uptake_RB'               , ...
        'ph_type1'                   , 'ph_type2'                  , 'ph_type3'                  , 'ph_type4'                  , ...
        'para_radiation'             , 'para_microenvironment'     , 'para_soil'                 , 'para_respiration'          , ...
        'para_soilCN'                ,  ...
        'para_soilCN_crop1'          , 'para_soilCN_crop2'         , 'para_soilCN_crop3'         ,'para_soilCN_crop4'          , ...
        'para_soilCN_crop11'         , 'para_soilCN_crop22'        , 'para_soilCN_crop33'        ,'para_soilCN_crop44'         , ...
        'para_leaf_crop1'            , 'para_leaf_crop2'           , 'para_leaf_crop3'           , 'para_leaf_crop4'           , ...
        'para_leaf_crop11'           , 'para_leaf_crop22'          , 'para_leaf_crop33'          , 'para_leaf_crop44'          , ...
        'para_canopy_crop1'          , 'para_canopy_crop2'         , 'para_canopy_crop3'         , 'para_canopy_crop4'         , ...
        'para_canopy_crop_fixed'     , ...
        'para_photosynthesisC3_crop1','para_photosynthesisC3_crop2','para_photosynthesisC3_crop3','para_photosynthesisC3_crop4', ...
        'para_photosynthesisC4_crop1','para_photosynthesisC4_crop2','para_photosynthesisC4_crop3','para_photosynthesisC4_crop4', ...
        'para_conductance_crop1'     ,'para_conductance_crop2'     ,'para_conductance_crop3'     ,'para_conductance_crop4');

if num_species == 1
    set(handles.parameter_but_species4,'Enable','on');
    set(handles.parameter_but_species2,'Enable','off');
    set(handles.parameter_but_species3,'Enable','off');
    set(handles.parameter_but_species4,'Enable','off');
elseif num_species == 2
    set(handles.parameter_but_species4,'Enable','on');
    set(handles.parameter_but_species2,'Enable','on');
    set(handles.parameter_but_species3,'Enable','off');
    set(handles.parameter_but_species4,'Enable','off');
elseif num_species == 3
    set(handles.parameter_but_species4,'Enable','on');
    set(handles.parameter_but_species2,'Enable','on');
    set(handles.parameter_but_species3,'Enable','on');
    set(handles.parameter_but_species4,'Enable','off');
elseif num_species == 4
    set(handles.parameter_but_species4,'Enable','on');
    set(handles.parameter_but_species2,'Enable','on');
    set(handles.parameter_but_species3,'Enable','on');
    set(handles.parameter_but_species4,'Enable','on');
end

set(handles.para_tab_soil,'Data',para_soil);
set(handles.para_tab_radiation,'Data',para_radiation);
%
set(handles.para_tab_respiration1,'Data',para_respiration);
set(handles.para_tab_microenvironment1,'Data',para_microenvironment);
%
set(handles.para_tab_leaf1,'Data',para_leaf_crop1);
set(handles.para_tab_leaf11,'Data',para_leaf_crop11);
set(handles.para_tab_canopy1,'Data',para_canopy_crop1);
set(handles.para_tab_canopy_fixed,'Data',para_canopy_crop_fixed);
set(handles.para_tab_leaf2,'Data',para_leaf_crop2);
set(handles.para_tab_leaf22,'Data',para_leaf_crop22);
set(handles.para_tab_canopy2,'Data',para_canopy_crop2);
set(handles.para_tab_leaf3,'Data',para_leaf_crop3);
set(handles.para_tab_leaf33,'Data',para_leaf_crop33);
set(handles.para_tab_canopy3,'Data',para_canopy_crop3);
set(handles.para_tab_leaf4,'Data',para_leaf_crop4);
set(handles.para_tab_leaf44,'Data',para_leaf_crop44);
set(handles.para_tab_canopy4,'Data',para_canopy_crop4);
%
set(handles.para_tab_conductance1,'Data',para_conductance_crop1);
set(handles.para_tab_conductance2,'Data',para_conductance_crop2);
set(handles.para_tab_conductance3,'Data',para_conductance_crop3);
set(handles.para_tab_conductance4,'Data',para_conductance_crop4);
%
if ph_type1 == 0
    set(handles.para_tab_photosynthesisC41,'Enable','on');
    set(handles.para_tab_photosynthesisC31,'Enable','off');
    set(handles.para_tab_photosynthesisC41,'Data',para_photosynthesisC4_crop1);
else
    set(handles.para_tab_photosynthesisC31,'Enable','on');
    set(handles.para_tab_photosynthesisC41,'Enable','off');
    set(handles.para_tab_photosynthesisC31,'Data',para_photosynthesisC3_crop1);
end
if ph_type2 == 0
    set(handles.para_tab_photosynthesisC42,'Enable','on');
    set(handles.para_tab_photosynthesisC32,'Enable','off');
    set(handles.para_tab_photosynthesisC42,'Data',para_photosynthesisC4_crop2);
else
    set(handles.para_tab_photosynthesisC32,'Enable','on');
    set(handles.para_tab_photosynthesisC42,'Enable','off');
    set(handles.para_tab_photosynthesisC32,'Data',para_photosynthesisC3_crop2);
end
if ph_type3 == 0
    set(handles.para_tab_photosynthesisC43,'Enable','on');
    set(handles.para_tab_photosynthesisC33,'Enable','off');
    set(handles.para_tab_photosynthesisC43,'Data',para_photosynthesisC4_crop3);
else
    set(handles.para_tab_photosynthesisC33,'Enable','on');
    set(handles.para_tab_photosynthesisC43,'Enable','off');
    set(handles.para_tab_photosynthesisC33,'Data',para_photosynthesisC3_crop3);
end 
if ph_type4 == 0
    set(handles.para_tab_photosynthesisC44,'Enable','on');
    set(handles.para_tab_photosynthesisC34,'Enable','off');
    set(handles.para_tab_photosynthesisC44,'Data',para_photosynthesisC4_crop4);
else
    set(handles.para_tab_photosynthesisC34,'Enable','on');
    set(handles.para_tab_photosynthesisC44,'Enable','off');
    set(handles.para_tab_photosynthesisC34,'Data',para_photosynthesisC3_crop4);
end

set(handles.para_tab_site_nutrient1,'Data',para_soilCN);
set(handles.para_tab_site_nutrient1,'Data',para_soilCN);
set(handles.para_tab_site_nutrient3,'Data',para_soilCN);
set(handles.para_tab_site_nutrient4,'Data',para_soilCN);


set(handles.para_tab_species_nutrient1,'Data',para_soilCN_crop1);
set(handles.para_tab_species_nutrient2,'Data',para_soilCN_crop2);
set(handles.para_tab_species_nutrient3,'Data',para_soilCN_crop3);
set(handles.para_tab_species_nutrient4,'Data',para_soilCN_crop4);

set(handles.para_tab_species_nutrient11,'Data',para_soilCN_crop11);
set(handles.para_tab_species_nutrient22,'Data',para_soilCN_crop22);
set(handles.para_tab_species_nutrient33,'Data',para_soilCN_crop33);
set(handles.para_tab_species_nutrient44,'Data',para_soilCN_crop44);

set(handles.par_but_Decom2,'Enable','off');
set(handles.par_but_Decom3,'Enable','off');
set(handles.par_but_Decom4,'Enable','off');
if Soil_nutrient
    set(handles.parameter_but_nutrient,'Enable','on');
    if N_Uptake_RB == 1
        set(handles.para_tab_species_nutrient1,'Enable','on');
        set(handles.para_tab_species_nutrient2,'Enable','on');
        set(handles.para_tab_species_nutrient3,'Enable','on');
        set(handles.para_tab_species_nutrient4,'Enable','on');
        
        set(handles.para_tab_species_nutrient11,'Enable','off');
        set(handles.para_tab_species_nutrient22,'Enable','off');
        set(handles.para_tab_species_nutrient33,'Enable','off');
        set(handles.para_tab_species_nutrient44,'Enable','off');
    elseif N_Uptake_RB ==0
        set(handles.para_tab_species_nutrient1,'Enable','off');
        set(handles.para_tab_species_nutrient2,'Enable','off');
        set(handles.para_tab_species_nutrient3,'Enable','off');
        set(handles.para_tab_species_nutrient4,'Enable','off');
        
        set(handles.para_tab_species_nutrient11,'Enable','on');
        set(handles.para_tab_species_nutrient22,'Enable','on');
        set(handles.para_tab_species_nutrient33,'Enable','on');
        set(handles.para_tab_species_nutrient44,'Enable','on');
    end
else
    set(handles.parameter_but_nutrient,'Enable','off');
end

set(handles.parameter_pan_soilradiation,'Visible','on');
set(handles.parameter_pan_respiration_microenvironment,'Visible','off');
set(handles.parameter_pan_leafcanopy1,'Visible','off');
set(handles.parameter_pan_leafcanopy2,'Visible','off');
set(handles.parameter_pan_leafcanopy3,'Visible','off');
set(handles.parameter_pan_leafcanopy4,'Visible','off');
set(handles.parameter_pan_conductance1,'Visible','off');
set(handles.parameter_pan_conductance2,'Visible','off');
set(handles.parameter_pan_conductance3,'Visible','off');
set(handles.parameter_pan_conductance4,'Visible','off');
set(handles.parameter_pan_photosynthesis1,'Visible','off');
set(handles.parameter_pan_photosynthesis2,'Visible','off');
set(handles.parameter_pan_photosynthesis3,'Visible','off');
set(handles.parameter_pan_photosynthesis4,'Visible','off');
set(handles.parameter_pan_nutrient1,'Visible','off');
set(handles.parameter_pan_nutrient2,'Visible','off');
set(handles.parameter_pan_nutrient3,'Visible','off');
set(handles.parameter_pan_nutrient4,'Visible','off');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes model_parameters wait for user response (see UIRESUME)
% uiwait(handles.model_parameters);


% --- Outputs from this function are returned to the command line.
function varargout = model_parameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in parameter_but_photosynthesis.
function tab_photosynthesis_Callback(hObject, eventdata, handles)
% hObject    handle to parameter_but_photosynthesis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parameter_but_photosynthesis
set(handles.parameter_pan_leafcanopy1,'Visible','off');
set(handles.pan2,'Visible','off');


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over parameter_but_leafcanopy.
function parameter_but_leafcanopy_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to parameter_but_leafcanopy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on parameter_but_leafcanopy and none of its controls.
function parameter_but_leafcanopy_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to parameter_but_leafcanopy (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on parameter_but_soilradiation and none of its controls.
function parameter_but_soilradiation_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to parameter_but_soilradiation (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

set(handles.parameter_pan_leafcanopy1,'Visible','off');
set(handles.parameter_pan_photosynthesis1,'Visible','off');
set(handles.parameter_pan_respiration_microenvironment,'Visible','off');
set(handles.parameter_pan_soilradiation,'Visible','on');

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over parameter_but_soilradiation.


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes when selected cell(s) is changed in para_tab_canopy1.
function para_tab_canopy1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to para_tab_canopy1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function parameter_but_soilradiation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameter_but_soilradiation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function parameter_but_leafcanopy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameter_but_leafcanopy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function parameter_but_photosynthesis_CreateFcn(hObject, eventdata, handles)


function parameter_but_respiration_mircoenvironment_CreateFcn(hObject, eventdata, handles)


function parameter_but_photosynthesis_DeleteFcn(hObject, eventdata, handles)


function parameter_but_OK_Callback(hObject, eventdata, handles)
para_soil                   = get(handles.para_tab_soil,'Data');
para_radiation              = get(handles.para_tab_radiation,'Data');
para_respiration            = get(handles.para_tab_respiration1,'Data');
para_microenvironment       = get(handles.para_tab_microenvironment1,'Data');
para_leaf_crop1             = get(handles.para_tab_leaf1,'Data');
para_leaf_crop11            = get(handles.para_tab_leaf11,'Data');
para_leaf_crop2             = get(handles.para_tab_leaf2,'Data');
para_leaf_crop22            = get(handles.para_tab_leaf22,'Data');
para_leaf_crop3             = get(handles.para_tab_leaf3,'Data');
para_leaf_crop33            = get(handles.para_tab_leaf33,'Data');
para_leaf_crop4             = get(handles.para_tab_leaf4,'Data');
para_leaf_crop44            = get(handles.para_tab_leaf44,'Data');
para_canopy_crop1           = get(handles.para_tab_canopy1,'Data');
para_canopy_crop_fixed      =  get(handles.para_tab_canopy_fixed,'Data');
para_canopy_crop2           = get(handles.para_tab_canopy2,'Data');
para_canopy_crop3           = get(handles.para_tab_canopy3,'Data');
para_canopy_crop4           = get(handles.para_tab_canopy4,'Data');
para_conductance_crop1      = get(handles.para_tab_conductance1,'Data');
para_conductance_crop2      = get(handles.para_tab_conductance2,'Data');
para_conductance_crop3      = get(handles.para_tab_conductance3,'Data');
para_conductance_crop4      = get(handles.para_tab_conductance4,'Data');
para_photosynthesisC3_crop1 = get(handles.para_tab_photosynthesisC31,'Data');
para_photosynthesisC3_crop2 = get(handles.para_tab_photosynthesisC32,'Data');
para_photosynthesisC3_crop3 = get(handles.para_tab_photosynthesisC33,'Data');
para_photosynthesisC3_crop4 = get(handles.para_tab_photosynthesisC34,'Data');
para_photosynthesisC4_crop1 = get(handles.para_tab_photosynthesisC41,'Data');
para_photosynthesisC4_crop2 = get(handles.para_tab_photosynthesisC42,'Data');
para_photosynthesisC4_crop3 = get(handles.para_tab_photosynthesisC43,'Data');
para_photosynthesisC4_crop4 = get(handles.para_tab_photosynthesisC44,'Data');
para_soilCN                 = get(handles.para_tab_site_nutrient1,'Data');
para_soilCN_crop1           = get(handles.para_tab_species_nutrient1,'Data');
para_soilCN_crop2           = get(handles.para_tab_species_nutrient2,'Data');
para_soilCN_crop3           = get(handles.para_tab_species_nutrient3,'Data');
para_soilCN_crop4           = get(handles.para_tab_species_nutrient4,'Data');
para_soilCN_crop11          = get(handles.para_tab_species_nutrient11,'Data');
para_soilCN_crop22          = get(handles.para_tab_species_nutrient22,'Data');
para_soilCN_crop33          = get(handles.para_tab_species_nutrient33,'Data');
para_soilCN_crop44          = get(handles.para_tab_species_nutrient44,'Data');

save './Temps/temp_variable.mat'...
    'para_soil' 'para_radiation' 'para_respiration' 'para_microenvironment'...
    'para_leaf_crop1' 'para_leaf_crop2' 'para_leaf_crop3' 'para_leaf_crop4' ...
    'para_leaf_crop11' 'para_leaf_crop22' 'para_leaf_crop33' 'para_leaf_crop44' ...
    'para_canopy_crop1' 'para_canopy_crop2' 'para_canopy_crop3' 'para_canopy_crop4' ...
    'para_canopy_crop_fixed'  ...
    'para_conductance_crop1' 'para_conductance_crop2' 'para_conductance_crop3' 'para_conductance_crop4' ...
    'para_photosynthesisC3_crop1' 'para_photosynthesisC3_crop2' 'para_photosynthesisC3_crop3' 'para_photosynthesisC3_crop4' ...
    'para_photosynthesisC4_crop1' 'para_photosynthesisC4_crop2' 'para_photosynthesisC4_crop3' 'para_photosynthesisC4_crop4' ...
    'para_soilCN' ...
    'para_soilCN_crop1' 'para_soilCN_crop2' 'para_soilCN_crop3' 'para_soilCN_crop4' ...
    'para_soilCN_crop11' 'para_soilCN_crop22' 'para_soilCN_crop33' 'para_soilCN_crop44' -append;
    
load ('./Temps/temporary.mat',...
    'working_name'  )
copyfile('./Temps/temp_variable.mat',working_name,'f');

load (  './Temps/temp_variable.mat',...
        'num_species', 'Soil_nutrient','ph_type1','ph_type2','ph_type3','ph_type4' ...
        );

if ph_type1 == 1
    isthere_photo1=para_photosynthesisC3_crop1;
else
    isthere_photo1=para_photosynthesisC4_crop1;
end
if ph_type2 == 1
    isthere_photo2=para_photosynthesisC3_crop2;
else
    isthere_photo2=para_photosynthesisC4_crop2;
end
if ph_type3 == 1
    isthere_photo3=para_photosynthesisC3_crop3;
else
    isthere_photo3=para_photosynthesisC4_crop3;
end
if ph_type4 == 1
    isthere_photo4=para_photosynthesisC3_crop4;
else
    isthere_photo4=para_photosynthesisC4_crop4;
end
if Soil_nutrient
    isthere_site_nutrient  = para_soilCN;
    isthere_crop_nutrient1 = para_soilCN_crop1;
    isthere_crop_nutrient2 = para_soilCN_crop2;
    isthere_crop_nutrient3 = para_soilCN_crop3;
    isthere_crop_nutrient4 = para_soilCN_crop4;
else
    isthere_site_nutrient  = 1;
    isthere_crop_nutrient1 = 1;
    isthere_crop_nutrient2 = 1;
    isthere_crop_nutrient3 = 1;
    isthere_crop_nutrient4 = 1;
end

load (  './Temps/temp_variable.mat',...
        'dat_decom_litter','dat_decom');
if Soil_nutrient
    if (sum(sum(dat_decom)) == 0 || sum(sum(dat_decom_litter)) == 0) 
        msgbox( 'Please go back to NUTRIENT & Species1 - Click on Add/Edit Input decomposition parameters',   ...
            'MLCan Error','Error');
        return
    end
end


if num_species == 1
    if (isempty(para_soil) || isempty(para_radiation) || isempty(para_respiration)|| isempty(para_microenvironment)||...
            isempty(para_leaf_crop1)|| isempty(para_canopy_crop1)|| isempty(para_conductance_crop1)|| isempty(isthere_photo1) || isempty(isthere_site_nutrient) || isempty(isthere_crop_nutrient1))
        msgbox('Please enter enough information','Error');
        return
    else
        close;
    end
    
elseif num_species == 2
    if (isempty(para_soil) || isempty(para_radiation) || isempty(para_respiration)|| isempty(para_microenvironment)||...
            isempty(para_leaf_crop1)|| isempty(para_canopy_crop1)|| isempty(para_conductance_crop1)|| isempty(isthere_photo1) || isempty(isthere_site_nutrient) || isempty(isthere_crop_nutrient1) || ...
            isempty(para_leaf_crop2)|| isempty(para_canopy_crop2)|| isempty(para_conductance_crop2)|| isempty(isthere_photo2) || isempty(isthere_crop_nutrient2))
        msgbox('Please enter enough information','Error');
        return
    else
        close;
    end
elseif num_species == 3
    if (isempty(para_soil) || isempty(para_radiation) || isempty(para_respiration)|| isempty(para_microenvironment)||...
            isempty(para_leaf_crop1)|| isempty(para_canopy_crop1)|| isempty(para_conductance_crop1)|| isempty(para_photosynthesisC3_crop1) || isempty(para_photosynthesisC4_crop1) || ...
            isempty(para_leaf_crop2)|| isempty(para_canopy_crop2)|| isempty(para_conductance_crop2)|| isempty(isthere_photo2) || isempty(isthere_crop_nutrient2) || ...
            isempty(para_leaf_crop3)|| isempty(para_canopy_crop3)|| isempty(para_conductance_crop3)|| isempty(isthere_photo3) || isempty(isthere_crop_nutrient3))
        msgbox('Please enter enough information','Error');
        return
    else
        close;
    end
elseif num_species == 4
    if (isempty(para_soil) || isempty(para_radiation) || isempty(para_respiration)|| isempty(para_microenvironment)||...
            isempty(para_leaf_crop1)|| isempty(para_canopy_crop1)|| isempty(para_conductance_crop1)|| isempty(para_photosynthesisC3_crop1) || isempty(para_photosynthesisC4_crop1) || ...
            isempty(para_leaf_crop2)|| isempty(para_canopy_crop2)|| isempty(para_conductance_crop2)|| isempty(isthere_photo2) || isempty(isthere_crop_nutrient2) || ...
            isempty(para_leaf_crop3)|| isempty(para_canopy_crop3)|| isempty(para_conductance_crop3)|| isempty(isthere_photo3) || isempty(isthere_crop_nutrient3) || ...
            isempty(para_leaf_crop4)|| isempty(para_canopy_crop4)|| isempty(para_conductance_crop4)|| isempty(isthere_photo4) || isempty(isthere_crop_nutrient4))
        msgbox('Please enter enough information','Error');
        return
    else
        close;
    end
end




function parameter_but_cancel_Callback(hObject, eventdata, handles)
close

function parameter_but_photosynthesis_ButtonDownFcn(hObject, eventdata, handles)


function axes1_CreateFcn(hObject, eventdata, handles)
paraiconsmall = imread('./users/icons/parameter_small.png'); 
image(paraiconsmall);
axis off
%imshow('./users/icons/parameter_small.png'); 


function para_tab_microenvironment1_CellEditCallback(hObject, eventdata, handles)


function para_tab_conductance_CellSelectionCallback(hObject, eventdata, handles)


function para_tab_respiration1_CellSelectionCallback(hObject, eventdata, handles)


function para_tab_photosynthesisC41_CellEditCallback(hObject, eventdata, handles)


function para_tab_soil_CellSelectionCallback(hObject, eventdata, handles)


function para_tab_radiation_CellEditCallback(hObject, eventdata, handles)


function para_tab_leaf1_CellEditCallback(hObject, eventdata, handles)


function para_tab_canopy1_CellEditCallback(hObject, eventdata, handles)


function parameter_pan_selctions_SelectionChangeFcn(hObject, eventdata, handles)
load ('./Temps/temp_variable.mat', 'current_species');
        
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'parameter_but_soilradiation'
        set(handles.parameter_pan_soilradiation,'Visible','on');
        set(handles.parameter_pan_respiration_microenvironment,'Visible','off');
        set(handles.parameter_pan_leafcanopy1,'Visible','off');
        set(handles.parameter_pan_leafcanopy2,'Visible','off');
        set(handles.parameter_pan_leafcanopy3,'Visible','off');
        set(handles.parameter_pan_leafcanopy4,'Visible','off');
        set(handles.parameter_pan_conductance1,'Visible','off');
        set(handles.parameter_pan_conductance2,'Visible','off');
        set(handles.parameter_pan_conductance3,'Visible','off');
        set(handles.parameter_pan_conductance4,'Visible','off');
        set(handles.parameter_pan_photosynthesis1,'Visible','off');
        set(handles.parameter_pan_photosynthesis2,'Visible','off');
        set(handles.parameter_pan_photosynthesis3,'Visible','off');
        set(handles.parameter_pan_photosynthesis4,'Visible','off');
        set(handles.parameter_pan_nutrient1,'Visible','off');
        set(handles.parameter_pan_nutrient2,'Visible','off');
        set(handles.parameter_pan_nutrient3,'Visible','off');
        set(handles.parameter_pan_nutrient4,'Visible','off');
    case 'parameter_but_respiration_mircoenvironment'
        set(handles.parameter_pan_soilradiation,'Visible','off');
        set(handles.parameter_pan_respiration_microenvironment,'Visible','on');
        set(handles.parameter_pan_leafcanopy1,'Visible','off');
        set(handles.parameter_pan_leafcanopy2,'Visible','off');
        set(handles.parameter_pan_leafcanopy3,'Visible','off');
        set(handles.parameter_pan_leafcanopy4,'Visible','off');
        set(handles.parameter_pan_conductance1,'Visible','off');
        set(handles.parameter_pan_conductance2,'Visible','off');
        set(handles.parameter_pan_conductance3,'Visible','off');
        set(handles.parameter_pan_conductance4,'Visible','off');
        set(handles.parameter_pan_photosynthesis1,'Visible','off');
        set(handles.parameter_pan_photosynthesis2,'Visible','off');
        set(handles.parameter_pan_photosynthesis3,'Visible','off');
        set(handles.parameter_pan_photosynthesis4,'Visible','off');
        set(handles.parameter_pan_nutrient1,'Visible','off');
        set(handles.parameter_pan_nutrient2,'Visible','off');
        set(handles.parameter_pan_nutrient3,'Visible','off');
        set(handles.parameter_pan_nutrient4,'Visible','off');
    case 'parameter_but_leafcanopy'
        set(handles.parameter_pan_soilradiation,'Visible','off');
        set(handles.parameter_pan_respiration_microenvironment,'Visible','off');
        if current_species == 1
            set(handles.parameter_pan_leafcanopy1,'Visible','on');
            set(handles.parameter_pan_leafcanopy2,'Visible','off');
            set(handles.parameter_pan_leafcanopy3,'Visible','off');
            set(handles.parameter_pan_leafcanopy4,'Visible','off');
        elseif current_species == 2
            set(handles.parameter_pan_leafcanopy1,'Visible','off');
            set(handles.parameter_pan_leafcanopy2,'Visible','on');
            set(handles.parameter_pan_leafcanopy3,'Visible','off');
            set(handles.parameter_pan_leafcanopy4,'Visible','off');
        elseif current_species == 3
            set(handles.parameter_pan_leafcanopy1,'Visible','off');
            set(handles.parameter_pan_leafcanopy2,'Visible','off');
            set(handles.parameter_pan_leafcanopy3,'Visible','on');
            set(handles.parameter_pan_leafcanopy4,'Visible','off');
        elseif current_species == 4
            set(handles.parameter_pan_leafcanopy1,'Visible','off');
            set(handles.parameter_pan_leafcanopy2,'Visible','off');
            set(handles.parameter_pan_leafcanopy3,'Visible','off');
            set(handles.parameter_pan_leafcanopy4,'Visible','on');
        end
        set(handles.parameter_pan_conductance1,'Visible','off');
        set(handles.parameter_pan_conductance2,'Visible','off');
        set(handles.parameter_pan_conductance3,'Visible','off');
        set(handles.parameter_pan_conductance4,'Visible','off');
        set(handles.parameter_pan_photosynthesis1,'Visible','off');
        set(handles.parameter_pan_photosynthesis2,'Visible','off');
        set(handles.parameter_pan_photosynthesis3,'Visible','off');
        set(handles.parameter_pan_photosynthesis4,'Visible','off');
        set(handles.parameter_pan_nutrient1,'Visible','off');
        set(handles.parameter_pan_nutrient2,'Visible','off');
        set(handles.parameter_pan_nutrient3,'Visible','off');
        set(handles.parameter_pan_nutrient4,'Visible','off');
    case 'parameter_but_conductance'
        set(handles.parameter_pan_soilradiation,'Visible','off');
        set(handles.parameter_pan_respiration_microenvironment,'Visible','off');
        set(handles.parameter_pan_leafcanopy1,'Visible','off');
        set(handles.parameter_pan_leafcanopy2,'Visible','off');
        set(handles.parameter_pan_leafcanopy3,'Visible','off');
        set(handles.parameter_pan_leafcanopy4,'Visible','off');
        if current_species == 1
            set(handles.parameter_pan_conductance1,'Visible','on');
            set(handles.parameter_pan_conductance2,'Visible','off');
            set(handles.parameter_pan_conductance3,'Visible','off');
            set(handles.parameter_pan_conductance4,'Visible','off');
        elseif current_species == 2
            set(handles.parameter_pan_conductance1,'Visible','off');
            set(handles.parameter_pan_conductance2,'Visible','on');
            set(handles.parameter_pan_conductance3,'Visible','off');
            set(handles.parameter_pan_conductance4,'Visible','off');
        elseif current_species == 3
            set(handles.parameter_pan_conductance1,'Visible','off');
            set(handles.parameter_pan_conductance2,'Visible','off');
            set(handles.parameter_pan_conductance3,'Visible','on');
            set(handles.parameter_pan_conductance4,'Visible','off');
        elseif current_species == 4
            set(handles.parameter_pan_conductance1,'Visible','off');
            set(handles.parameter_pan_conductance2,'Visible','off');
            set(handles.parameter_pan_conductance3,'Visible','off');
            set(handles.parameter_pan_conductance4,'Visible','on');
        end
        set(handles.parameter_pan_photosynthesis1,'Visible','off');
        set(handles.parameter_pan_photosynthesis2,'Visible','off');
        set(handles.parameter_pan_photosynthesis3,'Visible','off');
        set(handles.parameter_pan_photosynthesis4,'Visible','off');
        set(handles.parameter_pan_nutrient1,'Visible','off');
        set(handles.parameter_pan_nutrient2,'Visible','off');
        set(handles.parameter_pan_nutrient3,'Visible','off');
        set(handles.parameter_pan_nutrient4,'Visible','off');
    case 'parameter_but_photosynthesis'
        set(handles.parameter_pan_soilradiation,'Visible','off');
        set(handles.parameter_pan_respiration_microenvironment,'Visible','off');
        set(handles.parameter_pan_leafcanopy1,'Visible','off');
        set(handles.parameter_pan_leafcanopy2,'Visible','off');
        set(handles.parameter_pan_leafcanopy3,'Visible','off');
        set(handles.parameter_pan_leafcanopy4,'Visible','off');
        set(handles.parameter_pan_conductance1,'Visible','off');
        set(handles.parameter_pan_conductance2,'Visible','off');
        set(handles.parameter_pan_conductance3,'Visible','off');
        set(handles.parameter_pan_conductance4,'Visible','off');
        if current_species == 1
            set(handles.parameter_pan_photosynthesis1,'Visible','on');
            set(handles.parameter_pan_photosynthesis2,'Visible','off');
            set(handles.parameter_pan_photosynthesis3,'Visible','off');
            set(handles.parameter_pan_photosynthesis4,'Visible','off');
        elseif current_species == 2
            set(handles.parameter_pan_photosynthesis1,'Visible','off');
            set(handles.parameter_pan_photosynthesis2,'Visible','on');
            set(handles.parameter_pan_photosynthesis3,'Visible','off');
            set(handles.parameter_pan_photosynthesis4,'Visible','off');
        elseif current_species == 3
            set(handles.parameter_pan_photosynthesis1,'Visible','off');
            set(handles.parameter_pan_photosynthesis2,'Visible','off');
            set(handles.parameter_pan_photosynthesis3,'Visible','on');
            set(handles.parameter_pan_photosynthesis4,'Visible','off');
        elseif current_species == 4
            set(handles.parameter_pan_photosynthesis1,'Visible','off');
            set(handles.parameter_pan_photosynthesis2,'Visible','off');
            set(handles.parameter_pan_photosynthesis3,'Visible','off');
            set(handles.parameter_pan_photosynthesis4,'Visible','on');
        end
        set(handles.parameter_pan_nutrient1,'Visible','off');
        set(handles.parameter_pan_nutrient2,'Visible','off');
        set(handles.parameter_pan_nutrient3,'Visible','off');
        set(handles.parameter_pan_nutrient4,'Visible','off');
    case 'parameter_but_nutrient'
        set(handles.parameter_pan_soilradiation,'Visible','off');
        set(handles.parameter_pan_respiration_microenvironment,'Visible','off');
        set(handles.parameter_pan_leafcanopy1,'Visible','off');
        set(handles.parameter_pan_leafcanopy2,'Visible','off');
        set(handles.parameter_pan_leafcanopy3,'Visible','off');
        set(handles.parameter_pan_leafcanopy4,'Visible','off');
        set(handles.parameter_pan_conductance1,'Visible','off');
        set(handles.parameter_pan_conductance2,'Visible','off');
        set(handles.parameter_pan_conductance3,'Visible','off');
        set(handles.parameter_pan_conductance4,'Visible','off');
        set(handles.parameter_pan_photosynthesis1,'Visible','off');
        set(handles.parameter_pan_photosynthesis2,'Visible','off');
        set(handles.parameter_pan_photosynthesis3,'Visible','off');
        set(handles.parameter_pan_photosynthesis4,'Visible','off');
        if current_species == 1
            set(handles.parameter_pan_nutrient1,'Visible','on');
            set(handles.parameter_pan_nutrient2,'Visible','off');
            set(handles.parameter_pan_nutrient3,'Visible','off');
            set(handles.parameter_pan_nutrient4,'Visible','off');
            
            set(handles.para_tab_site_nutrient1,'Enable','on');
        elseif current_species == 2
            set(handles.parameter_pan_nutrient1,'Visible','off');
            set(handles.parameter_pan_nutrient2,'Visible','on');
            set(handles.parameter_pan_nutrient3,'Visible','off');
            set(handles.parameter_pan_nutrient4,'Visible','off');
            
            para_soilCN=get(handles.para_tab_site_nutrient1,'Data');
            set(handles.para_tab_site_nutrient2,'Data',para_soilCN);
            set(handles.para_tab_site_nutrient2,'Enable','off');
        elseif current_species == 3
            set(handles.parameter_pan_nutrient1,'Visible','off');
            set(handles.parameter_pan_nutrient2,'Visible','off');
            set(handles.parameter_pan_nutrient3,'Visible','on');
            set(handles.parameter_pan_nutrient4,'Visible','off');
            
            para_soilCN=get(handles.para_tab_site_nutrient1,'Data');
            set(handles.para_tab_site_nutrient3,'Data',para_soilCN);
            set(handles.para_tab_site_nutrient3,'Enable','off');
        elseif current_species == 4
            set(handles.parameter_pan_nutrient1,'Visible','off');
            set(handles.parameter_pan_nutrient2,'Visible','off');
            set(handles.parameter_pan_nutrient3,'Visible','off');
            set(handles.parameter_pan_nutrient4,'Visible','on');
            
            para_soilCN=get(handles.para_tab_site_nutrient1,'Data');
            set(handles.para_tab_site_nutrient4,'Data',para_soilCN);
            set(handles.para_tab_site_nutrient4,'Enable','off');
        end
end

% --- Executes on button press in togglebutton9.
function togglebutton9_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton9


% --- Executes on button press in parameter_but_species4.
function parameter_but_species1_Callback(hObject, eventdata, handles)
% hObject    handle to parameter_but_species4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parameter_but_species4

% --- Executes on button press in parameter_but_species2.
function parameter_but_species2_Callback(hObject, eventdata, handles)
% hObject    handle to parameter_but_species2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parameter_but_species2


% --- Executes on button press in parameter_but_species3.
function parameter_but_species3_Callback(hObject, eventdata, handles)
% hObject    handle to parameter_but_species3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parameter_but_species3


% --- Executes on button press in parameter_but_species4.
function parameter_but_species4_Callback(hObject, eventdata, handles)
% hObject    handle to parameter_but_species4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parameter_but_species4


% --- Executes when selected object is changed in parameter_pan_species.
function parameter_pan_species_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in parameter_pan_species 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'parameter_but_species1'
        set(handles.parameter_but_soilradiation,'Enable','on');
        set(handles.parameter_but_respiration_mircoenvironment,'Enable','on');
        current_species=1;
        save './Temps/temp_variable.mat'...
            'current_species' -append;
    case 'parameter_but_species2'
        set(handles.parameter_but_soilradiation,'Enable','off');
        set(handles.parameter_but_respiration_mircoenvironment,'Enable','off');
        current_species=2;
        save './Temps/temp_variable.mat'...
            'current_species' -append;
    case 'parameter_but_species3'
        set(handles.parameter_but_soilradiation,'Enable','off');
        set(handles.parameter_but_respiration_mircoenvironment,'Enable','off');
        current_species=3;
        save './Temps/temp_variable.mat'...
            'current_species' -append;
    case 'parameter_but_species4'
        set(handles.parameter_but_soilradiation,'Enable','off');
        set(handles.parameter_but_respiration_mircoenvironment,'Enable','off');
        current_species=4;
        save './Temps/temp_variable.mat'...
            'current_species' -append;
end


% --- Executes on button press in par_but_Decom3.
function par_but_Decom3_Callback(hObject, eventdata, handles)
% hObject    handle to par_but_Decom3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in par_but_Decom2.
function par_but_Decom2_Callback(hObject, eventdata, handles)
% hObject    handle to par_but_Decom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in par_but_Decom4.
function par_but_Decom4_Callback(hObject, eventdata, handles)
% hObject    handle to par_but_Decom4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in par_but_Decom1.
function par_but_Decom1_Callback(hObject, eventdata, handles)
% hObject    handle to par_but_Decom1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model_parameters_decom;
