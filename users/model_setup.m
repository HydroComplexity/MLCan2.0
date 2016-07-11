function varargout = model_setup(varargin)
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                        FUNCTION CODE INTERFACE                        %%
%%                              MODEL SETUP                              %%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
% This interface is used to call an interface for setting up and entering % 
% crop, geographic, number of layers information that will be assigned    %
% to the MLCan model                                                      %
%                                                                         %
%-------------------------------------------------------------------------%
% MODEL_SETUP M-file for model_setup.fig                                  %
%-------------------------------------------------------------------------%
%   Created by      : Phong Le                                            %
%   Date            : May 21, 2010                                        %
%   Last Modified   : May 23, 2010                                        %
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%% Begin initialization code - DO NOT EDIT                               %%
%-------------------------------------------------------------------------%
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @model_setup_OpeningFcn, ...
                   'gui_OutputFcn',  @model_setup_OutputFcn, ...
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
%-------------------------------------------------------------------------%
%% End initialization code - DO NOT EDIT                                 %%
%-------------------------------------------------------------------------%

% --- Executes just before model_setup is made visible.
function model_setup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to model_setup (see VARARGIN)

% Choose default command line output for model_setup
handles.output = hObject;

load  './Temps/temporary.mat' 'working_name';
load ('./Temps/temp_variable.mat',...
        'num_species'                                                   ,   ...
        'lat_face'      , 'long_face'   , 'elev_face'   ,  ...
        'DOY_start'     , 'DOY_end'     , ...
        'crop_name1'    , 'crop_name2'  , 'crop_name3'  , 'crop_name4'  ,   ...
        'ph_type1'      , 'ph_type2'    , 'ph_type3'    , 'ph_type4'    ,   ...
        'num_LAD1'      , 'num_LAD2'    , 'num_LAD3'    , 'num_LAD4'    ,   ...
        'num_root1'     , 'num_root2'   , 'num_root3'   , 'num_root4'   ,   ...
        'opt_root1'     , 'opt_root2'   , 'opt_root3'   , 'opt_root4'   );
%     'LAImin_face'   , ...
%         'num_can1'      , 'num_can2'    , 'num_can3'    , 'num_can4'    ,   ...

set(handles.edit_num_species,'String',num_species);
set(handles.txt_lat,'String',lat_face);
set(handles.txt_long,'String',long_face);
set(handles.txt_elevation,'String',elev_face);
set(handles.txt_specified_DOY_start,'String',DOY_start);
set(handles.txt_specified_DOY_end,'String',DOY_end);
% set(handles.txt_specified_DOY_start,'String',LAImin_face);

set(handles.txt_crop_name1,'String',crop_name1);
set(handles.txt_crop_name2,'String',crop_name2);
set(handles.txt_crop_name3,'String',crop_name3);
set(handles.txt_crop_name4,'String',crop_name4);

set(handles.radio_c3type1,'Value',ph_type1);
set(handles.radio_c3type2,'Value',ph_type2);
set(handles.radio_c3type3,'Value',ph_type3);
set(handles.radio_c3type4,'Value',ph_type4);
if ph_type1 == 0
    set(handles.radio_c4type1,'Value',1);
end
if ph_type2 == 0
    set(handles.radio_c4type2,'Value',1);
end
if ph_type3 == 0
    set(handles.radio_c4type3,'Value',1);
end
if ph_type4 == 0
    set(handles.radio_c4type4,'Value',1);
end

%set(handles.txt_canopy_layer1,'String',num_can1);
%set(handles.txt_canopy_layer2,'String',num_can2);
%set(handles.txt_canopy_layer3,'String',num_can3);
%set(handles.txt_canopy_layer4,'String',num_can4);

set(handles.txt_LAD_layer1,'String',num_LAD1);
set(handles.txt_LAD_layer2,'String',num_LAD2);
set(handles.txt_LAD_layer3,'String',num_LAD3);
set(handles.txt_LAD_layer4,'String',num_LAD4);

set(handles.txt_root_layer1,'String',num_root1);
set(handles.txt_root_layer2,'String',num_root2);
set(handles.txt_root_layer3,'String',num_root3);
set(handles.txt_root_layer4,'String',num_root4);

set(handles.model_setup,'Color',[0.941 0.941 0.941]);

if num_species == 1
    set(handles.push_species1,'enable','on');
    set(handles.push_species2,'enable','off');
    set(handles.push_species3,'enable','off');
    set(handles.push_species4,'enable','off');
end
if num_species == 2
    set(handles.push_species1,'enable','on');
    set(handles.push_species2,'enable','on');
    set(handles.push_species3,'enable','off');
    set(handles.push_species4,'enable','off');
end
if num_species == 3
    set(handles.push_species1,'enable','on');
    set(handles.push_species2,'enable','on');
    set(handles.push_species3,'enable','on');
    set(handles.push_species4,'enable','off');
end
if num_species == 4
    set(handles.push_species1,'enable','on');
    set(handles.push_species2,'enable','on');
    set(handles.push_species3,'enable','on');
    set(handles.push_species4,'enable','on');
end
%if opt_root1 == 1
%    set(handles.txt_root_layer1,'String',num_root1);
%else
%    set(handles.txt_root_layer1,'String',num_root_eqn1);
%end

% Update handles structure
guidata(hObject, handles);
%-------------------------------------------------------------------------%
% UIWAIT makes main_MLCan wait for user response (see UIRESUME)           %
% uiwait(handles.Main);                                                   %
%-------------------------------------------------------------------------%

% Dongkook: Fixing 2015a bugs
function uipanel2_ResizeFcn(hObject, eventdata, handles)


% --- Outputs from this function are returned to the command line.
function varargout = model_setup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in set_mod_ok.
function set_mod_ok_Callback(hObject, eventdata, handles)
%
%load './Temps/temp_variable.mat' 'canopy_systems' 'root_systems'
%
% Read string data
% num_species  = get(handles.edit_num_species,'String');
% lat_face     = get(handles.txt_lat,'String');
% long_face    = get(handles.txt_long,'String');
% elev_face    = get(handles.txt_elevation,'String');
% LAImin_face  = get(handles.txt_specified_DOY_start,'String');
% 
% crop_name1   = get(handles.txt_crop_name1,'String');
% crop_name2   = get(handles.txt_crop_name2,'String');
% crop_name3   = get(handles.txt_crop_name3,'String');
% crop_name4   = get(handles.txt_crop_name4,'String');
% 
% ph_type1     = get(handles.radio_c3type1,'Value');
% ph_type2     = get(handles.radio_c3type2,'Value');
% ph_type3     = get(handles.radio_c3type3,'Value');
% ph_type4     = get(handles.radio_c3type4,'Value');
% 
% num_can1     = get(handles.txt_canopy_layer1,'String');
% num_can2     = get(handles.txt_canopy_layer2,'String');
% num_can3     = get(handles.txt_canopy_layer3,'String');
% num_can4     = get(handles.txt_canopy_layer4,'String');
% 
% num_LAD1     = get(handles.txt_LAD_layer1,'String');
% num_LAD2     = get(handles.txt_LAD_layer2,'String');
% num_LAD3     = get(handles.txt_LAD_layer3,'String');
% num_LAD4     = get(handles.txt_LAD_layer4,'String');
% 
% num_root1    = get(handles.txt_root_layer1,'String');
% num_root2    = get(handles.txt_root_layer2,'String');
% num_root3    = get(handles.txt_root_layer3,'String');
% num_root4    = get(handles.txt_root_layer4,'String');

num_species_str = get(handles.edit_num_species,'String');
lat_str         = get(handles.txt_lat,'String');
long_str        = get(handles.txt_long,'String');
elev_str        = get(handles.txt_elevation,'String');
DOY_start       = get(handles.txt_specified_DOY_start,'String');
DOY_end         = get(handles.txt_specified_DOY_end,'String');
% LAImin_str     = get(handles.txt_specified_DOY_start,'String');

crop_name1      = get(handles.txt_crop_name1,'String');
crop_name2      = get(handles.txt_crop_name2,'String');
crop_name3      = get(handles.txt_crop_name3,'String');
crop_name4      = get(handles.txt_crop_name4,'String');

ph_type1        = get(handles.radio_c3type1,'Value');
ph_type2        = get(handles.radio_c3type2,'Value');
ph_type3        = get(handles.radio_c3type3,'Value');
ph_type4        = get(handles.radio_c3type4,'Value');

%num_can_str1    = get(handles.txt_canopy_layer1,'String');
%num_can_str2    = get(handles.txt_canopy_layer2,'String');
%num_can_str3    = get(handles.txt_canopy_layer3,'String');
%num_can_str4    = get(handles.txt_canopy_layer4,'String');

num_LAD_str1    = get(handles.txt_LAD_layer1,'String');
num_LAD_str2    = get(handles.txt_LAD_layer2,'String');
num_LAD_str3    = get(handles.txt_LAD_layer3,'String');
num_LAD_str4    = get(handles.txt_LAD_layer4,'String');

num_root_str1   = get(handles.txt_root_layer1,'String');
num_root_str2   = get(handles.txt_root_layer2,'String');
num_root_str3   = get(handles.txt_root_layer3,'String');
num_root_str4   = get(handles.txt_root_layer4,'String');
%
% % Convert string to number;
num_species     = str2num(num_species_str);
lat_face      	= str2num(lat_str);
long_face     	= str2num(long_str);
elev_face      	= str2num(elev_str);
DOY_start       = str2num(DOY_start);
DOY_end         = str2num(DOY_end);
% LAImin_face     = str2num(LAImin_str);

%num_can1        = str2num(num_can_str1);
%num_can2        = str2num(num_can_str2);
%num_can3        = str2num(num_can_str3);
%num_can4        = str2num(num_can_str4);

num_LAD1        = str2num(num_LAD_str1);
num_LAD2        = str2num(num_LAD_str2);
num_LAD3        = str2num(num_LAD_str3);
num_LAD4        = str2num(num_LAD_str4);

num_root1       = str2num(num_root_str1);
num_root2       = str2num(num_root_str2);
num_root3       = str2num(num_root_str3);
num_root4       = str2num(num_root_str4);
%
% Check empty entering information
if (isempty(lat_str) || isempty(long_str)|| ...
        isempty(num_LAD_str1)|| isempty(num_root_str1))
    msgbox('Please enter enough information ','MLCan Error: Empty information','error');
    % isempty(num_can_str1)||
    return
else % Check numeric entering information
    if (isempty(lat_face)   == 1 || isempty(long_face)  == 1 || isempty(elev_face) == 1 ||...
        isempty(num_LAD1)    == 1 || isempty(num_root1)   == 1)
    % isempty(num_can1)    == 1 || 
            msgbox('Numbers must be in numeric ','MLCan Error: Data type','error');
            return
    else
        if (mod(num_LAD1,1) ~= 0 || mod(num_root1,1) ~= 0 ||...
            num_LAD1 <= 0 || num_root1 <= 0)
            msgbox ('Layer numbers must be positive integer (> 0) ','MLCan Error: Data type','error');
            % mod(num_can1,1) ~= 0 ||
            % num_can1 <= 0 || 
            return
        else
            if (lat_face < -90 || lat_face > 90 || long_face < -180 || long_face > 180)
                msgbox ('Latitude must be in [-90 90] degree and Longitude in [-180 180] degree', 'MLCan Error: Geographic','error');
                return;
            else
                if (num_species == 1 && isempty(crop_name1))
                    msgbox('Please enter enough information ','MLCan Error: Empty information','error');
                    return
                else
                    if ((num_species == 2 && isempty(crop_name1)) || (num_species == 2 && isempty(crop_name2)))
                        msgbox('Please enter enough information ','MLCan Error: Empty information','error');
                        return
                    else
                        if ((num_species == 3 && isempty(crop_name1)) || (num_species == 3 && isempty(crop_name2)) ||...
                                (num_species == 3 && isempty(crop_name3)))
                            msgbox('Please enter enough information ','MLCan Error: Empty information','error');
                            return
                        else
                            if ((num_species == 4 && isempty(crop_name1)) || (num_species == 4 && isempty(crop_name2)) ||...
                                    (num_species == 4 && isempty(crop_name3)) || (num_species == 4 && isempty(crop_name4)))
                                msgbox('Please enter enough information ','MLCan Error: Empty information','error');
                                return
                            else
                                if (isempty(DOY_start) || isempty(DOY_end))
                                    msgbox('Please enter enough information ','MLCan Error: Empty information','error');
                                    return
                                else
%                             if (isempty(LAImin_face) == 1 || LAImin_face<0)
%                                 msgbox('Minimum LAI must be positive number (> 0) ','MLCan Error: Empty information','error');
%                                 return
%                             else
                        
                        %{
                if num_can > length(canopy_systems(:,1))
                    canopy_systems = [canopy_systems;zeros(num_can-length(canopy_systems(:,1)+1))];
                else
                    canopy_systems = canopy_systems(1:num_can,:);
                end
                %}

%                 canopy_systems  = zeros(num_can1, 3);
%                 root_systems    = zeros(num_root1,3);

                                save './Temps/temp_variable.mat'...
                                    'num_species'       'lat_face'      'long_face'     'elev_face'...
                                    'DOY_start'         'DOY_end'  ...
                                    'crop_name1'        'crop_name2'    'crop_name3'    'crop_name4'...
                                    'ph_type1'          'ph_type2'      'ph_type3'      'ph_type4'...
                                    'num_LAD1'          'num_LAD2'      'num_LAD3'      'num_LAD4'...
                                    'num_root1'         'num_root2'     'num_root3'     'num_root4'...
                                    -append;

%                                'LAImin_face' ...                
%                                     'num_can1'          'num_can2'      'num_can3'      'num_can4'...
                                load ('./Temps/temporary.mat',...
                                    'working_name'  )
                                copyfile('./Temps/temp_variable.mat',working_name,'f');
                                close;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end



% load ('./Temps/temp_variable.mat',...
%         'num_species'                                                   ,   ...
%         'lat_face'      , 'long_face'   , 'elev_face'   , 'LAImin_face' ,   ...
%         'crop_name1'    , 'crop_name2'  , 'crop_name3'  , 'crop_name4'  ,   ...
%         'ph_type1'      , 'ph_type2'    , 'ph_type3'    , 'ph_type4'    ,   ...
%         'num_can1'      , 'num_can2'    , 'num_can3'    , 'num_can4'    ,   ...
%         'num_LAD1'      , 'num_LAD2'    , 'num_LAD3'    , 'num_LAD4'    ,   ...
%         'num_root1'     , 'num_root2'   , 'num_root3'   , 'num_root4'   ,   ...
%         'opt_root1'     , 'opt_root2'   , 'opt_root3'   , 'opt_root4'   );

















% --- Executes on button press in set_mod_cancel.
function set_mod_cancel_Callback(hObject, eventdata, handles)
close



function txt_crop_name1_Callback(hObject, eventdata, handles)
% hObject    handle to txt_crop_name1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_crop_name1 as text
%        str2double(get(hObject,'String')) returns contents of txt_crop_name1 as a double


% --- Executes during object creation, after setting all properties.
function txt_crop_name1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_crop_name1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_lat_Callback(hObject, eventdata, handles)
% hObject    handle to txt_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_lat as text
%        str2double(get(hObject,'String')) returns contents of txt_lat as a double


% --- Executes during object creation, after setting all properties.
function txt_lat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_long_Callback(hObject, eventdata, handles)
% hObject    handle to txt_long (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_long as text
%        str2double(get(hObject,'String')) returns contents of txt_long as a double


% --- Executes during object creation, after setting all properties.
function txt_long_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_long (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_canopy_layer1_Callback(hObject, eventdata, handles)
% hObject    handle to txt_canopy_layer1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_canopy_layer1 as text
%        str2double(get(hObject,'String')) returns contents of txt_canopy_layer1 as a double


% --- Executes during object creation, after setting all properties.
function txt_canopy_layer1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_root_layer1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function txt_root_layer1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function axes1_CreateFcn(hObject, eventdata, handles)
modeliconsmall = imread('./users/icons/wheat_small.png');
image(modeliconsmall);
axis off
%imshow('./users/icons/wheat_small.png');


function axes2_CreateFcn(hObject, eventdata, handles)
guifig = imread('./users/GUI_Fig.jpg');
image(guifig);
axis off
%imshow('./users/GUI_Fig.jpg');
axis image

%-------------------------------------------------------------------------%
% Call root profile function for entering root profile information        %
%-------------------------------------------------------------------------%
function set_mod_but_root_Callback(hObject, eventdata, handles)
num_root1 = get(handles.txt_root_layer1,'String');
if (isempty(num_root1) == 1 || isempty(str2num(num_root1)) == 1)
    msgbox('Root number can not be a blank or string','MLCan error', 'error');
    return
else
    if (mod(str2num(num_root1),1) ~= 0 || str2num(num_root1) <= 0)
        msgbox ('Root layer number must be an positive integer','MLCan error', 'error');
        return
    else
        save './Temps/temp_variable.mat' 'num_root1' -append;
        setup_root_profile;
    end
end

%-------------------------------------------------------------------------%
% Call LAD profile function for entering LAD profile information          %
%-------------------------------------------------------------------------%
function set_mod_but_LAD1_Callback(hObject, eventdata, handles)
num_LAD1 = get(handles.txt_LAD_layer1,'String');
if (isempty(num_LAD1) == 1 || isempty(str2num(num_LAD1)) == 1)
    msgbox('LAD layer number can not be a blank or string','MLCan error', 'error');
    return
else
    if (mod(str2num(num_LAD1),1) ~= 0 || str2num(num_LAD1) <= 0)
        msgbox ('LAD layer number must be an positive integer','MLCan error','error');
        return
    else
        save './Temps/temp_variable.mat' 'num_LAD1' -append;
        setup_LAD_profile;
    end
end

function txt_LAD_layer1_Callback(hObject, eventdata, handles)


function txt_LAD_layer1_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txt_specified_DOY_start_Callback(hObject, eventdata, handles)


function txt_specified_DOY_start_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txt_elevation_Callback(hObject, eventdata, handles)


function txt_elevation_CreateFcn(hObject, eventdata, handles)

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



function edit_num_species_Callback(hObject, eventdata, handles)
% hObject    handle to edit_num_species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_num_species as text
%        str2double(get(hObject,'String')) returns contents of edit_num_species as a double


% --- Executes during object creation, after setting all properties.
function edit_num_species_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_num_species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel_species1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_species1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over txt_crop_name1.
function txt_crop_name1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to txt_crop_name1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over set_mod_name.
function set_mod_name_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to set_mod_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uipanel_species1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_species1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in push_species1.
function push_species1_Callback(hObject, eventdata, handles)
% hObject    handle to push_species1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Num_of_Species_str = get(handles.edit_num_species,'String');
if (isempty(Num_of_Species_str) == 1 || isempty(str2num(Num_of_Species_str)))
    msgbox('Number of species can not be a blank or string','MLCan error', 'error');
    return
end
set(handles.uipanel_species1,'Visible','on');
set(handles.uipanel_species2,'Visible','off');
set(handles.uipanel_species3,'Visible','off');
set(handles.uipanel_species4,'Visible','off');

% --- Executes on button press in push_species2.
function push_species2_Callback(hObject, eventdata, handles)
% hObject    handle to push_species2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Num_of_Species_str = get(handles.edit_num_species,'String');
if (isempty(Num_of_Species_str) == 1 || isempty(str2num(Num_of_Species_str)))
    msgbox('Number of species can not be a blank or string','MLCan error', 'error');
    return
end

%Num_canopy_layer1_str=get(handles.txt_canopy_layer1,'String');
Num_of_LAD_str = get(handles.txt_LAD_layer1,'String');
Num_of_root_str = get(handles.txt_root_layer1,'String');
if (isempty(Num_of_LAD_str) == 1    || isempty(str2num(Num_of_LAD_str))     || ...
        isempty(Num_of_root_str) == 1   || isempty(str2num(Num_of_root_str)))
    msgbox('Species 1''s # canopy, LAD, and root layers can not be blanks or strings','MLCan error', 'error');
    %isempty(Num_canopy_layer1_str) == 1 || isempty(str2num(Num_canopy_layer1_str)) || ...   
    return
else
%    set(handles.txt_canopy_layer2,'String',Num_canopy_layer1_str);
%    set(handles.txt_canopy_layer2,'Enable','off');
    set(handles.txt_LAD_layer2,'String',Num_of_LAD_str);
    set(handles.txt_LAD_layer2,'Enable','off');
    set(handles.txt_root_layer2,'String',Num_of_root_str);
    set(handles.txt_root_layer2,'Enable','off');
end

set(handles.uipanel_species1,'Visible','off');
set(handles.uipanel_species2,'Visible','on');
set(handles.uipanel_species3,'Visible','off');
set(handles.uipanel_species4,'Visible','off');

% --- Executes on button press in push_species3.
function push_species3_Callback(hObject, eventdata, handles)
% hObject    handle to push_species3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Num_of_Species_str = get(handles.edit_num_species,'String');
if (isempty(Num_of_Species_str) == 1 || isempty(str2num(Num_of_Species_str)))
    msgbox('Number of species can not be a blank or string','MLCan error', 'error');
    return
end

%Num_canopy_layer1_str=get(handles.txt_canopy_layer1,'String');
Num_of_LAD_str = get(handles.txt_LAD_layer1,'String');
Num_of_root_str = get(handles.txt_root_layer1,'String');
if (isempty(Num_of_LAD_str) == 1    || isempty(str2num(Num_of_LAD_str))     || ...
        isempty(Num_of_root_str) == 1   || isempty(str2num(Num_of_root_str)))
    msgbox('Species 1''s # canopy, LAD, and root layers can not be blanks or strings','MLCan error', 'error');
    
    %isempty(Num_canopy_layer1_str) == 1 || isempty(str2num(Num_canopy_layer1_str)) || ...
    return
else
%    set(handles.txt_canopy_layer3,'String',Num_canopy_layer1_str);
%    set(handles.txt_canopy_layer3,'Enable','off');
    set(handles.txt_LAD_layer3,'String',Num_of_LAD_str);
    set(handles.txt_LAD_layer3,'Enable','off');
    set(handles.txt_root_layer3,'String',Num_of_root_str);
    set(handles.txt_root_layer3,'Enable','off');
end

set(handles.uipanel_species1,'Visible','off');
set(handles.uipanel_species2,'Visible','off');
set(handles.uipanel_species3,'Visible','on');
set(handles.uipanel_species4,'Visible','off');

% --- Executes on button press in push_species4.
function push_species4_Callback(hObject, eventdata, handles)
% hObject    handle to push_species4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Num_of_Species_str = get(handles.edit_num_species,'String');
if (isempty(Num_of_Species_str) == 1 || isempty(str2num(Num_of_Species_str)))
    msgbox('Number of species can not be a blank or string','MLCan error', 'error');
    return
end

%Num_canopy_layer1_str=get(handles.txt_canopy_layer1,'String');
Num_of_LAD_str = get(handles.txt_LAD_layer1,'String');
Num_of_root_str = get(handles.txt_root_layer1,'String');
if (isempty(Num_of_LAD_str) == 1    || isempty(str2num(Num_of_LAD_str))     || ...
        isempty(Num_of_root_str) == 1   || isempty(str2num(Num_of_root_str)))
    msgbox('Species 1''s # canopy, LAD, and root layers can not be blanks or strings','MLCan error', 'error');
    
    % isempty(Num_canopy_layer1_str) == 1 || isempty(str2num(Num_canopy_layer1_str)) || ...
        
    return
else
%    set(handles.txt_canopy_layer4,'String',Num_canopy_layer1_str);
%    set(handles.txt_canopy_layer4,'Enable','off');
    set(handles.txt_LAD_layer4,'String',Num_of_LAD_str);
    set(handles.txt_LAD_layer4,'Enable','off');
    set(handles.txt_root_layer4,'String',Num_of_root_str);
    set(handles.txt_root_layer4,'Enable','off');
end

set(handles.uipanel_species1,'Visible','off');
set(handles.uipanel_species2,'Visible','off');
set(handles.uipanel_species3,'Visible','off');
set(handles.uipanel_species4,'Visible','on');


function txt_crop_name2_Callback(hObject, eventdata, handles)
% hObject    handle to txt_crop_name2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_crop_name2 as text
%        str2double(get(hObject,'String')) returns contents of txt_crop_name2 as a double


% --- Executes during object creation, after setting all properties.
function txt_crop_name2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_crop_name2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_crop_name3_Callback(hObject, eventdata, handles)
% hObject    handle to txt_crop_name3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_crop_name3 as text
%        str2double(get(hObject,'String')) returns contents of txt_crop_name3 as a double


% --- Executes during object creation, after setting all properties.
function txt_crop_name3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_crop_name3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_root_layer2_Callback(hObject, eventdata, handles)
% hObject    handle to txt_root_layer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_root_layer2 as text
%        str2double(get(hObject,'String')) returns contents of txt_root_layer2 as a double


% --- Executes during object creation, after setting all properties.
function txt_root_layer2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_root_layer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_mod_but_root2.
function set_mod_but_root2_Callback(hObject, eventdata, handles)
% hObject    handle to set_mod_but_root2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setup_root_profile2;


function txt_canopy_layer2_Callback(hObject, eventdata, handles)
% hObject    handle to txt_canopy_layer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_canopy_layer2 as text
%        str2double(get(hObject,'String')) returns contents of txt_canopy_layer2 as a double


% --- Executes during object creation, after setting all properties.
function txt_canopy_layer2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_canopy_layer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_LAD_layer2_Callback(hObject, eventdata, handles)
% hObject    handle to txt_LAD_layer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_LAD_layer2 as text
%        str2double(get(hObject,'String')) returns contents of txt_LAD_layer2 as a double


% --- Executes during object creation, after setting all properties.
function txt_LAD_layer2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_LAD_layer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_mod_but_LAD2.
function set_mod_but_LAD2_Callback(hObject, eventdata, handles)
% hObject    handle to set_mod_but_LAD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setup_LAD_profile2;


function txt_root_layer3_Callback(hObject, eventdata, handles)
% hObject    handle to txt_root_layer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_root_layer3 as text
%        str2double(get(hObject,'String')) returns contents of txt_root_layer3 as a double


% --- Executes during object creation, after setting all properties.
function txt_root_layer3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_root_layer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_mod_but_root3.
function set_mod_but_root3_Callback(hObject, eventdata, handles)
% hObject    handle to set_mod_but_root3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setup_root_profile3;


function txt_canopy_layer3_Callback(hObject, eventdata, handles)
% hObject    handle to txt_canopy_layer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_canopy_layer3 as text
%        str2double(get(hObject,'String')) returns contents of txt_canopy_layer3 as a double


% --- Executes during object creation, after setting all properties.
function txt_canopy_layer3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_canopy_layer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_LAD_layer3_Callback(hObject, eventdata, handles)
% hObject    handle to txt_LAD_layer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_LAD_layer3 as text
%        str2double(get(hObject,'String')) returns contents of txt_LAD_layer3 as a double


% --- Executes during object creation, after setting all properties.
function txt_LAD_layer3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_LAD_layer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_mod_but_LAD3.
function set_mod_but_LAD3_Callback(hObject, eventdata, handles)
% hObject    handle to set_mod_but_LAD3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setup_LAD_profile3;

function txt_crop_name4_Callback(hObject, eventdata, handles)
% hObject    handle to txt_crop_name4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_crop_name4 as text
%        str2double(get(hObject,'String')) returns contents of txt_crop_name4 as a double


% --- Executes during object creation, after setting all properties.
function txt_crop_name4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_crop_name4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_root_layer4_Callback(hObject, eventdata, handles)
% hObject    handle to txt_root_layer4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_root_layer4 as text
%        str2double(get(hObject,'String')) returns contents of txt_root_layer4 as a double


% --- Executes during object creation, after setting all properties.
function txt_root_layer4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_root_layer4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_mod_but_root4.
function set_mod_but_root4_Callback(hObject, eventdata, handles)
% hObject    handle to set_mod_but_root4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setup_root_profile4;


function txt_canopy_layer4_Callback(hObject, eventdata, handles)
% hObject    handle to txt_canopy_layer4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_canopy_layer4 as text
%        str2double(get(hObject,'String')) returns contents of txt_canopy_layer4 as a double


% --- Executes during object creation, after setting all properties.
function txt_canopy_layer4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_canopy_layer4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_LAD_layer4_Callback(hObject, eventdata, handles)
% hObject    handle to txt_LAD_layer4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_LAD_layer4 as text
%        str2double(get(hObject,'String')) returns contents of txt_LAD_layer4 as a double


% --- Executes during object creation, after setting all properties.
function txt_LAD_layer4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_LAD_layer4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_mod_but_LAD4.
function set_mod_but_LAD4_Callback(hObject, eventdata, handles)
% hObject    handle to set_mod_but_LAD4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setup_LAD_profile4;

% --- Executes on button press in pushbutton_set_num_of_species.
function pushbutton_set_num_of_species_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set_num_of_species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in push_set_num_of_species.
function push_set_num_of_species_Callback(hObject, eventdata, handles)
% hObject    handle to push_set_num_of_species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Num_of_Species_str = get(handles.edit_num_species,'String');
Num_of_Species = str2num(Num_of_Species_str);
if Num_of_Species == 1
    set(handles.push_species1,'enable','on');
    set(handles.push_species2,'enable','off');
    set(handles.push_species3,'enable','off');
    set(handles.push_species4,'enable','off');
end
if Num_of_Species == 2
    set(handles.push_species1,'enable','on');
    set(handles.push_species2,'enable','on');
    set(handles.push_species3,'enable','off');
    set(handles.push_species4,'enable','off');
end
if Num_of_Species == 3
    set(handles.push_species1,'enable','on');
    set(handles.push_species2,'enable','on');
    set(handles.push_species3,'enable','on');
    set(handles.push_species4,'enable','off');
end
if Num_of_Species == 4
    set(handles.push_species1,'enable','on');
    set(handles.push_species2,'enable','on');
    set(handles.push_species3,'enable','on');
    set(handles.push_species4,'enable','on');
end
integerTest=~mod(Num_of_Species,1);
if isempty(integerTest) == 1
        msgbox ('Number of species must be a integer between 1 and 4 ','MLCan Error','error');
        return
else
    if (Num_of_Species < 1 || Num_of_Species > 4 || integerTest ~= 1)
        msgbox ('Number of species must be a integer between 1 and 4 ','MLCan Error','error');
        return
    end
end


% --- Executes during object deletion, before destroying properties.
function uipanel1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to model_setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object deletion, before destroying properties.
function model_setup_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to model_setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function txt_specified_DOY_end_Callback(hObject, eventdata, handles)
% hObject    handle to txt_specified_DOY_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_specified_DOY_end as text
%        str2double(get(hObject,'String')) returns contents of txt_specified_DOY_end as a double


% --- Executes during object creation, after setting all properties.
function txt_specified_DOY_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_specified_DOY_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DOY_start_Callback(hObject, eventdata, handles)
% hObject    handle to txt_specified_DOY_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_specified_DOY_start as text
%        str2double(get(hObject,'String')) returns contents of txt_specified_DOY_start as a double


% --- Executes during object creation, after setting all properties.
function DOY_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_specified_DOY_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DOY_end_Callback(hObject, eventdata, handles)
% hObject    handle to txt_specified_DOY_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_specified_DOY_end as text
%        str2double(get(hObject,'String')) returns contents of txt_specified_DOY_end as a double


% --- Executes during object creation, after setting all properties.
function DOY_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_specified_DOY_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
