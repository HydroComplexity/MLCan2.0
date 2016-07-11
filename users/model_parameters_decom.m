function varargout = model_parameters_decom(varargin)
% MODEL_PARAMETERS_DECOM M-file for model_parameters_decom.fig
%      MODEL_PARAMETERS_DECOM, by itself, creates a new MODEL_PARAMETERS_DECOM or raises the existing
%      singleton*.
%
%      H = MODEL_PARAMETERS_DECOM returns the handle to a new MODEL_PARAMETERS_DECOM or the handle to
%      the existing singleton*.
%
%      MODEL_PARAMETERS_DECOM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL_PARAMETERS_DECOM.M with the given input arguments.
%
%      MODEL_PARAMETERS_DECOM('Property','Value',...) creates a new MODEL_PARAMETERS_DECOM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before model_parameters_decom_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to model_parameters_decom_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help model_parameters_decom

% Last Modified by GUIDE v2.5 20-Jan-2014 15:34:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @model_parameters_decom_OpeningFcn, ...
                   'gui_OutputFcn',  @model_parameters_decom_OutputFcn, ...
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


% --- Executes just before model_parameters_decom is made visible.
function model_parameters_decom_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to model_parameters_decom (see VARARGIN)

% Choose default command line output for model_parameters_decom
handles.output = hObject;

load './Temps/temp_variable.mat' 'dat_decom' 'dat_decom_litter' 'num_root1' 
dat_length = length(dat_decom(:,1));
if isstr(num_root1) == 1
    Dec_num = str2num(num_root1);
else
    Dec_num = num_root1;
end
if dat_length > Dec_num
    dat_decom = Dec_num(1:Dec_num,:);
else
	dat_decom_add = zeros(Dec_num-dat_length,4);
	dat_decom = [dat_decom; dat_decom_add];
end

set(handles.par_Decomposition_tab,'Data',dat_decom)
set(handles.par_Decomposition_tab_litter,'Data',dat_decom_litter)
par_Decomposition_tab_CellEditCallback(hObject, eventdata, handles)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes model_parameters_decom wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = model_parameters_decom_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function but_DEC_load_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile(...
    {'*.xls;*.xlsx'     , 'Microsoft Excel Files (*.xls,*.xlsx)';   ...
    '*.csv'             , 'CSV - comma delimited (*.csv)'       ;   ... 
    '*.txt'             , 'Text (Tab Delimited (*.txt)'         ;   ...
    '*.*'               , 'All Files (*.*)'                     },  ...
    'Load LAD from files');
    
if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
else
	load './Temps/temp_variable.mat' 'dat_decom' 'dat_decom_litter' 'num_root1'
    str_file = regexp(filename, '\.', 'split');
    if (strcmp(str_file(end),'xls') == 1 || strcmp(str_file(end),'xlsx') == 1)
        DEC_data_table = xlsread([pathname,filename]);      
    elseif strcmp(str_file(end),'csv') == 1
        DEC_data_table = csvread([pathname,filename], 1, 0);      
    elseif strcmp(str_file(end),'txt') == 1
        DEC_data_table = dlmread([pathname,filename], '\t', 1, 0);      
    else
        msgbox('The file extension is unrecognized, please choose another file','MLCan Error','Error');
        return
    end
    
    if isstr(num_root1)
        num_root1=str2num(num_root1);
    end
    if length(DEC_data_table(:,1)) ~= (num_root1+1)
        numstr=num2str(length(DEC_data_table(:,1)));
        notice = ['Number of parameters (',num_root1,') must be equal to data in imported file (',numstr,')'];
        msgbox(notice, 'MLCan Error: Different data length','Error');
        return        
    else
        if length(DEC_data_table(1,:)) ~= 4
            msgbox ('Incorrect format or input file''MLcan Error','Error')
            return
        else
            load './Temps/temp_variable.mat' 'litter_depth' 'dat_root1'
            if (sum(DEC_data_table(:,1) - [litter_depth; dat_root1(:,1)]) > 0.001 || sum(DEC_data_table(:,1) - [litter_depth; dat_root1(:,1)]) < - 0.001)
                msgbox ('Soil/litter depth distribution dose not match up with what is in MODEL SETUP - check the soil/litter depth','MLCan Error','Error')
                return
            else
                set(handles.par_Decomposition_tab,'Data',DEC_data_table(2:end,:))
                set(handles.par_Decomposition_tab_litter,'Data',DEC_data_table(1,:))
                par_Decomposition_tab_CellEditCallback(hObject, eventdata, handles);
            end
        end
    end
end





function but_DEC_cancel_Callback(hObject, eventdata, handles)
close

function but_DEC_ok_Callback(hObject, eventdata, handles)
dat_decom=get(handles.par_Decomposition_tab,'Data');
dat_decom_litter=get(handles.par_Decomposition_tab_litter,'Data');

save './Temps/temp_variable.mat' 'dat_decom' 'dat_decom_litter' -append
load ('./Temps/temporary.mat',...
    'working_name'  )
copyfile('./Temps/temp_variable.mat',working_name,'f');

close;


function par_Decomposition_tab_CellEditCallback(hObject, eventdata, handles)
dat_decom        = get(handles.par_Decomposition_tab,'Data');
soil_depth      = dat_decom(:,1);
MDR             = dat_decom(:,2);
LDR             = dat_decom(:,3);
HDR             = dat_decom(:,4);

plot(handles.axes_MDR_plot,MDR,soil_depth,'Color', 'b',...
        'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
        'MarkerFaceColor','g','MarkerSize',6);
axis(handles.axes_MDR_plot,[0 max(MDR)*1.2+10^(-10) min(soil_depth)/1.2 max(soil_depth)*1.1+10^(-10)]);
ylabel(handles.axes_MDR_plot,'\bf Soil Depth [m]');
legend(handles.axes_MDR_plot,'Microbial Death Rate');
grid(handles.axes_MDR_plot,'on');

plot(handles.axes_LDR_plot,LDR,soil_depth,'Color', 'r',...
            'LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b',...
            'MarkerFaceColor','m','MarkerSize',6);
axis(handles.axes_LDR_plot,[0 max(LDR)*1.2+10^(-10) min(soil_depth)/1.2 max(soil_depth)*1.1+10^(-10)]);
ylabel(handles.axes_LDR_plot,'\bf Soil Depth [m]');
legend(handles.axes_LDR_plot,'Litter Decomposition Rate');
grid(handles.axes_LDR_plot,'on');

plot(handles.axes_HDR_plot,HDR,soil_depth,'Color', 'b',...
        'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
        'MarkerFaceColor','g','MarkerSize',6);
axis(handles.axes_HDR_plot,[0 max(HDR)*1.2+10^(-10) min(soil_depth)/1.2 max(soil_depth)*1.1+10^(-10)]);
ylabel(handles.axes_HDR_plot,'\bf Soil Depth [m]');
legend(handles.axes_HDR_plot,'Humus Decomposition Rate');
grid(handles.axes_HDR_plot,'on');
