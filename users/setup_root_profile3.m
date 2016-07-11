function varargout = setup_root_profile3(varargin)
% SETUP_ROOT_PROFILE3 M-file for setup_root_profile3.fig
%      SETUP_ROOT_PROFILE3, by itself, creates a new SETUP_ROOT_PROFILE3 or raises the existing
%      singleton*.
%
%      H = SETUP_ROOT_PROFILE3 returns the handle to a new SETUP_ROOT_PROFILE3 or the handle to
%      the existing singleton*.
%
%      SETUP_ROOT_PROFILE3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETUP_ROOT_PROFILE3.M with the given input arguments.
%
%      SETUP_ROOT_PROFILE3('Property','Value',...) creates a new SETUP_ROOT_PROFILE3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before setup_root_profile3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to setup_root_profile3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help setup_root_profile3

% Last Modified by GUIDE v2.5 19-Jan-2014 02:12:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @setup_root_profile3_OpeningFcn, ...
                   'gui_OutputFcn',  @setup_root_profile3_OutputFcn, ...
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


% --- Executes just before setup_root_profile3 is made visible.
function setup_root_profile3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to setup_root_profile3 (see VARARGIN)

% Choose default command line output for setup_root_profile3
handles.output = hObject;

% load './Temps/temp_variable.mat' 'opt_root'
% if opt_root == 1
% 	set(handles.pan_root_obs,'Visible','on');
% 	set(handles.pan_root_eqn,'Visible','off');
% 	set(handles.but_root_load,'Enable','on');
%     set(handles.root_rad_obs,'Value',1);
% else
%     set(handles.pan_root_obs,'Visible','off');
%     set(handles.pan_root_eqn,'Visible','on');
%     set(handles.but_root_load,'Enable','off');
%     set(handles.root_rad_eqn,'Value',1);
% end    
load './Temps/temp_variable.mat' 'opt_root3'
if opt_root3 == 1
	set(handles.pan_root_obs,'Visible','on');
	set(handles.pan_root_eqn,'Visible','off');
	set(handles.but_root_load,'Enable','on');
    set(handles.root_rad_obs,'Value',1);
else
    set(handles.pan_root_obs,'Visible','off');
    set(handles.pan_root_eqn,'Visible','on');
    set(handles.but_root_load,'Enable','off');
    set(handles.root_rad_eqn,'Value',1);
end    

% Soil grid should not change----------------------------------------------
set(handles.setup_root_tab_eqn,'Enable','off')

% Load table and information for observed option---------------------------
    load './Temps/temp_variable.mat' 'num_root1' 'dat_root3'
    num_root3=num_root1;
    
    dat_length = length(dat_root3(:,1));
    %root_num3 = str2num(num_root3);

    if isstr(num_root3) == 1
        root_num3 = str2num(num_root3);
    else
        root_num3 = (num_root3);
    end
    
    if dat_length > root_num3
        dat_root3 = dat_root3(1:root_num3,:);
    elseif dat_length < root_num3
        dat_root_add = zeros(root_num3-dat_length,2);
        dat_root3 = [dat_root3; dat_root_add];
    end

    plot(handles.axes_root_plot,dat_root3,dat_root3,'Color', 'r',...
        'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
        'MarkerFaceColor','g','MarkerSize',6);
    
    set(handles.setup_root_tab,'Data',dat_root3)
    xlabel(handles.axes_root_plot,'\bf Root Length Density [m m^{-3}]');
    ylabel(handles.axes_root_plot,'\bf Root Depth [m]');
    legend(handles.axes_root_plot,'Root Profile');
    setup_root_tab_CellEditCallback(hObject, eventdata, handles)
% End of load information for observed option------------------------------

% Load table and information for equation option---------------------------
    load './Temps/temp_variable.mat' 'set_para_root3' 'set_root_func'
    para_root_eqn_name = set_para_root3(:,1);
    setappdata(0,'para_root_eqn_name',para_root_eqn_name);
    set(handles.setup_root_tab_para,'Data',set_para_root3);
    set(handles.setup_root_tab_eqn,'Data',set_root_func);
    grid(handles.axes_root_eqn,'on');
    setup_root_tab_eqn_CellEditCallback(hObject, eventdata, handles)
% End of load information for equation option------------------------------

% Load litter depth
    load './Temps/temp_variable.mat' 'litter_depth'
    set(handles.txt_litterDepth3,'String',litter_depth);   
    set(handles.txt_litterDepth3,'Enable','off');
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes setup_root_profile3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = setup_root_profile3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function but_root_load_Callback(hObject, eventdata, handles)

[filename, pathname] = uigetfile(...
    {'*.xls;*.xlsx'     , 'Microsoft Excel Files (*.xls,*.xlsx)';   ...
    '*.csv'             , 'CSV - comma delimited (*.csv)'       ;   ... 
    '*.txt'             , 'Text (Tab Delimited (*.txt)'         ;   ...
    '*.*'               , 'All Files (*.*)'                     },  ...
    'Load LAD from files');
    %'*.xls;*.xlsx'     , 'Microsoft Excel Files (*.xls,*.xlsx)';   ...

if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
else
	load './Temps/temp_variable.mat' 'num_root3'
    str_file = regexp(filename, '\.', 'split');
    if (strcmp(str_file(end),'xls') == 1 || strcmp(str_file(end),'xlsx') == 1)
        root_data_table = xlsread([pathname,filename]);      
    elseif strcmp(str_file(end),'csv') == 1
        root_data_table = csvread([pathname,filename], 1, 0);      
    elseif strcmp(str_file(end),'txt') == 1
        root_data_table = dlmread([pathname,filename], '\t', 1, 0);      
    else
        msgbox('The file extension is unrecognized, please choose another file','MLCan Error','Error');
        return
    end

    if length(root_data_table(:,1)) ~= str2num(num_root3)
        numstr=num2str(length(root_data_table(:,1)));
        notice = ['Number of root layer (',num_root3,') must be equal to data in imported file (',numstr,')'];
        msgbox(notice, 'MLCan Error: Different data length','Error');
        return        
    else
        if length(root_data_table(1,:)) ~= 2
            msgbox ('Incorrect format or input file','MLcan Error','Error')
            return
        else
            set(handles.setup_root_tab,'Data',root_data_table);   
            setup_root_tab_CellEditCallback(hObject, eventdata, handles);
        end
    end
end

function but_root_cancel_Callback(hObject, eventdata, handles)
close

function but_root_ok_Callback(hObject, eventdata, handles)
num_root_eqn    = num2str(getappdata(0,'nl_soil'));
opt_root3       = get(handles.root_rad_obs,'Value');
set_para_root3  = get(handles.setup_root_tab_para,'Data');
set_root_func   = get(handles.setup_root_tab_eqn,'Data');
if opt_root3 == 1
    dat_root3        = get(handles.setup_root_tab,'Data');
else
    zns_eqn         = getappdata(0,'zns');
    rootfr_eqn      = getappdata(0,'rootfr');
    dat_root3        = [zns_eqn, rootfr_eqn];
end
save './Temps/temp_variable.mat'...
        'dat_root3' 'opt_root3' 'set_para_root3' 'set_root_func' 'num_root_eqn' -append
close;


function setup_root_tab_CellEditCallback(hObject, eventdata, handles)
dat_root3    = get(handles.setup_root_tab,'Data');
root_depth  = dat_root3(:,1);
RLD         = dat_root3(:,2);

plot(handles.axes_root_plot,RLD,root_depth,'Color', 'r',...
        'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
        'MarkerFaceColor','g','MarkerSize',6);
axis(handles.axes_root_plot,[0 max(RLD)*1.2+0.1 min(root_depth)/1.2-0.2 max(root_depth)*1.1+0.1]);
xlabel(handles.axes_root_plot,'\bf Root Length Density [m m^{-3}]');
ylabel(handles.axes_root_plot,'\bf Root Depth [m]');
legend(handles.axes_root_plot,'Root Profile');
grid(handles.axes_root_plot,'on');


function set_pan_options_SelectionChangeFcn(hObject, eventdata,handles)

handles = guidata(hObject);

switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'root_rad_obs'
      set(handles.pan_root_obs,'Visible','on');
      set(handles.pan_root_eqn,'Visible','off');
      set(handles.but_root_load,'Enable','on');
    case 'root_rad_eqn'
      set(handles.pan_root_obs,'Visible','off');
      set(handles.pan_root_eqn,'Visible','on');
      set(handles.but_root_load,'Enable','off');
    otherwise
       % Code for when there is no match.
end
%updates the handles structure
guidata(hObject, handles);


function setup_root_tab_para_CellEditCallback(hObject, eventdata, handles)
para_root_eqn       = get(handles.setup_root_tab_para,'Data');
maxrootdepth        = cell2mat(para_root_eqn(1,2));
z50                 = cell2mat(para_root_eqn(2,2));
z95                 = cell2mat(para_root_eqn(3,2));
load './Temps/temp_variable.mat' 'set_para_root1'
Nrootcut            = cell2mat(set_para_root1(4,2));

if z50 > z95
    msgbox('z50 must be less than z95','MLCan Error','Warn');
    z50 = z95;    
    para_root_eqn   = {...
    ' Maximum depth'    , maxrootdepth  ,   ; ...
    ' z50'              , z50           ,   ; ...
    ' z95'              , z95           };
    
    set(handles.setup_root_tab_para,'Data',para_root_eqn);
    return
else
    zns = getappdata(0,'zns');
    dzs = getappdata(0,'dzs');
    zhs = getappdata(0,'zhs');
    [rootfr] = ROOTDIST_LOGISTIC(zns,zhs,dzs,z50,z95,maxrootdepth);
    % rootcut
    %    rootfr(1:Nrootcut,:) = 0;
    % root cut end
    
        rootfr_norm = rootfr./(dzs*1000);
        rootfr_norm = rootfr_norm./(sum(rootfr_norm));

    setappdata(0,'rootfr',rootfr);
    plot(handles.axes_root_eqn,rootfr_norm,zns,'Color', 'b',...
            'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
            'MarkerFaceColor','g','MarkerSize',6);
    axis(handles.axes_root_eqn,[0 max(rootfr_norm)+.05 min(zns)/1.2-0.2 max(zns)]);
    xlabel(handles.axes_root_eqn,'\bf Root Fraction');
    ylabel(handles.axes_root_eqn,'\bf z [m]');
    legend(handles.axes_root_eqn,'Root Profile');
    grid(handles.axes_root_eqn,'on');
end


function setup_root_tab_eqn_CellEditCallback(hObject, eventdata, handles)
load './Temps/temp_variable.mat' 'num_root1';
num_root3=num_root1;

%nl_soil=str2num(num_root3);
if isstr(num_root3) == 1
    nl_soil=str2num(num_root3);
else
    nl_soil=(num_root3);
end
tab_eqn         = get(handles.setup_root_tab_eqn,'Data');
scaleza         = cell2mat(tab_eqn(1,2));
scalezb         = cell2mat(tab_eqn(2,2));

%[zns, dzs, zhs, nl_soil] = SOILGRIDtest (dzs,depths);
[zns, dzs, zhs, nl_soil] = SOILGRIDtest (nl_soil,scaleza,scalezb);
setappdata(0,'zns',zns);
setappdata(0,'dzs',dzs);
setappdata(0,'zhs',zhs);
setappdata(0,'nl_soil',nl_soil);
setup_root_tab_para_CellEditCallback(hObject, eventdata, handles)


function [znode_all, dz_all, zlayer_all, nl_soil] = SOILGRIDtest (nl_soil,scaleza,scalezb)

%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                             FUNCTION CODE                             %%
%%                       CONSTRUCT THE SOIL GRID                         %%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   This function constructs the soil grid, with 0.1m grid spacing down   %
%   to 1m depth, and 0.5 meter spacing below that, down to the specified  %
%   maxdepth                                                              %
%                                                                         %
%   INPUTS:                                                               %
%       dzs 	= Layer spacings [m]                                      %
%       depths 	= Depth of soil column to grid for each layer             %
%                   spacing (dzs) [m] --> sum(depths) = total depth       %
%                   of soil column                                        %
%                                                                         %
%   OUTPUTS:                                                              %
%       znode   - node depths (center of layers)     [m]                  %
%       dznode  - layer thicknesses                  [m]                  %
%       zlayer  - depths of soil layer interfaces    [m]                  %
%       nl_soil - number of soil layers              [-]                  %
%                                                                         %
%-------------------------------------------------------------------------%
%   Created by  : Darren Drewry                                           %
%   Date        : January 01, 2008                                        %
%   Editted by  : Phong Le                                                %
%   Date        : December 23, 2009                                       %
%   Editted by  : Dongkook Woo                                            %
%   Date        : September 08, 2013                                      %
%% --------------------------------------------------------------------- %%  
%% 
    % Phong's edition
%     if (length(dzs) ~= length(depths))		% Check the consistency	  %
%         disp('*** ERROR: ''dzs'' and ''depths'' must have same length!');
%         return;
%     end
% %
% %    
%     % Construct Grid
%     zl = 0;
%     znode_all=[]; dz_all=[]; zlayer_all=[];
%     for ii = 1:length(dzs)
%         znode   = []; 
%         dz      = []; 
%         zlayer  = [];
%     %  
%         dzii    = dzs(ii);
%         depthii = depths(ii);
%     %    
%         znode   = linspace(dzii/2, depthii, (depthii-dzii)/dzii+1);
%         dz      = ones(1,length(znode))*dzii;
%         zlayer  = znode + dz/2;
%     %    
%         if (ii==1)
%             znode_all   = [znode_all, znode];
%             zlayer_all  = [zlayer_all, zlayer];
%         else
%             znode_all   = [znode_all, znode+zlayer_all(end)];
%             zlayer_all  = [zlayer_all, zlayer+zlayer_all(end)];
%         end
%         dz_all  = [dz_all, dz];    
%     end
% %    
%     nl_soil     = length(zlayer_all);
% %
% %    
%     % Make Column Vectors
%     znode_all   = znode_all(:);
%     dz_all      = dz_all(:);
%     zlayer_all  = zlayer_all(:);

% Dongkook's edition
%  The code was changed. It uses the same distribution  as Amenu (Common land model)

% ALLOCATE MEMORY
zsoi = nan(nl_soil,1);
dzsoi = nan(nl_soil,1);
zsoih = nan(nl_soil,1);

% Soil layer node depths, i.e., depth of layer center from surface [m]
for j = 1:nl_soil
    zsoi(j) = scaleza*(exp(scalezb*(j-0.5))-1);  %node depths
end

% Soil layer thicknesses
dzsoi(1)  = 0.5*(zsoi(1)+zsoi(2));
dzsoi(nl_soil)= zsoi(nl_soil)-zsoi(nl_soil-1);
for j = 2:nl_soil-1
    dzsoi(j)= 0.5*(zsoi(j+1)-zsoi(j-1));
end

% Soil layer interface depths from the surface [m]
zsoih(nl_soil) = zsoi(nl_soil) + 0.5*dzsoi(nl_soil);
for j = 1:nl_soil-1
    zsoih(j)= 0.5*(zsoi(j)+zsoi(j+1));
end

% Make Column Vectors
znode_all = zsoi(:);
dz_all = dzsoi(:);
zlayer_all = zsoih(:);
        
%
%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%% <<<<<<<<<<<<<<<<<<<<<<<<< END OF FUNCTION >>>>>>>>>>>>>>>>>>>>>>>>>>>>%%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%

function [rootfr] = ROOTDIST_LOGISTIC(zn,zh,dz,z50,z95,maxrootdepth)
% 
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                             FUNCTION CODE                             %%
%%                           ROOTDIST_LOGISTIC                           %%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   This function calculate the logistic root fraction distribution       %                                                  %
%   following the the methodology of Schenk and Jackson(2002). Root       %
%   fraction distribution will be used to calculate the root hydraulic    %
%   conductivity in soil layers.                                          %
%   See Amenu and Kumar, 2008 - Eqn 16, pp. 60                            %
%                                                                         %
%   INPUTS:                                                               %
%       zn = node depths [m]                                              %
%       dz = layer thicknesses [m]                                        %
%       Z50 = depth to which 50% of root biomass exists                   %
%       Z95 = depth to which 95% of root biomass exists                   %
%                                                                         %
%-------------------------------------------------------------------------%
%   Created by  : Darren Drewry                                           %
%   Editted by  : Phong Le, Dongkook Woo                                  %
%   Date        : December 23, 2009                                       %
%   Editted by  : Dongkook Woo                                            %
%   Date        : September 08, 2013                                      %
%% --------------------------------------------------------------------- %%  
%% 
    cc     = 1.27875 / (log10(z50) - log10(z95));                           % see Root_distribution_logistic.pdf and Amenu & Kumar, 2008.
%
    rootfr = -dz*cc/z50 .* (zn/z50).^(cc-1) .* (1+(zn/z50).^cc).^-2;
%
% Cut off roots at maximum depth
    % Phong's edition
    %mind   = find(zn>maxrootdepth, 1, 'first');
    %rootfr(mind:end) = 0;
    % Dongkook's edition
    mind = find(zh>maxrootdepth, 1, 'first');
    rootfr(mind+1:end) = 0;

 %
    rootfr = rootfr / sum(rootfr);
%
%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%% <<<<<<<<<<<<<<<<<<<<<<<<< END OF FUNCTION >>>>>>>>>>>>>>>>>>>>>>>>>>>>%%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%



function txt_litterDepth3_Callback(hObject, eventdata, handles)
% hObject    handle to txt_litterDepth3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_litterDepth3 as text
%        str2double(get(hObject,'String')) returns contents of txt_litterDepth3 as a double


% --- Executes during object creation, after setting all properties.
function txt_litterDepth3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_litterDepth3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
