function varargout = model_forcings(varargin)
% MODEL_FORCINGS M-file for model_forcings.fig
%      MODEL_FORCINGS, by itself, creates a new MODEL_FORCINGS or raises
%      the existing
%      singleton*.
%
%      H = MODEL_FORCINGS returns the handle to a new MODEL_FORCINGS or the handle to
%      the existing singleton*.
%
%      MODEL_FORCINGS('CALLBACK',hObject,eventData,handles,...) calls the
%      local
%      function named CALLBACK in MODEL_FORCINGS.M with the given input arguments.
%
%      MODEL_FORCINGS('Property','Value',...) creates a new MODEL_FORCINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before model_forcings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to model_forcings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help model_forcings

% Last Modified by GUIDE v2.5 19-Jan-2014 15:06:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @model_forcings_OpeningFcn, ...
                   'gui_OutputFcn',  @model_forcings_OutputFcn, ...
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


% --- Executes just before model_forcings is made visible.
function model_forcings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to model_forcings (see VARARGIN)

% Choose default command line output for model_forcings
handles.output = hObject;

% Load table and information for observed option---------------------------
load './Temps/temp_variable.mat'...
    'Soil_heat' 'root_name' 'nutrient_name1' 'nutrient_name2' 'Soil_C_pool1' 'root_init' 'root_init_litter'   ...
    'nutrient_int' 'nutrient_int_litter' 'num_root1' 'dat_root1' 'Soil_nutrient'   ...
    'working_forcings'
dat_length = length(root_init(:,1));
if isstr(num_root1) == 1
    root_num = str2num(num_root1);
else
    root_num = num_root1;
end
if dat_length > root_num
    root_init = root_init(1:root_num,:);
elseif dat_length < root_num
    root_init_add = zeros(root_num-dat_length,3);
    root_init = [root_init; root_init_add];
end
root_init = [dat_root1(:,1),root_init(:,2:3)];

if Soil_nutrient == 1
    dat_length_nutrient = length(nutrient_int(:,1));
    if dat_length_nutrient > root_num
        nutrient_int = nutrient_int(1:root_num,:);
    elseif dat_length_nutrient < root_num
        nutrient_int_add = zeros(root_num-dat_length_nutrient,size(nutrient_int,2));
        if Soil_C_pool1 == 1 % 1 Soil organic carbon pool
            nutrient_int_add2 = zeros(root_num,5);
        else % 2 Soil organic carbon pool
            nutrient_int_add2 = zeros(root_num,7);
        end
        nutrient_int = [nutrient_int; nutrient_int_add];
        nutrient_int = [nutrient_int nutrient_int_add2];
    end
    nutrient_int = [dat_root1(:,1),nutrient_int(:,2:end)];
end

button_load_picture(hObject, eventdata, handles);
set(handles.forc_gbut_options,'SelectionChangeFcn',@forc_gbut_options_SelectionChangeFcn);
set(handles.forc_pan_forcing,'Visible','on');
set(handles.forc_pan_condition,'Visible','off');
set(handles.forc_pan_nutrient,'Visible','off');
set(handles.forc_but_load,'Enable','off');
set(handles.forc_but_load_nutrient,'Enable','off');

set(handles.forc_table_root,'ColumnName',root_name,'Data',root_init);
set(handles.forc_table_root,'ColumnEditable',logical([0 1 Soil_heat]));
set(handles.forc_table_root_litter,'ColumnName',root_name,'Data',root_init_litter);
set(handles.forc_table_root_litter,'ColumnEditable',logical([0 1 Soil_heat]));

if Soil_nutrient == 1
    set(handles.forc_tabbut_nutrient,'Enable','on');
else
    set(handles.forc_tabbut_nutrient,'Enable','off');
end


% 1 Soil organic carbon pool
if Soil_C_pool1 == 1
    set(handles.forc_table_biogeochemical,'ColumnName',nutrient_name1,'Data',nutrient_int);
    set(handles.forc_table_biogeochemical,'ColumnEditable',logical([0 1 1 1 1 1 1 1]));
    set(handles.forc_table_biogeochemical_litter,'ColumnName',nutrient_name1,'Data',nutrient_int_litter);
    set(handles.forc_table_biogeochemical_litter,'ColumnEditable',logical([0 1 1 1 1 1 1 1]));
    set(handles.forc_but_litterC,'Enable','off');
    set(handles.forc_but_litterC,'Enable','off');
    set(handles.forc_but_humusC,'Enable','off');
%    set(handles.forc_but_litterCdr,'Enable','off');
%    set(handles.forc_but_litterCdr,'Visible','off');
    set(handles.forc_but_litterCN,'Enable','off');
    set(handles.forc_but_litterCN,'Enable','off');
%    set(handles.forc_but_humusCdr,'Enable','off');    
    
    set(handles.forc_but_organicC,'Enable','on');
    set(handles.forc_but_organicC,'Enable','on');
%    set(handles.forc_but_organicCdr,'Enable','on');
%    set(handles.forc_but_organicCdr,'Visible','on');
    set(handles.forc_but_organicCN,'Enable','on');
    set(handles.forc_but_organicCN,'Enable','on');
    
else
    % 2 Soil organic carbon pool
    set(handles.forc_table_biogeochemical,'ColumnName',nutrient_name2,'Data',nutrient_int);
    set(handles.forc_table_biogeochemical,'ColumnEditable',logical([0 1 1 1 1 1 1 1]));    
    set(handles.forc_table_biogeochemical_litter,'ColumnName',nutrient_name2,'Data',nutrient_int_litter);
    set(handles.forc_table_biogeochemical_litter,'ColumnEditable',logical([0 1 1 1 1 1 1 1]));  
    set(handles.forc_but_litterC,'Enable','on');
    set(handles.forc_but_litterC,'Enable','on');
    set(handles.forc_but_humusC,'Enable','on');
%    set(handles.forc_but_litterCdr,'Enable','on');
%    set(handles.forc_but_litterCdr,'Visible','on');
    set(handles.forc_but_litterCN,'Enable','on');
    set(handles.forc_but_litterCN,'Enable','on');
%    set(handles.forc_but_humusCdr,'Enable','on'); 
    
    set(handles.forc_but_organicC,'Enable','off');
    set(handles.forc_but_organicC,'Enable','off');
%    set(handles.forc_but_organicCdr,'Enable','off');
%    set(handles.forc_but_organicCdr,'Visible','off');
    set(handles.forc_but_organicCN,'Enable','off');
    set(handles.forc_but_organicCN,'Enable','off');
end
if Soil_nutrient == 1
    forc_but_nitrate_Callback(hObject, eventdata, handles);
    forc_but_microbialC_Callback(hObject, eventdata, handles);
end

% Forcings
if isempty(working_forcings)
    rad_new_forcings(hObject, eventdata, handles)
    set(handles.forc_rad_new_forcings,'Value',1);
    set(handles.forc_rad_creat_import,'Value',1);
    rad_creat_import(hObject, eventdata, handles);
else
    % Check whether the data has been moved
    % Especially Forcing
    if -exist(working_forcings)
        % Load information from model_forcings file
        load (working_forcings,...
            'starttime', 'endtime'  , 'start_str'   , 'end_str',...
            'days_step', 'mins_step', 'hours_step');
        
        rad_load_forcings(hObject, eventdata, handles)
        set(handles.forc_rad_load_forcings,'Value',1);
        fullpath_forcings = [pwd,strrep(working_forcings,'./','/')];
        set(handles.forc_txt_openforcings,'String',fullpath_forcings);
        
        %
        % Save information to working temporary variable file
        save './Temps/temp_variable.mat'...
            'starttime' 'endtime'   'start_str'     'end_str'...
            'days_step' 'mins_step' 'hours_step'    'working_forcings' -append;
        %
        % Display information from forcing files to the screen
        set(handles.txt_start_time,'String',start_str);
        set(handles.txt_end_time,'String',end_str);
        set(handles.txt_start_time,'Enable','inactive');
        set(handles.txt_end_time,'Enable','inactive');
        set(handles.forc_txt_step_days,'Enable','inactive');
        set(handles.forc_txt_step_hours,'Enable','inactive');
        set(handles.forc_txt_step_mins,'Enable','inactive');
        set(handles.forc_txt_step_days,'String',num2str(days_step));
        set(handles.forc_txt_step_hours,'String',num2str(hours_step));
        set(handles.forc_txt_step_mins,'String',num2str(mins_step));
    else
        uiwait(msgbox( 'Forcing data has been moved to another place. please relocate the forcing data','Error'));
    end
    

end
forc_table_root_CellEditCallback(hObject, eventdata, handles);

if isempty(working_forcings)
    set(handles.forc_but_view,'Enable','off');
else
    set(handles.forc_but_view,'Enable','on');
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes model_forcings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = model_forcings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function but_start_time_Callback(hObject, eventdata, handles)
try
    load './Temps/temp_variable.mat' 'starttime'
    [starttime, year_s, month_s, day_s, time_s] = uigetdate(now);
catch exception
    [starttime, year_s, month_s, day_s, time_s] = uigetdate(now);
end

if (isempty(year_s) || isempty(month_s) || isempty(day_s) || isempty(time_s))
    return
else
    start_str  = [num2str(month_s), ' - ', num2str(day_s), ' - ',num2str(year_s), '  |  ', num2str(time_s)];
    set(handles.txt_start_time,'String',start_str);
end
save './Temps/temp_variable.mat' 'starttime' 'start_str' 'year_s' 'month_s' 'day_s' 'time_s' -append;

function but_end_time_Callback(hObject, eventdata, handles)
try
    load './Temps/temp_variable.mat' 'endtime'
    [endtime, year_e, month_e, day_e, time_e] = uigetdate(endtime);
catch exception
    [endtime, year_e, month_e, day_e, time_e] = uigetdate(now);
end

if (isempty(year_e) || isempty(month_e) || isempty(day_e) || isempty(time_e))
    return
else
    end_str  = [num2str(month_e), ' - ', num2str(day_e), ' - ',num2str(year_e), '  |  ', num2str(time_e)];
    set(handles.txt_end_time,'String',end_str);
end
save './Temps/temp_variable.mat' 'endtime' 'end_str' 'year_e' 'month_e' 'day_e' 'time_e' -append;


function txt_start_time_Callback(hObject, eventdata, handles)


function txt_start_time_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [out,year_c,month_c,day_c,time_c] = uigetdate(varargin)
% UIGETDATE  date selection dialog box
%    T = UIGETDATE(D) displays a dialog box in form of a calendar 
%    
%    UIGETDATE expects serial date number or standard MATLAB Date 
%    format (see DATESTR) as input data und returns serial date number 
%    for the selected date and time.
%
%    UIGETDATE by itself uses the current date and time as input data
%
% Example:
%         t = datestr( uigetdate('16-Aug-1974 03:00') )
% 
% See also datevec, datestr, datenum

%   version: v1.0
%   author:  Elmar Tarajan [MCommander@gmx.de]

if nargin == 0
   varargin{1} = now;
end% if

if ~ishandle(varargin{1})
   %
   datvec = datevec(varargin{1});
   %
   scr = get(0,'ScreenSize');
   h.units = 'pixels';
   h.parent = figure(h,'menubar','none', ...
            'numbertitle','off', ...
            'resize','off', ...
            'handlevisibility','on', ...
            'visible','off', ...            
            'WindowStyle','modal', ...
            'Tag','uigetdate', ...
            'position',[ (scr(3:4)- [197 199])/2 197 199 ]);
   %
   pos = [5 5 0 0];
   uicontrol(h,'style','edit','position',pos+[0 0 104 26])
   uicontrol('style','slider','units','pixels','position',pos+[3 2 100 20], ...
             'sliderstep',[.0005 .0005],'min',-10,'max',10,'value',0, ...
             'callback','uigetdate(gcbo,''time'')','UserData',0)
   %
   h.style           = 'edit';
   h.fontweight      = 'bold';
   h.foregroundcolor = [.2 .2 .2];
   uicontrol(h,'enable','on','position',pos+[ 17 2 73 20],'Tag','time', ...
               'String',sprintf('%02d:%02d',datvec(4:5)))
   %
   % textbanners
   tmp = [2 20 101 4 ; 2 2 101 3 ; 2 2 3 22 ; 17 2 2 22 ; 88 2 2 22 ; 101 2 2 22 ];
   for i=1:6 ; uicontrol(h,'style','text','position',pos+tmp(i,:)) ; end% for
   %
   uicontrol(h,'style','edit','position',pos+[105 0 84 26],'visible','on')   
   uicontrol(h,'style','pushbutton','position',pos+[108 2 78 21],'Tag','ok', ...
               'CData',repmat(repmat([0.3:0.01:1 1:-0.01:0.3],18,1),[1 1 3]), ...
               'string','ok','Callback','uigetdate(gcbo,''ok'')')
   %
   pos = [5 32 0 0];
   uicontrol(h,'style','edit','position',pos+[0 0 189 136],'enable','inactive','Tag','cday', ...
      'UserData',datvec(3))   
   h.style      = 'pushbutton';
   h.fontweight = 'normal';
   for i=95:-19:0
      for j=0:26:156
         uicontrol(h,'position',pos+[j+3 i+2 27 20],'Enable','off', ...
                     'foregroundcolor',[.2 .2 .2],'Tag','day', ...
                     'callback','uigetdate(gcbo,''day'')')
      end% for
   end% for
   %
   tmp = {'Mon' 'Tue' 'Wed' 'Thu' 'Fri' 'Sat' 'Sun'};
   for j=0:6
      uicontrol(h,'style','text','position',pos+[j*26+4 119 25 13],'string',tmp{j+1}, ...
                  'backgroundcolor',[0.4 0.4 0.4],'foregroundcolor',[.9 .9 .9])         
   end% for
   %
   pos = [5 169 0 0];
   uicontrol(h,'style','edit','position',pos+[0 0 189 26])
   h.style = 'slider';
   uicontrol(h,'position',pos+[3 2 100 20],'sliderstep',[0.00025 1], ...
               'min',-2000,'max',2000,'Value',datvec(2), ...
               'callback','uigetdate(gcbo,''months'')')
   uicontrol(h,'position',pos+[112 2 74 20],'sliderstep',[0.00025 1], ...
               'min',0,'max',4000,'value',datvec(1), ...
               'callback','uigetdate(gcbo,''year'')')
   %
   h.style           = 'edit';
   h.enable          = 'inactive';
   h.fontweight      = 'bold';
   h.foregroundcolor = [.2 .2 .2];
   tmp = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun'... 
          'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
   uicontrol(h,'position',pos+[ 17 2 73 20],'Tag','months','String',tmp{datvec(2)},'Userdata',tmp)
   uicontrol(h,'position',pos+[126 2 47 20],'Tag','year','String',num2str(datvec(1)))
   %
   % textbanners
   h.style = 'text';
   tmp = [2 20 185 4 ; 2 2 185 3 ; 2 2 3 22 ; 17 2 2 22 ; 88 2 2 22 ; ...
      101 2 13 22 ; 126 2 2 22 ; 171 2 2 22 ; 184 2 3 22];
   for i=1:9
      uicontrol(h,'position',pos+tmp(i,:))
   end% for
   %
   set(h.parent,'visible','on')
   setday(varargin{1})
   %
   set(findobj(gcf,'string',num2str(datvec(3))),'CData',geticon)
   %
   uiwait
   try
      year_c    = get(findobj(gcf,'Tag','year'),'String');
      month_c   = get(findobj(gcf,'Tag','months'),'String');
      day_c     = get(findobj(gcf,'Tag','cday'),'UserData');
      time_c    = get(findobj(gcf,'Tag','time'),'String');
      
      out = datenum([num2str( ...
               get(findobj(gcf,'Tag','cday'),'UserData')) '-' ...
               get(findobj(gcf,'Tag','months'),'String') '-' ...
               get(findobj(gcf,'Tag','year'),'String') ' ' ...
               get(findobj(gcf,'Tag','time'),'String') ':00']);
      delete(findobj(0,'Tag','uigetdate'))
   catch
      out = [];
      %closereq
   end% try
   
   return
end% if

switch varargin{2}
   case 'months'
      h = findobj(gcbf,'Tag','months');
      months = get(h,'UserData');
      set(h,'String',months{mod(get(gcbo,'Value')-1,12)+1})
      set(findobj(gcbf,'Tag','ok'),'Enable','off')      
      %
   case 'year'
      set(findobj(gcbf,'Tag','year'),'String',get(gcbo,'Value'))
      set(findobj(gcbf,'Tag','ok'),'Enable','off')
      %
   case 'day'
      h= findobj(gcf,'Tag','day');
      set(h,'CData',[])

      set(varargin{1},'CData',geticon)
      set(findobj(gcbf,'Tag','cday'),'Userdata',get(varargin{1},'String'))
      set(findobj(gcbf,'Tag','ok'),'Enable','on')
      try ; uicontrol(h(3)) ; end% try
      return
      %
   case 'time'
      try
         if toc<0.1
            step = get(gcbo,'UserData');
            set(gcbo,'UserData',step+1)
            step = floor(step*sign(get(gcbo,'value'))/2);
         else
            set(gcbo,'UserData',1)
            step = sign(get(gcbo,'value'));
            set(gcbo,'value',0)
         end% if
         %
         handles.time = findobj(gcbf,'Tag','time');
         time = sum(sscanf(get(handles.time,'String'),'%d:%d').*[60;1]);
         time = time+step;
         if time<0
            time = 1439;
         elseif time>1439
            time = 0;
         end% if
         time = sprintf('%02.f:%02.f',floor(time/60),(time/60-floor(time/60))*60);
         set(handles.time,'String',time)
         %
         tic
         return
      catch
         tic
      end% try
      drawnow
      %
   case 'ok'
      uiresume
      return
      %
end% switch
setday(['1-' get(findobj(gcbf,'Tag','months'),'String') '-' ...
             get(findobj(gcbf,'Tag','year'),'String')])
  %
  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setday(datvec)
%-------------------------------------------------------------------------------
datvec = datevec(datvec);
datvec(3) = 1;
%
day = [7 1 2 3 4 5 6];
day = day(weekday(datestr(datvec)));
%
monthend = eomday(datvec(1),datvec(2));
%
ind = [zeros(1,42-monthend-day+1) monthend:-1:1 zeros(1,day-1)];
%
enable = repmat({'on'},42,1);
enable(ind==0) = {'off'};
%
count = strrep(strrep(cellstr(num2str(ind')),' 0',''),' ','');
%
h = findobj(gcf,'Tag','day');
set(h,{'String'},count,{'Enable'},enable,'backgroundcolor',[0.7 0.7 0.7],'CData',[])
set(h(ind~=0),'backgroundcolor',[.925 .922 .9002]);
set(h(ind~=0&repmat([1 1 0 0 0 0 0],1,6)),'backgroundcolor',[1 .8 .8])
  %
  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function icon = geticon
%-------------------------------------------------------------------------------
tmp = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 ;
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 ; ...
       1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 ; ...
       1 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 ; ...
       1 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 ; ...
       0 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 ; ...
       0 0 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 ; ...
       0 0 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 ; ...
       0 0 0 1 1 0 0 0 0 0 1 1 1 1 1 1 1 ; ...
       1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 ; ...
       1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 ; ...
       1 1 1 0 0 0 0 1 1 1 1 1 1 1 1 1 1 ; ...
       1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 ; ...
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ; ...
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ];
tmp(tmp==1)=NaN;
tmp(tmp==0)=1;
icon(:,:,1) = tmp;
tmp(tmp==1)=0.25;
icon(:,:,2) = tmp;
tmp(tmp==.25)=0;
icon(:,:,3) = tmp;


function forc_rad_new_forcings_Callback(hObject, eventdata, handles)


function forc_rad_load_forcings_Callback(hObject, eventdata, handles)


function forc_txt_openforcings_Callback(hObject, eventdata, handles)


function forc_txt_openforcings_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function forc_rad_new_forcings_CreateFcn(hObject, eventdata, handles)


function forc_rad_load_forcings_CreateFcn(hObject, eventdata, handles)


function forc_gbut_options_SelectionChangeFcn(hObject, eventdata)
 
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'forc_rad_new_forcings'
        rad_new_forcings(hObject, eventdata, handles);
        
        creat_blank=get(handles.forc_rad_creat_blank,'Value');
        if creat_blank
            rad_creat_blank(hObject, eventdata, handles);
        else
            rad_creat_import(hObject, eventdata, handles);
        end
        
    case 'forc_rad_load_forcings'
        rad_load_forcings(hObject, eventdata, handles);
    otherwise
        % Code for when there is no match.
end
%updates the handles structure
guidata(hObject, handles);


% --- Executes when selected object is changed in forc_gbut_options_create.
function forc_gbut_options_create_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in forc_gbut_options_create 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'forc_rad_creat_blank'
      rad_creat_blank(hObject, eventdata, handles);
      
    case 'forc_rad_creat_import'
      rad_creat_import(hObject, eventdata, handles);
    otherwise
       % Code for when there is no match.
end
%updates the handles structure
guidata(hObject, handles);


function rad_creat_blank(hObject, eventdata, handles)
set(handles.forc_but_create_mat_file,'Enable','on');
set(handles.forc_but_createimport_mat_file,'Enable','off');
%
set(handles.but_start_time,'Enable','on');
set(handles.txt_start_time,'Enable','on');
%
set(handles.but_end_time,'Enable','on');
set(handles.txt_end_time,'Enable','on');
%
set(handles.forc_txt_step_mins,'Enable','on');


function rad_creat_import(hObject, eventdata, handles)
set(handles.forc_but_create_mat_file,'Enable','off');
set(handles.forc_but_createimport_mat_file,'Enable','on');
%
set(handles.but_start_time,'Enable','off');
set(handles.txt_start_time,'Enable','off');
%
set(handles.but_end_time,'Enable','off');
set(handles.txt_end_time,'Enable','off');
%
set(handles.forc_txt_step_mins,'Enable','off');


function rad_new_forcings(hObject, eventdata, handles)
    set(handles.forc_txt_openforcings,'Enable','off');
    set(handles.forc_but_openforcings,'Enable','off');
    %
    set(handles.forc_txt_newforcings,'Enable','on');
    set(handles.forc_but_newforcings,'Enable','on');
    %
    set(handles.but_start_time,'Enable','on');
    set(handles.txt_start_time,'Enable','inactive');
    %
    set(handles.but_end_time,'Enable','on');
    set(handles.txt_end_time,'Enable','inactive');
    %
    set(handles.forc_txt_step_days,'Enable','off');
    set(handles.forc_txt_step_hours,'Enable','off');
    set(handles.forc_txt_step_mins,'Enable','on');
    %
    set(handles.forc_but_create_mat_file,'Enable','on');
    set(handles.forc_but_createimport_mat_file,'Enable','on');
    %
    set(handles.forc_rad_creat_import,'Enable','on');
    set(handles.forc_rad_creat_blank,'Enable','on');
      
    
function rad_load_forcings(hObject, eventdata, handles)
	set(handles.forc_txt_openforcings,'Enable','on');
    set(handles.forc_but_openforcings,'Enable','on');
    %
    set(handles.forc_txt_newforcings,'Enable','off');
    set(handles.forc_but_newforcings,'Enable','off');
    %
    set(handles.but_start_time,'Enable','off');
    set(handles.txt_start_time,'Enable','off');
    %
    set(handles.but_end_time,'Enable','off');
    set(handles.txt_end_time,'Enable','off');
    %
    set(handles.forc_txt_step_days,'Enable','off');
    set(handles.forc_txt_step_hours,'Enable','off');
    set(handles.forc_txt_step_mins,'Enable','off');
    %
    set(handles.forc_but_create_mat_file,'Enable','off');
	set(handles.forc_but_createimport_mat_file,'Enable','off');
    %
    set(handles.forc_rad_creat_import,'Enable','off');
    set(handles.forc_rad_creat_blank,'Enable','off');
    
    
function txt_end_time_Callback(hObject, eventdata, handles)


function txt_end_time_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function button_load_picture(hObject, eventdata, handles)
RGB_calendar=imread('./users/icons/calendar.png');
set(handles.but_start_time, 'Cdata', RGB_calendar);
set(handles.but_end_time, 'Cdata', RGB_calendar);
RGB_open=imread('./users/icons/open.png');
set(handles.forc_but_openforcings, 'Cdata', RGB_open);
set(handles.forc_but_newforcings, 'Cdata', RGB_open);


function but_start_time_CreateFcn(hObject, eventdata, handles)

function forc_but_newforcings_Callback(hObject, eventdata, handles)
[filename, pathname] = uiputfile('*.mat', 'Create forcing files');
if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
else
    fullpath_forcings   = [pathname, filename];
    workpath            = getappdata(0,'workpath');
    rem_path_forcings   = strrep(fullpath_forcings,workpath,'');
    if strcmp(rem_path_forcings,fullpath_forcings) == 1
        Msgbox('Do not create Forcings file outside the working folder','MLCan Error');
        return
    else
        working_forcings = ['.',rem_path_forcings];
        forcings_show = [pathname,filename];
        set(handles.forc_txt_newforcings,'String',forcings_show);
        setappdata(0,'working_forcings',working_forcings);
    end
end

function forc_but_openforcings_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.mat', 'Open forcing files');
if pathname == 0 %if the user pressed cancelled, then we exit this callback
    set(handles.forc_but_view,'Enable','off');
    return
else
	fullpath_forcings   = [pathname, filename];
    workpath            = getappdata(0,'workpath');
    rem_path_forcings   = strrep(fullpath_forcings,workpath,'');
    if strcmp(rem_path_forcings,fullpath_forcings) == 1
        Msgbox('Please copy forcings inside the working folder before loading','MLCan Error');
        return
    else
        working_forcings = ['.',rem_path_forcings];
        forcings_show = [pathname,filename];
        set(handles.forc_txt_openforcings,'String',forcings_show);
        save './Temps/temp_variable.mat' 'working_forcings' 'fullpath_forcings' -append;
        %
        % Load information from model_forcings file
        load (working_forcings,...
            'starttime', 'endtime'  , 'start_str'   , 'end_str',...
            'days_step', 'mins_step', 'hours_step');
        %
        % Save information to working temporary variable file
        save './Temps/temp_variable.mat'...
            'starttime' 'endtime'   'start_str'     'end_str'...
            'days_step' 'mins_step' 'hours_step'    'working_forcings' 'fullpath_forcings' -append;
        %
        % Display information from forcing files to the screen
        set(handles.txt_start_time,'String',start_str);
        set(handles.txt_end_time,'String',end_str);
        set(handles.txt_start_time,'Enable','inactive');
        set(handles.txt_end_time,'Enable','inactive');
        set(handles.forc_txt_step_days,'Enable','inactive');
        set(handles.forc_txt_step_hours,'Enable','inactive');
        set(handles.forc_txt_step_mins,'Enable','inactive');    
        set(handles.forc_txt_step_days,'String',num2str(days_step));
        set(handles.forc_txt_step_hours,'String',num2str(hours_step));
        set(handles.forc_txt_step_mins,'String',num2str(mins_step));
        
        set(handles.forc_but_view,'Enable','on');
        Message_forcings = ['Forcing file ', working_forcings, ' was loaded'];
        msgbox(Message_forcings,'Message','warn');
    end
end


function forc_but_create_mat_file_Callback(hObject, eventdata, handles)
button = questdlg('Create a new blank forcing file?',...
'Continue Operation','Yes','No','No');
if strcmp(button,'Yes')
    working_forcings = getappdata(0,'working_forcings');
    load './Temps/temp_variable.mat' 'lat_face' 'long_face' 'elev_face'...
        'starttime' 'endtime' 'start_str' 'end_str';
    days_str    = get(handles.forc_txt_step_days,'String');
    hours_str   = get(handles.forc_txt_step_hours,'String');
    mins_str    = get(handles.forc_txt_step_mins,'String');

    days_step   = str2num(days_str);
    hours_step  = str2num(hours_str);
    mins_step   = str2num(mins_str);

    if days_step > endtime-starttime
        msgbox('Time step is greater than simulation time');
        return
    elseif (hours_step < 0 || hours_step>24)
        msgbox('In correct hour');
        return
    elseif (mins_step<0 || mins_step>60)
        msgbox('In correct minute');
        return
    else
        time_step = days_step + hours_step/24 + mins_step/1440;
    end
        check_integer = (endtime - starttime)/time_step;
    if mod(abs(check_integer-round(check_integer)),1) > 1e-5;
        assignin('base','endtime',endtime);
        assignin('base','starttime',starttime);
        assignin('base','time_step',time_step);
        msgbox('Check time information','MLCan Error');
        return
    else
        num_step = round(check_integer)+1;
    end

    if starttime >= endtime % check starting and ending time
        msgbox('Starting time can not be sooner than ending time','MLCan Error','Error');
        return
    else
        id                = [1:1:num_step]';
        LAT               = lat_face;
        LONG              = long_face;
        ELEV              = elev_face;
        ea_crop           = ones(num_step,0);
        %LAI_crop     = ones(num_step,0);
        %LWin_crop    = ones(num_step,0);
        LAI_species1_crop = ones(num_step,0);
        LAI_species2_crop = ones(num_step,0);
        LAI_species3_crop = ones(num_step,0);
        LAI_species4_crop = ones(num_step,0);
        Pa_crop           = ones(num_step,0);
        PPT_crop          = ones(num_step,0);
        Rg_crop           = ones(num_step,0);
        Ta_crop           = ones(num_step,0);
        U_crop            = ones(num_step,0);
        ustar_crop        = ones(num_step,0);
        VPD_crop          = ones(num_step,0);
        RH_crop           = ones(num_step,0);
        
        % Assign value of year.
        time_simu = starttime;
        for i = 1:num_step
            time_str  = datevec(time_simu);
            year_crop(i)    = time_str(1);
            decyear_crop(i) = year_crop(1) + (i-1)*time_step/yearday(year_crop(i));
            day_crop(i)     = floor(date2doy(time_simu));
            doy_crop(i)     = date2doy(time_simu);
            month_crop(i)   = month(doy_crop(i));
            hour_crop(i)    = time_str(4)+time_str(5)/60;
            time_simu       = time_simu + time_step;
        end
        year_crop       = year_crop';
        decyear_crop   	= decyear_crop';
        day_crop        = day_crop';
        doy_crop        = doy_crop';
        month_crop      = month_crop';
        hour_crop      	= hour_crop';

        ZEN_crop        = ZEN_calculation(lat_face,long_face,hour_crop,doy_crop);

%         save (working_forcings,...
%             'starttime'     , 'endtime'     , 'start_str'   , 'end_str'     ,                                   ...
%             'LAT'           , 'LONG'        , 'ELEV'        ,                                                   ...
%             'days_step'     , 'hours_step'  , 'mins_step'   , 'num_step'    , 'time_step'   , 'id'          ,   ...
%             'ea_crop'       , 'LAI_crop'    , 'LWin_crop'   , 'Pa_crop'     , 'PPT_crop'    , 'Rg_crop'     ,   ...
%             'Ta_crop'       , 'U_crop'      , 'ustar_crop'  , 'VPD_crop'    , 'ZEN_crop'    , 'year_crop'   ,   ...          
%             'decyear_crop'  , 'doy_crop'    , 'decdoy_crop' , 'hour_crop'   );  

        save (working_forcings,...
            'starttime'         , 'endtime'           , 'start_str'         , 'end_str'           ,                        ...
            'days_step'         , 'hours_step'        , 'mins_step'         , 'num_step'          , 'time_step'         ,  ...
            'id'                , 'year_crop'         , 'doy_crop'          , 'hour_crop'         , 'day_crop'          ,  ...
            'month_crop'        , 'LAI_species1_crop' , 'LAI_species2_crop' , 'LAI_species3_crop' , 'LAI_species4_crop' ,  ...
            'Ta_crop'           , 'U_crop'            , 'Rg_crop'           , 'PPT_crop'          , 'VPD_crop'          , ...
            'Pa_crop'           , 'ea_crop'           , 'RH_crop'           , 'ZEN_crop'          , 'ustar_crop');

        Message_forcings = ['Forcing file ', working_forcings, ' was created'];
        msgbox(Message_forcings,'Message','warn');
        
        set(handles.forc_but_view,'Enable','on');
    end
elseif strcmp(button,'No')
    return
end


function forc_but_createimport_mat_file_Callback(hObject, eventdata, handles)
% [filename, pathname] = uigetfile(...
%     {'*.xls;*.xlsx'     , 'Microsoft Excel Files (*.xls,*.xlsx)';   ...
%     '*.csv'             , 'CSV - comma delimited (*.csv)'       ;   ... 
%     '*.txt'             , 'Text (Tab Delimited (*.txt)'         ;   ...
%     '*.*'               , 'All Files (*.*)'                     },  ...
%     'Import data from files');
[filename, pathname] = uigetfile(...
    {'*.xls;*.xlsx'     , 'Microsoft Excel Files (*.xls,*.xlsx)';   ...
    '*.*'               , 'All Files (*.*)'                     },  ...
    'Import data from files');

if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
else     
%     working_forcings = getappdata(0,'working_forcings');
%     load './Temps/temp_variable.mat' 'lat_face' 'long_face' 'elev_face'...
%         'starttime' 'endtime' 'start_str' 'end_str';
% 
%     days_str    = get(handles.forc_txt_step_days,'String');
%     hours_str   = get(handles.forc_txt_step_hours,'String');
%     mins_str    = get(handles.forc_txt_step_mins,'String');
% 
%     days_step   = str2num(days_str);
%     hours_step  = str2num(hours_str);
%     mins_step   = str2num(mins_str);
% 
%     if days_step > endtime-starttime
%         msgbox('Time step is greater than simulation time');
%         return
%     elseif (hours_step < 0 || hours_step>24)
%         msgbox('In correct hour');
%         return
%     elseif (mins_step<0 || mins_step>60)
%         msgbox('In correct minute');
%         return
%     else
%         time_step = days_step + hours_step/24 + mins_step/1440;
%     end
%         check_integer = (endtime - starttime)/time_step;
%     if mod(abs(check_integer-round(check_integer)),1) > 1e-5;
%         msgbox('Check time information','MLCan Error');
%         return
%     else
%         num_step = round(check_integer)+1;
%         assignin('base','num_step',num_step);
%     end
%     
%     str_file = regexp(filename, '\.', 'split');
%     if (strcmp(str_file(end),'xls') == 1 || strcmp(str_file(end),'xlsx') == 1)
%         data_table = xlsread([pathname,filename]);      
%     elseif strcmp(str_file(end),'csv') == 1
%         data_table = csvread([pathname,filename], 1, 0);      
%     elseif strcmp(str_file(end),'txt') == 1
%         data_table = dlmread([pathname,filename], '\t', 1, 0);      
%     else
%         msgbox('The file extension is unrecognized, please choose another file','MLCan Error','Error');
%         return
%     end
%     
%     if starttime >= endtime % check starting and ending time
%         msgbox('Starting time can not be sooner than ending time','Error');
%         return
%     else
%         if length(data_table(:,1)) ~= num_step
%             numstr=num2str(length(data_table(:,1)));
%             notice = ['Number of timestep (',num2str(num_step),') must be equal to the number of rows in data file (',numstr,')'];
%             msgbox(notice, 'MLCan Error: Different data length','Error');
%             return        
%         else
%             id           = [1:1:num_step]';
%             LAT          = lat_face;
%             LONG         = long_face;
%             ELEV         = elev_face;
%             LAI_crop     = data_table(:,1);
%             LWin_crop    = data_table(:,2);
%             PPT_crop     = data_table(:,3);
%             Pa_crop      = data_table(:,4);            
%             Rg_crop      = data_table(:,5);
% 
%             Ta_crop      = data_table(:,6);
%             U_crop       = data_table(:,7);
%             VPD_crop     = data_table(:,8);
% 
%             ea_crop      = data_table(:,9);
%             ustar_crop   = data_table(:,10);
%             % Assign value of year.
%             time_simu = starttime;
%             for i = 1:num_step
%                 time_str  = datevec(time_simu);
%                 year_crop(i)    = time_str(1);
%                 decyear_crop(i) = year_crop(1) + (i-1)*time_step/yearday(year_crop(i));
%                 doy_crop(i)     = floor(date2doy(time_simu));
%                 decdoy_crop(i)  = date2doy(time_simu);
%                 hour_crop(i)    = time_str(4)+time_str(5)/60;
%                 time_simu       = time_simu + time_step;
%             end
%             year_crop       = year_crop';
%             decyear_crop   	= decyear_crop';
%             doy_crop        = doy_crop';
%             decdoy_crop     = decdoy_crop';
%             hour_crop      	= hour_crop';
% 
%             ZEN_crop        = ZEN_calculation(LAT, LONG, doy_crop,hour_crop);
% 
%             save (working_forcings,...
%                 'starttime'     , 'endtime'     , 'start_str'   , 'end_str'     ,                                   ...
%                 'LAT'           , 'LONG'        , 'ELEV'        ,                                                   ...
%                 'days_step'     , 'hours_step'  , 'mins_step'   , 'num_step'    , 'time_step'   , 'id'          ,   ...
%                 'ea_crop'       , 'LAI_crop'    , 'LWin_crop'   , 'Pa_crop'     , 'PPT_crop'    , 'Rg_crop'     ,   ...
%                 'Ta_crop'       , 'U_crop'      , 'ustar_crop'  , 'VPD_crop'    , 'ZEN_crop'    , 'year_crop'   ,   ...
%                 'decyear_crop'  , 'doy_crop'    , 'decdoy_crop' , 'hour_crop'   );  
% 
%             Message_forcings = ['Forcing file ', working_forcings, ' was created and imported'];
%             msgbox(Message_forcings,'MLCan message','warn');
%         end
%     end


    working_forcings = getappdata(0,'working_forcings');
    load './Temps/temp_variable.mat' 'lat_face' 'long_face' 'elev_face';

    str_file = regexp(filename, '\.', 'split');
    MsgBOX=msgbox('Processing, please wait','MLCan Notice','warn');
    if (strcmp(str_file(end),'xls') == 1 || strcmp(str_file(end),'xlsx') == 1)
        data_table = xlsread([pathname,filename],'Canopy');      
    elseif strcmp(str_file(end),'csv') == 1
        data_table = csvread([pathname,filename], 1, 0);      
    elseif strcmp(str_file(end),'txt') == 1
        data_table = dlmread([pathname,filename], '\t', 1, 0);      
    else
        close(MsgBOX);
        msgbox('The file extension is unrecognized, please choose another file','MLCan Error','Error');
        return
    end
    close(MsgBOX);

    if (strcmp(str_file(end),'xls') == 1 || strcmp(str_file(end),'xlsx') == 1)
        data_table_RB = xlsread([pathname,filename],'Nutrient_RootBiomass');
    else
        msgbox('The file extension is unrecognized, please choose another file','MLCan Error','Error');
        return
    end
    if (strcmp(str_file(end),'xls') == 1 || strcmp(str_file(end),'xlsx') == 1)
        data_table_DEM = xlsread([pathname,filename],'Nutrient_DEM');
    else
        msgbox('The file extension is unrecognized, please choose another file','MLCan Error','Error');
        return
    end
    
    if size(data_table,2) ~= 17
        msgbox('The imported file (Canopy) does not include all the data needed','MLCan Error','Error');
        return
    elseif size(data_table_RB,2) ~= 12
        msgbox('The imported file (Nutrient_RootBiomass) does not include all the data needed','MLCan Error','Error');
        return
    elseif size(data_table_DEM,2) ~= 20
        msgbox('The imported file (Nutrient_DEM) does not include all the data needed','MLCan Error','Error');
        return
    else
        
        
        id                  = data_table(:,1);
        year_crop           = data_table(:,2);
        month_crop          = data_table(:,3);
        day_crop            = data_table(:,4);
        doy_crop            = data_table(:,5); 
        
        hour_crop           = data_table(:,6);        
        Ta_crop             = data_table(:,7);
        Pa_crop             = data_table(:,8);
        Rg_crop             = data_table(:,9);
        PPT_crop            = data_table(:,10);

        VPD_crop            = data_table(:,11);
        U_crop              = data_table(:,12);
        ustar_crop          = data_table(:,13);
        
        LAI_species1_crop   = data_table(:,14);
        LAI_species2_crop   = data_table(:,15);
        LAI_species3_crop   = data_table(:,16);
        LAI_species4_crop   = data_table(:,17);
        
        [ea_crop,RH_crop]   = Ta_VPD_calculation(Ta_crop,VPD_crop);
        %zen_crop          = data_table(:,11);
        [ZEN_crop] = ZEN_calculation(lat_face,long_face,hour_crop,doy_crop);
        % 
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        BMla_species1_crop    = data_table_RB(:,1);
        BMla_species2_crop    = data_table_RB(:,2);
        BMla_species3_crop    = data_table_RB(:,3);
        BMla_species4_crop    = data_table_RB(:,4);
        
        CNratio_species1_crop = data_table_RB(:,5);
        CNratio_species2_crop = data_table_RB(:,6);
        CNratio_species3_crop = data_table_RB(:,7);
        CNratio_species4_crop = data_table_RB(:,8);
        
        LMI_species1_crop     = data_table_RB(:,9);
        LMI_species2_crop     = data_table_RB(:,10);
        LMI_species3_crop     = data_table_RB(:,11);
        LMI_species4_crop     = data_table_RB(:,12);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        BMlaDEM_species1_crop  = data_table_DEM(:,1);
        BMlaDEM_species2_crop  = data_table_DEM(:,2);
        BMlaDEM_species3_crop  = data_table_DEM(:,3);
        BMlaDEM_species4_crop  = data_table_DEM(:,4);
        
        BMrbDEM_species1_crop  = data_table_DEM(:,5);
        BMrbDEM_species2_crop  = data_table_DEM(:,6);
        BMrbDEM_species3_crop  = data_table_DEM(:,7);
        BMrbDEM_species4_crop  = data_table_DEM(:,8);
        
        CNratioA_species1_crop  = data_table_DEM(:,9);
        CNratioA_species2_crop  = data_table_DEM(:,10);
        CNratioA_species3_crop  = data_table_DEM(:,11);
        CNratioA_species4_crop  = data_table_DEM(:,12);       

        CNratioB_species1_crop  = data_table_DEM(:,13);
        CNratioB_species2_crop  = data_table_DEM(:,14);
        CNratioB_species3_crop  = data_table_DEM(:,15);
        CNratioB_species4_crop  = data_table_DEM(:,16);       
        
        NDEM_species1_crop      = data_table_DEM(:,17);
        NDEM_species2_crop      = data_table_DEM(:,18);
        NDEM_species3_crop      = data_table_DEM(:,19);
        NDEM_species4_crop      = data_table_DEM(:,20);
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


        % Start time
%         isleap = @(x) (mod(x,4)==0 & mod(x,100)~=0) | mod(x,400) == 0;
%         leapyrs = (isleap(year_crop(1)));
%         if leapyrs == 1
%             day_start_here = day(day_crop(1));
%         else
%             if day_crop(1) > 59
%                 day_start_here = day(day_crop(1)+1);
%             else
%                 day_start_here = day(day_crop(1));
%             end
%         end
        day_start_here = (day_crop(1));
        
        starttime=datenum(year_crop(1),month_crop(1),day_start_here);
        month_st=datestr(starttime);
        month_str=month_st(4:6);
        hour_str=num2str(floor(hour_crop(1)));
        if str2num(hour_str) < 10
            start_str=[month_str, ' - ' num2str(day_start_here) ' - ' num2str(year_crop(1))...
                ' | ' '0' hour_str ':' num2str(mod(hour_crop(1),1)*60)];
        elseif str2num(hour_str) >= 10
            start_str=[month_str, ' - ' num2str(day_start_here) ' - ' num2str(year_crop(1))...
                ' | ' hour_str ':' num2str(mod(hour_crop(1),1)*60)];
        end
        if start_str(end-1) == ':'
            if num2str(mod(hour_crop(end),1)*60) == 0
                start_str(end-1:end+1)=[':0' num2str(mod(hour_crop(end),1)*60)];                                
            else
                start_str(end:end+1)=[num2str(mod(hour_crop(end),1)*60)];                
            end
        end
        
        % End time
%         isleap = @(x) (mod(x,4)==0 & mod(x,100)~=0) | mod(x,400) == 0;
%         leapyrs = (isleap(year_crop(end)));
%         if leapyrs == 1
%             day_end_here = day(day_crop(end));
%         else
%             if day_crop(end) > 59
%                 day_end_here = day(day_crop(end)+1);
%             else
%                 day_end_here = day(day_crop(end));
%             end
%         end
        day_end_here = (day_crop(end));
        
        endtime=datenum(year_crop(end),month_crop(end),day_end_here);
        month_st=datestr(endtime);
        month_str=month_st(4:6);
        hour_str=num2str(floor(hour_crop(end)));
        if str2num(hour_str) < 10
            end_str=[month_str, ' - ' num2str(day_end_here) ' - ' num2str(year_crop(end))...
                ' | ' '0' hour_str ':' num2str(mod(hour_crop(end),1)*60)];
        elseif str2num(hour_str) >= 10
            end_str=[month_str, ' - ' num2str(day_end_here) ' - ' num2str(year_crop(end))...
                ' | ' hour_str ':' num2str(mod(hour_crop(end),1)*60)];
        end
        if end_str(end-1) == ':'
            if (mod(hour_crop(end),1)*60) == 0
                end_str(end-1:end+1)=[':0' num2str(mod(hour_crop(end),1)*60)];                
            else
                end_str(end-1:end+1)=[':' num2str(mod(hour_crop(end),1)*60)];
            end
        end
        
        % time step
        days_step_test=[doy_crop(2)-doy_crop(1) doy_crop(3)-doy_crop(2) doy_crop(4)-doy_crop(3)];
        days_step=floor(mode(days_step_test));
        hours_step_test=[hour_crop(2)-hour_crop(1) hour_crop(3)-hour_crop(2) hour_crop(4)-hour_crop(3)];
        hours_step=floor(mode(hours_step_test));
        mins_step_test=[hour_crop(2)-hour_crop(1) hour_crop(3)-hour_crop(2) hour_crop(4)-hour_crop(3)];
        mins_step=mod((mode(mins_step_test))*60,60);
        
        num_step=id(end);
        time_step=mode(days_step_test);
        
        % Display information from forcing files to the screen
        set(handles.txt_start_time,'String',start_str);
        set(handles.txt_end_time,'String',end_str);
        set(handles.txt_start_time,'Enable','inactive');
        set(handles.txt_end_time,'Enable','inactive');
        set(handles.forc_txt_step_days,'Enable','inactive');
        set(handles.forc_txt_step_hours,'Enable','inactive');
        set(handles.forc_txt_step_mins,'Enable','inactive');
        set(handles.forc_txt_step_days,'String',num2str(days_step));
        set(handles.forc_txt_step_hours,'String',num2str(hours_step));
        set(handles.forc_txt_step_mins,'String',num2str(mins_step));
        
        if isempty(working_forcings)
            msgbox('Not enough information, check the data path','MLCan Error: File format');
            return
        else
            save (working_forcings,...
                'starttime'         , 'endtime'           , 'start_str'         , 'end_str'           ,                        ...
                'days_step'         , 'hours_step'        , 'mins_step'         , 'num_step'          , 'time_step'         ,  ...
                'id'                , 'year_crop'         , 'doy_crop'          , 'hour_crop'         , 'day_crop'          ,  ...
                'month_crop'        , 'LAI_species1_crop' , 'LAI_species2_crop' , 'LAI_species3_crop' , 'LAI_species4_crop' ,  ...
                'Ta_crop'           , 'U_crop'            , 'Rg_crop'           , 'PPT_crop'          , 'VPD_crop'          ,  ...
                'Pa_crop'           , 'ea_crop'           , 'RH_crop'           , 'ZEN_crop'          , 'ustar_crop'        ,  ...
                'BMla_species1_crop', 'BMla_species2_crop', 'BMla_species3_crop', 'BMla_species4_crop', ...
                'CNratio_species1_crop', 'CNratio_species2_crop', 'CNratio_species3_crop', 'CNratio_species4_crop', ... 
                'LMI_species1_crop' , 'LMI_species2_crop' , 'LMI_species3_crop' , 'LMI_species4_crop' , ...
                'BMlaDEM_species1_crop', 'BMlaDEM_species2_crop','BMlaDEM_species3_crop','BMlaDEM_species4_crop',...
                'BMrbDEM_species1_crop','BMrbDEM_species2_crop','BMrbDEM_species3_crop','BMrbDEM_species4_crop', ...
                'CNratioA_species1_crop','CNratioA_species2_crop','CNratioA_species3_crop','CNratioA_species4_crop', ...
                'CNratioB_species1_crop','CNratioB_species2_crop','CNratioB_species3_crop','CNratioB_species4_crop', ...
                'NDEM_species1_crop','NDEM_species2_crop','NDEM_species3_crop','NDEM_species4_crop');
              
            Message_forcings = ['Forcing file ', working_forcings, ' was created and imported'];
            set(handles.forc_but_view,'Enable','on');
            
            workpath            = getappdata(0,'workpath');
            fullpath_forcings   = [workpath working_forcings(2:end)];
            save './Temps/temp_variable.mat' 'working_forcings' 'fullpath_forcings' -append;
            
            msgbox(Message_forcings,'MLCan message','warn');
        end
    end
end


function forc_txt_newforcings_Callback(hObject, eventdata, handles)


function forc_txt_newforcings_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function forc_but_ok_Callback(hObject, eventdata, handles)

root_init           =  get(handles.forc_table_root,'Data');
root_init_litter    =  get(handles.forc_table_root_litter,'Data');
nutrient_int =  get(handles.forc_table_biogeochemical,'Data');
nutrient_int_litter =  get(handles.forc_table_biogeochemical_litter,'Data');
if sum((root_init(:,2))<=0) > 0
    msgbox( 'Soil moisture should not be <= 0 - Please go back to ICs FOR CANOPY','MLCan Error','Error');
    return
else
     load './Temps/temp_variable.mat' 'Soil_nutrient'
    if ((sum(sum(nutrient_int(:,2:end)))==0 && Soil_nutrient == 1)||(sum(sum(nutrient_int_litter(:,2:end)))==0 && Soil_nutrient == 1))
        msgbox( 'Please set initial conditions for nutrient','MLCan Error','Error');
        return
    else
        save './Temps/temp_variable.mat' 'root_init' 'root_init_litter' 'nutrient_int' 'nutrient_int_litter' -append
        
        load ('./Temps/temporary.mat',...
            'working_name'  )
        copyfile('./Temps/temp_variable.mat',working_name,'f');
        close
    end
end
function forc_but_cancel_Callback(hObject, eventdata, handles)
close

function forc_txt_step_days_Callback(hObject, eventdata, handles)


function forc_txt_step_days_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function forc_txt_step_mins_Callback(hObject, eventdata, handles)


function forc_txt_step_mins_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function forc_txt_step_hours_Callback(hObject, eventdata, handles)


function forc_txt_step_hours_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function axes1_CreateFcn(hObject, eventdata, handles)
forceiconsmall = imread('./users/icons/forcings_small.png');
image(forceiconsmall);
axis off
%imshow('./users/icons/forcings_small.png');


function forc_tabbut_forcing_N_Callback(hObject, eventdata, handles)
set(handles.forc_pan_condition,'Visible','off');
set(handles.forc_pan_forcing,'Visible','on');
set(handles.forc_pan_nutrient,'Visible','off');
set(handles.forc_pan_condition,'Visible','off');
set(handles.forc_but_load,'Enable','off');
set(handles.forc_but_load_nutrient,'Enable','off');
set(handles.forc_gbut_options_create,'Visible','on');
load './Temps/temp_variable.mat'...
    'working_forcings'
if isempty(working_forcings)
    set(handles.forc_but_view,'Enable','off');
else
    set(handles.forc_but_view,'Enable','on');
end

function forc_tabbut_condition_Callback(hObject, eventdata, handles)
set(handles.forc_pan_condition,'Visible','on');
set(handles.forc_pan_forcing,'Visible','off');
set(handles.forc_pan_nutrient,'Visible','off');
set(handles.forc_but_load,'Enable','on');
set(handles.forc_but_load,'Visible','on');
set(handles.forc_but_load_nutrient,'Enable','off');
set(handles.forc_but_load_nutrient,'Visible','off');
set(handles.forc_gbut_options_create,'Visible','off');
set(handles.forc_but_view,'Enable','off');

% --- Executes on button press in forc_tabbut_nutrient.
function forc_tabbut_nutrient_Callback(hObject, eventdata, handles)
set(handles.forc_pan_condition,'Visible','off');
set(handles.forc_pan_forcing,'Visible','off');
set(handles.forc_pan_nutrient,'Visible','on');
set(handles.forc_but_load,'Enable','off');
set(handles.forc_but_load,'Visible','off');
set(handles.forc_but_load_nutrient,'Enable','on');
set(handles.forc_but_load_nutrient,'Visible','on');
set(handles.forc_gbut_options_create,'Visible','off');
set(handles.forc_but_view,'Enable','off');

function forc_table_canopy_CellEditCallback(hObject, eventdata, handles)


%function forc_but_view_mat_file_Callback(hObject, eventdata, handles)


function forc_table_root_CellEditCallback(hObject, eventdata, handles)
root_ic = get(handles.forc_table_root,'Data');
soil_depth = root_ic(:,1);
soil_moisture = root_ic(:,2);
soil_temperature = root_ic(:,3);

if get(handles.forc_table_root,'ColumnEditable') == logical([1 1 0])
    plot(handles.forc_axes_ic,soil_moisture,-soil_depth,'Color', 'b',...
            'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
            'MarkerFaceColor','g','MarkerSize',6);
    axis(handles.forc_axes_ic,[min(soil_moisture)/1.2 max(soil_moisture)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
    ylabel(handles.forc_axes_ic,'Soil depth [m]');
    %xlabel(handles.forc_axes_ic,'Soil moisture [-]');
    legend(handles.forc_axes_ic,'Soil moisture [-]');
    grid on;
   
else
    plot(handles.forc_axes_ic,soil_moisture,-soil_depth,'Color', 'b',...
            'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
            'MarkerFaceColor','g','MarkerSize',6);
    axis(handles.forc_axes_ic,[min(soil_moisture)/1.2 max(soil_moisture)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
    ylabel(handles.forc_axes_ic,'Soil depth [m]');
    %xlabel(handles.forc_axes_ic,'Soil moisture [-]');
    legend(handles.forc_axes_ic,'Soil moisture [-]');
    grid on;
    
    plot(handles.forc_axes_ic2,soil_temperature,-soil_depth,'Color', 'r',...
            'LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b',...
            'MarkerFaceColor','m','MarkerSize',6);
    axis(handles.forc_axes_ic2,[min(soil_temperature)/1.2 max(soil_temperature)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
    ylabel(handles.forc_axes_ic2,'Soil depth [m]');
    %xlabel(handles.forc_axes_ic2,'Temperature [C]');
    legend(handles.forc_axes_ic2,'Soil temperature [C]');
    grid on;    
        
    %[AX,H1,H2] = plotyy(handles.forc_axes_ic,soil_depth,soil_moisture,soil_depth,soil_temperature,'plot');
    %axis(AX(1),[0 max(soil_depth)+1 min(soil_moisture)/1.2 max(soil_moisture)*1.1 + 0.1]);
    %axis(AX(2),[0 max(soil_depth)+1 min(soil_temperature)/1.2 max(soil_temperature)*1.1+1]);
    %xlabel('Root depth [m]');
    %set(AX(1),'ycolor','b');
    %set(AX(2),'ycolor','r');
    %set(H1, 'Color', 'b','LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6);
    %set(H2, 'Color', 'r','LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','m','MarkerSize',6);
    %legend('Soil moisture [-]','Soil temperature [C]');
    %grid on;
end

% --- Executes when entered data in editable cell(s) in forc_table_biogeochemical.
function forc_table_biogeochemical_CellEditCallback(hObject, eventdata, handles)
load './Temps/temp_variable.mat' 'Soil_C_pool1'
nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
soil_depth = nutrient_ic(:,1);
if Soil_C_pool1 == 1 % 1 soil C pool
    MicrobialC = nutrient_ic(:,3);
%    MicrobialDR = nutrient_ic(:,8);
else % 2 soil C pools
    MicrobialC = nutrient_ic(:,4);
%    MicrobialDR = nutrient_ic(:,10);
end
   
plot(handles.forc_axes_ic3,MicrobialC,-soil_depth,'Color', 'b',...
    'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',6);
axis(handles.forc_axes_ic3,[min(MicrobialC)/1.2 max(MicrobialC)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
ylabel(handles.forc_axes_ic3,'Soil depth [m]');
%xlabel(handles.forc_axes_ic3,'Microbial Carbon [gC/m^{-3}]');
legend(handles.forc_axes_ic3,'Microbial Carbon [gC/m^{-3}]');
grid on;

% plot(handles.forc_axes_ic4,MicrobialDR,-soil_depth,'Color', 'r',...
%     'LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b',...
%     'MarkerFaceColor','m','MarkerSize',6);
% axis(handles.forc_axes_ic4,[min(MicrobialDR)/1.2 max(MicrobialDR)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
% ylabel(handles.forc_axes_ic4,'Soil depth [m]');
% %xlabel(handles.forc_axes_ic4,'Microbial death rate [/t]');
% legend(handles.forc_axes_ic4,'Microbial death rate [/t]');
% grid on;



% % % 
% function [Zen_deg] = ZEN_calculation(LAT,LONG,hour_crop,doy_crop)
% doy_crop=doys;
% hour_crop=hours;
% sin_del = 0.39785*sin((278.97+0.9856*doy_crop+1.9165*sin((356.6+0.9856*doy_crop)*pi/180))*pi/180);
% f = (279.575 + 0.9856*doy_crop)*pi/180;
% ET = (-104.7*sin(f)+596.2*sin(2*f)+4.3*sin(3*f)-12.7*sin(4*f)-429.3*cos(f)-2.0*cos(2*f)+19.3*cos(3*f))/3600;
% cos_Zen = sin(LAT*pi/180)*sin_del + cos(LAT*pi/180)*(1-sin_del.^2).^(1/2).*cos(15*(hour_crop-12+(-LONG-75)/15+ET)*pi/180);
% Zen_cal = acos(cos_Zen);
% Zen_deg = Zen_cal*180/pi;

%function [zen, dec] = zenith_angle(LAT, LONG, hours, doys)
function [zen] = ZEN_calculation(LAT,LONG,hours,doys)
%
%   Compute the zenith angle for a particular site and range of times
%   
%   Taken From Campbell and Norman, Eqns 11.1 - 11.5
%   (An Introduction to Environmental Biophysics, p168-170)
%   Written By : Darren Drewry
%   Modified By: Dongkook Woo
%
%   INPUTS:
%       LAT - latitude of site [degrees]
%       hours - decimal hours (0..24)
%       doys - integer day of year
%       LONG   - Longitude [degrees]
%
%  OUTPUTS:
%       zen  - Zenith angle [degrees]
%       dec  - Declination angle [degrees]
%

d2r = pi/180;   % Conversion factor from degrees to radians
r2d = 180/pi;   % Conversion factor from radians to degrees

% Dongkook Woo - Edit
% doy should start with 1
doys=floor(round((doys+1).*10^10)./10^10);
% Dongkook Woo - Edit End

% Solar Declination Angle [rad] - ranges from -23.45 degrees to + 23.45 degrees
d1 = d2r*(356.6 + 0.9856*doys); % [rad]
d2 = sin(d1);                   % [rad]
d3 = 1.9165 * d2;               % [deg]
d4 = 278.97 + 0.9856*doys + d3; % [deg]
d5 = sin( d2r*d4 );             % [rad]
d6 = 0.39785 * d5;              % [deg]
dec = asin( d6 );

ff = d2r*(279.575 + 0.9856*doys);    % [rad]
ET = (-104.7*sin(ff) + 596.2*sin(2*ff) + 4.3*sin(3*ff) - 12.7*sin(4*ff) - 429.3*cos(ff) - 2.0*cos(2*ff) + 19.3*cos(3*ff)) / 3600;

% Compute the longitude correction 1/15 grades
% Dongkook Woo - Edit
%LC=LONG*(1/15);
n=ceil(LONG./15);
LC=(15*n-LONG)*(1/15);
% Dongkook Woo - Edit End

t0 = 12 - LC - ET;              %[hours] 

LAT = d2r*LAT;  % [rad]
zen = acos( sin(LAT)*sin(dec) + cos(LAT)*cos(dec).*cos(d2r*15*(hours-t0)) );    %[Zenith Angle]

zen = zen * r2d;  % [degrees]
dec = dec * r2d;  % [degrees]



function [ea,RH] = Ta_VPD_calculation(Ta,VPD)
% This function computes the atmospheric vapor pressure and also the
% relative humidity using atmospheric temperature and water vapor deficit

% Taken from An Introduction to Environmental Biophysics, p41
a = 0.611;  % [kPa]
%b = 17.27;
%c = 237.3;
b = 17.502;
c = 240.97; % [C]
esat = a*exp((b*Ta)./(c+Ta));
ea = esat - VPD; % [kPa]
RH = ea./esat; % [kPa]

function forc_but_load_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile(...
    {'*.xls;*.xlsx'     , 'Microsoft Excel Files (*.xls,*.xlsx)';   ...
    '*.csv'             , 'CSV - comma delimited (*.csv)'       ;   ... 
    '*.txt'             , 'Text (Tab Delimited (*.txt)'         ;   ...
    '*.*'               , 'All Files (*.*)'                     },  ...
    'Import data from files');
    %'*.xls;*.xlsx'     , 'Microsoft Excel Files (*.xls,*.xlsx)';   ...
    
if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
else
	str_file = regexp(filename, '\.', 'split');
    if (strcmp(str_file(end),'xls') == 1 || strcmp(str_file(end),'xlsx') == 1)
        condition_excel = xlsread([pathname,filename]);      
    elseif strcmp(str_file(end),'csv') == 1
        condition_excel = csvread([pathname,filename], 1, 0);      
    elseif strcmp(str_file(end),'txt') == 1
        condition_excel = dlmread([pathname,filename], '\t', 1, 0);      
    else
        msgbox('The file extension is unrecognized, please choose another file','MLCan Error','Error');
        return
    end
    
    load './Temps/temp_variable.mat' 'Soil_heat' 'root_init' 'root_init_litter' 'litter_depth' 'dat_root1'
    if (length(condition_excel(1,:)) ~= 2 && length(condition_excel(1,:)) ~= 3)
        msgbox ('Incorrect format or input file for initial condition')
       	return        
    elseif (size(dat_root1,1) ~= size(condition_excel,1)-1)
        msgbox ('Number of soil layer dose not match up with what is in MODEL SETUP - check the input excel data','MLCan Error','Error')
        return
    elseif (sum(condition_excel(:,1) - [litter_depth; dat_root1(:,1)]) > 0.001 || sum(condition_excel(:,1) - [litter_depth; dat_root1(:,1)]) < - 0.001)
        msgbox ('Soil/litter depth distribution dose not match up with what is in MODEL SETUP - check the soil/litter depth','MLCan Error','Error')
        return
    else
        if (length(condition_excel(1,:)) == 2)     
            if Soil_heat == 1
                msgbox('Not enough information, check the input file or column 3 in excel file','MLCan Error: File format');
                return
            else
                condition_excel = [condition_excel(:,1:2),[root_init_litter(:,3); root_init(:,3)]];
            end
        else
            if Soil_heat == 0
                msgbox('Soil heat model is turned off, the last column is ignored');
                condition_excel = [condition_excel(:,1:2),zeros(size(condition_excel,1),1)];
            end
        end
    end
    set(handles.forc_table_root,'Data',condition_excel(2:end,:));   
    set(handles.forc_table_root_litter,'Data',condition_excel(1,:));
    forc_table_root_CellEditCallback(hObject, eventdata, handles);
end

% --- Executes on button press in forc_but_load_nutrient.
function forc_but_load_nutrient_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile(...
    {'*.xls;*.xlsx'     , 'Microsoft Excel Files (*.xls,*.xlsx)';   ...
    '*.csv'             , 'CSV - comma delimited (*.csv)'       ;   ... 
    '*.txt'             , 'Text (Tab Delimited (*.txt)'         ;   ...
    '*.*'               , 'All Files (*.*)'                     },  ...
    'Import data from files');
    %'*.xls;*.xlsx'     , 'Microsoft Excel Files (*.xls,*.xlsx)';   ...
    
if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
else
	str_file = regexp(filename, '\.', 'split');
    if (strcmp(str_file(end),'xls') == 1 || strcmp(str_file(end),'xlsx') == 1)
        condition_excel = xlsread([pathname,filename]);      
    elseif strcmp(str_file(end),'csv') == 1
        condition_excel = csvread([pathname,filename], 1, 0);      
    elseif strcmp(str_file(end),'txt') == 1
        condition_excel = dlmread([pathname,filename], '\t', 1, 0);      
    else
        msgbox('The file extension is unrecognized, please choose another file','MLCan Error','Error');
        return
    end
    
    load './Temps/temp_variable.mat' 'Soil_C_pool1' 'root_init'
    if (length(condition_excel(1,:)) ~= 6 && length(condition_excel(1,:)) ~= 7)
        msgbox ('Incorrect format or input file for initial condition')
       	return
    else
        if (length(condition_excel(1,:)) == 6)     
            if Soil_C_pool1 == 0
                msgbox('Not enough information, check the option or excel file','MLCan Error: File format');
                return
            else
                %condition_excel = [condition_excel(:,1:2),root_init(:,3)];
            end
        else
            if Soil_C_pool1 == 1
                msgbox('Too much information, check the option or excel file','MLCan Error: File format');
                return
            else
                %condition_excel = [condition_excel(:,1:2),root_init(:,3)];
            end
        end
        load './Temps/temp_variable.mat' 'litter_depth' 'dat_root1'
        if (sum(condition_excel(:,1) - [litter_depth; dat_root1(:,1)]) > 0.001 || sum(condition_excel(:,1) - [litter_depth; dat_root1(:,1)]) < - 0.001)
            msgbox ('Soil/litter depth distribution dose not match up with what is in MODEL SETUP - check the soil/litter depth','MLCan Error','Error')
            return
        end
    end
    set(handles.forc_table_biogeochemical,'Data',condition_excel(2:end,:));   
    set(handles.forc_table_biogeochemical_litter,'Data',condition_excel(1,:));
    forc_table_biogeochemical_CellEditCallback(hObject, eventdata, handles);
end

function forc_but_newforcings_CreateFcn(hObject, eventdata, handles)


function forc_tabbut_forcing_N_CreateFcn(hObject, eventdata, handles)


%function uipanel9_SelectionChangeFcn(hObject, eventdata, handles)


function uipanel9_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in forc_but_humusCdr.
% function forc_but_humusCdr_Callback(hObject, eventdata, handles)
% % hObject    handle to forc_but_humusCdr (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
% soil_depth = nutrient_ic(:,1);
% HumusCdr = nutrient_ic(:,9);
%   
% plot(handles.forc_axes_ic4,HumusCdr,-soil_depth,'Color', 'r',...
%     'LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b',...
%     'MarkerFaceColor','m','MarkerSize',6);
% axis(handles.forc_axes_ic4,[min(HumusCdr)/1.2 max(HumusCdr)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
% ylabel(handles.forc_axes_ic4,'Soil depth [m]');
% %xlabel(handles.forc_axes_ic4,'Microbial death rate [/t]');
% legend(handles.forc_axes_ic4,'Soil humus decomposition rate [gC/m^{3}/t]');
% grid on;

% % --- Executes on button press in forc_but_microbialCdr.
% function forc_but_microbialCdr_Callback(hObject, eventdata, handles)
% % hObject    handle to forc_but_microbialCdr (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% load './Temps/temp_variable.mat' 'Soil_C_pool1'
% nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
% soil_depth = nutrient_ic(:,1);
% if Soil_C_pool1 == 1 % 1 soil C pool
%     MicrobialDR = nutrient_ic(:,8);
% else % 2 soil C pools
%     MicrobialDR = nutrient_ic(:,10);
% end
%   
% plot(handles.forc_axes_ic4,MicrobialDR,-soil_depth,'Color', 'r',...
%     'LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b',...
%     'MarkerFaceColor','m','MarkerSize',6);
% axis(handles.forc_axes_ic4,[min(MicrobialDR)/1.2 max(MicrobialDR)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
% ylabel(handles.forc_axes_ic4,'Soil depth [m]');
% %xlabel(handles.forc_axes_ic4,'Microbial death rate [/t]');
% legend(handles.forc_axes_ic4,'Microbial death rate [/t]');
% grid on;

% --- Executes on button press in forc_but_litterCdr.
% function forc_but_litterCdr_Callback(hObject, eventdata, handles)
% % hObject    handle to forc_but_litterCdr (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
% soil_depth = nutrient_ic(:,1);
% LitterCdr = nutrient_ic(:,8);
%   
% plot(handles.forc_axes_ic4,LitterCdr,-soil_depth,'Color', 'r',...
%     'LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b',...
%     'MarkerFaceColor','m','MarkerSize',6);
% axis(handles.forc_axes_ic4,[min(LitterCdr)/1.2 max(LitterCdr)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
% ylabel(handles.forc_axes_ic4,'Soil depth [m]');
% %xlabel(handles.forc_axes_ic4,'Microbial death rate [/t]');
% legend(handles.forc_axes_ic4,'Soil litter decomposition rate [gC/m^{3}/t]');
% grid on;

% --- Executes on button press in forc_but_litterCN.
function forc_but_litterCN_Callback(hObject, eventdata, handles)
% hObject    handle to forc_but_litterCN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load './Temps/temp_variable.mat' 'Soil_C_pool1'
nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
soil_depth = nutrient_ic(:,1);
if Soil_C_pool1 == 1 % 1 soil C pool
    LitterCN = nutrient_ic(:,6);
else % 2 soil C pools
    LitterCN = nutrient_ic(:,7);
end

plot(handles.forc_axes_ic4,LitterCN,-soil_depth,'Color', 'r',...
    'LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b',...
    'MarkerFaceColor','m','MarkerSize',6);
axis(handles.forc_axes_ic4,[min(LitterCN)/1.2 max(LitterCN)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
ylabel(handles.forc_axes_ic4,'Soil depth [m]');
%xlabel(handles.forc_axes_ic4,'Microbial death rate [/t]');
legend(handles.forc_axes_ic4,'Soil litter C:N ratio [/t]');
grid on;

% --- Executes on button press in forc_but_ammonium.
function forc_but_ammonium_Callback(hObject, eventdata, handles)
% hObject    handle to forc_but_ammonium (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load './Temps/temp_variable.mat' 'Soil_C_pool1'
nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
soil_depth = nutrient_ic(:,1);
if Soil_C_pool1 == 1 % 1 soil C pool
    Ammonium = nutrient_ic(:,4);
else % 2 soil C pools
    Ammonium = nutrient_ic(:,5);
end
   
plot(handles.forc_axes_ic4,Ammonium,-soil_depth,'Color', 'r',...
    'LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b',...
    'MarkerFaceColor','m','MarkerSize',6);
axis(handles.forc_axes_ic4,[min(Ammonium)/1.2 max(Ammonium)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
ylabel(handles.forc_axes_ic4,'Soil depth [m]');
%xlabel(handles.forc_axes_ic3,'Microbial Carbon [gC/m^{-3}]');
legend(handles.forc_axes_ic4,'Soil ammonium-N [gN/m^{-3}]');
grid on;

% --- Executes on button press in forc_but_nitrate.
function forc_but_nitrate_Callback(hObject, eventdata, handles)
% hObject    handle to forc_but_nitrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load './Temps/temp_variable.mat' 'Soil_C_pool1'
nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
soil_depth = nutrient_ic(:,1);
if Soil_C_pool1 == 1 % 1 soil C pool
    Nitrate = nutrient_ic(:,5);
else % 2 soil C pools
    Nitrate = nutrient_ic(:,6);
end
   
plot(handles.forc_axes_ic4,Nitrate,-soil_depth,'Color', 'r',...
    'LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b',...
    'MarkerFaceColor','m','MarkerSize',6);
axis(handles.forc_axes_ic4,[min(Nitrate)/1.2 max(Nitrate)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
ylabel(handles.forc_axes_ic4,'Soil depth [m]');
%xlabel(handles.forc_axes_ic3,'Microbial Carbon [gC/m^{-3}]');
legend(handles.forc_axes_ic4,'Soil nitrate-N [gN/m^{-3}]');
grid on;


% --- Executes on button press in forc_but_litterC.
function forc_but_litterC_Callback(hObject, eventdata, handles)
% hObject    handle to forc_but_litterC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
soil_depth = nutrient_ic(:,1);
LitterC = nutrient_ic(:,2);
   
plot(handles.forc_axes_ic3,LitterC,-soil_depth,'Color', 'b',...
    'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',6);
axis(handles.forc_axes_ic3,[min(LitterC)/1.2 max(LitterC)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
ylabel(handles.forc_axes_ic3,'Soil depth [m]');
%xlabel(handles.forc_axes_ic3,'Soil organic carbon [gC/m^{-3}]');
legend(handles.forc_axes_ic3,'Soil litter carbon [gC/m^{3}]');
grid on;

% --- Executes on button press in forc_but_humusC.
function forc_but_humusC_Callback(hObject, eventdata, handles)
% hObject    handle to forc_but_humusC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
soil_depth = nutrient_ic(:,1);
HumusC = nutrient_ic(:,3);
   
plot(handles.forc_axes_ic3,HumusC,-soil_depth,'Color', 'b',...
    'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',6);
axis(handles.forc_axes_ic3,[min(HumusC)/1.2 max(HumusC)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
ylabel(handles.forc_axes_ic3,'Soil depth [m]');
%xlabel(handles.forc_axes_ic3,'Soil organic carbon [gC/m^{-3}]');
legend(handles.forc_axes_ic3,'Soil humus carbon [gC/m^{3}]');
grid on;

% --- Executes on button press in forc_but_microbialC.
function forc_but_microbialC_Callback(hObject, eventdata, handles)
% hObject    handle to forc_but_microbialC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load './Temps/temp_variable.mat' 'Soil_C_pool1'
nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
soil_depth = nutrient_ic(:,1);
if Soil_C_pool1 == 1 % 1 soil C pool
    MicrobialC = nutrient_ic(:,3);
else % 2 soil C pools
    MicrobialC = nutrient_ic(:,4);
end
   
plot(handles.forc_axes_ic3,MicrobialC,-soil_depth,'Color', 'b',...
    'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',6);
axis(handles.forc_axes_ic3,[min(MicrobialC)/1.2 max(MicrobialC)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
ylabel(handles.forc_axes_ic3,'Soil depth [m]');
%xlabel(handles.forc_axes_ic3,'Microbial Carbon [gC/m^{-3}]');
legend(handles.forc_axes_ic3,'Soil microbial carbon [gC/m^{-3}]');
grid on;

% --- Executes on button press in forc_but_organicC.
function forc_but_organicC_Callback(hObject, eventdata, handles)
% hObject    handle to forc_but_organicC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
soil_depth = nutrient_ic(:,1);
OrganicC = nutrient_ic(:,2);
   
plot(handles.forc_axes_ic3,OrganicC,-soil_depth,'Color', 'b',...
    'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',6);
axis(handles.forc_axes_ic3,[min(OrganicC)/1.2 max(OrganicC)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
ylabel(handles.forc_axes_ic3,'Soil depth [m]');
%xlabel(handles.forc_axes_ic3,'Soil organic carbon [gC/m^{-3}]');
legend(handles.forc_axes_ic3,'Soil organic carbon [gC/m^{3}]');
grid on;


% --- Executes on button press in forc_but_organicCdr.
% function forc_but_organicCdr_Callback(hObject, eventdata, handles)
% % hObject    handle to forc_but_organicCdr (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
% soil_depth = nutrient_ic(:,1);
% OrganicCdr = nutrient_ic(:,7);
%    
% plot(handles.forc_axes_ic4,OrganicCdr,-soil_depth,'Color', 'r',...
%     'LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b',...
%     'MarkerFaceColor','m','MarkerSize',6);
% axis(handles.forc_axes_ic4,[min(OrganicCdr)/1.2 max(OrganicCdr)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
% ylabel(handles.forc_axes_ic4,'Soil depth [m]');
% %xlabel(handles.forc_axes_ic4,'Soil organic carbon [gC/m^{-3}]');
% legend(handles.forc_axes_ic4,'Soil organic carbon decomposition rate [m^{3}/gC/t]');
% grid on;


% --- Executes on button press in forc_but_organicCN.
function forc_but_organicCN_Callback(hObject, eventdata, handles)
% hObject    handle to forc_but_organicCN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nutrient_ic = get(handles.forc_table_biogeochemical,'Data');
soil_depth = nutrient_ic(:,1);
OrganicCN = nutrient_ic(:,6);
   
plot(handles.forc_axes_ic4,OrganicCN,-soil_depth,'Color', 'r',...
    'LineStyle','--','LineWidth',2,'Marker','o','MarkerEdgeColor','b',...
    'MarkerFaceColor','m','MarkerSize',6);
axis(handles.forc_axes_ic4,[min(OrganicCN)/1.2 max(OrganicCN)*1.1+10^(-10) min(-soil_depth)*1.2 0]);
ylabel(handles.forc_axes_ic4,'Soil depth [m]');
%xlabel(handles.forc_axes_ic4,'Soil organic carbon [gC/m^{-3}]');
legend(handles.forc_axes_ic4,'Soil organic matter C:N ratio [-]');
grid on;


% --- Executes during object creation, after setting all properties.
function forc_gbut_options_create_CreateFcn(hObject, eventdata, handles)
% hObject    handle to forc_gbut_options_create (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function forc_gbut_options_create_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to forc_gbut_options_create (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in forc_but_view.
function forc_but_view_Callback(hObject, eventdata, handles)
% hObject    handle to forc_but_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Forcing_view;


% --- Executes on button press in pushbutton39.
function pushbutton39_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton40.
function pushbutton40_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel9.
function uipanel9_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel9 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'forc_tabbut_forcing'
        set(handles.forc_pan_condition,'Visible','off');
        set(handles.forc_pan_forcing,'Visible','on');
        set(handles.forc_pan_nutrient,'Visible','off');
        set(handles.forc_pan_condition,'Visible','off');
        set(handles.forc_but_load,'Enable','off');
        set(handles.forc_but_load_nutrient,'Enable','off');
        set(handles.forc_gbut_options_create,'Visible','on');
        load './Temps/temp_variable.mat'...
            'working_forcings'
        if isempty(working_forcings)
            set(handles.forc_but_view,'Enable','off');
        else
            set(handles.forc_but_view,'Enable','on');
        end
    otherwise
end
%updates the handles structure
guidata(hObject, handles);
