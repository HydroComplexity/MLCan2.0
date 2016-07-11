function varargout = Forcing_view(varargin)
% FORCING_VIEW MATLAB code for Forcing_view.fig
%      FORCING_VIEW, by itself, creates a new FORCING_VIEW or raises the existing
%      singleton*.
%
%      H = FORCING_VIEW returns the handle to a new FORCING_VIEW or the handle to
%      the existing singleton*.
%
%      FORCING_VIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FORCING_VIEW.M with the given input arguments.
%
%      FORCING_VIEW('Property','Value',...) creates a new FORCING_VIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Forcing_view_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Forcing_view_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Forcing_view

% Last Modified by GUIDE v2.5 19-Jan-2014 13:22:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Forcing_view_OpeningFcn, ...
                   'gui_OutputFcn',  @Forcing_view_OutputFcn, ...
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



% --- Executes just before Forcing_view is made visible.
function Forcing_view_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Forcing_view (see VARARGIN)

% Choose default command line output for Forcing_view
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Forcing_view wait for user response (see UIRESUME)
% uiwait(handles.figure1);

load './Temps/temp_variable.mat' 'working_forcings' 'num_species' 'Soil_nutrient' 'N_Uptake_RB'
load(working_forcings);

if Soil_nutrient == 1
    if N_Uptake_RB == 1
        dat_forcing=[year_crop, month_crop, day_crop, doy_crop, hour_crop,  ...
            Ta_crop, Pa_crop, Rg_crop, PPT_crop, VPD_crop, U_crop, ustar_crop, ...
            LAI_species1_crop, LAI_species2_crop, LAI_species3_crop, LAI_species4_crop,  ...
            BMla_species1_crop, BMla_species2_crop, BMla_species3_crop, BMla_species4_crop, ...
            CNratio_species1_crop, CNratio_species2_crop, CNratio_species3_crop, CNratio_species4_crop, ... 
            LMI_species1_crop, LMI_species2_crop, LMI_species3_crop, LMI_species4_crop];        
        
        set(handles.forcing_tab2,'Data',dat_forcing);
        set(handles.forcing_tab,'Visible','off');
        set(handles.forcing_tab2,'Visible','on');
        set(handles.forcing_tab3,'Visible','off');
        
        forcing_str = {...
            'Air temperature [''C]'                 , 'Barometric pressure [kPa]'             , 'Global solar radition [W/m2]'          , 'Precipitation [mm]'                    , ...
            'Vapor pressure deficit [kPa]'          , 'Wind speed [m/2]'                      , 'Friction Velocity [m/s]'               , ....
            'Leaf Area Index for species 1 [m2/m2]' , 'Leaf Area Index for species 2 [m2/m2]' , 'Leaf Area Index for species 3 [m2/m2]' , 'Leaf Area Index for species 4 [m2/m2]' ,...
            'Aboveground biomass drop as litter for species 1 [gC/m2]','Aboveground biomass drop as litter for species 2 [gC/m2]','Aboveground biomass drop as litter for species 3 [gC/m2]','Aboveground biomass drop as litter for species 4 [gC/m2]',...
            'Carbon to Nitrogen ratio for species 1 [-]','Carbon to Nitrogen ratio for species 2 [-]','Carbon to Nitrogen ratio for species 3 [-]','Carbon to Nitrogen ratio for species 4 [-]',...
            'Root biomass for species 1 [g/m2]','Root biomass for species 2 [g/m2]','Root biomass for species 3 [g/m2]','Root biomass for species 4 [g/m2]'};
        
    elseif N_Uptake_RB == 0
        dat_forcing=[year_crop, month_crop, day_crop, doy_crop, hour_crop,  ...
            Ta_crop, Pa_crop, Rg_crop, PPT_crop, VPD_crop, U_crop, ustar_crop, ...
            LAI_species1_crop, LAI_species2_crop, LAI_species3_crop, LAI_species4_crop, ...
            BMlaDEM_species1_crop, BMlaDEM_species2_crop, BMlaDEM_species3_crop, BMlaDEM_species4_crop,...
            BMrbDEM_species1_crop, BMrbDEM_species2_crop, BMrbDEM_species3_crop, BMrbDEM_species4_crop, ...
            CNratioA_species1_crop, CNratioA_species2_crop, CNratioA_species3_crop, CNratioA_species4_crop, ...
            CNratioB_species1_crop, CNratioB_species2_crop, CNratioB_species3_crop, CNratioB_species4_crop, ...
            NDEM_species1_crop, NDEM_species2_crop, NDEM_species3_crop, NDEM_species4_crop];
        
        set(handles.forcing_tab3,'Data',dat_forcing);
        set(handles.forcing_tab,'Visible','off');
        set(handles.forcing_tab2,'Visible','off');
        set(handles.forcing_tab3,'Visible','on');
        
        forcing_str = {...
            'Air temperature [''C]'                 , 'Barometric pressure [kPa]'             , 'Global solar radition [W/m2]'          , 'Precipitation [mm]'                    , ...
            'Vapor pressure deficit [kPa]'          , 'Wind speed [m/2]'                      , 'Friction Velocity [m/s]'               , ....
            'Leaf Area Index for species 1 [m2/m2]' , 'Leaf Area Index for species 2 [m2/m2]' , 'Leaf Area Index for species 3 [m2/m2]' , 'Leaf Area Index for species 4 [m2/m2]' ,...
            'Aboveground biomass drop as litter for species 1 [gC/m2]','Aboveground biomass drop as litter for species 2 [gC/m2]','Aboveground biomass drop as litter for species 3 [gC/m2]','Aboveground biomass drop as litter for species 4 [gC/m2]'...
            'Belowground biomass drop as litter for species 1 [gC/m2]','Belowground biomass drop as litter for species 2 [gC/m2]','Belowground biomass drop as litter for species 3 [gC/m2]','Belowground biomass drop as litter for species 4 [gC/m2]'...
            'Aboveground carbon to nitrogen ratio for species 1 [-]','Aboveground carbon to nitrogen ratio for species 2 [-]','Aboveground carbon to nitrogen ratio for species 3 [-]','Aboveground carbon to nitrogen ratio for species 4 [-]'...
            'Belowground carbon to nitrogen ratio for species 1 [-]','Belowground carbon to nitrogen ratio for species 2 [-]','Belowground carbon to nitrogen ratio for species 3 [-]','Belowground carbon to nitrogen ratio for species 4 [-]'...
            'Nitrogen demand for species 1 [gN/m2]','Nitrogen demand for species 2 [gN/m2]','Nitrogen demand for species 3 [gN/m2]','Nitrogen demand for species 4 [gN/m2]'...
            };
    end
elseif Soil_nutrient == 0
    dat_forcing=[year_crop, month_crop, day_crop, doy_crop, hour_crop,  ...
        Ta_crop, Pa_crop, Rg_crop, PPT_crop, VPD_crop, U_crop, ustar_crop, ...
        LAI_species1_crop, LAI_species2_crop, LAI_species3_crop, LAI_species4_crop];
    
    set(handles.forcing_tab,'Data',dat_forcing);
    set(handles.forcing_tab,'Visible','on');
    set(handles.forcing_tab2,'Visible','off');
    set(handles.forcing_tab3,'Visible','off');
    
    forcing_str = {...
        'Air temperature [''C]'                 , 'Barometric pressure [kPa]'             , 'Global solar radition [W/m2]'          , 'Precipitation [mm]'                    , ...
        'Vapor pressure deficit [kPa]'          , 'Wind speed [m/2]'                      , 'Friction Velocity [m/s]'               , ....
        'Leaf Area Index for species 1 [m2/m2]' , 'Leaf Area Index for species 2 [m2/m2]' , 'Leaf Area Index for species 3 [m2/m2]' , 'Leaf Area Index for species 4 [m2/m2]' ,...
        };
end




           
set(handles.forcing_lst_str,'String',forcing_str);
but_plot_Callback(hObject, eventdata, handles)

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = Forcing_view_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in but_close.
function but_close_Callback(hObject, eventdata, handles)
% hObject    handle to but_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close

% --- Executes on button press in but_ok.
function but_ok_Callback(hObject, eventdata, handles)
% hObject    handle to but_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load './Temps/temp_variable.mat' 'Soil_nutrient' 'N_Uptake_RB'
%
if Soil_nutrient == 1
    if N_Uptake_RB == 1
        dat_forcing       = get(handles.forcing_tab2,'Data');
    elseif N_Uptake_RB == 0
        dat_forcing       = get(handles.forcing_tab3,'Data');
    end
elseif Soil_nutrient == 0
    dat_forcing       = get(handles.forcing_tab,'Data');
end



if Soil_nutrient == 1
    if N_Uptake_RB == 1
        year_crop         = dat_forcing(:,1);
        month_crop        = dat_forcing(:,2);
        day_crop          = dat_forcing(:,3);
        doy_crop          = dat_forcing(:,4);
        hour_crop         = dat_forcing(:,5);
        
        Ta_crop           = dat_forcing(:,6);
        Pa_crop           = dat_forcing(:,7);
        Rg_crop           = dat_forcing(:,8);
        PPT_crop          = dat_forcing(:,9);
        VPD_crop          = dat_forcing(:,10);
        
        U_crop            = dat_forcing(:,11);
        ustar_crop        = dat_forcing(:,12);
        
        LAI_species1_crop = dat_forcing(:,13);
        LAI_species2_crop = dat_forcing(:,14);
        LAI_species3_crop = dat_forcing(:,15);
        LAI_species4_crop = dat_forcing(:,16);
        
        BMla_species1_crop    = dat_forcing(:,17);
        BMla_species2_crop    = dat_forcing(:,18);
        BMla_species3_crop    = dat_forcing(:,19);
        BMla_species4_crop    = dat_forcing(:,20);
        
        CNratio_species1_crop = dat_forcing(:,21);
        CNratio_species2_crop = dat_forcing(:,22);
        CNratio_species3_crop = dat_forcing(:,23);
        CNratio_species4_crop = dat_forcing(:,24);
        
        LMI_species1_crop     = dat_forcing(:,25);
        LMI_species2_crop     = dat_forcing(:,26);
        LMI_species3_crop     = dat_forcing(:,27);
        LMI_species4_crop     = dat_forcing(:,28);
        
        load './Temps/temp_variable.mat' 'working_forcings'
        save(working_forcings, ...
            'year_crop', 'doy_crop', 'hour_crop', 'day_crop', 'month_crop', ...
            'LAI_species1_crop', 'LAI_species2_crop', 'LAI_species3_crop', 'LAI_species4_crop', 'Ta_crop', ...
            'U_crop', 'Rg_crop', 'PPT_crop', 'VPD_crop', 'Pa_crop', 'ustar_crop', ...
            'BMla_species1_crop', 'BMla_species2_crop', 'BMla_species3_crop', 'BMla_species4_crop', ...
            'CNratio_species1_crop', 'CNratio_species2_crop', 'CNratio_species3_crop', 'CNratio_species4_crop', ...
            'LMI_species1_crop' , 'LMI_species2_crop' , 'LMI_species3_crop' , 'LMI_species4_crop', '-append');
    elseif N_Uptake_RB == 0
        year_crop         = dat_forcing(:,1);
        month_crop        = dat_forcing(:,2);
        day_crop          = dat_forcing(:,3);
        doy_crop          = dat_forcing(:,4);
        hour_crop         = dat_forcing(:,5);
        
        Ta_crop           = dat_forcing(:,6);
        Pa_crop           = dat_forcing(:,7);
        Rg_crop           = dat_forcing(:,8);
        PPT_crop          = dat_forcing(:,9);
        VPD_crop          = dat_forcing(:,10);
        
        U_crop            = dat_forcing(:,11);
        ustar_crop        = dat_forcing(:,12);
        
        LAI_species1_crop = dat_forcing(:,13);
        LAI_species2_crop = dat_forcing(:,14);
        LAI_species3_crop = dat_forcing(:,15);
        LAI_species4_crop = dat_forcing(:,16);
        
        BMlaDEM_species1_crop  = dat_forcing(:,17);
        BMlaDEM_species2_crop  = dat_forcing(:,18);
        BMlaDEM_species3_crop  = dat_forcing(:,19);
        BMlaDEM_species4_crop  = dat_forcing(:,20);
        
        BMrbDEM_species1_crop  = dat_forcing(:,21);
        BMrbDEM_species2_crop  = dat_forcing(:,22);
        BMrbDEM_species3_crop  = dat_forcing(:,23);
        BMrbDEM_species4_crop  = dat_forcing(:,24);
        
        CNratioA_species1_crop  = dat_forcing(:,25);
        CNratioA_species2_crop  = dat_forcing(:,26);
        CNratioA_species3_crop  = dat_forcing(:,27);
        CNratioA_species4_crop  = dat_forcing(:,28);       

        CNratioB_species1_crop  = dat_forcing(:,29);
        CNratioB_species2_crop  = dat_forcing(:,30);
        CNratioB_species3_crop  = dat_forcing(:,31);
        CNratioB_species4_crop  = dat_forcing(:,32);       
        
        NDEM_species1_crop      = dat_forcing(:,33);
        NDEM_species2_crop      = dat_forcing(:,34);
        NDEM_species3_crop      = dat_forcing(:,35);
        NDEM_species4_crop      = dat_forcing(:,36);
        
        load './Temps/temp_variable.mat' 'working_forcings'
        save(working_forcings, ...
            'year_crop', 'doy_crop', 'hour_crop', 'day_crop', 'month_crop', ...
            'LAI_species1_crop', 'LAI_species2_crop', 'LAI_species3_crop', 'LAI_species4_crop', 'Ta_crop', ...
            'U_crop', 'Rg_crop', 'PPT_crop', 'VPD_crop', 'Pa_crop', 'ustar_crop', ...
            'BMlaDEM_species1_crop', 'BMlaDEM_species2_crop','BMlaDEM_species3_crop','BMlaDEM_species4_crop',...
            'BMrbDEM_species1_crop','BMrbDEM_species2_crop','BMrbDEM_species3_crop','BMrbDEM_species4_crop', ...
            'CNratioA_species1_crop','CNratioA_species2_crop','CNratioA_species3_crop','CNratioA_species4_crop', ...
            'CNratioB_species1_crop','CNratioB_species2_crop','CNratioB_species3_crop','CNratioB_species4_crop', ...
            'NDEM_species1_crop','NDEM_species2_crop','NDEM_species3_crop','NDEM_species4_crop' ,'-append');
                        
    end
elseif Soil_nutrient == 0
    year_crop         = dat_forcing(:,1);
    month_crop        = dat_forcing(:,2);
    day_crop          = dat_forcing(:,3);
    doy_crop          = dat_forcing(:,4);
    hour_crop         = dat_forcing(:,5);
    
    Ta_crop           = dat_forcing(:,6);
    Pa_crop           = dat_forcing(:,7);
    Rg_crop           = dat_forcing(:,8);
    PPT_crop          = dat_forcing(:,9);
    VPD_crop          = dat_forcing(:,10);
    
    U_crop            = dat_forcing(:,11);
    ustar_crop        = dat_forcing(:,12);
    
    LAI_species1_crop = dat_forcing(:,13);
    LAI_species2_crop = dat_forcing(:,14);
    LAI_species3_crop = dat_forcing(:,15);
    LAI_species4_crop = dat_forcing(:,16);
    
    load './Temps/temp_variable.mat' 'working_forcings'
    save(working_forcings, ...
        'year_crop', 'doy_crop', 'hour_crop', 'day_crop', 'month_crop', ...
        'LAI_species1_crop', 'LAI_species2_crop', 'LAI_species3_crop', 'LAI_species4_crop', 'Ta_crop', ...
        'U_crop', 'Rg_crop', 'PPT_crop', 'VPD_crop', 'Pa_crop', 'ustar_crop','-append');
end

close;


% --- Executes on selection change in forcing_lst_str.
function forcing_lst_str_Callback(hObject, eventdata, handles)
% hObject    handle to forcing_lst_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns forcing_lst_str contents as cell array
%        contents{get(hObject,'Value')} returns selected item from forcing_lst_str


% --- Executes during object creation, after setting all properties.
function forcing_lst_str_CreateFcn(hObject, eventdata, handles)
% hObject    handle to forcing_lst_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in forcing_chk_grid.
function forcing_chk_grid_Callback(hObject, eventdata, handles)
% hObject    handle to forcing_chk_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of forcing_chk_grid


% --- Executes on button press in forcing_chk_hold.
function forcing_chk_hold_Callback(hObject, eventdata, handles)
% hObject    handle to forcing_chk_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of forcing_chk_hold


% --- Executes on button press in but_plot.
function but_plot_Callback(hObject, eventdata, handles)
% hObject    handle to but_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
load './Temps/temp_variable.mat' 'working_forcings' 'num_species' 'DOY_start' 'DOY_end' 'Soil_nutrient' 'N_Uptake_RB'
load(working_forcings);

%
if Soil_nutrient == 1
    if N_Uptake_RB == 1
        dat_forcing       = get(handles.forcing_tab2,'Data');
    elseif N_Uptake_RB == 0
        dat_forcing       = get(handles.forcing_tab3,'Data');
    end
elseif Soil_nutrient == 0
    dat_forcing       = get(handles.forcing_tab,'Data');
end

doy_crop          = dat_forcing(:,4);
doys=[DOY_start:DOY_end];
inds = find(ismember(floor(doy_crop), doys));

%

index_forcing_p = get(handles.forcing_lst_str,'Value');
Name_forcing_p  = get(handles.forcing_lst_str,'String');
Plot_variable   = dat_forcing(:,index_forcing_p+5);
Plot_name       = Name_forcing_p{index_forcing_p};
x_value         = [1:size(Plot_variable,1)]';

Ind_Nan         = isnan(Plot_variable);

forcing_hold        = get(handles.forcing_chk_hold,'Value');
if forcing_hold == 1
    hold on
else
    hold off
end

if isnan(Plot_variable(Ind_Nan)) == 1    
    plot(x_value(Ind_Nan),zeros(size(x_value(Ind_Nan),1),1),'MarkerFaceColor',[1 0 0],...
        'MarkerEdgeColor',[0 0 0],'Marker','square','LineStyle','none');
    hold on
end


% plot(handles.axes_forcing_plot,x_value,Plot_variable,'Color', 'b',...
%         'LineStyle','--','LineWidth',2,'Marker','s','MarkerEdgeColor','k',...
%         'MarkerFaceColor','g','MarkerSize',6);
plot(handles.axes_forcing_plot,x_value,Plot_variable,'Color', 'b',...
    'LineStyle','-','LineWidth',1);
if isnan(Plot_variable(Ind_Nan)) == 1
   % axis(handles.axes_forcing_plot,[0 size(Plot_variable,1) 0 max(Plot_variable)*1.1+0.0000000000001]);
    legend(handles.axes_forcing_plot,'Not A Number, needs fixing');
    if index_forcing_p == 7
        legend(handles.axes_forcing_plot,'Not A Number will be fixed while simulation');
    end
else
    axis(handles.axes_forcing_plot,[0 size(Plot_variable,1) min(Plot_variable)/1.1 max(Plot_variable)*1.1+0.0000000000001]);    
    legend(handles.axes_forcing_plot,Plot_name);
end
xlabel('Time step');
ylabel(handles.axes_forcing_plot,Plot_name);

% Indicator of simulation period
hold on
plot(handles.axes_forcing_plot,inds(1)*ones(1,10000),(-10000/2:10000/2-1),'Color', 'r',...
    'LineStyle','-','LineWidth',1);
plot(handles.axes_forcing_plot,inds(end)*ones(1,10000),(-10000/2:10000/2-1),'Color', 'r',...
    'LineStyle','--','LineWidth',1);
%


forcing_grid = get(handles.forcing_chk_grid,'Value');
if forcing_grid == 1
    grid on
else
    grid off
end


% --- Executes on button press in but_TimeStepCheck.
function but_TimeStepCheck_Callback(hObject, eventdata, handles)
% hObject    handle to but_TimeStepCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MSG1=msgbox('This check only works for less than a year simulation','MLCan Infomation','warn');
uiwait(MSG1);

%load './Temps/temp_variable.mat' 'working_forcings'
%load(working_forcings,'starttime','endtime')
load './Temps/temp_variable.mat' 'Soil_nutrient' 'N_Uptake_RB'
%
if Soil_nutrient == 1
    if N_Uptake_RB == 1
        dat_forcing       = get(handles.forcing_tab2,'Data');
    elseif N_Uptake_RB == 0
        dat_forcing       = get(handles.forcing_tab3,'Data');
    end
elseif Soil_nutrient == 0
    dat_forcing       = get(handles.forcing_tab,'Data');
end

year_crop         = dat_forcing(:,1);    
month_crop        = dat_forcing(:,2);
day_crop          = dat_forcing(:,3);
doy_crop          = dat_forcing(:,4);
hour_crop         = dat_forcing(:,5);


% isleap = @(x) (mod(x,4)==0 & mod(x,100)~=0) | mod(x,400) == 0;
% leapyrs = (isleap(year_crop(1)));
% if leapyrs == 1
%     day_start_here = day(day_crop(1));
% else
%     if day_crop(1) > 59
%         day_start_here = day(day_crop(1)+1);
%     else
%         day_start_here = day(day_crop(1));
%     end
% end
day_start_here = (day_crop(1));

% isleap = @(x) (mod(x,4)==0 & mod(x,100)~=0) | mod(x,400) == 0;
% leapyrs = (isleap(year_crop(end)));
% if leapyrs == 1
%     day_end_here = day(day_crop(end));
% else
%     if day_crop(end) > 59
%         day_end_here = day(day_crop(end)+1);
%     else
%         day_end_here = day(day_crop(end));
%     end
% end
day_end_here = (day_crop(end));

% Find leap year
[uni_year,ind_uni_year]=unique(year_crop);
% for i=1:size(uni_year,1)
%     if mod(uni_year(i),400) == 0
%         is_leap_year(i)=uni_year(i);
%     elseif mod(uni_year(i),100) == 0
%         not_leap_year(i)=uni_year(i);
%     elseif mod(uni_year(i),4) == 0
%         is_leap_year(i)=uni_year(i);
%     else
%         not_leap_year(i)=uni_year(i);
%     end
% end
% leap_year=unique(is_leap_year);
% ind_zero1=(leap_year==0);
% leap_year(ind_zero1)=[];
% N_leap_year=unique(not_leap_year);
% ind_zero2=(N_leap_year==0);
% N_leap_year(ind_zero2)=[];

% Set base step!
days_step_test=[doy_crop(2)-doy_crop(1) doy_crop(3)-doy_crop(2) doy_crop(4)-doy_crop(3)];
days_step=floor(mode(days_step_test));
hours_step_test=[hour_crop(2)-hour_crop(1) hour_crop(3)-hour_crop(2) hour_crop(4)-hour_crop(3)];
hours_step=floor(mode(hours_step_test));
mins_step_test=[hour_crop(2)-hour_crop(1) hour_crop(3)-hour_crop(2) hour_crop(4)-hour_crop(3)];
mins_step=mod((mode(mins_step_test))*60,60);

base_time_step_hr=days_step*24+hours_step+mins_step/60;

% Set start and end time
hour_start_here=hour_crop(1)-base_time_step_hr;
if hour_start_here < 0
    hour_start_here=24-base_time_step_hr;
end
hour_end_here=hour_crop(end);
starttime=datenum(year_crop(1),month_crop(1),day_start_here,hour_start_here,0,0);
endtime=datenum(year_crop(end),month_crop(end),day_end_here,hour_end_here,0,0);

Start_Simul = starttime;
End_Simul   = endtime+0.98;%+(base_time_step_hr*1/24);

% Expected Time step
Ex_t_step= [Start_Simul:base_time_step_hr*1/24:End_Simul];

if isempty(Ex_t_step) == 1
    msgbox('1st or/and last line''s time step is/are wrong. Please check','MLCan Error: File format','warn');
    return
end


% Y M D H MN S
% Check Time step
YMDHMS       = datevec(Ex_t_step);
Year_Should  = YMDHMS(:,1);
Month_Should = YMDHMS(:,2);
Day_Should   = YMDHMS(:,3);

for i=1:size(YMDHMS,1)
    Hour_Should(i)  = mod((hour_crop(1)+base_time_step_hr*(i-1)),24);
end

%for i=1:size(ind_uni_year,1)
%    if i == 1
%        for j=1:ind_uni_year(i)
        for j=1:size(YMDHMS,1)
            DOY_Should(j)=(Day_Should(1)+Hour_Should(1)/24)+base_time_step_hr/24*(j-1)-1-base_time_step_hr/24;
            if DOY_Should(j) < 0 && DOY_Should(j) > -0.00001 
                DOY_Should(j) = 0;
            end
        end
%    else
%        for j=ind_uni_year(i-1)+1:ind_uni_year(i)
%            DOY_Should(j)=1+base_time_step_hr/24*(j-(ind_uni_year(i-1)+1))-1;
%        end
%    end
%end

if DOY_Should(1) < 0
    Message_forcings_view = ['DOY has a problem. General cases for this issue are that the first day of the forcing hour should start from 0.5, not 0.0 and/or the last day of the forcing hour should end with 0.0. Please check!'];
    msgbox(Message_forcings_view,'MLCan Error: File format','warn');
    return
end

if (size(year_crop,1)  ~= size(Year_Should,1))
    if (size(year_crop,1)-size(Year_Should,1) == 48)
        Message_forcings_view = ['Number of time step has a problem. It should be ', num2str(size(Year_Should,1)),...
            ' not ' , num2str(size(year_crop,1)), '. A general case for this issue is that the first day of simulation should start 0.0, not 0.5. Please check!'];
        msgbox(Message_forcings_view,'MLCan Error: File format','warn');
        return
    else
        Message_forcings_view = ['Number of time step has a problem. It should be ', num2str(size(Year_Should,1)),...
            ' not ' , num2str(size(year_crop,1)), '. Please check!'];
        msgbox(Message_forcings_view,'MLCan Error: File format','warn');
        return
    end
end

ind_diff_year  = (year_crop  ~= Year_Should);
ind_diff_month = (month_crop ~= Month_Should);
ind_diff_day   = (day_crop   ~= Day_Should);


% Find difference
for i=1:size(YMDHMS,1)    
    ind_diff_doy(i)  = (abs(doy_crop(i)-DOY_Should(i))>0.001);
    ind_diff_hour(i) = (abs(hour_crop(i)-Hour_Should(i))>0.001);
    
    if ind_diff_year(i) == 1
        Message_forcings_view = ['Row ' num2str(i),' variable Year has a problem. It should be ', num2str(Year_Should(i)), '. Please check!'];
        msgbox(Message_forcings_view,'MLCan Error: File format','warn');
        return
    else
        if ind_diff_month(i) == 1
            Message_forcings_view = ['Row ' num2str(i),' variable Month has a problem. It should be ', num2str(Month_Should(i)), '. Please check!'];
            msgbox(Message_forcings_view,'MLCan Error: File format','warn');    
            return
        else
            if ind_diff_day(i) == 1
                Message_forcings_view = ['Row ' num2str(i),' variable Day has a problem. It should be ', num2str(Day_Should(i)), '. Please check!'];
                msgbox(Message_forcings_view,'MLCan Error: File format','warn');
                return
            else
                if ind_diff_doy(i) == 1
                    Message_forcings_view = ['Row ' num2str(i),' variable DOY has a problem. It should be ', num2str(DOY_Should(i)), '. Please check!'];
                    msgbox(Message_forcings_view,'MLCan Error: File format','warn');
                    return
                else
                    if ind_diff_hour(i) == 1
                        Message_forcings_view = ['Row ' num2str(i),' variable Hour has a problem. It should be ', num2str(Hour_Should(i)), '. Please check!'];
                        msgbox(Message_forcings_view,'MLCan Error: File format','warn');
                        return
                    end
                end
            end
            
        end
    end
end
msgbox('Time Steps are Correct!','MLCan message','warn')


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
