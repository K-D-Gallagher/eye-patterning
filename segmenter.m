function varargout = segmenter(varargin)
% segmenter MATLAB code for segmenter.fig
%      segmenter, by itself, creates a new segmenter or raises the existing
%      singleton*.
%
%      H = segmenter returns the handle to a new segmenter or the handle to
%      the existing singleton*.
%
%      segmenter('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in segmenter.M with the given input arguments.
%
%      segmenter('Property','Value',...) creates a new segmenter or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the segmenter before segmenter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segmenter_OpeningFcn via varargin.
%
%      *See segmenter Options on GUIDE's Tools menu.  Choose "segmenter allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segmenter

% Last Modified by GUIDE v2.5 31-Mar-2020 10:50:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segmenter_OpeningFcn, ...
                   'gui_OutputFcn',  @segmenter_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization code - also loads data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before segmenter is made visible.
function segmenter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segmenter (see VARARGIN)

% Choose default command line output for segmenter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%  load in raw images from workspace
handles.raw_images1 = evalin('base','raw_images1');
handles.raw_images2 = evalin('base','raw_images2');
handles.raw = handles.raw_images1;

handles.cell_area = evalin('base','cell_area');
handles.cell_centroid = evalin('base','cell_centroid');

% find maximum number of time points and store in handles.max
handles.max = size(handles.raw,3);

% feed the max number of images to the 'num_images' function
set(handles.num_images,'String',handles.max);

% set current time point to 1, for loading the data
handles.current = 1;

% load in watershed from workspace
handles.L = evalin('base','L');

% load in 'tracks', 'tracking display', and 'tracking_error_display'
handles.tracks = evalin('base','tracks');
handles.display_tracks = evalin('base','display_tracks');
ise = evalin( 'base', 'exist(''birth'',''var'') == 1' );
if ise
    handles.tracking_error = evalin('base','birth_death_display');
    handles.birth = evalin('base','birth');
    handles.death = evalin('base','death');
else
    handles.tracking_error = zeros(size(handles.L,1),size(handles.L,2),3,size(handles.L,3));
    handles.birth = cell(size(handles.L,3));
    handles.death = cell(size(handles.L,3));
end

% set initialization values of toggle check boxes
set(handles.toggle_raw,'Value',1)
set(handles.toggle_watershed,'Value',1)
set(handles.toggle_tracking_cues, 'Value',0)
set(handles.toggle_all_tracks, 'Value', 0)
set(handles.show_me_where, 'Value', 0)

%%%
% initialize undo helper - we'll have a maximum of 5 undo's
%%%
% create a counter to record which save number we're on
handles.undo_counter = 0;
% create a vector to record which time point each undo belong to
handles.undo_time_point = nan(5,1);
% create matrix to store back up copies of images to implement undo
% function.  first two dimensions record the actual image.  dimension 3
% records which undo copy we're on
handles.undo_helper = zeros(size(handles.L,1),size(handles.L,2),5);

% store max size
handles.max_x_size = size(handles.L,1);
handles.max_y_size = size(handles.L,2);

% initialize zoom & x and y axes limits

handles.zoom = 1;
if size(handles.L,1) > size(handles.L,2)
    handles.x_limit = [0.5,size(handles.L,1)+0.5]; 
    handles.y_limit = handles.x_limit;
else
    handles.y_limit = [0.5,size(handles.L,2)+0.5];
    handles.x_limit = handles.y_limit;
end

% store baseline x and y axes limits - by baseline, I mean the default zoom
% configuration
handles.baseline_x_limit = handles.x_limit;
handles.baseline_y_limit = handles.y_limit;

% initialize tracking cutoff
set(handles.cutoff,'String',num2str(15));
handles.cutoff_value = str2double(get(handles.cutoff,'String'));

% initialize birth/death count
handles.num_born = length(handles.birth{handles.current});
handles.num_dead = length(handles.death{handles.current});

% track links - (time, variables). variables are: cellID, centroid_X,
% centroid_Y, time, then all 4 variables repeated again for next cell.
ise = evalin( 'base', 'exist(''track_links'',''var'') == 1' );
if ise
    handles.track_links = evalin('base','track_links');
else
    handles.track_links = zeros(1,8);
end


update_arrows(handles);
guidata(hObject,handles)
display_layers(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% toggle raw layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in raw1.
function raw1_Callback(hObject, eventdata, handles)
% hObject    handle to raw1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of raw1

handles.raw = handles.raw_images1;

set(handles.raw2,'value',0)

display_layers(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes on button press in raw2.
function raw2_Callback(hObject, eventdata, handles)
% hObject    handle to raw2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of raw2

handles.raw = handles.raw_images2;

set(handles.raw1,'value',0)

display_layers(hObject, eventdata, handles);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resize GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set(findobj(handles.figure1),'Units','Normalized')
guidata(hObject,handles)


% --- Outputs from this function are returned to the command line.
function varargout = segmenter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Foward and back button display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Handles the display of the forward and back buttons.
function update_arrows(handles)
set(handles.current_image,'String', handles.current);
set(handles.birth_count,'String',handles.num_born);
set(handles.death_count,'String',handles.num_dead);
if (handles.current < handles.max) 
    set(handles.forward, 'Enable', 'On');
else
    set(handles.forward, 'Enable', 'Off');
end
if (handles.current > 1)
    set(handles.back, 'Enable', 'On');
else
    set(handles.back, 'Enable', 'Off');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Toggle on/off different layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Toggle Raw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in toggle_raw.
function toggle_raw_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_raw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_raw

% store state of toggle as a handle
current_state = get(hObject,'Value');
handles.toggle_raw.Value = current_state;

display_layers(hObject, eventdata, handles);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Toggle Watershed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in toggle_watershed.
function toggle_watershed_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_watershed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_watershed

% store state of toggle as a handle
current_state = get(hObject,'Value');
handles.toggle_watershed.Value = current_state;

display_layers(hObject, eventdata, handles)
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Toggle tracking cues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in toggle_tracking_cues.
function toggle_tracking_cues_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_tracking_cues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_tracking_cues

% store state of toggle as a handle
current_state = get(hObject,'Value');
handles.toggle_tracking_cues.Value = current_state;

% turn off all tracks layer so both can't be displayed at the same time
set(handles.toggle_all_tracks,'Value',0);

display_layers(hObject, eventdata, handles)
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Toggle all tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in toggle_all_tracks.
function toggle_all_tracks_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_all_tracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_all_tracks

% store state of toggle as a handle
current_state = get(hObject,'Value');
handles.toggle_all_tracks.Value = current_state;

% turn of tracking cues layer so both can't be displayed at the same time
set(handles.toggle_tracking_cues,'Value',0);

display_layers(hObject, eventdata, handles)
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINK TRACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select cell 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in cell1.
function cell1_Callback(hObject, eventdata, handles)
% hObject    handle to cell1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% record mouseclick position
% [x,y] = getpts(handles.axes1);
[x,y] = ginputax(gca,1);

% record current time
curr_time = handles.current;

ID = handles.L(round(y),round(x),curr_time)

handles.cell1_ID = ID;
handles.cell1_time = curr_time;

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select cell 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in cell2.
function cell2_Callback(hObject, eventdata, handles)
% hObject    handle to cell2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% record mouseclick position
% [x,y] = getpts(gca);
[x,y] = ginputax(gca,1);

% record current time
curr_time = handles.current;

ID = handles.L(round(y),round(x),curr_time)

handles.cell2_ID = ID;
handles.cell2_time = curr_time;

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% designate link
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in create_link.
function create_link_Callback(hObject, eventdata, handles)
% hObject    handle to create_link (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currlength = size(handles.track_links,1);
currlength = currlength + 1;

handles.track_links(currlength,1) = handles.cell1_ID;
handles.track_links(currlength,2) = handles.cell1_time;
handles.track_links(currlength,3) = handles.cell_centroid{handles.cell1_time}{handles.cell1_ID}(1);
handles.track_links(currlength,4) = handles.cell_centroid{handles.cell1_time}{handles.cell1_ID}(2);
handles.track_links(currlength,5) = handles.cell2_ID;
handles.track_links(currlength,6) = handles.cell2_time;
handles.track_links(currlength,7) = handles.cell_centroid{handles.cell2_time}{handles.cell2_ID}(1);
handles.track_links(currlength,8) = handles.cell_centroid{handles.cell2_time}{handles.cell2_ID}(2);

guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT, ADD, DELETE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit - this launches imfreehand to allow user to draw line that can then
% be selected for insertion or deletion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in edit
function edit_Callback(hObject, eventdata, handles)
% hObject    handle to edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% initialze imfreehand
hFH = imfreehand(gca,'Closed',false);

% extract positions and separate out x-coord and y-coord
xy = hFH.getPosition;
y1 = xy(:,1);
x1 = xy(:,2);

handles.x = x1;
handles.y = y1;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non zero element code from:
% https://www.mathworks.com/matlabcentral/answers/104614-grouping-continuous-nonzero-in-a-row-vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nonzero x elements
try
    wrapX       = [0, x1, 0] ;
catch
    wrapX       = [0, x1', 0] ;
end
tempX       = diff( wrapX > 0 ) ;
blockStartX = find( tempX == 1 ) + 1 ;
blockEndX   = find( tempX == -1 ) ;
blocksX     = arrayfun( @(bId) wrapX(blockStartX(bId):blockEndX(bId)), ...
    1:numel(blockStartX), 'UniformOutput', false ) ;
x1 = blocksX{1};

% nonzero y elements
try
    wrapY       = [0, y1, 0] ;
catch
    wrapY       = [0, y1', 0] ;
end
tempY       = diff( wrapY > 0 ) ;
blockStartY = find( tempY == 1 ) + 1 ;
blockEndY   = find( tempY == -1 ) ;
blocksY     = arrayfun( @(bId) wrapY(blockStartY(bId):blockEndY(bId)), ...
    1:numel(blockStartY), 'UniformOutput', false ) ;
y1 = blocksY{1};

% save as handles
handles.x = x1;
handles.y = y1;

% interpolate mouse positions form imfreehand
temp = interpolateFreeHand(x1,y1,handles.L(:,:,handles.current));

% store mask of interpolated mouse position in handles.temp
handles.temp = temp;

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function add_Callback(hObject, eventdata, handles)
% hObject    handle to add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%%%
% save Undo info - include condition for reaching max number of back up
% images, which requires redefining save number (i.e. undo 5 -> undo 4,
% undo 4 -> undo 3, etc)
%%%
if handles.undo_counter < 6 %&& handles.undo_counter > 0
    handles.undo_counter = handles.undo_counter + 1;
    handles.undo_time_point(handles.undo_counter) = handles.current;
    handles.undo_helper(:,:,handles.undo_counter) = handles.L(:,:,handles.current);
elseif handles.undo_counter > 0
    % redefine previous back up numbers, bumping them all up one number to
    % make room for the current back up image
    handles.undo_helper(:,:,1) = handles.undo_helper(:,:,2);
    handles.undo_helper(:,:,2) = handles.undo_helper(:,:,3);
    handles.undo_helper(:,:,3) = handles.undo_helper(:,:,4);
    handles.undo_helper(:,:,4) = handles.undo_helper(:,:,5);
    
    handles.undo_time_point(1) = handles.undo_time_point(2);
    handles.undo_time_point(2) = handles.undo_time_point(3);
    handles.undo_time_point(3) = handles.undo_time_point(4);
    handles.undo_time_point(4) = handles.undo_time_point(5);
    
    % back up current image and record time point to which it belongs
    handles.undo_helper(:,:,5) = handles.L(:,:,handles.current);
    handles.undo_time_point(5) = handles.current;
else % if we're at 0th copy of undo
     return
end

%%%
% make change
%%%

% pull out binary mask of imfreehand coordinates from temp handler
temp = handles.temp;

% pull out current segmentation mask from L handler
current_image = handles.L(:,:,handles.current);

% burn in modifications and re-watershed
merge = addCellsWatershed(temp,current_image);

% update L handle to reflect modified segmentation
handles.L(:,:,handles.current) = merge;

% clear handles.temp
handles.temp = 0;

% update handles
display_layers(hObject, eventdata, handles)
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function delete_Callback(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%%%
% save Undo info
%%%
if handles.undo_counter < 6 %&& handles.undo_counter > 0
    handles.undo_counter = handles.undo_counter + 1;
    handles.undo_time_point(handles.undo_counter) = handles.current;
    handles.undo_helper(:,:,handles.undo_counter) = handles.L(:,:,handles.current);
elseif handles.undo_counter > 0
    % redefine previous back up numbers, bumping them all up one number to
    % make room for the current back up image
    handles.undo_helper(:,:,1) = handles.undo_helper(:,:,2);
    handles.undo_helper(:,:,2) = handles.undo_helper(:,:,3);
    handles.undo_helper(:,:,3) = handles.undo_helper(:,:,4);
    handles.undo_helper(:,:,4) = handles.undo_helper(:,:,5);
    
    handles.undo_time_point(1) = handles.undo_time_point(2);
    handles.undo_time_point(2) = handles.undo_time_point(3);
    handles.undo_time_point(3) = handles.undo_time_point(4);
    handles.undo_time_point(4) = handles.undo_time_point(5);
    
    % back up current image and record time point to which it belongs
    handles.undo_helper(:,:,5) = handles.L(:,:,handles.current);
    handles.undo_time_point(5) = handles.current;
else % if we're at 0th copy of undo
     return
end

%%%
% make change
%%%

% pull out binary mask of imfreehand coordinates from temp handler
temp = handles.temp;

% pull out current segmentation mask from L handler
current_image = handles.L(:,:,handles.current);

% burn in modifications and re-watershed
merge = deleteCellsWatershed(temp,current_image);

% update L handle to reflect modified segmentation
handles.L(:,:,handles.current) = merge;

% clear handles.temp
handles.temp = 0;

% update handles
display_layers(hObject, eventdata, handles)
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Undo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in undo.
function undo_Callback(hObject, eventdata, handles)
% hObject    handle to undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.undo_counter ~= 0
    handles.L(:,:,handles.undo_time_point(handles.undo_counter)) = ...
        handles.undo_helper(:,:,handles.undo_counter);
    handles.undo_counter = handles.undo_counter - 1;
end

% update gui info and axes
display_layers(hObject, eventdata, handles)
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT_IMAGE, FOWARD, BACK, EXPORT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Current image - edit box (both displays current image and allows user to
% type different image number and travel there)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function current_image_Callback(hObject, eventdata, handles)
% hObject    handle to current_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of current_image as text
%        str2double(get(hObject,'String')) returns contents of
%        current_image as a double

% find number associated with current_image callback
new_num = str2num(handles.current_image.String);

% if there isn't a number associated with current_image callback, set
% handles.current_image equal to handles.current
if  isempty(new_num)
    set(handles.current_image,'String',handles.current);
    display_layers(hObject, eventdata, handles)
    return
end

% if the handle.current_image callback falls within the range of loaded
% images,then update handles.current == handles.current_image
if  new_num > 0 && new_num <= handles.max
    handles.current = new_num;
    display_layers(hObject, eventdata, handles)
else 
    set(handles.current_image,'String',handles.current);
    display_layers(hObject, eventdata, handles)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Current image - CreateFcn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function current_image_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Back
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
% hObject    handle to back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.current > 1)
    curr_time = handles.current - 1;
    handles.current = curr_time;
    handles.num_born = length(handles.birth{curr_time});
    handles.num_dead = length(handles.death{curr_time});
    
    display_layers(hObject, eventdata, handles)
    guidata(hObject, handles);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in forward.
function forward_Callback(hObject, eventdata, handles)
% hObject    handle to forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.current < handles.max)
    curr_time = handles.current + 1;
    handles.current = curr_time;
    handles.num_born = length(handles.birth{curr_time});
    handles.num_dead = length(handles.death{curr_time});
    handles.test = 1;
    
    display_layers(hObject, eventdata, handles)
    guidata(hObject, handles);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Birth/death stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function birth_count_Callback(hObject, eventdata, handles)
% hObject    handle to birth_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of birth_count as text
%        str2double(get(hObject,'String')) returns contents of birth_count as a double

function death_count_Callback(hObject, eventdata, handles)
% hObject    handle to death_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of death_count as text
%        str2double(get(hObject,'String')) returns contents of death_count as a double

% --- Executes on button press in show_me_where.
function show_me_where_Callback(hObject, eventdata, handles)
% hObject    handle to show_me_where (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_me_where

current_state = get(hObject,'Value');
handles.show_me.Value = current_state;

display_layers(hObject, eventdata, handles)
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export and refresh tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save to workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in save_workspace.
function save_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to save_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','L_GUI',handles.L);
assignin('base','tracks_GUI',handles.tracks);
assignin('base','display_tracks_GUI',handles.display_tracks);
assignin('base','track_links',handles.track_links);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in write_file.
function write_file_Callback(hObject, eventdata, handles)
% hObject    handle to write_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i=1:handles.max
    name = strcat(strcat('edited_watershed-', num2str(i)),'.tiff');
    imwrite(handles.L(:,:,i),name)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cutoff text box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoff as text
%        str2double(get(hObject,'String')) returns contents of cutoff as a double

handles.cutoff_value = str2double(get(handles.cutoff,'String'));

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refresh tracks
% This is kinda clunky code.  Rather than refering to functions outside of
% the GUI, I embedded the 'tracks.cells' and 'cellStats' functions within
% the GUI so I could use the counter for these functions in order to update
% a display that lets you know their progress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in refresh_tracks.
function refresh_tracks_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_tracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% read in image to display while everything updates
%construction_sign = imread('construction_sign.jpg');
update_image = uint8(zeros(size(handles.L,1),size(handles.L,2),3));
%left_right_diff = (size(handles.L,1) - size(construction_sign,1))/2;
%top_bottom_diff = (size(handles.L,2) - size(construction_sign,2))/2;
%update_image(left_right_diff:left_right_diff+size(construction_sign,1)-1,top_bottom_diff:top_bottom_diff+size(construction_sign,2)-1,1:3) = construction_sign;


% print message for timepoint 1
imshow(update_image);
text((size(handles.L,1)/3),(size(handles.L,2)/8),100,...
    'Preparing',...
    'Color', 'white', 'FontSize',20);
drawnow

% Prealocate cell arrays for storing centroid and PixelIdxList
handles.cell_area = cell(handles.max,1);
handles.cell_centroid = cell(handles.max,1);
for ii = 1:handles.max
    handles.cell_area{ii,1} = cell(length(unique(handles.L(:,:,ii))),1);
    handles.cell_centroid{ii,1} = cell(length(unique(handles.L(:,:,ii))),1);
end

% Set cutoff for tracking
cutoff = handles.cutoff_value;

% Prealocate cell array for storing hungarian tracking maps
cellPairs = cell(size(handles.L,3),1);

% Update cell area and centroid for time point 1 outside of for loop
stats = regionprops(handles.L(:,:,1),'centroid','PixelIdxList');
for iii = 1:size(stats,1)
    handles.cell_centroid{1,1}{iii,1} = stats(iii).Centroid;
    [x_coord, y_coord] = ind2sub(size(handles.L,1),stats(iii).PixelIdxList);
    handles.cell_area{1,1}{iii,1} = [x_coord,y_coord];
end

% Record cell area and centroid info of first time point into variables for
% the tracking script below
S0 = stats;
R0 = vertcat(S0.Centroid);

% do remainder of time points in loop and concurrently calculate cell
% area, centroid, and tracks
for t = 2:handles.max
    
    % regionprops
    stats = regionprops(handles.L(:,:,t),'PixelIdxList','Centroid');
    R = vertcat(stats.Centroid);
    
    %%%
    % record cell area and centroid in cell arrays
    %%%
    for iii = 1:size(stats,1)
        handles.cell_centroid{t,1}{iii,1} = stats(iii).Centroid;
        [x_coord, y_coord] = ind2sub(size(handles.L,1),stats(iii).PixelIdxList);
        handles.cell_area{t,1}{iii,1} = [x_coord,y_coord];
    end
    
    %%%
    % tracking
    %%%
    % Pre-filtering step based upon distance.
    D = pdist2(R,R0);
    [currentCells,oldCells] = find(D < cutoff);
    
    Ov = inf*ones(size(D));
    for ii = 1:length(oldCells)
        sharedArea = sum(ismembc(stats(currentCells(ii)).PixelIdxList,S0(oldCells(ii)).PixelIdxList));
        
        Ov(currentCells(ii),oldCells(ii)) = ...
            (length(stats(currentCells(ii)).PixelIdxList) + length(S0(oldCells(ii)).PixelIdxList))/sharedArea - 2;
    end
    
    [ M ] = track.munkres(Ov);
    matchedCells = find(M);
    cellPairs{t} = zeros(length(stats),1);
    cellPairs{t}(matchedCells) = M(matchedCells);
    
    S0 = stats;
    R0 = R;
    
    %%%
    % print a message
    %%%
    imshow(update_image);
    text((size(handles.L,1)/3),(size(handles.L,2)/8),100,...
        sprintf('Tracking timepoint %s',num2str(t)),...
        'Color', 'white', 'FontSize',20);
    drawnow

    
end

handles.cellPairs = cellPairs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manually link tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 2:size(handles.track_links,1)
    
    %%%
    % print a message
    %%%
    imshow(update_image);
    text((size(handles.L,1)/3),(size(handles.L,2)/8),100,...
        sprintf('Manually linking track %s',num2str(j)),...
        'Color', 'white', 'FontSize',20);
    drawnow
    
    % check which time point is larger
    if handles.track_links(j,2) > handles.track_links(j,6)
        test1 = 0;
    elseif handles.track_links(j,2) < handles.track_links(j,6)
        test1 = 1;
    end
    
    % pull out centroid components to ID cells
    if test1
        cell1_x = handles.track_links(j,3);
        cell1_y = handles.track_links(j,4);
        cell2_x = handles.track_links(j,7);
        cell2_y = handles.track_links(j,8);
        cell1_time = handles.track_links(j,2);
        cell2_time = handles.track_links(j,6);
    else
        cell1_x = handles.track_links(j,7);
        cell1_y = handles.track_links(j,8);
        cell2_x = handles.track_links(j,3);
        cell2_y = handles.track_links(j,4);
        cell1_time = handles.track_links(j,6);
        cell2_time = handles.track_links(j,2);
    end
    
    handles.cell1_x = cell1_x;
    handles.cell1_y = cell1_y;
    handles.cell2_x = cell2_x;
    handles.cell2_y = cell2_y;
    handles.cell1_time = cell1_time;
    handles.cell2_time = cell2_time;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find cell IDs based on centroid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % cell1
    for jj = 1:length(handles.cell_centroid{cell1_time})-1
        centpair = [cell1_x cell1_y];
        compareQuant = handles.cell_centroid{cell1_time}{jj};
        handles.centpair = centpair;
        handles.compareQuant = compareQuant;
        handles.jj = jj;
        guidata(hObject, handles);
        if compareQuant == centpair
            cell1_ID = jj;
        end
    end
    
    % cell2
    for jj = 1:length(handles.cell_centroid{cell2_time})-1
        centpair = [cell2_x cell2_y];
        compareQuant = handles.cell_centroid{cell2_time}{jj};
        if compareQuant == centpair
            cell2_ID = jj;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    % update track entries
    %%%%%%%%%%%%%%%%%%%%%%
    
    cellPairs{cell2_time}(cell2_ID) = cell1_ID;
    
    % clear all other cells that are linked to cell1_ID
    tempCellPairs = zeros(length(cellPairs{cell2_time}));
    for e = 1:length(cellPairs{cell2_time})
        if e ~= cell2_ID
            if cellPairs{cell2_time}(e) == cell1_ID
                cellPairs{cell2_time}(e) = 0;
            end
        end
    end
    
end

handles.cellPairs = cellPairs;


% display image for final steps
imshow(update_image);
text((size(handles.L,1)/3),(size(handles.L,2)/8),100,...
    'Finalizing', 'Color', 'white', 'FontSize',20);
drawnow

% generate new tracking matrix
[handles.tracks, handles.display_tracks] = ...
    generateTracksGUI(handles.L, handles.cell_area, cellPairs);

% generate new tracking error display
[handles.tally, handles.tracking_error] = ...
    trackingErrorDisplay(handles.L,handles.tracks,handles.cell_area);


guidata(hObject, handles);

update_arrows(handles);
guidata(hObject, handles);
display_layers(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zoom & Pan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zoom in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in zoom_in.
function zoom_in_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% read in baseline zoom configuration
baseline_x_limit = handles.baseline_x_limit;
baseline_y_limit = handles.baseline_y_limit;

% read in current axes limits
current_x_limit = handles.x_limit;
current_y_limit = handles.y_limit;

% read in current zoom state
current_zoom = handles.zoom;

if current_zoom < 5
    
    % read in current zoom configuration and update as function of callback
    handles.zoom = current_zoom + 1;
    
    % find new axes limit range
    new_x_limits = current_x_limit .* 0.9;
    new_y_limits = current_y_limit .* 0.9;
    
    % find upper and lower limits of new axes range
    x_lower = new_x_limits(1);
    x_upper = new_x_limits(2);
    y_lower = new_y_limits(1);
    y_upper = new_y_limits(2);
    
    % find total axes limit range
    x_range = x_upper - x_lower;
    y_range = y_upper - y_lower;
    
    % find offset from boundaries of total image size
    x_displacement = current_x_limit(2) - x_range;
    y_displacement = current_y_limit(2) - y_range;
    
    % determine center position for new axes limits
    new_zoom_x = [x_displacement/2 + current_x_limit(1), current_x_limit(2) - x_displacement/2];
    new_zoom_y = [y_displacement/2 + current_y_limit(1), current_y_limit(2) - y_displacement/2];
    
    % update handles
    handles.x_limit = new_zoom_x;
    handles.y_limit = new_zoom_y;
    
else
    
    handles.x_limit = current_x_limit;
    handles.y_limit = current_y_limit;
    
end

% update GUI info
guidata(hObject, handles);
display_layers(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zoom out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in zoom_out.
function zoom_out_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% read in baseline zoom configuration
baseline_x_limit = handles.baseline_x_limit;
baseline_y_limit = handles.baseline_y_limit;

% read in current axes limits
old_x_limits = handles.x_limit;
old_y_limits = handles.y_limit;

% read in current zoom state
current_zoom = handles.zoom;

% update zoom state based on callback
if current_zoom ~= 0
    
    handles.zoom = current_zoom - 1;
    
    % find new axes limit range
    new_x_limits = old_x_limits .* (1/0.9);
    new_y_limits = old_y_limits .* (1/0.9);
    
    % find upper and lower limits of new axes range
    x_lower = new_x_limits(1);
    x_upper = new_x_limits(2);
    y_lower = new_y_limits(1);
    y_upper = new_y_limits(2);
    
    % find total axes limit range
    x_range = x_upper - x_lower;
    y_range = y_upper - y_lower;
    
    % find old axes limit range
    x_range_old = old_x_limits(2) - old_x_limits(1);
    y_range_old = old_y_limits(2) - old_y_limits(1);
    
    % find offset from boundaries of total image size
    x_displacement = x_range - x_range_old;
    y_displacement = y_range - y_range_old;
    
    % determine center position for new axes limits
    new_zoom_x = [old_x_limits(1) - x_displacement/2 , old_x_limits(2) + x_displacement/2];
    new_zoom_y = [old_y_limits(1) - y_displacement/2, old_y_limits(2) + y_displacement/2];
    
    % update handles
    handles.x_limit = new_zoom_x;
    handles.y_limit = new_zoom_y;
    
else
    
    handles.zoom = 1;
    handles.x_limit = baseline_x_limit;
    handles.y_limit = baseline_y_limit;
    
end

guidata(hObject, handles);
display_layers(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pan.
function pan_Callback(hObject, eventdata, handles)
% hObject    handle to pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pan toggle

new_x_limit = get(gca,'XLim');
new_y_limit = get(gca,'YLim');

handles.x_limit = new_x_limit;
handles.y_limit = new_y_limit;

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function display_layers(hObject, eventdata, handles)
% hObject    handle to watershed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cur = handles.current;
x_resolution = size(handles.L,1);
y_resolution = size(handles.L,2);

% read in zoom state and axes limits
x_limit = handles.x_limit;
y_limit = handles.y_limit;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAW, WATERSHED, TRACKING CUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if handles.toggle_raw.Value && ...
   handles.toggle_watershed.Value && ...
   handles.toggle_tracking_cues.Value && ... 
   not(handles.toggle_all_tracks.Value)

    %%%
    % Raw + watershed
    %%%
    %find watershed outline and initialize new rgb image
    [row,col] = find(not(handles.L(:,:,cur)));
    rgbImage = cat(3, handles.raw(:,:,cur), handles.raw(:,:,cur), handles.raw(:,:,cur));
    
    % display watershed in red
    for i = 1:length(row)
        rgbImage(row(i), col(i), 1) = 255;
        rgbImage(row(i), col(i), 2) = 0;
        rgbImage(row(i), col(i), 3) = 0;
    end
    
    % set transparency of tracking error cues and display over
    % raw/watershed
    transparency = 1 - im2double(handles.raw(:,:,cur)) - 0.4;
    imshow(rgbImage);
    xlim(x_limit);
    ylim(y_limit);
    hold on
    overlay = imshow(handles.tracking_error(:,:,:,cur));
    set(overlay,'AlphaData',transparency);
    
    
    % Store GUI info
    update_arrows(handles);
    guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAW, WATERSHED, ALL TRACKS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif handles.toggle_raw.Value && ... 
       handles.toggle_watershed.Value && ...
       not(handles.toggle_tracking_cues.Value) && ...
       handles.toggle_all_tracks.Value
   
    % find watershed pixels
    [row,col] = find(not(handles.L(:,:,cur)));
    
    % load in raw data as RGB
    rgbImage = cat(3, handles.raw(:,:,cur), handles.raw(:,:,cur), handles.raw(:,:,cur));
    
    % burn watershed into raw image in red
    for i = 1:length(row)
        rgbImage(row(i), col(i), 1) = 255;
        rgbImage(row(i), col(i), 2) = 0;
        rgbImage(row(i), col(i), 3) = 0;
    end
    
    % set transparency for overlay
    transparency = 1 - im2double(handles.raw(:,:,cur)) - 0.4;
    %transparency = (ones(x_resolution,y_resolution)/4);
    
    % display images including overlay with transparency
    imshow(rgbImage);
    xlim(x_limit);
    ylim(y_limit);
    hold on
    h = imshow(handles.display_tracks(:,:,:,cur));
    set(h, 'AlphaData',transparency);
    
    % update GUI
    update_arrows(handles);
    guidata(hObject, handles);
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAW, TRACKING CUES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif handles.toggle_raw.Value && ...
        not(handles.toggle_watershed.Value) && ...
        handles.toggle_tracking_cues.Value && ...
        not(handles.toggle_all_tracks.Value)
    
    % load in raw image
    rgbImage = cat(3, handles.raw(:,:,cur), handles.raw(:,:,cur), handles.raw(:,:,cur));
    
    % set transparency for overlaying tracking cues
    transparency = 1 - im2double(handles.raw(:,:,cur)) - 0.4;
    
    % display images and overlay tracking cues with transparency
    imshow(rgbImage);
    xlim(x_limit);
    ylim(y_limit);
    hold on
    overlay = imshow(handles.tracking_error(:,:,:,cur));
    set(overlay,'AlphaData',transparency);
    
    % update GUI data
    update_arrows(handles);
    guidata(hObject, handles);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAW, ALL TRACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif handles.toggle_raw.Value && ...
       not(handles.toggle_watershed.Value) && ...
       not(handles.toggle_tracking_cues.Value) && ... 
       handles.toggle_all_tracks.Value
   
   % load in raw image
   rgbImage = cat(3, handles.raw(:,:,cur), handles.raw(:,:,cur), handles.raw(:,:,cur));
   
   % set transparency for overlay
   transparency = 1 - im2double(handles.raw(:,:,cur)) - 0.4;
   
   % display images including overlay with transparency
   imshow(rgbImage);
   xlim(x_limit);
   ylim(y_limit);
   hold on
   h = imshow(handles.display_tracks(:,:,:,cur));
   set(h, 'AlphaData',transparency);
   
   % update GUI data
   update_arrows(handles);
   guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WATERSHED, TRACKING CUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif not(handles.toggle_raw.Value) && ...
       handles.toggle_watershed.Value && ...
       handles.toggle_tracking_cues.Value && ... 
       not(handles.toggle_all_tracks.Value)
   
   % create stand alone rgb image of watershed shown in red
   rgbImage = zeros(size(handles.raw(:,:,cur),1),size(handles.raw(:,:,cur),2),3);
   [row,col] = find(not(handles.L(:,:,cur)));
   for i = 1:length(row)
       rgbImage(row(i), col(i), 1) = 255;
       rgbImage(row(i), col(i), 2) = 0;
       rgbImage(row(i), col(i), 3) = 0;
   end
   
   % set transparency for overlaying tracking cues
   transparency = (not(isnan(handles.tracking_error(:,:,:,cur))));
   transparency = transparency(:,:,1);
   
   % display images and overlay tracking cues with transparency
   imshow(rgbImage);
   xlim(x_limit);
   ylim(y_limit);
   hold on
   overlay = imshow(handles.tracking_error(:,:,:,cur));
   set(overlay,'AlphaData',transparency);
   
   % update GUI data
   update_arrows(handles);
   guidata(hObject, handles);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WATERSHED, ALL TRACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif not(handles.toggle_raw.Value) && ...
       handles.toggle_watershed.Value && ...
       not(handles.toggle_tracking_cues.Value) && ... 
       handles.toggle_all_tracks.Value
   
   % create stand alone rgb image of watershed shown in red
   rgbImage = handles.display_tracks(:,:,:,cur);
   [row,col] = find(not(handles.L(:,:,cur)));
   for i = 1:length(row)
       rgbImage(row(i), col(i), 1) = 255;
       rgbImage(row(i), col(i), 2) = 0;
       rgbImage(row(i), col(i), 3) = 0;
   end
   
   imshow(rgbImage);
   xlim(x_limit);
   ylim(y_limit);
   hold on
   
   % update GUI data
   update_arrows(handles);
   guidata(hObject, handles);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRACKING CUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif not(handles.toggle_raw.Value) && ...
       not(handles.toggle_watershed.Value) && ...
       handles.toggle_tracking_cues.Value && ... 
       not(handles.toggle_all_tracks.Value)
   
   imshow(handles.tracking_error(:,:,:,cur));
   xlim(x_limit);
   ylim(y_limit);
   hold on
   update_arrows(handles);
   guidata(hObject,handles);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALL TRACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif not(handles.toggle_raw.Value) && ...
       not(handles.toggle_watershed.Value) && ...
       not(handles.toggle_tracking_cues.Value) && ... 
       handles.toggle_all_tracks.Value
   
   imshow(handles.display_tracks(:,:,:,cur));
   xlim(x_limit);
   ylim(y_limit);
   hold on
   update_arrows(handles);
   guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAW & WATERSHED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif handles.toggle_raw.Value && ...
       handles.toggle_watershed.Value && ...
       not(handles.toggle_tracking_cues.Value) && ...
       not(handles.toggle_all_tracks.Value)
    
    % find watershed outline and initialize new rgb image
    [row,col] = find(not(handles.L(:,:,cur)));
    rgbImage = cat(3, handles.raw(:,:,cur), handles.raw(:,:,cur), handles.raw(:,:,cur));
    
    % display watershed in red
    for i = 1:length(row)
        rgbImage(row(i), col(i), 1) = 255;
        rgbImage(row(i), col(i), 2) = 0;
        rgbImage(row(i), col(i), 3) = 0;
    end
    
    % plot image
    imshow(rgbImage);
    xlim(x_limit);
    ylim(y_limit);
    update_arrows(handles);
    guidata(hObject, handles);
    hold on
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WATERSHED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif handles.toggle_watershed.Value && ... 
       not(handles.toggle_raw.Value) && ...
       not(handles.toggle_tracking_cues.Value) && ... 
       not(handles.toggle_all_tracks.Value)
   
    rgbImage = zeros(size(handles.raw(:,:,cur),1),size(handles.raw(:,:,cur),2),3);
    [row,col] = find(not(handles.L(:,:,cur)));
    for i = 1:length(row)
        rgbImage(row(i), col(i), 1) = 255;
        rgbImage(row(i), col(i), 2) = 0;
        rgbImage(row(i), col(i), 3) = 0;
    end
    imshow(rgbImage);
    xlim(x_limit);
    ylim(y_limit);
    hold on
    update_arrows(handles);
    guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif handles.toggle_raw.Value && ... 
       not(handles.toggle_watershed.Value) && ...
       not(handles.toggle_tracking_cues.Value) && ... 
       not(handles.toggle_all_tracks.Value)
   
    rgbImage = cat(3, handles.raw(:,:,cur), handles.raw(:,:,cur), handles.raw(:,:,cur));
    imshow(rgbImage)
    xlim(x_limit);
    ylim(y_limit);
    hold on
    update_arrows(handles);
    guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTHING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif not(handles.toggle_watershed.Value) && ... 
       not(handles.toggle_raw.Value) && ...
       not(handles.toggle_tracking_cues.Value) && ... 
       not(handles.toggle_all_tracks.Value)
   
    rgbImage = zeros(size(handles.raw(:,:,cur),1),size(handles.raw(:,:,cur),2),3);
    imshow(rgbImage)
    hold on
    update_arrows(handles);    
    guidata(hObject, handles);
end

if handles.show_me_where.Value
    
    curr_birth = handles.birth{handles.current};
    curr_death = handles.death{handles.current};
    
    cell_centroid = handles.cell_centroid;
    
    if curr_birth > 0
    for j = 1:length(curr_birth)
        plot(cell_centroid{cur}{curr_birth(j)}(1),cell_centroid{cur}{curr_birth(j)}(2),'*c','MarkerSize',20);
    end
    end
    
    if curr_death > 0
    for j = 1:length(curr_death)
        plot(cell_centroid{cur}{curr_death(j)}(1),cell_centroid{cur}{curr_death(j)}(2),'*r','MarkerSize',20);
    end
    end
    
end
hold off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keyboard shortcuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

 % determine the key that was pressed 
 keyPressed = eventdata.Key;
 
 % edit shortcut
 if strcmpi(keyPressed,'e')
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % freehand
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % initialze imfreehand
     hFH = imfreehand(gca,'Closed',false);
     
     % extract positions and separate out x-coord and y-coord
     xy = hFH.getPosition;
     y1 = xy(:,1);
     x1 = xy(:,2);
     
     handles.x_original = x1;
     handles.y_original = y1;
     guidata(hObject, handles);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % non zero element code from:
     % https://www.mathworks.com/matlabcentral/answers/104614-grouping-continuous-nonzero-in-a-row-vector
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     %%%%%%%%%%%%%%%%%%%%
     % nonzero x elements
     %%%%%%%%%%%%%%%%%%%%
     try
         wrapX       = [0, x1', 0] ;
     catch
         wrapX       = [0, x1, 0] ;
     end
     tempX       = diff( wrapX > 1 ) ;
     blockStartX = find( tempX == 1 ) + 1;
     blockEndX   = find( tempX == -1 );
     
     %%%%%%%%%%%%%%%%%%%%
     % nonzero y elements
     %%%%%%%%%%%%%%%%%%%%
     try
         wrapY       = [0, y1', 0] ;
     catch
         wrapY       = [0, y1, 0] ;
     end
     tempY       = diff( wrapY > 1 ) ;
     blockStartY = find( tempY == 1 ) + 1;
     blockEndY   = find( tempY == -1 );
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % find larger starting value
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if blockStartX > blockStartY
         totalStart = blockStartX;
     else
         totalStart = blockStartY;
     end
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%
     % find smaller ending value
     %%%%%%%%%%%%%%%%%%%%%%%%%%%
     if blockEndX < blockEndY
         totalEnd = blockEndX;
     else
         totalEnd = blockEndY;
     end
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % determine whether we're at the edge of the image and need to add
     % another indice or not to ensure watershed executes properly along
     % boundary
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     startX = false;
     endX = false;
     startY = false;
     endY = false;
     
     if blockEndX ~= (length(x1) + 1)
         if blockEndX ~= 2
             endX = true;
         end
     end
     
     if blockEndY ~= (length(y1) + 1)
         if blockEndY ~= 2
             endY = true;
         end
     end
     
     if blockStartX ~= (length(x1) + 1)
         if blockStartX ~= 2
             startX = true;
         end
     end
     
     if blockStartY ~= (length(y1) + 1)
         if blockStartY ~= 2
             startY = true;
         end
     end
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % pull out nonzero indices
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     blocksX     = arrayfun( @(bId) wrapX(totalStart(bId):totalEnd(bId)), ...
         1:numel(totalStart), 'UniformOutput', false ) ;
     
     blocksY     = arrayfun( @(bId) wrapY(totalStart(bId):totalEnd(bId)), ...
         1:numel(totalStart), 'UniformOutput', false ) ;
     
     % store temp version before checking if we're at boundary
     temp_x1 = blocksX{1};
     temp_y1 = blocksY{1};
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % if at boundary, add 1 to proper dimension so watershed executes
     % properly
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     if startX         
         try
             temp_x1 = [1, temp_x1];
             temp_y1 = [temp_y1(1), temp_y1];
         catch
             temp_x1 = [1, temp_x1'];
             temp_y1 = [temp_y1(1), temp_y1'];
         end
     elseif endX
         try
             temp_x1 = [temp_x1, 1];
             temp_y1 = [temp_y1,temp_y1(end)];
         catch
             temp_x1 = [temp_x', 1];
             temp_y1 = [temp_y1',temp_y1(end)];
         end
     elseif startY
         try
             temp_x1 = [temp_x1(1),temp_x1];
             temp_y1 = [1, temp_y1];
         catch
             temp_x1 = [temp_x1(1), temp_x1'];
             temp_y1 = [1, temp_y1'];
         end
     elseif endY
         try
             temp_x1 = [temp_x1,temp_x1(end)];
             temp_y1 = [temp_y1,1];
         catch
             temp_x1 = [temp_x1',temp_x1(end)];
             temp_y1 = [temp_y1',1];
         end
     end
     x1 = temp_x1;
     y1 = temp_y1;
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % save as handles & interpolate
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if startX || startY || endX || endY
         extra_dilation = true;
     else
         extra_dilation = false;
     end
     handles.extra_dilation = extra_dilation;
     handles.startX = startX;
     handles.endX = endX;
     handles.startY = startY;
     handles.endY = endY;
     
     handles.x = x1;
     handles.y = y1;
     guidata(hObject, handles);
     
     % interpolate mouse positions form imfreehand
     if extra_dilation
         temp = interpolateFreeHand_extraDilation(x1,y1,handles.L(:,:,handles.current),...
             handles.extra_dilation,handles.startX,handles.endX,handles.startY,handles.endY);
     else
         temp = interpolateFreeHand(x1,y1,handles.L(:,:,handles.current));
     end
     
     
     % store mask of interpolated mouse position in handles.temp
     handles.temp = temp;
     
     guidata(hObject, handles);
 
 % add shortcut
 elseif strcmpi(keyPressed,'a')
     % set focus to the button
     uicontrol(handles.add);
     % call the callback
     add_Callback(handles.add,[],handles);
     
 % delete shortcut
 elseif strcmpi(keyPressed,'d')
     % set focus to the button
     uicontrol(handles.delete);
     % call the callback
     delete_Callback(handles.delete,[],handles);
     
 % raw layer
 elseif strcmpi(keyPressed,'1')
     
     % grap current state
     current_state = get(handles.toggle_raw,'Value');
     
     % flip current state
     if current_state == 1
         set(handles.toggle_raw,'Value',0);
     else
         set(handles.toggle_raw,'Value',1);
     end
     
     % update GUI info and axes
     display_layers(hObject, eventdata, handles);
     guidata(hObject, handles);
     
 % watershed layer
 elseif strcmpi(keyPressed,'2')
     
     % grap current state
     current_state = get(handles.toggle_watershed,'Value');
     
     % flip current state
     if current_state == 1
         set(handles.toggle_watershed,'Value',0);
     else
         set(handles.toggle_watershed,'Value',1);
     end
     
     % update GUI info and axes
     display_layers(hObject, eventdata, handles);
     guidata(hObject, handles);
     
 % tracking cues layer
 elseif strcmpi(keyPressed,'3')
     
     % grap current state
     current_state = get(handles.toggle_tracking_cues,'Value');
     
     % flip current state
     if current_state == 1
         set(handles.toggle_tracking_cues,'Value',0);
     else
         set(handles.toggle_tracking_cues,'Value',1);
     end
     
     % turn off all tracks layer so both can't be displayed at the same
     % time
     set(handles.toggle_all_tracks,'Value',0);
     
     % update GUI info and axes
     display_layers(hObject, eventdata, handles);
     guidata(hObject, handles);
     
 % all tracks layer
 elseif strcmpi(keyPressed,'4')
     
     % grap current state
     current_state = get(handles.toggle_all_tracks,'Value');
     
     % flip current state
     if current_state == 1
         set(handles.toggle_all_tracks,'Value',0);
     else
         set(handles.toggle_all_tracks,'Value',1);
     end
     
     % turn off tracking cues layer so both can't be displayed at the same
     % time
     set(handles.toggle_tracking_cues,'Value',0);
     
     % update GUI info and axes
     display_layers(hObject, eventdata, handles);
     guidata(hObject, handles);
     
 % forward time
 elseif strcmpi(keyPressed,'l')
     % set focus to the button
     uicontrol(handles.forward);
     % call the callback
     forward_Callback(handles.forward,[],handles);
     
 % backward time
 elseif strcmpi(keyPressed,'k')
     % set focus to the button
     uicontrol(handles.back);
     % call the callback
     back_Callback(handles.back,[],handles);
 
 % zoom in
 elseif strcmpi(keyPressed,'w')
     % set focus to the button
     uicontrol(handles.zoom_in);
     % call the callback
     zoom_in_Callback(handles.zoom_in,[],handles);
     
 % zoom out
 elseif strcmpi(keyPressed,'q')
     % set focus to the button
     uicontrol(handles.zoom_out);
     % call the callback
     zoom_out_Callback(handles.zoom_out,[],handles);
    
 % pan
 elseif strcmpi(keyPressed,'z')
     
     % toggle pan state
     pan toggle
     
     % Enabling 'pan' activates WindowListener's that will block inputs
     % from WindowKeyPressFnc, so we need to disable these.  The following
     % code does exactly this.  I learned how to do this from the following
     % post: http://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan
     
     % Set focus to figure1
     fig = gcf;
     
     % get proper handles
     hManager = uigetmodemanager(fig);
     
     % the following logic loop makes the code backwards compatible with
     % older version of matlab.  previous, WindowListenerHandles had string
     % options (i.e. try...) for whether it was enable or not.  In post
     % 2014.b versions, it has boolian options (i.e. catch...)
     try
         set(hManager.WindowListenerHandles, 'Enable', 'off');  % HG1
     catch
         [hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
     end
     
     % enable the WindowsKeyPressFcn (which this is all embedded w/in!)
     set(handles.figure1, 'WindowKeyPressFcn', ...
         @(hObject,eventdata)kgGUI('figure1_WindowKeyPressFcn', ...
         hObject,eventdata,guidata(hObject)));
     
     % store new axes limits
     new_x_limit = get(gca,'XLim');
     new_y_limit = get(gca,'YLim');
     
     % update handles
     handles.x_limit = new_x_limit;
     handles.y_limit = new_y_limit;
     
     % update GUI info
     guidata(hObject, handles);
     
      % show_me_where
 elseif strcmpi(keyPressed,'t')
     
     % grap current state
     current_state = get(handles.show_me_where,'Value');
     
     % flip current state
     if current_state == 1
         set(handles.show_me_where,'Value',0);
     else
         set(handles.show_me_where,'Value',1);
     end
     
     % update GUI info and axes
     display_layers(hObject, eventdata, handles);
     guidata(hObject, handles);
     
 end



% --- Executes during object creation, after setting all properties.
function cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function birth_count_CreateFcn(hObject, eventdata, handles)
% hObject    handle to birth_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function death_count_CreateFcn(hObject, eventdata, handles)
% hObject    handle to death_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in debug.
function debug_Callback(hObject, eventdata, handles)
% hObject    handle to debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard
