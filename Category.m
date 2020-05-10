function varargout = Category(varargin)
% CATEGORY MATLAB code for Category.fig
%      CATEGORY, by itself, creates a new CATEGORY or raises the existing
%      singleton*.
%
%      H = CATEGORY returns the handle to a new CATEGORY or the handle to
%      the existing singleton*.
%
%      CATEGORY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CATEGORY.M with the given input arguments.
%
%      CATEGORY('Property','Value',...) creates a new CATEGORY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Category_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Category_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Category

% Last Modified by GUIDE v2.5 07-Nov-2015 23:54:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Category_OpeningFcn, ...
                   'gui_OutputFcn',  @Category_OutputFcn, ...
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
end


% --- Executes just before Category is made visible.
function Category_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Category (see VARARGIN)

% Choose default command line output for Category
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Category wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
end
function varargout = Category_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = 'culture.xlsx';
x=xlsread(filename);
x(:,2:10) = [];
set(handles.uitable2,'data',x);

xmean=mean(x);
xnew = (x-xmean);

axes(handles.axes1);
plot(x, 'o');
title('Original Data');

covariancematrix=cov(xnew);
[V,D] = eig(covariancematrix);
D=diag(D);
maxeigval=V(:,find(D==max(D)));

finaldata=maxeigval'*xnew';
set(handles.uitable1,'data',finaldata);
axes(handles.axes6);
stem(finaldata, 'DisplayName', 'finaldata');
title('PCA 1D output ');
axes(handles.axes7);
hold on
title('Final Classification')
for i=1:size(finaldata,2)
    if  finaldata(i)>=0
        plot(x(i),'o') 
        plot(x(i),'r*')
    else
        plot(x(i),'o')
        plot(x(i),'g*')
    end
end
hold off
set(handles.edit1, 'String', sprintf('Classification of "region" using PCA\n\nMean variable: %0.1f\ncovariancematrix %0.1f\nMaximum Eigenvectors %0.1f', xmean, covariancematrix, maxeigval));
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
end
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
