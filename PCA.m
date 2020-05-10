function varargout = PCA(varargin)
% PCA MATLAB code for PCA.fig
%      PCA, by itself, creates a new PCA or raises the existing
%      singleton*.
%
%      H = PCA returns the handle to a new PCA or the handle to
%      the existing singleton*.
%
%      PCA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCA.M with the given input arguments.
%
%      PCA('Property','Value',...) creates a new PCA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PCA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PCA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PCA

% Last Modified by GUIDE v2.5 08-Nov-2015 23:10:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PCA_OpeningFcn, ...
                   'gui_OutputFcn',  @PCA_OutputFcn, ...
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


% --- Executes just before PCA is made visible.
function PCA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PCA (see VARARGIN)

% Choose default command line output for PCA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PCA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
end
function varargout = PCA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
end

function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = 'EM_Algo.xlsx';
Y = xlsread(file);
set(handles.uitable5, 'data', Y);
filename = 'culture.xlsx';
X = xlsread(filename);
X(:,[1]) = [];
[coeff1] = pca(X,'algorithm','als');
set(handles.uitable1,'data',coeff1);
[pc,W,data_mean,xr,evals,percentVar] = ppca(X,3);
set(handles.uitable3, 'data', xr);
set(handles.uitable2,'data',pc);
axes(handles.axes1);
plot3(pc(1,:),pc(2,:),pc(3,:),'.');
title('{\bf PCA}');
xlabel(['PC 1 (',num2str(round(percentVar(1)*10)/10),'%)',]);
ylabel(['PC 2 (',num2str(round(percentVar(2)*10)/10),'%)',]);
zlabel(['PC 3 (',num2str(round(percentVar(3)*10)/10),'%)',]);
I = imread('OUTPUT1.jpg');
axes(handles.axes2);
imshow(I);
J = imread('OUTPUT2.jpg');
axes(handles.axes3);
imshow(J);
C = textscanu('OUTPUT.txt', 'UTF-16', 9, 13, 'waitbar');
u = uicontrol(handles.edit3);
set(u , 'String', C);
end

function [pc,W,data_mean,xr,evals,percentVar]=ppca(data,k)
 
  [C,ss,M,X,Ye]=ppca_mv(data',k,0,0);
  xr=Ye';
  W=C';
  data_mean=M';
  pc=X';
  
  for i=1:size(data,1)  
   v(i)=var(data(i,~isnan(data(i,:)))); 
  end
  total_variance=sum(v(~isnan(v)));
  
  evals=nan(1,k);
  for i=1:k 
    data_recon = (pinv(W(i,:))*pc(i,:)); % without mean correction (does not change the variance)
    evals(i)=sum(var(data_recon'));
  end
  
  percentVar=evals./total_variance*100;
  
%    cumsumVarPC=nan(1,k);  
%   for i=1:k 
%     data_recon = (pinv(W(1:i,:))*pc(1:i,:)) + repmat(data_mean,1,size(data,2));
%     cumsumVarPC(i)=sum(var(data_recon')); 
%   end
%   cumsumVarPC
end  

function [C, ss, M, X,Ye] = ppca_mv(Ye,d,dia,plo)

threshold = 1e-5;  

if plo 
    set(gcf,'Double','on'); 
end

[N,D] = size(Ye);
    
Obs   = ~isnan(Ye);
hidden = find(~Obs);
missing = length(hidden);

if missing
  for i=1:D
      M(i) = mean(Ye(find(Obs(:,i)),i)); 
  end
else
    M = mean(Ye);
end
Ye = Ye - repmat(M,N,1);

if missing   
    Ye(hidden)=0;
end

r     = randperm(N); 
C     = Ye(r(1:d),:)';     
C     = randn(size(C));
CtC   = C'*C;
X     = Ye * C * inv(CtC);
recon = X*C'; recon(hidden) = 0;
ss    = sum(sum((recon-Ye).^2)) / ( (N*D)-missing);

count = 1; 
old   = Inf;

while count     
   
    Sx = inv( eye(d) + CtC/ss );   
    ss_old = ss;
    if missing 
        proj = X*C'; 
        Ye(hidden) = proj(hidden); 
    end  
    X = Ye*C*Sx/ss;          
    
    SumXtX = X'*X;                              
    C      = (Ye'*X)  / (SumXtX + N*Sx );    
    CtC    = C'*C;
    ss     = ( sum(sum( (C*X'-Ye').^2 )) + N*sum(sum(CtC.*Sx)) + missing*ss_old ) /(N*D); 
    
    objective = N*(D*log(ss) +trace(Sx)-log(det(Sx)) ) +trace(SumXtX) -missing*log(ss_old);           
    rel_ch    = abs( 1 - objective / old );
    old       = objective;
    
    count = count + 1;
    if ( rel_ch < threshold) && (count > 5)
        count = 0;
    end
end  

C = orth(C);
[vecs,vals] = eig(cov(Ye*C));
[vals,ord] = sort(diag(vals));
ord = flipud(ord);
vecs = vecs(:,ord);

C = C*vecs;
X = Ye*C;

Ye = Ye + repmat(M,N,1);
end

function C = textscanu(filename, encoding, del_sym, eol_sym, wb)
h = [];
if nargin < 2
    if ispc
        encoding = 'UTF-16LE';
    else
        encoding = 'UTF-16BE';
    end
end
if ispc
    % Windows defaults
    if strcmp(encoding(end-1:end),'BE')
        byte_order = 'b';
    else
        byte_order = 'l';
    end
else
    % Unix defaults
    if strcmp(encoding(end-1:end),'LE')
        byte_order = 'l';
    else
        byte_order = 'b';
    end
end
if nargin < 3
    del_sym = 9; % column delimitator symbol (TAB=9)
end
if nargin < 4
    if ispc
        eol_sym = 13; % end of line symbol (CR=13, LF=10)
        eol_len = 2; % LF & CR
    else
        eol_sym = 10;
        eol_len = 1; % LF
    end
end
if nargin >= 4
    eol_len = 2;
    t = ispc;
    if (t == 1 && eol_sym ~= 13) || ~t
        eol_len = 1;
    end
end
if nargin > 4
    if strcmp(wb, 'waitbar') == 1;
        h = waitbar(0,''); % display waitbar
        set(h,'name','textscanu')
    end
end
warning off MATLAB:iofun:UnsupportedEncoding;

% read input
fid = fopen(filename, 'r', byte_order, encoding);
S = fscanf(fid, '%c');
fclose(fid);

% remove trailing end-of-line delimitators
while abs(S(end)) == eol_sym || abs(S(end)) == 10
    S = S(1:end-1);
end

% remove Byte Order Marker
BOM = unicode2native(S(1));
if BOM == 26
    S = S(2:end);
end

% add an end of line mark at the end of the file
S = [S char(eol_sym)];

% locates column delimitators and end of lines
del = find(abs(S) == del_sym); 
eol = find(abs(S) == eol_sym);

% get number of rows and columns in input
row = numel(eol);
col = 1 + numel(del) / row;
C = cell(row,col); 

if col - fix(col) ~= 0
    if ishandle(h)
        close(h)
    end
    error(['Error: The file doesn''t have the same number'...
        'of columns per row or row-ends are malformed.'])
end

m = 1;
n = 1;
sos = 1;

if col == 1
    for r = 1:row
        if ishandle(h)
            waitbar( r/row, h, [num2str(r), '/', num2str(row)...
                ' file rows processed'] )
        end
        eos = eol(n) - 1;
        C(r,col) = {S(sos:eos)};
        n = n + 1;
        sos = eos + eol_len + 1;
    end
else
    for r = 1:row
        if ishandle(h)
            waitbar( r/row, h, [num2str(r), '/', num2str(row)...
                ' file rows processed'] )
        end
        for c = 1:col-1
            eos = del(m) - 1;
            C(r,c) = {S(sos:eos)};
            sos = eos + 2;
            m = m + 1;
        end
        % last string in the row
        sos = eos + 2;
        eos = eol(n) - 1;
        C(r,col) = {S(sos:eos)};
        n = n + 1;
        sos = eos + eol_len + 1;
    end
end
if ishandle(h)
    close(h)
end
end
