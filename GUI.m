function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 09-Apr-2023 15:58:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function zszf_Callback(hObject, eventdata, handles)
% hObject    handle to zszf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zszf as text
%        str2double(get(hObject,'String')) returns contents of zszf as a double


% --- Executes during object creation, after setting all properties.
function zszf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zszf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zspl_Callback(hObject, eventdata, handles)
% hObject    handle to zspl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zspl as text
%        str2double(get(hObject,'String')) returns contents of zspl as a double


% --- Executes during object creation, after setting all properties.
function zspl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zspl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spzq_Callback(hObject, eventdata, handles)
% hObject    handle to spzq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spzq as text
%        str2double(get(hObject,'String')) returns contents of spzq as a double


% --- Executes during object creation, after setting all properties.
function spzq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spzq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spdk_Callback(hObject, eventdata, handles)
% hObject    handle to spdk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spdk as text
%        str2double(get(hObject,'String')) returns contents of spdk as a double


% --- Executes during object creation, after setting all properties.
function spdk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spdk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function qspl_Callback(hObject, eventdata, handles)
% hObject    handle to qspl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of qspl as text
%        str2double(get(hObject,'String')) returns contents of qspl as a double


% --- Executes during object creation, after setting all properties.
function qspl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to qspl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zsl_Callback(hObject, eventdata, handles)
% hObject    handle to zsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zsl as text
%        str2double(get(hObject,'String')) returns contents of zsl as a double


% --- Executes during object creation, after setting all properties.
function zsl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sjxs_Callback(hObject, eventdata, handles)
% hObject    handle to sjxs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sjxs as text
%        str2double(get(hObject,'String')) returns contents of sjxs as a double


% --- Executes during object creation, after setting all properties.
function sjxs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sjxs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fsxs_Callback(hObject, eventdata, handles)
% hObject    handle to fsxs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fsxs as text
%        str2double(get(hObject,'String')) returns contents of fsxs as a double


% --- Executes during object creation, after setting all properties.
function fsxs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fsxs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function juli_Callback(hObject, eventdata, handles)
% hObject    handle to juli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of juli as text
%        str2double(get(hObject,'String')) returns contents of juli as a double


% --- Executes during object creation, after setting all properties.
function juli_CreateFcn(hObject, eventdata, handles)
% hObject    handle to juli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;

Z=str2num(get(handles.juli,'string'));
n=str2num(get(handles.zsl,'string'));
B=str2num(get(handles.spdk,'string'));
zf=str2num(get(handles.zszf,'string'));
pl=str2num(get(handles.zspl,'string'));
f0=str2num(get(handles.qspl,'string'));
r=str2num(get(handles.fsxs,'string'));
alpha=str2num(get(handles.sjxs,'string'));

c=3e8;
tau=2*Z*n/c;

T=tau;                % 信号持续时间
gamma=B/T;                    % 调频率
ratio=2000;                  % 过采样率
Fs = ratio*B;               % 采样频率
dt = 1/Fs;                  % 采样间隔
N = ceil(T/dt);             % 采样点数
t = (0:N-1)/N*T;      % 时间轴
f = (0:N-1)/N*Fs;     % 频率轴

fb=gamma*tau;E0=5;
%f0是激光源初始频率；gamma是可调谐激光器调谐速度，单位Hz/s；
%fb是拍频大小，E0是光波振幅
dt=t+tau;
Er=E0*exp(1i*(2*pi*f0*t+pi*gamma*t.^2+phi(t,zf,pl)));
dEr=E0*exp(1i*(2*pi*f0*dt+pi*gamma*dt.^2+phi(dt,zf,pl)));
%phi是在t时刻光源随机波动的光相位，Er是本振参考光的光场
%r=0.8;alpha=0.3/1000;
R=r*exp(-alpha*tau*c/n);
Es=sqrt(R)*dEr;         %Es为测试光光场
It=E0^2*(1+R+2*sqrt(R)*cos(2*pi*(f0*tau+fb*t+0.5*gamma*tau^2+phi(t,zf,pl)-phi(dt,zf,pl))));

axes(handles.axes1);
plot(t*1e6,f0*tau+fb*t+0.5*gamma*tau^2+phi(t,zf,pl)-phi(dt,zf,pl)),xlabel('时间 - [\mus]'),ylabel('相位 - [rad]'),title('拍频信号瞬时相位');

%利用希尔伯特变换，将拍频信号I(t)转换到复指数形式
It=2*sqrt(R)*E0^2*exp(1i*2*pi*(f0*tau+fb*t+0.5*gamma*tau.^2+phi(t,zf,pl)-phi(dt,zf,pl)));
%It1=It.*conj(se(t));   
It1=2*sqrt(R)*E0^2*exp(1i*2*pi*(f0*tau+fb*t+0.5*gamma*tau.^2-phi(dt,zf,pl)));

It2=ifft(fft(It1).*exp(1i*pi*f.^2/gamma));
%It2=2*sqrt(R)*E0^2*exp(j*2*pi*(f0*tau+fb*t))*s(t);
s=ifft(fft(conj(se(t,zf,pl))).*exp(1i*pi*f.^2/gamma));                     %信号卷积去斜滤波器
It3=It2.*conj(s);

axes(handles.axes2)
plot(t*1e6,phase(It3),'b');
hold on
plot(t*1e6,f0*tau+fb*t+0.5*gamma*tau^2+phi(t,zf,pl)-phi(dt,zf,pl),'r');
xlabel('时间 - [\mus]'),ylabel('相位 - [rad]'),title('拍频信号I(t)和I_{3}(t)的相位');
legend('I_{3}(t)相位','I(t)相位');

axes(handles.axes3);
plot(t*1e6,2*pi*f0*dt+pi*gamma*dt.^2+phi(dt,zf,pl),'r');
hold on
plot(t*1e6,2*pi*f0*t+pi*gamma*t.^2+phi(t,zf,pl),'b');
xlabel('时间 - [\mus]'),ylabel('MHz'),title('LO和测试光频率曲线');
legend('LO','测试光');

axes(handles.axes4);
Itideal=E0^2*(1+R+2*sqrt(R)*cos(2*pi*(f0*tau+fb*t+0.5*gamma*tau^2)));
spec=fft(Itideal);
spec=fftshift(spec);
plot(t/n*c,abs(spec),'b');
xlabel('距离 - [m]');ylabel('振幅');title('理想条件下的散射峰谱');

axes(handles.axes5);
spec=fft(It);
spec=fftshift(spec);
grid on
plot(t/n*c,abs(spec))
xlabel('距离 - [m]');ylabel('振幅');title('原始拍频信号散射峰谱');

axes(handles.axes6);
spectrum=fft(It3);
spec=spectrum.*exp(1i*pi*f.^2/gamma);
spec=fftshift(spec);
grid on
plot(t/n*c,abs(spec))
xlabel('距离 - [m]');ylabel('振幅');title('去斜滤波器补偿后的散射峰谱');

axes(handles.axes7);
spec=fft(It);
spec=fftshift(spec);
plot(t/n*c,abs(spec),'b');
hold on
spectrum=fft(It3);

spec=spectrum.*exp(1i*pi*f.^2/gamma);
spec=fftshift(spec);
grid on
plot(t/n*c,abs(spec),'r');
xlabel('距离 - [m]');ylabel('振幅');title('去斜滤波器前后散射峰谱对比');
legend('原始频谱','去斜滤波器补偿后频谱');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    delete(allchild(handles.axes1));
    delete(allchild(handles.axes2));
    delete(allchild(handles.axes3));
    delete(allchild(handles.axes4));
    delete(allchild(handles.axes5));
    delete(allchild(handles.axes6));
    delete(allchild(handles.axes7));
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
