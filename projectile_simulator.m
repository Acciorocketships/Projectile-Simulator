function varargout = projectile_simulator(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @projectile_simulator_OpeningFcn, ...
                   'gui_OutputFcn',  @projectile_simulator_OutputFcn, ...
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

function projectile_simulator_OpeningFcn(hObject, eventdata, handles, ...
                                                             varargin)
global paused;
paused = false;
handles.output = hObject;
guidata(hObject, handles);

function varargout = projectile_simulator_OutputFcn(hObject, eventdata, ...
                                                                handles) 
varargout{1} = handles.output;

function v0_Callback(hObject, eventdata, handles)

function v0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y0_Callback(hObject, eventdata, handles)

function y0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function angle0_Callback(hObject, eventdata, handles)

function angle0_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dt_Callback(hObject, eventdata, handles)

function dt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xmax_Callback(hObject, eventdata, handles)

function xmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rho_Callback(hObject, eventdata, handles)

function rho_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Cd_Callback(hObject, eventdata, handles)

function Cd_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bouncecoeff_Callback(hObject, eventdata, handles)

function bouncecoeff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mass_Callback(hObject, eventdata, handles)

function mass_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A_Callback(hObject, eventdata, handles)

function A_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxbounce_Callback(hObject, eventdata, handles)

function maxbounce_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function windmag_Callback(hObject, eventdata, handles)

function windmag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function winddirection_Callback(hObject, eventdata, handles)

function winddirection_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function random_Callback(hObject, eventdata, handles)

function func_Callback(hObject, eventdata, handles)

function func_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function autocalculate_Callback(hObject, eventdata, handles)


%% Executes when Simulate is pressed
function sim_Callback(hObject, eventdata, handles)
% Initialize (every time simulate is pressed)
global paused
paused = false;
set(handles.pause,'String','Pause');

% Get values from text boxes
v0 = str2double(get(handles.v0,'String'));
y0 = str2double(get(handles.y0,'String'));
angle0 = str2double(get(handles.angle0,'String'));
xmax = str2double(get(handles.xmax,'String'));
dt = str2double(get(handles.dt,'String'));
bouncecoeff = str2double(get(handles.bouncecoeff,'String'));
mass = str2double(get(handles.mass,'String'));
A = str2double(get(handles.A,'String'));
rho = str2double(get(handles.rho,'String'));
Cd = str2double(get(handles.Cd,'String'));
D = 1/2*rho*A*Cd;
random = get(handles.random,'Value');
maxbounce = str2double(get(handles.maxbounce,'String'));
windmag = str2double(get(handles.windmag,'String'));
winddirection = str2double(get(handles.winddirection,'String'));

% Generate random landscape (or user-specified landscape)
if random
    c = rand(1,11);
    c = c - 0.5;
    c(1) = (c(1)*3+1.5) * 15;
    c(2) = (c(2)+.375) * 0.7;
    c(3) = (c(3)+1.5) * 30;
    c(4) = 2*pi/((c(4)+1.5) * 200);
    c(5) = c(5)*2*pi/c(4);
    c(6) = (c(6)+1.5) * 15;
    c(7) = 2*pi/((c(7)+1.5) * 70);
    c(8) = c(8)*2*pi/c(7);
    c(9) = (c(9)+1.5) * 5;
    c(10) = 2*pi/((c(10)+1.5) * 30);
    c(11) = c(11)*2*pi/c(10);
    groundfunc = @(x) c(1) + x*c(2) + c(3)*sin(c(4)*x-c(5)) + ...
                      c(6)*sin(c(7)*x-c(8)) + c(9)*sin(c(10)*x-c(11));
else
    groundfunc = ['@(x) ' get(handles.func, 'String')];
    groundfunc = str2func(groundfunc);
end

% Call function to calculate trajectories
ground = 1:dt:xmax;
[x,y,vx,vy] = euler_projectile_2d(groundfunc,xmax,dt,y0,v0,angle0,D,...
                        mass,bouncecoeff,maxbounce,winddirection,windmag);

% Animate the trajectory
for i = 1:length(x)
    while paused
        pause(0.1);
    end
    plot(x(1:i),y(1:i),'r--')
    hold on;
    scatter(x(i),y(i),'k');
    plot(ground,groundfunc(ground),'b');
    hold off;
    axis([0 xmax 0 xmax]);
    pause(0.01);
    set(handles.x, 'String', int16(x(i)));
    set(handles.y, 'String', int16(y(i)));
    set(handles.vx, 'String', int16(vx(i)));
    set(handles.vy, 'String', int16(vy(i)));
end
return


%% Executes when "pause" is pressed
function pause_Callback(hObject, eventdata, handles)
% Toggle the value of "paused" and the string displayed
global paused;
if ~paused
    paused = true;
    set(handles.pause,'String','Resume');
else
    paused = false;
    set(handles.pause,'String','Pause');
end


%% Calculates the trajectory
function [x,y,vx,vy] = euler_projectile_2d(groundfunc,xmax,dt,y0,v0,...
                angle0,D,mass,bouncecoeff,maxbounce,winddirection,windmag)
% Preallocation
x = zeros(1,length(0:dt:xmax));
y = zeros(size(x));
vx = zeros(size(x));
vy = zeros(size(y));

% Initial Values
v(1) = v0;
x(1) = 0;
y(1) = y0;
k = 2;
bounce = 0;

% Constants
g = 9.81;

% Initial Calculation
vx(1) = v(1) * cosd(angle0);
vy(1) = v(1) * sind(angle0);
windx = windmag * cosd(winddirection);
windy = windmag * sind(winddirection);

% Calculate Trajectory
while (x(k-1) < xmax) && (x(k-1) >= 0) && (y(k-1) >= 0) && ...
                                          (bounce <= maxbounce)
    % Euler's Fomula
    vx(k) = vx(k-1) + (-D * (vx(k-1)-windx) * ...
                            abs(vx(k-1)-windx)/mass) * dt;
    vy(k) = vy(k-1) + (-D * (vy(k-1)-windy) * ...
                            abs(vy(k-1)-windy)/mass - g) * dt;
    x(k) = x(k-1) + vx(k)*dt;
    y(k) = y(k-1) + vy(k)*dt;
    % In case of bounce
    if y(k) < groundfunc(x(k))
        if k == 2
            y(1) = groundfunc(x(2));
            y(2) = y(1) + vy(k)*dt;
        else
        bounce = bounce + 1;
        y(k) = groundfunc(x(k));
        % Magnitude of the velocity
        mag = (vx(k)^2+vy(k)^2)^0.5;
        % Derivative of ground function
        dgrndfunc = (groundfunc(x(k))-groundfunc(x(k-1)))/(x(k)-x(k-1));
        % Angle of separation is calculated using dot product
        sepangle = acos((vx(k)+vy(k)*dgrndfunc)/((dgrndfunc^2+1)^0.5*mag));
        % Angle of the ground
        grndangle = atan(dgrndfunc);
        % Calculating the deflected angle
        deflectangle = grndangle + sepangle;
        % Calculating components
        vy(k) = bouncecoeff * mag * sin(deflectangle);
        vx(k) = bouncecoeff * mag * cos(deflectangle);
        end
    end
    k = k + 1;
end
% Delete excess preallocated zeros
x(k:end) = [];
y(k:end) = [];
vx(k:end) = [];
vy(k:end) = [];
return
