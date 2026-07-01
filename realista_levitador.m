% ============================================================
% Simulación 2D con k-Wave
% Dos fuentes lineales de 40 kHz en aire
%
% Fuentes:
%   (x,y) = (55,0) mm
%   (x,y) = (55,34) mm
%
% Longitud de cada fuente: 10 mm (eje x)
%
% Región superior e inferior:
% impedancia acústica 1000x mayor que aire
%
% Requiere:
%   https://www.k-wave.org
% ============================================================

clear;
clc;
close all;

%% ============================================================
% PARÁMETROS FÍSICOS
% =============================================================

f0 = 40e3;                 % Frecuencia [Hz]

c_air   = 343;             % Velocidad sonido aire [m/s]
rho_air = 1.21;            % Densidad aire [kg/m3]

% Medio rígido aproximado
Zratio = 20;

rho_wall = rho_air * Zratio;
c_wall   = c_air;


%% ============================================================
% MALLA
% =============================================================

dx = 0.25e-3;              % Paso espacial [m]
dy = dx;

Nx = 114e-3/dx;
Ny = 56e-3/dy;

kgrid = kWaveGrid(Nx, dx, Ny, dy);

x_vec = (0:Nx-1)*dx;
y_vec = (0:Ny-1)*dy;

%% ============================================================
% MAPAS DEL MEDIO
% =============================================================

medium.sound_speed = c_air * ones(Nx,Ny);
medium.density     = rho_air * ones(Nx,Ny);

%% ============================================================
% PAREDES DE ALTA IMPEDANCIA
% =============================================================

% Convertimos posiciones físicas a índices

y_bottom = round(34e-3/dy);
x_bottom = round(20e-3/dy);
x_bottom2 = round(78e-3/dy);

% Región rígida inferior
medium.sound_speed(:,1) = c_wall;
medium.density(:,1)     = rho_wall;

% Región rígida superior
medium.sound_speed(1:x_bottom2,y_bottom:end) = c_wall;
medium.density(1:x_bottom2,y_bottom:end)     = rho_wall;

medium.sound_speed(1:x_bottom,:) = c_wall;
medium.density(1:x_bottom,:)     = rho_wall;

% figure
% surface(medium.sound_speed,'linestyle','none')
% figure
% surface(medium.density,'linestyle','none')

%% ============================================================
% FUENTES
% =============================================================

source.p_mask = zeros(Nx,Ny);

% Posiciones físicas
x0_mm = 55;

y1_mm = 0;
y2_mm = 34;

L_mm = 10;

% Conversión a índices

x0 = round((x0_mm*1e-3)/dx);

y1 = round((y1_mm*1e-3)/dy)+1;
y2 = round((y2_mm*1e-3)/dy);

Lpts = round((L_mm*1e-3)/dx);

xline = (x0-floor(Lpts/2)):(x0+floor(Lpts/2));

% Fuente inferior
source.p_mask(xline,y1) = 1;

% Fuente superior
source.p_mask(xline,y2) = 1;

%% ============================================================
% SEÑAL TEMPORAL
% =============================================================

CFL = 0.2;

kgrid.makeTime(c_air, CFL);

source_mag = 10;

source.p = source_mag * ...
    sin(2*pi*f0*kgrid.t_array);

%% ============================================================
% SENSOR
% =============================================================

sensor.mask = ones(Nx,Ny);

%% ============================================================
% SIMULACIÓN
% =============================================================

input_args = {...
    'PMLSize', 20,...
    'PMLInside', false,...
    'DisplayMask', source.p_mask,...
    'PlotPML', false,...
    'DataCast', 'single'};

sensor_data = kspaceFirstOrder2D(...
    kgrid,...
    medium,...
    source,...
    sensor,...
    input_args{:});

%% ============================================================
% CAMPO RMS
% =============================================================

p_rms = reshape(rms(sensor_data,2),Nx,Ny);

p_rms = p_rms ./ max(p_rms(:));

PdB = 20*log10(p_rms + 1e-12);

%% ============================================================
%% REPRESENTACIÓN
% =============================================================

figure;

imagesc(x_vec*1e3, y_vec*1e3, PdB.');

set(gca,'YDir','normal');

xlabel('x [mm]');
ylabel('y [mm]');

title('Campo acústico 2D - k-Wave');

axis image;

colormap(gray);

colorbar;

caxis([-40 0]);

hold on;

% % Dibujar transductores
% plot(y1_mm*ones(size(xline)), xline*dx*1e3,...
%     'w','LineWidth',3);
% 
% plot(y2_mm*ones(size(xline)), xline*dx*1e3,...
%     'w','LineWidth',3);
% 
% hold on

surface([0 20; 0 20],[0 0; 56 56],[10 10; 10 10],'FaceColor','k')
surface([20 78; 20 78],[34 34; 56 56],[10 10; 10 10],'FaceColor','k')
set(gca,'XDir','reverse');
axis equal

%% ============================================================
% INFORMACIÓN
% =============================================================

disp('Simulación completada');

lambda = c_air/f0;

fprintf('Longitud de onda = %.2f mm\n', lambda*1e3);
fprintf('Puntos por longitud de onda = %.2f\n', lambda/dx);