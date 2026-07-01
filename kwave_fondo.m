clear; clc; close all;

% Parametros fisicos
f0      = 40e3;
c_air   = 343;
rho_air = 1.21;

% Medio rigido aproximado (impedancia 20x mayor que aire)
Zratio   = 20;
rho_wall = rho_air * Zratio;
c_wall   = c_air;

% Malla
dx = 0.25e-3;
dy = dx;

Nx = 114e-3/dx;
Ny = 56e-3/dy;

kgrid = kWaveGrid(Nx, dx, Ny, dy);

x_vec = (0:Nx-1)*dx;
y_vec = (0:Ny-1)*dy;

% Medio: aire por defecto
medium.sound_speed = c_air * ones(Nx,Ny);
medium.density     = rho_air * ones(Nx,Ny);

% Paredes de alta impedancia
y_bottom  = round(34e-3/dy);
x_bottom  = round(20e-3/dy);
x_bottom2 = round(78e-3/dy);

medium.sound_speed(:,1)                      = c_wall;
medium.density(:,1)                          = rho_wall;
medium.sound_speed(1:x_bottom2,y_bottom:end) = c_wall;
medium.density(1:x_bottom2,y_bottom:end)     = rho_wall;
medium.sound_speed(1:x_bottom,:)             = c_wall;
medium.density(1:x_bottom,:)                 = rho_wall;

% Fuentes: dos lineas de 10 mm enfrentadas a y = 0 mm e y = 34 mm
source.p_mask = zeros(Nx,Ny);

x0_mm = 60;  y1_mm = 0;  y2_mm = 34;  L_mm = 10;

x0   = round((x0_mm*1e-3)/dx);
y1   = round((y1_mm*1e-3)/dy) + 1;
y2   = round((y2_mm*1e-3)/dy);
Lpts = round((L_mm*1e-3)/dx);

xline = (x0-floor(Lpts/2)):(x0+floor(Lpts/2));

source.p_mask(xline, y1) = 1;
source.p_mask(xline, y2) = 1;

% Senal temporal
CFL = 0.2;
kgrid.makeTime(c_air, CFL);

source.p = 10 * sin(2*pi*f0*kgrid.t_array);

% Sensor
sensor.mask = ones(Nx,Ny);

% Simulacion
input_args = {...
    'PMLSize',     20,             ...
    'PMLInside',   false,          ...
    'DisplayMask', source.p_mask,  ...
    'PlotPML',     false,          ...
    'DataCast',    'single'};

sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% Campo RMS normalizado en dB
p_rms = reshape(rms(sensor_data, 2), Nx, Ny);
p_rms = p_rms ./ max(p_rms(:));
PdB   = 20*log10(p_rms + 1e-12);

% Figura
figure;
imagesc(x_vec*1e3, y_vec*1e3, PdB.');
set(gca, 'YDir', 'normal');
xlabel('x [mm]');
ylabel('y [mm]');
title('Campo acustico 2D - k-Wave');
axis image;
colormap(gray);
colorbar;
caxis([-40 0]);
hold on;
surface([0 20; 0 20],   [0 0; 56 56],   [10 10; 10 10], 'FaceColor', 'k');
surface([20 78; 20 78], [34 34; 56 56], [10 10; 10 10], 'FaceColor', 'k');
set(gca, 'XDir', 'reverse');
axis equal;
