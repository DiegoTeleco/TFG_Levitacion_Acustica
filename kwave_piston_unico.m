clearvars; close all; clc;

% Parametros fisicos
c0     = 340;
rho0   = 1.2;
f0     = 40e3;
lambda = c0 / f0;
a      = 5e-3;

% Malla
ppw = 12;
dz  = lambda / ppw;
dr  = dz;
Lz  = 0.08;
Lr  = 0.06;
Nz  = round(Lz / dz);
Nr  = round(Lr / dr);

ir_max    = round(a / dr);
ir_centro = round(Nr/2);

kgrid = kWaveGrid(Nr, dr, Nz, dz);

% Medio
medium.sound_speed = c0;
medium.density     = rho0;
medium.alpha_coeff = 0;
medium.alpha_power = 0;

% Tiempo: 50 ciclos en total, se graban los ultimos 10
kgrid.makeTime(medium.sound_speed);
dt = kgrid.dt;

T0            = 1 / f0;
nCyclesTotal  = 50;
nCyclesRecord = 10;
Nt = round(nCyclesTotal * T0 / dt);
kgrid.setTime(Nt, dt);
record_start_index = Nt - round(nCyclesRecord * T0 / dt) + 1;

% Fuente: piston circular centrado en r = 0, z = 0
source.p_mask = zeros(Nr, Nz);
source.p_mask(ir_centro-ir_max+1 : ir_centro+ir_max, 1) = 1;

source_mag = 10;

% Sensor
sensor.mask               = ones(Nr, Nz);
sensor.record             = {'p_rms', 'p_final'};
sensor.record_start_index = record_start_index;

% Argumentos de simulacion
input_args = { ...
    'PMLInside', false, ...
    'PMLSize',   20,    ...
    'PlotSim',   false, ...
    'DataCast',  'single'};

% Simulacion 1: seno -> parte imaginaria del campo complejo
source.p      = source_mag * sin(2*pi*f0 * kgrid.t_array);
source.p      = filterTimeSeries(kgrid, medium, source.p);
source.p_mode = 'additive';

sensor_sin = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% Simulacion 2: coseno -> parte real del campo complejo
source.p      = source_mag * cos(2*pi*f0 * kgrid.t_array);
source.p      = filterTimeSeries(kgrid, medium, source.p);
source.p_mode = 'additive';

sensor_cos = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% Postprocesado
PML_size = 20;
iz_max   = Nz - PML_size - 5;

p_real_2D = reshape(sensor_cos.p_final, Nr, Nz);
p_imag_2D = reshape(sensor_sin.p_final, Nr, Nz);
p_abs_2D  = sqrt(2) * reshape(sensor_sin.p_rms, Nr, Nz);

% Eliminar zona de la fuente y recortar PML
p_real_2D(source.p_mask == 1) = 0;
p_imag_2D(source.p_mask == 1) = 0;
p_abs_2D(source.p_mask  == 1) = 0;

p_real_2D = p_real_2D(:, 1:iz_max);
p_imag_2D = p_imag_2D(:, 1:iz_max);
p_abs_2D  = p_abs_2D(:,  1:iz_max);

% Normalizar
P_real_n = p_real_2D / max(abs(p_real_2D(:)));
P_imag_n = p_imag_2D / max(abs(p_imag_2D(:)));
P_abs_n  = p_abs_2D  / max(p_abs_2D(:));

% Ejes
z_vec = (0:iz_max-1) * dz;
r_vec = kgrid.x_vec;
z_lim = [0      z_vec(end)];
r_lim = [-15e-3 15e-3];

% Figura: tres paneles (a) real, (b) imaginario, (c) absoluto
fig = figure('Units', 'centimeters', 'Position', [2 2 18 7]);

titles = {'(a)', '(b)', '(c)'};
datos  = {P_real_n, P_imag_n, P_abs_n};
clims  = {[-1 1], [-1 1], [0 0.6]};

for ii = 1:3
    subplot(1, 3, ii);
    imagesc(r_vec, z_vec, datos{ii}');
    set(gca, 'YDir', 'normal');
    colormap(gray);
    clim(clims{ii});
    xlabel('$r$ (m)', 'Interpreter', 'latex', 'FontSize', 10);
    if ii == 1
        ylabel('$z$ (m)', 'Interpreter', 'latex', 'FontSize', 10);
    else
        ylabel('');
        set(gca, 'YTickLabel', []);
    end
    xlim(r_lim);
    ylim(z_lim);
    title(titles{ii}, 'Interpreter', 'latex', 'FontSize', 10, ...
        'FontWeight', 'normal');
    axis square;
end

set(fig, 'Color', 'w');
exportgraphics(fig, 'kwave_piston_unico_2D.pdf', 'ContentType', 'vector');
