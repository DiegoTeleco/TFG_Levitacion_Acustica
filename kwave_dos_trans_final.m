clearvars; close all; clc;

% Parametros fisicos
c0     = 340;
rho0   = 1.2;
f0     = 40e3;
lambda = c0 / f0;
k      = 2*pi / lambda;
a      = 5e-3;
zs     = 4 * c0 / (2*f0);

fprintf('lambda = %.2f mm\n', lambda*1e3);
fprintf('a      = %.1f mm\n', a*1e3);
fprintf('2*zs   = %.2f mm\n', 2*zs*1e3);

% Malla
ppw = 40;
dz  = lambda / ppw;
dr  = dz;

margen = 0.01;
Lz     = 2*zs + 2*margen;
Lr     = 0.04;

Nz = round(Lz / dz);
Nr = round(Lr / dr);

ir_max    = round(a / dr);
ir_centro = round(Nr / 2);

kgrid = kWaveGrid(Nr, dr, Nz, dz);

[~, iz_left]  = min(abs(kgrid.y_vec - (-zs)));
[~, iz_right] = min(abs(kgrid.y_vec - ( zs)));

% Medio
medium.sound_speed = c0;
medium.density     = rho0;
medium.alpha_coeff = 0;
medium.alpha_power = 0;

% Tiempo: 60 ciclos en total, se graban los ultimos 10
kgrid.makeTime(medium.sound_speed);
dt = kgrid.dt;

T0            = 1 / f0;
nCyclesTotal  = 60;
nCyclesRecord = 10;
Nt = round(nCyclesTotal * T0 / dt);
kgrid.setTime(Nt, dt);
record_start_index = Nt - round(nCyclesRecord * T0 / dt) + 1;

% Fuentes: dos pistones circulares enfrentados centrados en r = 0
source.p_mask = zeros(Nr, Nz);
source.p_mask(ir_centro-ir_max+1 : ir_centro+ir_max, iz_left)  = 1;
source.p_mask(ir_centro-ir_max+1 : ir_centro+ir_max, iz_right) = 1;

source_mag    = 10;
source.p      = source_mag * sin(2*pi*f0 * kgrid.t_array);
source.p      = filterTimeSeries(kgrid, medium, source.p);
source.p_mode = 'additive';

% Sensor
sensor.mask               = ones(Nr, Nz);
sensor.record             = {'p_rms', 'p'};
sensor.record_start_index = record_start_index;

% Simulacion
input_args = { ...
    'PMLInside', false, ...
    'PMLSize',   20,    ...
    'PlotSim',   false, ...
    'DataCast',  'single'};

sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% Postprocesado: amplitud normalizada
p_rms  = reshape(sensor_data.p_rms, Nr, Nz);
p_amp  = sqrt(2) * p_rms;
p_amp(source.p_mask == 1) = 0;
p_norm = p_amp / max(p_amp(:));

% Campo complejo por proyeccion fasorial sobre los ciclos grabados
omega  = 2*pi*f0;
Nt_rec = size(sensor_data.p, 2);
t_abs  = ((record_start_index-1):(record_start_index-1+Nt_rec-1)) * dt;

cos_ref = cos(omega * t_abs);
sin_ref = sin(omega * t_abs);

P_re = reshape((2/Nt_rec) * (sensor_data.p * cos_ref'), Nr, Nz);
P_im = reshape((2/Nt_rec) * (sensor_data.p * sin_ref'), Nr, Nz);
p_complex = P_re - 1i*P_im;

% Suavizado gaussiano para la figura 2D (sin Image Processing Toolbox)
sigma_px = 2.5;
hw       = ceil(3 * sigma_px);
gx       = -hw:hw;
gauss1d  = exp(-gx.^2 / (2*sigma_px^2));
gauss1d  = gauss1d / sum(gauss1d);
gauss2d  = gauss1d' * gauss1d;
p_smooth = conv2(p_norm, gauss2d, 'same');

% Ejes
z_axis = kgrid.y_vec * 1e3;
r_axis = kgrid.x_vec * 1e3;

z_left_mm  = kgrid.y_vec(iz_left)  * 1e3;
z_right_mm = kgrid.y_vec(iz_right) * 1e3;

% Nodos teoricos de presion: z_n = (2n+1)*lambda/4
n_nodos = -10:10;
z_nodos = (2*n_nodos + 1) * lambda/4 * 1e3;
z_nodos = z_nodos(abs(z_nodos) < zs*1e3);

p_axis = p_norm(ir_centro, :);

% =========================================================================
% Figura 1: Campo de presiones 2D
% =========================================================================
fig1 = figure('Units', 'centimeters', 'Position', [2 2 18 8]);

imagesc(z_axis, r_axis, p_smooth);
set(gca, 'YDir', 'normal');
colormap(gray);
cb = colorbar;
cb.Label.String      = 'Amplitud normalizada';
cb.Label.Interpreter = 'latex';
xlabel('$z$ (mm)', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$r$ (mm)', 'Interpreter', 'latex', 'FontSize', 11);
clim([0 1]);
ylim([-15 15]);
xlim([z_left_mm - 5, z_right_mm + 5]);

hold on;
plot([z_left_mm  z_left_mm],  [-a*1e3  a*1e3], 'w-', 'LineWidth', 3);
plot([z_right_mm z_right_mm], [-a*1e3  a*1e3], 'w-', 'LineWidth', 3);
hold off;

exportgraphics(fig1, 'kwave_dos_trans_presiones_2D.pdf', 'ContentType', 'vector');

% =========================================================================
% Figura 2: Corte axial en r = 0
% =========================================================================
fig2 = figure('Units', 'centimeters', 'Position', [2 2 16 7]);

plot(z_axis, p_axis, 'k-', 'LineWidth', 1.4);
hold on;
for ii = 1:length(z_nodos)
    hn = xline(z_nodos(ii), 'r:', 'LineWidth', 0.9);
    hn.HandleVisibility = 'off';
end
xline(z_left_mm,  'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(z_right_mm, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
hold off;

xlabel('$z$ (mm)', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Amplitud normalizada', 'Interpreter', 'latex', 'FontSize', 11);
xlim([z_left_mm - 5, z_right_mm + 5]);
ylim([0 1.15]);
grid on; box on;

exportgraphics(fig2, 'kwave_dos_trans_presiones_axial.pdf', 'ContentType', 'vector');

% =========================================================================
% Figura 3: Potencial de Gor'kov y fuerza de radiacion axial
% =========================================================================
p_amp_clean = abs(p_complex);

% Enmascarar campo cercano de los transductores
mask_iz = false(1, Nz);
mask_iz(max(1,iz_left-3):min(Nz,iz_left+3))   = true;
mask_iz(max(1,iz_right-3):min(Nz,iz_right+3)) = true;
p_amp_clean(:, mask_iz) = NaN;

% Potencial de Gor'kov: U ~ <p^2>
U_field = p_amp_clean.^2;

% Fuerza axial: Fz = -dU/dz por diferencias finitas centradas
U_right  = U_field(:, 3:end);
U_left   = U_field(:, 1:end-2);
dUdz_full = zeros(Nr, Nz);
dUdz_full(:, 2:end-1) = (U_right - U_left) / (2*dz);
Fz = -dUdz_full;

Fz_axis = Fz(ir_centro, :);
U_axis  = U_field(ir_centro, :);
p_norm_ax = p_axis(:)';

U_max  = max(U_axis, [], 'omitnan');
U_norm = U_axis / U_max;
Fz_norm = Fz_axis / max(abs(Fz_axis));

z_row = reshape(z_axis, 1, []);

fig3 = figure('Units', 'centimeters', 'Position', [2 2 18 8]);

col_U  = [0.2157  0.4941  0.7216];
col_Fz = [0.8353  0.3686  0.0000];
col_p  = [0.9294  0.6941  0.1255];

hold on;
hU  = plot(z_row, U_norm,    '-', 'Color', col_U,  'LineWidth', 1.6);
hFz = plot(z_row, Fz_norm,   '-', 'Color', col_Fz, 'LineWidth', 1.6);
hp  = plot(z_row, p_norm_ax, '-', 'Color', col_p,  'LineWidth', 1.6);

yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
for ii = 1:length(z_nodos)
    hn = xline(z_nodos(ii), 'r:', 'LineWidth', 0.9);
    hn.HandleVisibility = 'off';
end
xline(z_left_mm,  'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(z_right_mm, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
hold off;

xlabel('$z$ (mm)', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$U,\; F_z,\; |p|$ normalizadas', 'Interpreter', 'latex', 'FontSize', 11);
xlim([z_left_mm+1, z_right_mm-1]);
ylim([-1.5 1.5]);
grid on; box on;

legend([hU, hFz, hp], ...
    {'$U_\mathrm{norm}$', '$F_{z\,\mathrm{norm}}$', '$|p|_\mathrm{norm}$'}, ...
    'Interpreter', 'latex', 'FontSize', 10, 'Location', 'southeast');

exportgraphics(fig3, 'kwave_dos_trans_fuerzas.pdf', 'ContentType', 'vector');
