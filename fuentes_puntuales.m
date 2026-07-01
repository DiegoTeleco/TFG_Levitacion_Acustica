clear all, close all
x=0;%-0.015:0.0002:0.015;
y=-0.0295:0.0002:0.0295;

[X,Y]=meshgrid(x,y);

freq=40000;
c=340;
k=2*pi*freq/c;

x1=0;
y1=-0.03:0.0001:-0.005;

for ii=1:length(y1)

r1=sqrt((X-x1).^2+(Y-y1(ii)).^2);
r2=sqrt((X+x1).^2+(Y+y1(ii)).^2);

p1p(ii,:)=exp(1i*k*r1);%./sqrt(r1);
p2p(ii,:)=exp(1i*k*r2);%./sqrt(r2);
p1(ii,:)=exp(1i*k*r1)./sqrt(r1);
p2(ii,:)=exp(1i*k*r2)./sqrt(r2);

end

figure
subplot(1,2,2)
imagesc(y1,y,real(p1p.'+p2p.'))
colormap gray
hold on
axis equal
xlim([min(y1) max(y1)]);
plot(y1,y1,'--b',y1,-y1,'--b','linewidth',2)
xlabel('$z_s$ (m)','interpreter','latex')
ylabel('$z$ (m)','interpreter','latex')
for i=4
plot(-[i*340/(2*freq) i*340/(2*freq)],[min(y) max(y)],'r','linewidth',2)
end
title('(b)','interpreter','latex')
%colorbar
subplot(1,2,1)
imagesc(y1,y,real(p1.'+p2.'))
colormap gray
axis equal
xlabel('$z_s$ (m)','interpreter','latex')
ylabel('$z$ (m)','interpreter','latex')
setfigpaper([17,.6],10,'arial')
title('(a)','interpreter','latex')
%colorbar
xlim([min(y1) max(y1)]);
clim([-20 20])

setfigpaper([17,0.45],10,'arial')


clear all

for kk=1:2
x=0;%-0.015:0.0002:0.015;
y=-0.0295:0.0002:0.0295;

[X,Y]=meshgrid(x,y);

freq=40000;
c=340;
k=2*pi*freq/c;

x1=0;
if kk==1
y1=(4)*340/(2*40000);
else 
y1=(4.5)*340/(2*40000);
end 

r1=sqrt((X-x1).^2+(Y-y1).^2);
r2=sqrt((X+x1).^2+(Y+y1).^2);

P1p(:)=exp(1i*k*r1);%./sqrt(r1);
P2p(:)=exp(1i*k*r2);%./sqrt(r2);
P1(:)=exp(1i*k*r1)./sqrt(r1);
P2(:)=exp(1i*k*r2)./sqrt(r2);

figure
%plot(y,(abs(real(P1p+P2p))/(max(real(P1p+P2p)))),'--k',y,(abs(real(P1+P2)/max(real(P1+P2)))),'-r','linewidth',2)
plot(y,(abs(real(P1p+P2p))/(max(real(P1p+P2p)))),'--k',y,2*(abs(real(P1+P2)/max(real(P1+P2)))),'-r','linewidth',2)
ylabel('$|p|$','interpreter','latex')
xlabel('$z$ (m)','interpreter','latex')
hold on
plot(y1,0,'sr',-y1, 0,'sr','markersize',10)
ylim([0 0.25])

setfigpaper([17,0.35],10,'arial')
end

clear all

x=-0.035:0.0001:0.035;
y=-0.0495:0.0001:0.0495;

[X,Y]=meshgrid(x,y);

freq=40000;
c=340;
k=2*pi*freq/c;
for kk=1:1
figure

if kk==2
y1=[(-4-1/2):2:(-4-1/2)].*340/(2*freq);
else 
y1=[(-4):2:(-4)].*340/(2*freq);
end
x1=0;
for ii=1:1

r1=sqrt((X-x1).^2+(Y-y1(ii)).^2);
r2=sqrt((X+x1).^2+(Y+y1(ii)).^2);

p1(:,:)=exp(1i*k*r1);%./sqrt(r1);
p2(:,:)=exp(1i*k*r2);%./sqrt(r2);
p=(p1+p2);

subplot(1,3,1)
imagesc(x,y,real(p1+p2))
title('(a)','interpreter','latex')
colormap gray
axis equal
xlabel('$x$ (m)','interpreter','latex')
ylabel('$z$ (m)','interpreter','latex')
ylim([-0.05 0.05])
clim([-20 20])
subplot(1,3,2)
imagesc(x,y,imag(p1+p2))
title('(b)','interpreter','latex')
colormap gray
axis equal
xlabel('$x$ (m)','interpreter','latex')
ylabel('$z$ (m)','interpreter','latex')
ylim([-0.05 0.05])
clim([-20 20])
subplot(1,3,3)
imagesc(x,y,abs(p1+p2))
title('(c)','interpreter','latex')
colormap gray
axis equal
xlabel('$x$ (m)','interpreter','latex')
ylabel('$z$ (m)','interpreter','latex')
ylim([-0.05 0.05])
clim([0 20])
setfigpaper([17,0.35],10,'arial')
end

%% Fuerza de radiación

% Velocidad acústica
dx = x(2)-x(1);
dy = y(2)-y(1);

[dpdx,dpdy] = gradient(p,dx,dy);

omega=2*pi*40000;
rho0=1.21;
c=340;
R=0.001;

vx = -(1/(1j*omega*rho0))*dpdx;
vy = -(1/(1j*omega*rho0))*dpdy;

% POTENCIAL DE GOR'KOV
U = (4*pi*R^3)*(abs(p).^2/(2*rho0*c^2) ...
    - (3*rho0/4)*( abs(vx).^2 + abs(vy).^2 ));

% FUERZA DE RADIACIÓN
% F = -grad(U)

[dUdx,dUdy] = gradient(U,dx,dy);

Fx = -dUdx;
Fy = -dUdy;


%% =======================================================
% FIGURA 3: fuerza de radiación
% ========================================================
figure
subplot(1,3,1)
imagesc(x,y,U);
hold on
ylim([-0.005 0.005])
title('(a)','interpreter','latex')
set(gca,'YDir','normal');
xlabel('$x$ (m)','interpreter','latex');
ylabel('$z$ (m)','interpreter','latex');
colormap gray
subplot(1,3,2)
imagesc(x,y,abs(sqrt(Fx.^2+Fy.^2)))
hold on
ylim([-0.005 0.005])
title('(b)','interpreter','latex')
set(gca,'YDir','normal');
xlabel('$x$ (m)','interpreter','latex');
ylabel('$z$ (m)','interpreter','latex');
colormap gray
%colormap hot


subplot(1,3,3)
imagesc(x,y,real(p));
title('(c)','interpreter','latex')
set(gca,'YDir','normal');
xlabel('$x$ (m)','interpreter','latex');
ylabel('$z$ (m)','interpreter','latex');
colormap gray
ylim([-0.005 0.005])
xlim([-0.01 0.01])
hold on;


% Máscaras
skip=10;

X1=X(1:skip:end,(end-1)/2-3);
Y1=Y(1:skip:end,(end-1)/2-3);
Fx1=Fx(1:skip:end,(end-1)/2-3);
Fy1=Fy(1:skip:end,(end-1)/2-3);


idx_pos = Fy1 > 0;
%idx_neg = Fy1 < 0;

% Flechas Fy > 0  (rojo)

quiver( ...
    X1(idx_pos),...
    Y1(idx_pos),...
    Fx1(idx_pos),...
    Fy1(idx_pos),...
    0.25,...
    'r');

% Flechas Fy < 0 (azul)

X1=X(1:skip:end,(end-1)/2+3);
Y1=Y(1:skip:end,(end-1)/2+3);
Fx1=Fx(1:skip:end,(end-1)/2+3);
Fy1=Fy(1:skip:end,(end-1)/2+3);

idx_neg = Fy1 < 0;

quiver( ...
    X1(idx_neg),...
    Y1(idx_neg),...
    Fx1(idx_neg),...
    Fy1(idx_neg),...
    0.25,...
    'b');



xlabel('$x$ (mm)','interpreter','latex');
ylabel('$y$ (mm)','interpreter','latex');
setfigpaper([17,0.35],10,'arial')

end

%%

Pot=U(1:end,(end-1)/2);
FF=Fy(1:end,(end-1)/2);
PP=abs(p(1:end,(end-1)/2));
figure
plot(y,Pot/max(Pot(400:end-400)),y,FF/max(FF(400:end-400)),y,PP/max(PP(400:end-400)),'linewidth',2)
xlabel('$z$ (mm)','interpreter','latex');
ylabel('$U$, $F_z$, $|p|$ normalizadas','interpreter','latex');
setfigpaper([17,0.35],10,'arial')
xlim([-0.017 0.017])
%ylim([-1 1])

