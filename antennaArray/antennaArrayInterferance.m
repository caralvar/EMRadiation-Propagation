%% Proyecto 1 Aplicaciones de radiacion y propagacion
% Simulaciones de campos de arreglos de antenas conformados por lineas
% infinitas de corriente a lo largo de z.
%
% Prof: Miguel Ruphuy
% Estudiantes:  Carlos Alvarado Salazar B50339
%               Sebasstián Mejía García B54260

%% Opciones del programa
%Lineas de corriente frecuencia y fase:
%La primera de izquierda a derecha corresponde a la que tiene posicion mas
%negativa en y.
phase = 0;
I = [10*exp(-1j*deg2rad(0)) 10*exp(-1j*deg2rad(0)) 10*exp(-1j*deg2rad(0))];
f = 10e9; %Frecuencia de la fuente 
step = 0.0005; %Resolucion del espacio
%fin Opciones.

%% Definicion de constantes y variables.
eo = 8.854187817e-12; %Permitividad del espacio libre
uo = 1.256637061e-6; %Permiabilidad del espacio libre
w = 2*pi*f; %Frecuenca angular
lambda = physconst('LightSpeed')/f; %Longitud de onda
k = 2*pi/lambda; %Numero de onda

%% Establecimiento de las coordenadas y matrices 
%Los vectores X y Y conforman el plano deseado.
[X,Y] = meshgrid(0:step:5.2*lambda,-5.2*lambda:step:5.2*lambda);

%Definicion de los vectores de matrices que almacenaran Ez Hx y Hy para cada
%antena
Hx = zeros(size(X,1)-1,size(X,2)-1,length(I));
Hy = Hx;
Ez = Hx(1:end-1,1:end-1,:);

%Definicion de las matrices que almacenaran Ez Hx y Hy totales

HxTot = Hx(:,:,1);
HyTot = HxTot;
EzTot = Ez(:,:,1);

%% Obtención de los campos para cada línea de corriente.
% Se utiliza el indice 'in' donde el indice menor hace referencia a la
% linea de corriente mas hacia -y 

for in = 1:length(I)
    %Distancia desde la linea hasta el pto de observacion
    rho = sqrt(X.^2+(Y-(in-ceil(length(I)/2))*0.5*lambda).^2); 
    %Funcion hankel de segundo orden de tipo cero evaluada en rho*k
    H2 = besselh(0,2,k*rho);
    %Obtencion de las derivadas necesarias.
    %Primera derivada respecto a x ^ y
    dxH2 = diff(H2,1,2)./step;
    dyH2 = diff(H2,1,1)./step;
    %Segunda derivada respecto a x ^ y
    dx2H2 = diff(dxH2,1,2)./step;
    dy2H2 = diff(dyH2,1,1)./step;
    %Recorte de filas y columnas para tener dimensiones constantes
    dxH2 = dxH2(1:end-1,:);
    dyH2 = dyH2(:,1:end-1);
    dx2H2 = dx2H2(1:end-2,:);
    dy2H2 = dy2H2(:,1:end-2);
    % Campo electrico Ez 
    Ez(:,:,in) = (I(in)/(4*uo*eo*w))*(dx2H2+dy2H2);
    % Campo magnetico Hy 
    Hy(:,:,in) = (I(in)*1j/(4*uo))*(dxH2);
    % Campo magnetico Hx 
    Hx(:,:,in) = (I(in)*-1j/(4*uo))*(dyH2);
end

%% Suma total de los campos

for in = 1:length(I)
    EzTot = EzTot + Ez(:,:,in);
    HyTot = HyTot + Hy(:,:,in);
    HxTot = HxTot + Hx(:,:,in);
end

%% Graficas solicitadas
%Linea a 30 grados para obtener la direccion experimentalmente
xx = 0:0.01:0.15;
yy = xx.*(sqrt(3)/3);
zz = 1.9e11*ones(size(xx));

%Creacion de las figuras
%Campo Ez
figure;
mesh(X(1,1:end-2),Y(1:end-2,1),real(EzTot));
zlim([-2e11 2e11]);
xlim([0 5*lambda]);
ylim([-5*lambda 5*lambda]);
colorbar;
c = colorbar;
c.Label.String = 'E_z (V/m)';
caxis([-2e11 2e11]);
title('Intensidad de campo eléctrico E_z');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
view(0,90);
hold on;
plot3(xx,yy,zz,'k--');

%Campo Hy
zz = 4.96e8*ones(size(xx));
figure;
mesh(X(1,1:end-1),Y(1:end-1,1),real(HyTot));
zlim([-5e8 5e8]);
xlim([0 5*lambda]);
ylim([-5*lambda 5*lambda]);
c = colorbar;
c.Label.String = 'H_y (A/m)';
caxis([-5e8 5e8]);
title('Intensidad de campo magnético H_y');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
view(0,90);
hold on;
plot3(xx,yy,zz,'k--');
