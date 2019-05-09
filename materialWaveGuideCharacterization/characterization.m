%% Proyecto 2 Apl. Radiación y propagación
% Carlos Alvarado Salazar B50339
% Prof. Miguel Ruphuy 

% Determinacion de la permitividad y permiabilidad de un material en un
% rango de frecuencia especifico utilizando los parametros s11 y s21
% obtenidos en un analizador de dos puertos. El metodo utilizado esta
% basado en las relaciones Kramers-Kronig desarrollado por Vasundara VV y
% Ruyen Ro. 

%% Determinacion de constantes

e0 = 8.854187817e-12; %Permitividad del espacio libre
mu0 = 1.256637061e-6; %Permiabilidad del espacio libre
eta0 = 376.73031340781; % Impedancia del espacio libre
d0 = 10e-3; % Grosor del material
d1 = 37.5e-3; % Distancia del puerto al material

%% Informacion obtenida de los parametros s

% El codigo requiere que dos 1rchivos .mat contengan los resultados de los
% parametros s. Se deben llamar sp11.mat y sp21.mat En la primer columna
% deben contener la frecuencia, en la segunda la parte real del parametro y
% en la tercera la parte imaginaria. Deben tener el mismo numero de
% columnas.
load("Material1S11.mat");
PS11 = Mat1S11;
load("Material1S21.mat");
PS21 = Mat1s21;

% El dato de la frecuencia puede estar en GHz, MHz o Hz entonces se multiplica
% por 10^9, 10^6 o 1 respectivamente.
Hz = 1;
GHz = 1e9;
MHz = 1e6;
PS11(:,1) = PS11(:,1).*GHz;
PS21(:,1) = PS21(:,1).*GHz;


f = PS11(:,1); %Vector con valores de la frecuencia
w = 2*pi*f;  %Frecuencia angular
%Constante de propagacion del espacio libre.
beta0 = w/physconst('LightSpeed');
%Angulo electrico de desfase provocado por la distancia de los puertos al
%material
theta = beta0*d1;

%Parametros s complejos:
s11 = (PS11(:,2)+1j*PS11(:,3)).*exp(2j*theta);
s21 = (PS21(:,2)+1j*PS21(:,3)).*exp(2j*theta);


%% Calculos

K = ((s11.^2 - s21.^2)+1)./(2.*s11); %factor K

Refl = zeros(size(K)); % Coeficiente de reflexion
for i = 1:size(Refl)    % La seleccion de si suma o resta la raiz depende
    Refl(i) = K(i) + sqrt(K(i)^2-1); % de si la magnitud de Refl es menor que 1
    if abs(Refl(i)) > 1
       Refl(i) = K(i) - sqrt(K(i)^2 -1);
    end
end

% Coeficiente de transmision
Trans =  ((s11 + s21) - Refl)./(1 - (s11 + s21).*Refl);

% Impedancia del material
eta = eta0*((1 + Refl)./(1 - Refl));

% Constante de atenuacion
alpha = (-1*log(abs(Trans)))/d0;

%% Calculo de la integral para determinar el desfase de beta (m)
Integral = zeros(length(PS11)); %Vector que guarda el resultado de la integral
betakk = zeros(size(Integral(:,1))); %Vector que guarda el resultado de beta

for i = 1:length(Integral) %Se prepara la funcion de la integral
    Integral(:,i) = (w.*alpha./beta0)./(w.^2-w(i).^2);
    %Como se presentan discontinuidades infinitas se promedian utilizando
    %el promedio movil (moving average)
    Integral(:,i) = medfilt1(Integral(:,i),6);
    %Se calcula la integral utilizando una aproximacion trapezoidal y se
    %calcula el beta a partir de la relacion de kramers-kronig
    betakk(i) = beta0(i)*(1+(2/pi)*(trapz(Integral(:,i),w)));
end

% Primero se calcula el beta para m=1, m=-1 y m=0
beta = zeros(length(Trans(:,1)),3);
for i = -1:1
    beta(:,i+2) = (-angle(Trans)+2*i*pi)/d0;
end
   % beta_final = (-angle(Trans)+2*0*pi)/d0;

%Luego se obtiene el beta final segun el beta obtenido 
%a partir de la relacion de Kramers-Kronig
%el valor de m del beta cambia de 0 a 1 en 4,785 GHz. 
beta_final = zeros(size(beta(:,1)));
for i = 1:length(beta_final(:,1))
    beta_final(i) = beta(i,2);
%     if f(i)>1.719e9
%         beta_final(i) = beta(i,1);
%      end
%     
%     if f(i)>2.7585e9
%         beta_final(i) = beta(i,3);
%      end
end

%Finalmente se calcula la permitividad y permeabilidad
%Permitividad relativa
e_r = (alpha+1j*beta_final)./(1j.*w.*eta*e0);
%Permeabilidad relativa
mu_r = (eta.*(alpha+1j*beta_final))./(1j.*w*mu0);


%% Figuras
figure;
plot(f,beta(:,1),'--m');
hold on;
plot(f,beta(:,2),'b','LineWidth',1);
plot(f,beta(:,3),'--r');
plot(f,betakk,'-.k','LineWidth',1);
grid;
title('Coeficiente de propagación');
xlabel('Frecuencia (Hz)');
ylabel('\beta','FontSize',14);
%xlim([3e9 5e9]);
legend('m = -1','m = 0','m = 1','\beta_{K-K}');

figure;
plot(f,real(e_r),'b','LineWidth',1.2);
hold on;
plot(f,imag(e_r),'--r','LineWidth',1.2);
grid;
title('Permitividad relativa del material');
xlabel('Frecuencia (Hz)');
ylabel('\epsilon_r','FontSize',14);
%xlim([3e9 5e9]);
%ylim([-3 14]);
legend('Real','Imaginaria','Location','northwest');

figure;
plot(f,real(mu_r),'b','LineWidth',1.2);
hold on;
plot(f,imag(mu_r),'--r','LineWidth',1.2);
grid;
title('Permeabilidad relativa del material');
xlabel('Frecuencia (Hz)');
ylabel('\mu_r','FontSize',14);
%xlim([3e9 5e9]);
%ylim([-1 2]);
legend('Real','Imaginaria','Location','northwest');








