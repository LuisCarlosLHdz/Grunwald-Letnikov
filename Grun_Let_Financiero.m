function [ x_Sol,y_Sol,z_Sol] = Grun_Let_Financiero( s,x_Ini,y_Ini,z_Ini,C_ini,h,alpha )
%Grun_Let_Financiero 
%Inputs:
%     s: Numero de soluciones
%     x_Ini: Condicion Inicial en x
%     y_Ini: Condicion Inicial en y
%     z_Ini: Condicion Inicial en z
%     C_Ini: Condicion Inicial en C
%     h: Paso de Integracion
%     alpha: El orden fraccionario de integracion
% Outputs:
%     X_Sol: Solucion en x
%     Y_Sol: Solucion en y
%     Z_Sol: Solucion en z
%     C: Parametro de solucion C
%--------------------------------------------------
%--------------------------------------------------
%Vectores
x_Sol = zeros(1,s+1); %Vector de la señal de solucion en x
y_Sol = zeros(1,s+1); %Vector de la señal de solucion en x
z_Sol = zeros(1,s+1); %Vector de la señal de solucion en x
C = zeros(1,s+1); %Vector del parametro C
S = zeros(1,s+1); %Vector del parametro S
%--------------------------------------------------
%Constantes
a = 1;
b = 0.1;
c = 1;
%--------------------------------------------------
%Condiciones Iniciales
x_Sol(1,1) = x_Ini;
y_Sol(1,1) = y_Ini;
z_Sol(1,1) = z_Ini;
C(1,1) = C_ini;
%--------------------------------------------------
%Solucion de Grünwald-Letnikov
for k=1:1:s
    %--------------------------------------------------
    for j=1:1:k
        %Cálculo para C_{j}
        C(1,j+1) = (1-((1 + alpha)/j))*C(1,j);
        S(1,(k-(j-1)))= C(1,j+1);
    end
    %--------------------------------------------------
    Multx = x_Sol.*S;
    Multy = y_Sol.*S;
    Multz = z_Sol.*S;
    Sumx = 0;
    Sumy = 0;
    Sumz = 0;
    %--------------------------------------------------
    for j=1:1:k
        Sumx = Sumx + Multx(1,j);
        Sumy = Sumy + Multy(1,j);
        Sumz = Sumz + Multz(1,j);
    end
    %--------------------------------------------------
    %Solución para el Sistema Financiero
    [x,y,z] = SisFinanciero(x_Sol(1,k),y_Sol(1,k),z_Sol(1,k),a,b,c);
    x_Sol(1,k+1) = x*(h^(alpha)) - Sumx;
    y_Sol(1,k+1) = y*(h^(alpha)) - Sumy;
    z_Sol(1,k+1) = z*(h^(alpha)) - Sumz;
end
%--------------------------------------------------
end
