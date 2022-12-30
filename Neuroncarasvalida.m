
clc
clear all
close all
%Imagen1
img_1=imread('alex1.jpg');

img_1_bin=imbinarize(img_1(:,:,1));
img_1_canny=edge(img_1(:,:,1),'Canny');
figure(1)
image(img_1)
dimension=300;
img_1_resize=imresize(img_1_canny,[dimension dimension]);
figure(3)
image(img_1_resize);
imshow(img_1_resize);
cara=img_1_resize;
  

[ nf nc] = size(cara);
x(1,:) = cara(1,:);
for k = 2:nf
    x = [ x  cara(k,:) ];
end
[ nxf nxc ] = size(x);
nx = nxf;

% ne = nxc;
% ns = 4;
load pesoscaras1;    % Carga nm v w bias
if(bias == 1)
   ne = ne + 1;
   x = [ x ones(nx,1) ];   
end

in = x';
m = v'*in;
% n = 1.0./(1+exp(-m));    % Sigmoidea 1
n = 2.0./(1+exp(-m)) - 1; % sigmoidea 2
% n = exp(-m.^2);         % Gaussiana
%SEGUNDA CAPA
o=v1'*n;
% n1 = 1.0./(1+exp(-o));    % Sigmoidea 1
n1 = 2.0./(1+exp(-o)) - 1; % sigmoidea 2
% n1 = exp(-o.^2);         % Gaussiana

out = w'*n1;
y = out;
[maxy k] = max(y);
Y_ediv=0;
[y_clase np]=size(y);
%------------------
for a=1:y_clase
    Y_ediv=exp(y(a,1))+Y_ediv;
end
y_porcentaje=[];
total=0;
cuenta=[];
for a =1:y_clase
    y_porcentaje(a,1)=exp(y(a,1))/Y_ediv;
    total=total+y_porcentaje(a,1);
    cuenta(a,1)=a;
end
y_porcentaje=y_porcentaje*100;
for i=1:2 %cantidad de datos comparados

disp("imagen "+i+": "+ y_porcentaje(i,1))
end

valor_maximo=max(y_porcentaje)
%----------------------


if(k == 1)
    disp('La Cara es Roy');
elseif(k == 2)
    disp('La Cara es Alex');
end

