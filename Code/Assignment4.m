%ELEC 4700 
%Assignment 4
%Tariq Aboushaer
%101064544

close all;
clear;
clc;
set(0, 'DefaultFigureWindowStyle', 'docked')

R1 = 1;
Cap = 0.25;
R2 =  2;
L = 0.2;
R3 = 10;
alpha = 100;
R4 = 0.1;
Ro = 1000;
Cn = 0;
iter = 1000;
t = 1; 
delta = t/iter;
Std = 0.05;

[G, C] = GModel(R1, Cap, R2, L, R3, alpha, R4, Ro, Cn);
F = GenF(1,0,1);

%% Question 1

display('Question 1');

display('G matrix');
display(G);
display('C matrix');
display(C);


fig_dc = figure;
hold on;
Vi = [];
Vo = [];
for Vin=-10:0.1:10
    F = GenF(Vin, 0,1);

    e = G\F;
    Vi = [Vi e(1)];
    Vo = [Vo e(5)];
end
title('DC Sweep (TA 101064544)');
xlabel('Vin (V)');
ylabel('Node Voltage (V)');
plot(-10:0.1:10, Vi);
plot(-10:0.1:10, Vo);
legend('Vi', 'VO');


fig_ac = figure;
Vo = [];
F = GenF(1, 0,1); 
for w=1E0:1:1E4
    e = (G+2*pi*w*1j*C)\F;
    Vo = [Vo 20*log10(abs(e(5)/F(8)))];
end
semilogx(1E0:1:1E4, Vo);
hold on;
title('AC Sweep (TA 101064544)');
xlabel('f (Hz)');
ylabel('Gain (dB)');
    

fig_cap = figure;
hold on;
Vo = [];
F = GenF(1, 0, 1); 
Std = 0.05;
CDis = Std.*randn(50000,1) + Cap;
w = pi;
for index=1:50000
    C(1,1) = CDis(index);
    C(1,2) = -CDis(index);
    C(2,1) = -CDis(index);
    C(2,2) = CDis(index);
    e = (G+2*pi*w*1j*C)\F;
    Vo = [Vo 20*log10(abs(e(5)/F(8)))];
end
title('Capacitor Sweep (TA 101064544)');
xlabel('Gain (dB)');
histogram(Vo);

%% Question 2

display('Question 2');

display('G matrix');
display(G);
display('C matrix');
display(C);


Vin = zeros(1,iter);
Vin(0.03*iter:iter) = 1;
F = GenF(Vin, zeros(1,iter), iter);
VList = transient(C, G, F, iter, delta);

Vout = VList(5,:,:);
Vout = Vout(1,:);
figure();
hold on;
plot(linspace(0,t,iter), Vout);
plot(linspace(0,t,iter), Vin);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Step Function Transient Response (TA 101064544)');
legend('Vout', 'Vin');


FF = abs(fftshift(fft(Vout)));
figure();
hold on;
plot(((1:length(FF))/iter)-0.5,20*log10(FF));
FF = abs(fftshift(fft(Vin)));
plot(((1:length(FF))/iter)-0.5,20*log10(FF));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
legend('Vout','Vin');
title('Step Function Frequency Response (TA 101064544)');


Vin = sin(linspace(0,1,iter)*2*pi*1/0.03);
F = GenF(Vin, zeros(1,iter), iter);
VList = transient(C, G, F, iter, delta);

Vout = VList(5,:,:);
Vout = Vout(1,:);
figure();
hold on;
plot(linspace(0,t,iter), Vout);
plot(linspace(0,t,iter), Vin);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Sine Function Transient Response (TA 101064544)');
legend('Vout', 'Vin');


figure();
hold on;
FF = abs(fftshift(fft(Vout)));
plot(((1:length(FF))/iter)-0.5,20*log10(FF));
FF = abs(fftshift(fft(Vin)));
plot(((1:length(FF))/iter)-0.5,20*log10(FF));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
legend('Vout','Vin');
title('Sine Function Frequency Response (TA 101064544)');


Std = 0.03;
Mean = 0.06;
Vin = gaussmf(linspace(0,1,iter),[Std Mean]);
F = GenF(Vin, zeros(1,iter), iter);
VList = transient(C, G, F, iter, delta);

Vout = VList(5,:,:);
Vout = Vout(1,:);
figure();
hold on;
plot(linspace(0,t,iter), Vout);
plot(linspace(0,t,iter), Vin);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Gaussian Function Transient Response (TA 101064544)');
legend('Vout', 'Vin');



figure();
hold on;
FF = abs(fftshift(fft(Vout)));
plot(((1:length(FF))/iter)-0.5,20*log10(FF));
FF = abs(fftshift(fft(Vin)));
plot(((1:length(FF))/iter)-0.5,20*log10(FF));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
legend('Vout','Vin');
title('Gaussian Function Frequency Response (TA 101064544)');

%% Question 3

display('Question 3');

display('G matrix');
display(G);
display('Updated C matrix');
display(C);


Std = 0.03;
Mean = 0.06;
Vin = gaussmf(linspace(0,1,iter),[Std Mean]);
In = 0.001*rand(iter,1);
F = GenF(Vin, In, iter);
VList = transient(C, G, F, iter, delta);

Vout = VList(5,:,:);
Vout = Vout(1,:);
figure();
hold on;
plot(linspace(0,t,iter), Vout);
plot(linspace(0,t,iter), Vin);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Transient Response with Noise std=0.03 (TA 101064544)');
legend('Vout', 'Vin');


figure();
hold on;
FF = abs(fftshift(fft(Vout)));
plot(((1:length(FF))/iter)-0.5,20*log10(FF));
FF = abs(fftshift(fft(Vin)));
plot(((1:length(FF))/iter)-0.5,20*log10(FF));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
legend('Vout','Vin');
title('Gaussian Function Frequency Response with Noise (TA 101064544)');


figure();
hold on;
FF = abs(fftshift(fft(Vout)));
plot(((1:length(FF))/iter)-0.5,20*log10(FF));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
title('Gaussian Function Frequency Response with Various Bandwidth (TA 101064544)');

[G, C] = GModel(R1, Cap, R2, L, R3, alpha, R4, Ro, 0.001);
VList = transient(C, G, F, iter, delta);
Vout = VList(5,:,:);
Vout = Vout(1,:);
FF = abs(fftshift(fft(Vout)));
plot(((1:length(FF))/iter)-0.5,20*log10(FF));

[G, C] = GModel(R1, Cap, R2, L, R3, alpha, R4, Ro, 0.01);
VList = transient(C, G, F, iter, delta);
Vout = VList(5,:,:);
Vout = Vout(1,:);
FF = abs(fftshift(fft(Vout)));
plot(((1:length(FF))/iter)-0.5,20*log10(FF));

[G, C] = GModel(R1, Cap, R2, L, R3, alpha, R4, Ro, 0.0000001);
VList = transient(C, G, F, iter, delta);
Vout = VList(5,:,:);
Vout = Vout(1,:);
FF = abs(fftshift(fft(Vout)));
plot(((1:length(FF))/iter)-0.5,20*log10(FF));

legend('Cn = 0.00001', 'Cn = 0.001', 'Cn = 0.01', 'Cn = 0.0000001');

Cn = 0.00001;
[G, C] = GModel(R1, Cap, R2, L, R3, alpha, R4, Ro, Cn);

Std = 0.03;
Mean = 0.06;

iter = 1000;
delta = t/iter;
Vin = gaussmf(linspace(0,1,iter),[Std Mean]);
In = 0.001*rand(iter,1);
F = GenF(Vin, In, iter);
VList = transient(C, G, F, iter, delta);
Vout = VList(5,:,:);
Vout = Vout(1,:);

figure();
hold on;
plot(linspace(0,t,iter), Vout);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Gaussian Function Transient Response with Various Time Steps (TA 101064544)');

iter = 10000;
delta = t/iter;
Vin = gaussmf(linspace(0,1,iter),[Std Mean]);
In = 0.001*rand(iter,1);
F = GenF(Vin, In, iter);
VList = transient(C, G, F, iter, delta);
Vout = VList(5,:,:);
Vout = Vout(1,:);

plot(linspace(0,t,iter), Vout);

iter = 100;
delta = t/iter;
Vin = gaussmf(linspace(0,1,iter),[Std Mean]);
In = 0.001*rand(iter,1);
F = GenF(Vin, In, iter);
VList = transient(C, G, F, iter, delta);
Vout = VList(5,:,:);
Vout = Vout(1,:);

plot(linspace(0,t,iter), Vout);

legend('1000 Steps', '10000 Steps', '100 Steps');



%% Functions

function [G, C] = GModel(R1, Cap, R2, L, R3, alpha, R4, RO, Cn)
G = [
    1.0000   -1.0000         0         0         0         0         0    1.0000    ;
   -1.0000    1.5000         0         0         0    1.0000         0         0    ;
         0         0    0.1000         0         0   -1.0000         0         0    ;
         0         0         0   10.0000  -10.0000         0    1.0000         0    ;
         0         0         0  -10.0000   10.0010         0         0         0    ;
         0    1.0000   -1.0000         0         0         0         0         0    ;
         0         0  -10.0000    1.0000         0         0         0         0    ;
    1.0000         0         0         0         0         0         0         0    ;
    ];
    
C = [
       Cap      -Cap         0         0         0         0         0         0    ;
      -Cap       Cap         0         0         0         0         0         0    ;
         0         0        Cn         0         0         0         0         0    ;
         0         0         0         0         0         0         0         0    ;
         0         0         0         0         0         0         0         0    ;
         0         0         0         0         0        -L         0         0    ;
         0         0         0         0         0         0         0         0    ;
         0         0         0         0         0         0         0         0    ;
    ];


end

function F = GenF(Vin, In, iterations)
F = zeros(8,1,iterations);
for i=1:iterations
    F(3,1,i) = -In(i);
    F(8,1,i) = Vin(i);
end
end

function VList = transient(C, G, F, iterations, delta)
VList = zeros(8,1, iterations);

for i=2:iterations
    e = C/delta + G;
    VList(:,:,i) = e\(C*VList(:,:,i-1)/delta + F(:,:,i));
end
end

