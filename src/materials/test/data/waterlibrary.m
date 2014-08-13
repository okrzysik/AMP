% Test of AMP WaterLibrary.cc
% This script is used to generate data for the "extra" tests in
% testWaterLibrary.cc.  Specifically, the goal is to test the water
% properties in extended regions of the correlations (i.e. supercritical,
% two-phase, and vapor regions)
% This script internally works in British units (these are the units
% used in the correlations), but the output is converted to SI units
% so the values can be used and tested directly without conversion.

%% Correlation data

% This data is taken from:
%  M. P. Paulsen, et. al, \
%  RETRAN-3D, \
%  Electric Power Research Institute, Technical Document NP-7450, Volume 1, \
%  September 1998

cf1 = [6.97088786e+01, 3.33752999e+01, 2.31824073e+00, 1.84059951e-01, -5.24550228e-03, 2.87800703e-03, 1.75365232e-03, -4.33485962e-04, 3.32569928e-05];
cf2 = [8.40861880e+05, 3.63741321e+05, -4.63450667e+05, 1.13030634e+05, -4.35021730e+02, -3.89898819e+03, 6.69739943e+02, -4.73072638e+01, 1.26512506e+00];
cf3 = [9.06003044e+02, -1.42681352e+01, 1.52223326e+00, -6.97399296e-01, 1.74309166e-01, -2.31971770e-02, 1.69401915e-03, -6.45477171e-05, 1.00300310e-06];

cg1 = [1.10583687e+03, 1.43694377e+01, 8.01828862e-01, 1.61723291e-02, -1.50114750e-03, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, -1.23767556e-05, 3.00477330e-06, -2.06239073e-07];
cg2 = [-2.23426500e+06, 1.23124763e+06, -1.97884787e+05, 1.85998804e+01, -2.76570132e+00, 1.03603388e+03, -2.14342313e+02, 1.69050776e+01, -4.86432213e-01];
cg3 = [9.05997825e+02, 5.56195754e+00, 3.43418961e+00, -6.40639063e-01, 5.91857948e-02, -2.72537857e-03, 5.00633694e-05];

ct1 = [3.27627555e+01, 9.76361700e-01, 1.85722603e-04, -4.68267433e-07;
       3.36088021e-03, -5.59528176e-05, 1.61859599e-07, -1.18020438e-10];
ct2 = [6.39080121e+02, -3.05521724e+00, 8.71323187e-03, -6.26940368e-06, -9.84470000e-18;
       -4.30285724e-01, 2.67330344e-03, -5.19838047e-06, 3.03782556e-09, 3.30970404e-13;
       1.17452458e-04, -6.83920099e-07, 1.16801177e-09, -4.26007418e-13, -2.73208716e-16;
       -1.47397729e-08, 8.01885817e-11, -1.16490136e-13, 4.55940029e-18, 5.82551114e-20;
       7.10432734e-13, -3.64950063e-15, 4.45738757e-18, 1.67839872e-21, -3.75680409e-24];
ct3 = [-1.17910086e+04, 2.82927434e+01, -2.67818156e-02, 1.21874275e-05, -2.09203315e-09;
       1.25616091e+02, -3.33344885e-01, 3.32690127e-04, -1.47789033e-07, 2.46325837e-11;
       -1.08371337e-01, 2.92817773e-04, -2.97243646e-07, 1.34263911e-10, -2.27558572e-14;
       3.27807185e-05, -8.97095936e-08, 9.24624831e-11, -4.24915552e-14, 7.33831675e-18;
       -3.42556493e-09, 9.52769245e-12, -1.00140904e-14, 4.70391440e-18, -8.31504474e-22];
ct4 = [3.79525685e+03, -6.34703101e+00, 2.86722833e-03, 5.95359981e-09, 4.79820744e-11;
       -3.91008624e+00, 1.22274782e-02, -1.40466470e-05, 7.50567946e-09, -1.60869365e-12;
       3.41050016e-05, 7.01090011e-10, -1.03020187e-10, 5.73109933e-15, 3.72079545e-17;
       1.52737754e-07, -5.35686632e-10, 6.82322598e-13, -3.66809614e-16, 6.94600462e-20;
       -1.43717975e-11, 5.00673134e-14, -6.36551955e-17, 3.47371135e-20, -6.84230608e-24];


cn0 = [1.19569017e-09, 3.75918044e-11, -2.44737963e-13;
       1.61732587e-13, -2.13832839e-14, 9.30548445e-17;
       7.49270857e-17, 4.21761414e-18, -1.15125164e-20];
cn1 = [-4.11796175e+00, -3.81129454e-04, 4.30826594e-06, -9.16012013e-09, 8.01792467e-12;
       -4.81606702e-06, 7.74478673e-08, -6.98846761e-10, 1.91672052e-12, -1.76028859e-15;
       -1.82062504e-09, 1.44078593e-11, -2.08217075e-14, -3.60362511e-17, 7.40712432e-20];
cn2 = [-1.40308618e+03, 1.80259476e+00, -2.09727921e-04;
       3.81719502e-01, -5.39444475e-04, 1.85520370e-07;
       -6.44950116e-05, 8.43763766e-08, -2.71375500e-11;
       7.82381786e-09, -1.05383465e-11, 3.62959076e-15];
cn3 = [2.64811680e+11, -2.65241662e+11, 2.07124646e+13, 1.78718981e+13;
       -8.21923200e+08, 8.21923200e+08, -6.23137536e+10, -5.64731689e+10;
       8.44992000e+05, -8.44992000e+05, 6.23324160e+07, 5.93464320e+07;
       -2.88000000e+02, 2.88000000e+02, -2.07360000e+04, -2.07360000e+04];

cp = [-2.33680867e+00, -2.68897816e-04, 1.72494526e-08];
cx = [2.58413286e+02, 1.98390114e-02, -5.81896505e-06, 7.60177892e-10];
ct = [8.71957734e-03, -1.70527032e-06, 1.08279812e-10];
cj = [1.36221661e+00, -1.49851697e-04, 2.73875216e-08, -2.91620586e-12];
d = -2.325680393600000e-09;

pcrit = 3208.2;
PaToPsi = 1.45037738e-4;
JKgToBtuLbm = 4.29922614e-4;
F2K = @(T)((5/9)*(T-32)+273.15);
Ft3ToM3 = 6.24279605761446e-2;

%% Saturated liquid enthalpy as a function of pressure

% p = 800 -- Eq. 4a
p = 800;
lp = log(p);
h = 0.0;
for i = 0:8
    h = h + cf1(i+1) * lp^i;
end
fprintf('Saturated liquid enthalpy at p=%12.8e is %12.8e \n',p/PaToPsi,h/JKgToBtuLbm);

% p = 2000 -- Eq. 4b
p = 2000;
lp = log(p);
h = 0.0;
for i = 0:8
    h = h + cf2(i+1) * lp^i;
end
hf2000 = h; % Store for later 2-phase calculations
fprintf('Saturated liquid enthalpy at p=%12.8e is %12.8e \n',p/PaToPsi,h/JKgToBtuLbm);

% p = 3000 -- Eq. 4c
p = 3000;
h = 0.0;
for i = 0:8
    h = h + cf3(i+1) * ((pcrit-p)^0.41)^i;
end
fprintf('Saturated liquid enthalpy at p=%12.8e is %12.8e \n',p/PaToPsi,h/JKgToBtuLbm);

%% Saturated vapor enthalpy

% p = 800 -- Eq. 5a
p = 800;
lp = log(p);
h = 0.0;
for i = 0:11
    h = h + cg1(i+1) * lp^i;
end
fprintf('Saturated vapor enthalpy at p=%12.8e is %12.8e \n',p/PaToPsi,h/JKgToBtuLbm);

% p = 2000 -- Eq. 5b
p = 2000;
lp = log(p);
h = 0.0;
for i = 0:8
    h = h + cg2(i+1) * lp^i;
end
hg2000 = h;
fprintf('Saturated vapor enthalpy at p=%12.8e is %12.8e \n',p/PaToPsi,h/JKgToBtuLbm);

% p = 3000 -- Eq. 5c
p = 3000;
h = 0.0;
for i = 0:6
    h = h + cg3(i+1) * ((pcrit-p)^0.41)^i;
end
fprintf('Saturated vapor enthalpy at p=%12.8e is %12.8e \n',p/PaToPsi,h/JKgToBtuLbm);

%% Temperature as a function of p,h

% p=4000, h=400 -- Eq. 6b
p = 4000;
h = 400;
T = 0.0;
for i = 0:4
    for j = 0:4
        T = T + ct2(i+1,j+1) * p^i * h^j;
    end
end
fprintf('Temperature at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,F2K(T));

% p=4000, h=1400 -- Eq. 6d
p = 4000;
h = 1400;
T = 0.0;
for i = 0:4
    for j = 0:4
        T = T + ct4(i+1,j+1) * p^i * h^j;
    end
end
fprintf('Temperature at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,F2K(T));

% p=2000, h=400 -- Eq. 6a
p = 2000;
h = 400;
T = 0.0;
for i = 0:1
    for j = 0:3
        T = T + ct1(i+1,j+1) * p^i * h^j;
    end
end
fprintf('Temperature at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,F2K(T));

% p=2000, h=800 -- Eq. 6a at Hf(2000)
p = 2000;
h = 800;
T = 0.0;
for i = 0:1
    for j = 0:3
        T = T + ct1(i+1,j+1) * p^i * hf2000^j;
    end
end

fprintf('Temperature at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,F2K(T));
% p=2000, h=1400 -- Eq. 6c
p = 2000;
h = 1400;
T = 0.0;
for i = 0:4
    for j = 0:4
        T = T + ct3(i+1,j+1) * p^i * h^j;
    end
end
fprintf('Temperature at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,F2K(T));

%% Specific volume as a function of p, h

% p=4000, h=200 -- Eq. 8a
p = 4000;
h = 200;
V1 = 0.0;
for i = 0:2
    for j = 2:4
        V1 = V1 + cn0(i+1,j-1) * p^i * (250-h)^j;
    end
end
V2 = 0.0;
for i = 0:2
    for j = 0:4
        V2 = V2 + cn1(i+1,j+1) * p^i * h^j;
    end
end
V = V1 + exp(V2);
fprintf('Specific volume at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,V*Ft3ToM3);
       
% p=4000, h=600 -- Eq. 8b
p = 4000;
h = 600;
V = 0.0;
for i = 0:2
    for j = 0:4
        V = V + cn1(i+1,j+1) * p^i * h^j;
    end
end
V = exp(V);
fprintf('Specific volume at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,V*Ft3ToM3);

% p=4000, h=950 -- Eq. 8d
p = 4000;
h = 950;
C = zeros(4,1);
for i = 0:2
    C(1) = C(1) + cp(i+1)*p^i;
end
C(1) = exp(C(1));
for i = -1:2
    C(2) = C(2) + cx(i+2)*p^i;
end
for i = 0:2
    C(3) = C(3) + ct(i+1)*C(1)*p^i;
end
for i = -1:2
    C(4) = C(4) + cj(i+2)*p^i;
end
V = 0.0;
for i = 1:4
    for j = 1:4
        V = V + cn3(i,j)*C(j)*h^(i-1)*d;
    end
end
fprintf('Specific volume at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,V*Ft3ToM3);

% p=4000, h=1400 -- Eq. 8c
p = 4000;
h = 1400;
V = 0.0;
for i = -1:2
    for j = 0:2
        V = V + cn2(i+2,j+1) * p^i * h^j;
    end
end
fprintf('Specific volume at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,V*Ft3ToM3);

% p=2000, h=200 -- Eq. 8a
p = 2000;
h = 200;
V1 = 0.0;
for i = 0:2
    for j = 2:4
        V1 = V1 + cn0(i+1,j-1) * p^i * (250-h)^j;
    end
end
V2 = 0.0;
for i = 0:2
    for j = 0:4
        V2 = V2 + cn1(i+1,j+1) * p^i * h^j;
    end
end
V = V1 + exp(V2);
fprintf('Specific volume at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,V*Ft3ToM3);

% p=2000, h=500 -- Eq. 8b
p = 2000;
h = 500;
V = 0.0;
for i = 0:2
    for j = 0:4
        V = V + cn1(i+1,j+1) * p^i * h^j;
    end
end
V = exp(V);
fprintf('Specific volume at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,V*Ft3ToM3);

% p=2000, h=800 -- Weighted average of 8b with hf(2000) and 8c with
% hg(2000)
p = 2000;
h = 800;
X = (h-hf2000)/(hg2000-hf2000);
Vb = 0.0;
for i = 0:2
    for j = 0:4
        Vb = Vb + cn1(i+1,j+1) * p^i * hf2000^j;
    end
end
Vb = exp(Vb);
Vc = 0.0;
for i = -1:2
    for j = 0:2
        Vc = Vc + cn2(i+2,j+1) * p^i * hg2000^j;
    end
end
V = Vb + X*(Vc-Vb);
fprintf('Specific volume at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,V*Ft3ToM3);

% p=2000, h=1400 -- Eq. 8c
p = 2000;
h = 1400;
V = 0.0;
for i = -1:2
    for j = 0:2
        V = V + cn2(i+2,j+1) * p^i * h^j;
    end
end
fprintf('Specific volume at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,V*Ft3ToM3);

% p=3000, h=1040 -- Eq. 8d
p = 3000;
h = 1040;
C = zeros(4,1);
for i = 0:2
    C(1) = C(1) + cp(i+1)*p^i;
end
C(1) = exp(C(1));
for i = -1:2
    C(2) = C(2) + cx(i+2)*p^i;
end
for i = 0:2
    C(3) = C(3) + ct(i+1)*C(1)*p^i;
end
for i = -1:2
    C(4) = C(4) + cj(i+2)*p^i;
end
V = 0.0;
for i = 1:4
    for j = 1:4
        V = V + cn3(i,j)*C(j)*h^(i-1)*d;
    end
end
fprintf('Specific volume at h=%12.8e, p=%12.8e is %12.8e \n',h/JKgToBtuLbm,p/PaToPsi,V*Ft3ToM3);
