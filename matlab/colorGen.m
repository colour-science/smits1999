function [All,A] = colorGen(minLambda,maxLambda,bins)
% Args: minLambda, maxLambda, numBins
% computes a matrix where each column is a spectrum (white, cyan, magenta, 
% yellow, red, green, blue).

global A;
global x;
global N;
global An;

N = bins;
A = Init(minLambda,maxLambda);
An = null(A);

white = [1 1 1]';
cyan = [0 1 1]';
magenta = [1 0 1]';
yellow = [1 1 0]';
red = [1 0 0]';
green = [0 1 0]';
blue = [0 0 1]';

options(14) = 10000;

x = pinv(A) *white;
whiteS = x + null(A) * constr('colorFun',zeros(N-3,1),options);

x = pinv(A) *cyan;
cyanS = x + null(A) * constr('colorFun',zeros(N-3,1),options);

x = pinv(A) *magenta;
magentaS = x + null(A) * constr('colorFun',zeros(N-3,1),options);

x = pinv(A) *yellow;
yellowS = x + null(A) * constr('colorFun',zeros(N-3,1),options);

x = pinv(A) *red;
redS = x + null(A) * constr('colorFun',zeros(N-3,1),options);

x = pinv(A) *green;
greenS = x + null(A) * constr('colorFun',zeros(N-3,1),options);

x = pinv(A) *blue;
blueS = x + null(A) * constr('colorFun',zeros(N-3,1),options);

All = [whiteS cyanS magentaS yellowS redS greenS blueS];



%
%%%
%

function A = Init(minLambda, maxLambda)
global N;

XYZtoRGB = InitXYZtoRGB;

vals = linspace(minLambda,maxLambda,N+1);	
sum = zeros(3,1);
for i = 1:N
   A(:,i) = zeros(3,1);
   for j = vals(i):(vals(i+1) -1)   %don't add the endpoints twice...
       sum = sum + GetXYZ(j);       %NOTE: this cause problems for non-integer
       A(:,i) = A(:,i) + GetRGB(j); %endpoints....
   end;
end;
sum
return;


%
%%%
%


function rgb = GetRGB(lambda) 
global XYZtoRGB
rgb = XYZtoRGB * GetXYZ(lambda);
return;


%
%%%
%



function xyz = GetXYZ(lambda)

riXCurve = [ 0.00013 0.00023 0.00041 0.00074 0.00137 0.00223 0.00424 0.00765 0.01431 0.02319 0.04351 0.07763 0.13438 0.21477 0.2839 0.3285 0.34828 0.34806 0.3362 0.3187 0.2908 0.2511 0.19536 0.1421 0.09564 0.05795 0.03201 0.0147 0.0049 0.0024 0.0093 0.0291 0.06327 0.1096 0.1655 0.22575 0.2904 0.3597 0.43345 0.51205 0.5945 0.6784 0.7621 0.8425 0.9163 0.9786 1.0263 1.0567 1.0622 1.0456 1.0026 0.9384 0.85445 0.7514 0.6424 0.5419 0.4479 0.3608 0.2835 0.2187 0.1649 0.1212 0.0874 0.0636 0.04677 0.0329 0.0227 0.01584 0.01136 0.00811 0.00579 0.00411 0.00289 0.00205 0.00144 0.001 0.00069 0.00048 0.00033 0.00023 0.00017 0.00012 8e-05 6e-05 4.1e-05 2.9e-05 2e-05 1.4e-05 1e-05];

riYCurve = [
 0 0 1e-05 2e-05 4e-05 6e-05 0.00012 0.00022 0.0004 0.00064 0.0012 0.00218 0.004 0.0073 0.0116 0.01684 0.023 0.0298 0.038 0.048 0.06 0.0739 0.09098 0.1126 0.13902 0.1693 0.20802 0.2586 0.323 0.4073 0.503 0.6082 0.71 0.7932 0.862 0.91485 0.954 0.9803 0.99495 1 0.995 0.9786 0.952 0.9154 0.87 0.8163 0.757 0.6949 0.631 0.5668 0.503 0.4412 0.381 0.321 0.265 0.217 0.175 0.1382 0.107 0.0816 0.061 0.04458 0.032 0.0232 0.017 0.01192 0.00821 0.00573 0.0041 0.00293 0.00209 0.00105 0.00105 0.00074 0.00052 0.00036 0.00025 0.00017 0.00012 8e-05 6e-05 4e-05 3e-05 2e-05 1.4e-05 1e-05 7e-06 5e-06 3e-06];

riZCurve = [ 0.00061 0.00108 0.00195 0.00349 0.00645 0.01055 0.02005 0.03621 0.06785  0.1102 0.2074 0.3713 0.6456 1.03905 1.3856 1.62296 1.74706 1.7826 1.77211 1.7441 1.6692 1.5281 1.28764 1.0419 0.81295 0.6162 0.46518 0.3533 0.272 0.2123 0.1582 0.1117 0.07825 0.05725 0.04216 0.02984 0.0203 0.0134 0.00875 0.00575 0.0039 0.00275 0.0021 0.0018 0.00165 0.0014 0.0011 0.001 0.0008 0.0006 0.00034 0.00024 0.00019 0.0001 5e-05 3e-05 2e-05 1e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

ls = (360:5:800);

xyz = [interp1(ls,riXCurve,lambda)
       interp1(ls,riYCurve,lambda) 
       interp1(ls,riZCurve,lambda)];
return;

%
%%%
%

function XYZtoRGB = InitXYZtoRGB()
global XYZtoRGB;
%rx = .625; ry = .340;  %sony data
%gx = .280; gy = .595;
%bx = .155; by = .070;
%rx = .6; ry = .350;  %less saturated version
%gx = .29; gy = .58;
%bx = .175; by = .090;
%rx = .628; ry = .346;   %CRT from Rogers
%gx = .268; gy = .588;
%bx = .150; by = .070;
%rx = .67; ry = .33;   %NTSC from Rogers
%gx = .21; gy = .71;
%bx = .140; by = .080;
rx = .64; ry = .33;   %data from W3C http://www.w3.org/Graphics/Color/sRGB 
gx = .3;  gy = .6;
bx = .15; by = .06;
wx = .3333; wy = .3333; % equal energy white
%wx = .283; wy = .298;  % sony white
wY = 106.8;  %should be the area under the XYZ curves, picked to make white = 1

Yovery = wY / wy;
D =  (rx * (gy - by) + gx * (by - ry) + bx * (ry - gy));
Cr = (Yovery  * (wx * (gy - by) - wy * (gx - bx) + gx * by - bx * gy) / D);
Cg = (Yovery  * (wx * (by - ry) - wy * (bx - rx) + bx * ry - rx * by) / D);
Cb = (Yovery  * (wx * (ry - gy) - wy * (rx - gx) + rx * gy - gx * ry) / D);

toXYZ = [rx * Cr  gx * Cg   bx * Cb 
         ry * Cr  gy * Cg   by * Cb
         (1.0 - (rx + ry)) * Cr    (1.0 - (gx + gy)) * Cg   (1.0 - (bx + by)) * Cb];

XYZtoRGB = inv(toXYZ);
