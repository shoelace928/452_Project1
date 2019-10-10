%% AERO452 Project 1
%% Collaborators: Lacey Davis and Ankit Maurya
% October 24, 2019 

close all; clear all; clc; 

%% Constants: 
mu_e = 398600 ; %km3/s2
TOL = 10^-8 ; 
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;

%% Target (A) Initial State

% Astra 1F as of October 3, 2019 at 3:16:47.9 UT
tle = [1 23842 96021 19275.63666510 .00000119 0 0 0 9994 ;2 23842 0.0474 31.7281 .0004336 154.1450 98.6995 1.00266434 85773] ;
[incA1, epochA1, RAANA1, eccA1, argA1, MeA1, nA1] = tle_convert(tle) ;
[rA1,vA1] = TLE_State(RAANA1,argA1,MeA1,nA1,incA1,eccA1) ;

%% Chaser (B) Intitial State 
eccB1 = 0 ;
incB1 = incA1 ; %save that delta v!!! 
delrB1 = [100; 0; 0]; %km, not close enough yet, drel_xo 
rB1 = rA1 + delrB1 ;   %km, initial position of chaser
vB1 = [0; 0; sqrt(mu_e/norm(rB1))] ;   %km/s, circular orbit, update to math new rB1


%% Potentially Useful Functions

%Julian Date
     function [ JD ] = julian(m,d,y,tf)
%julian is a function that will convert UT into Julian Date
%   inputs are month, day, year, and UT time fraction

Jo = 367*y - floor((7*(y+floor((m+9)/12)))/4) + floor((275*m)/9) + d + 1721013.5 ;

JD = Jo + (tf/24) ;

      end

%TLE to COEs 
     function [inc, epoch, RAAN, ecc, arg, Me, n] = tle_convert(tle)
        inc = tle(2,3) * (pi/180) ;   %radians, inclination
        epoch = tle(1,4) ;    %year and day fraction
        RAAN = tle(2,4) * (pi/180) ;  %radians, right ascension of ascending node
        ecc = tle(2,5)/10e6 ;  %eccentricity, divide by factors of 10 to move decimal to front
        arg = tle(2,6) * (pi/180) ;   %radians, argument of periapse
        Me = tle(2,7) * (pi/180) ;    %radians, mean anomaly at epoch
        n = tle(2,8) ;    %mean motion at epoch 
     end 

%COEs to State
function [R,V] = TLE_State(RAAN,OMEGA,ME,MM,INC,ecc)

%Given: RAAN,OMEGA = AoP, ME (Mean Anomaly), Mean Motion = MM, Inc, ecc 
muearth = 398600;
tol = 10^-8;
MM = MM *1/(60^2*24)*2*pi;
a = (muearth/(MM^2))^(1/3);
h = sqrt(a*muearth*(1-ecc^2));
 

        if ME < pi % finding correct guess for E
            E = ME + ecc/2;
        elseif ME > pi 
            E = ME - ecc/2;
        end
        
        ratio = 1; % sets initial ratio
        n = 1; % sets initial iteration value
        
        while ratio(end) > tol % since ratio values are being stored, end us used
            ratio(n) = (E(n) - ecc*sin(E(n)) - ME)/(1 - ecc*cos(E(n))); % ratio calc
            E(n+1) = E(n) - ratio(n); % Calculating E 
            n = n+1; % increasing iteration value
        end
        
   tantheta2 = (sqrt((1+ecc)/(1-ecc))) * tan((E(end)/2)); % goes into true anomaly calc
   TA = atan(tantheta2) * 2; % calculates true anomaly (radians)
   %TA = rad2deg(TA1);
   
R_peri = h^2/muearth*(1/(1+ecc*cos(TA)))*[cos(TA) sin(TA) 0]; % in the p hat direction 
V_peri = muearth/h * [-sin(TA) ecc+cos(TA) 0]; % in the q hat direction

QxX = [-sin(RAAN)*cos(INC)*sin(OMEGA)+cos(RAAN)*cos(OMEGA)...
-sin(RAAN)*cos(INC)*cos(OMEGA)-cos(RAAN)*sin(OMEGA) sin(RAAN)...
*sin(INC); cos(RAAN)*cos(INC)*sin(OMEGA)+sin(RAAN)*cos(OMEGA)...
cos(RAAN)*cos(INC)*cos(OMEGA)-sin(RAAN)*sin(OMEGA) -cos(RAAN)...
*sin(INC); sin(INC)*sin(OMEGA) sin(INC)*cos(OMEGA) cos(INC)];

R = QxX*transpose(R_peri);
V = QxX*transpose(V_peri);

end
