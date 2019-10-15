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

nA1 = nA1/(24*3600) ; %rev/s

e_vecA1 = (1/mu_e)*(rA1.*(norm(vA1)^2-mu_e/norm(rA1)) - (norm(rA1)*dot(rA1/norm(rA1),vA1)).*vA1);
thetaA1 = acos(dot(e_vecA1/eccA1,rA1/norm(rA1)));  %radians, position of A at t=0     

%% Chaser (B) Intitial State 
eccB1 = 0 ;
incB1 = incA1 ; %save that delta v!!! 
thetaB1 = thetaA1 ; 
drelrB1 = [100; 0; 0]; %km, not close enough yet, drelxo so B is orbiting above A   
rB1 = rA1 + drelrB1 ;   %km, initial position of chaser
vB1 = -sqrt(mu_e/norm(rB1)).*[-sin(thetaB1); eccB1+cos(thetaB1); 0] ;
 
h_AB = cross(rA1, vA1) ;
    %comoving frame, LVLH    
i_A = rA1/norm(rA1) ;
k_A = h_AB/norm(h_AB) ; 
j_A = cross(k_A,i_A) ;

Qxx = [i_A'; j_A'; k_A'] ;

av_A = h_AB/(norm(rA1)^2) ;    %angular velocity of A

drelv_B1 = (vB1 - vA1 - cross(av_A,drelrB1)) ;
drelvB1 = Qxx*drelv_B1 ;

%% Football: 20-40km relative position
drelvB2 = [(nA1*40*2); 0; 0] ; 
drelrB2 = [0; 40; 0] ;
tot_delta2 = drelvB2*2 ; 

%% Vbar - Station Keeping 1km 
drelvB3 = [0 ; 0; 0 ]; 
drelrB3 = [0; 1; 0] ;

[tot_delta3] = VbarStationkeeping(15*60, nA1, drelrB3, drelrB2, drelvB2, drelvB3) ;

%% Vbar - Station Keeping 300m 
drelvB4 = [0; 0; 0] ;
drelrB4 = [0; .3; 0] ;

[tot_delta4] = VbarStationkeeping(15*60, nA1, drelrB4, drelrB3, drelvB3, drelvB4) ;

%% Coelliptical 20 m 
drelvB5 = [0; ((-3/2)*nA1*.02); 0] ; 
drelrB5 = [0; .020; 0] ; 

[tot_delta5] = VbarStationkeeping(15*60, nA1, drelrB5, drelrB4, drelvB4, drelvB5) ;
%% Hop to 5 m 

%% Vbar Approach to 0m 

tot_deltav = tot_delta2 + tot_delta3 + tot_delta4 + tot_delta5 ;
%% Potentially Useful Functions


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

%TLEs to State
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

%Relative Kinematics 
function [r_rel, v_rel] = relrva (r_A, v_A, r_B, v_B, mu_earth) 
%Relative Position 
h_AB = cross(r_A, v_A) ;
    %comoving frame, LVLH    
i_A = r_A/norm(r_A) ;
k_A = h_AB/norm(h_AB) ; 
j_A = cross(k_A,i_A) ;

Qxx = [i_A'; j_A'; k_A'] ;

av_A = h_AB/(norm(r_A)^2) ;    %angular velocity of A
aa_A = -2*(dot(v_A,r_A)/(norm(r_A)^2)).*av_A ;    %angular acceleration of A

abs_A = (-mu_earth/(norm(r_A)^3))*r_A ;     %absolute acceleration
abs_B = (-mu_earth/(norm(r_B)^3))*r_B ;   %absolute acceleration

r_relAB = r_B - r_A ;
r_rel = Qxx*r_relAB ;

v_relAB = (v_B - v_A - cross(av_A,r_relAB))*1000 ;
v_rel = Qxx*v_relAB ;

% a_relAB = (abs_B - abs_A - cross(aa_A,r_relAB) - cross(av_A,(cross(av_A,r_relAB))) - (2*(cross(av_A,(v_relAB/1000)))))*1000 ;
% a_rel = Qxx*a_relAB ;

end 

%Two Impulse Relative 
function [deltav] = VbarStationkeeping(t,n,deltarfinal,deltarinital,vbeforeburn1,vafterburn2)
%Outcome of maneuver is v-bar stationkeeping. Get onto orbit by using two
%impulse maneuver

%Inputs: maneuver time (t) in sec, mean motion,relative position after
%maneuver, relative position before maneuver, velocity of chaser before
%first impulse, velocity of chaser after second impulse

%Outputs: delta v for first impulse (vector form), delta v for second 
%impulse (vector form), total delta v for entire maneuver (scalar)

phirr = [4-3*cos(n*t),0,0;...
         6*(sin(n*t)-n*t),1,0;...
         0,0,cos(n*t)];
phirv = [(1/n)*sin(n*t),(2/n)*(1-cos(n*t)),0;...
         (2/n)*(cos(n*t)-1),(1/n)*(4*sin(n*t)-3*n*t),0;...
         0,0,(1/n)*sin(n*t)];
phivr = [3*n*sin(n*t),0,0;...
         6*n*(cos(n*t)-1),0,0;...
         0,0,-n*sin(n*t)];
phivv = [cos(n*t),2*sin(n*t),0;...
         -2*sin(n*t),4*cos(n*t)-3,0;...
         0,0,cos(n*t)];
     
vafterburn1 = inv(phirv)*(deltarfinal-phirr*deltarinital);
v0 = vafterburn1 - vbeforeburn1;

vbeforeburn2 = phivr*deltarinital+phivv*vafterburn1;
vf = vafterburn2 - vbeforeburn2;

deltav0 = norm(v0); %total delta v required for first impulse
deltavf = norm(vf); %total delta v required for second impulse
deltav = deltav0+deltavf; %total delta v required for maneuver
end

function [xx,yy,zz] = earth_sphere(varargin)
%EARTH_SPHERE Generate an earth-sized sphere.
%   [X,Y,Z] = EARTH_SPHERE(N) generates three (N+1)-by-(N+1)
%   matrices so that SURFACE(X,Y,Z) produces a sphere equal to 
%   the radius of the earth in kilometers. The continents will be
%   displayed.
%
%   [X,Y,Z] = EARTH_SPHERE uses N = 50.
%
%   EARTH_SPHERE(N) and just EARTH_SPHERE graph the earth as a 
%   SURFACE and do not return anything.
%
%   EARTH_SPHERE(N,'mile') graphs the earth with miles as the unit rather
%   than kilometers. Other valid inputs are 'ft' 'm' 'nm' 'miles' and 'AU'
%   for feet, meters, nautical miles, miles, and astronomical units
%   respectively.
%
%   EARTH_SPHERE(AX,...) plots into AX instead of GCA.
% 
%  Examples: 
%    earth_sphere('nm') produces an earth-sized sphere in nautical miles
%
%    earth_sphere(10,'AU') produces 10 point mesh of the Earth in
%    astronomical units
%
%    h1 = gca;
%    earth_sphere(h1,'mile')
%    hold on
%    plot3(x,y,z)
%      produces the Earth in miles on axis h1 and plots a trajectory from
%      variables x, y, and z
%   Clay M. Thompson 4-24-1991, CBM 8-21-92.
%   Will Campbell, 3-30-2010
%   Copyright 1984-2010 The MathWorks, Inc. 
%% Input Handling
[cax,args,nargs] = axescheck(varargin{:}); % Parse possible Axes input
error(nargchk(0,2,nargs)); % Ensure there are a valid number of inputs
% Handle remaining inputs.
% Should have 0 or 1 string input, 0 or 1 numeric input
j = 0;
k = 0;
n = 50; % default value
units = 'km'; % default value
for i = 1:nargs
    if ischar(args{i})
        units = args{i};
        j = j+1;
    elseif isnumeric(args{i})
        n = args{i};
        k = k+1;
    end
end
if j > 1 || k > 1
    error('Invalid input types')
end
%% Calculations
% Scale factors
Scale = {'km' 'm'  'mile'            'miles'           'nm'              'au'                 'ft';
         1    1000 0.621371192237334 0.621371192237334 0.539956803455724 6.6845871226706e-009 3280.839895};
% Identify which scale to use
try
    myscale = 6378.1363*Scale{2,strcmpi(Scale(1,:),units)};
catch %#ok<*CTCH>
    error('Invalid units requested. Please use m, km, ft, mile, miles, nm, or AU')
end
     
% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;
x = myscale*cosphi*cos(theta);
y = myscale*cosphi*sintheta;
z = myscale*sin(phi)*ones(1,n+1);
%% Plotting
if nargout == 0
    cax = newplot(cax);
    % Load and define topographic data
    load('topo.mat','topo','topomap1');
    % Rotate data to be consistent with the Earth-Centered-Earth-Fixed
    % coordinate conventions. X axis goes through the prime meridian.
    % http://en.wikipedia.org/wiki/Geodetic_system#Earth_Centred_Earth_Fixed_.28ECEF_or_ECF.29_coordinates
    %
    % Note that if you plot orbit trajectories in the Earth-Centered-
    % Inertial, the orientation of the contintents will be misleading.
    topo2 = [topo(:,181:360) topo(:,1:180)]; %#ok<NODEF>
    
    % Define surface settings
    props.FaceColor= 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = topo2;
    % Create the sphere with Earth topography and adjust colormap
    surface(x,y,z,props,'parent',cax)
    colormap(topomap1)
% Replace the calls to surface and colormap with these lines if you do 
% not want the Earth's topography displayed.
%     surf(x,y,z,'parent',cax)
%     shading flat
%     colormap gray
    
    % Refine figure
    axis equal
    xlabel(['X [' units ']'])
    ylabel(['Y [' units ']'])
    zlabel(['Z [' units ']'])
    view(127.5,30)
else
    xx = x; yy = y; zz = z;
end
end

