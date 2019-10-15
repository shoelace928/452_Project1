function [v0,vf,deltav] = VbarStationkeeping(t,n,deltarfinal,deltarinital,vbeforeburn1,vafterburn2)
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

