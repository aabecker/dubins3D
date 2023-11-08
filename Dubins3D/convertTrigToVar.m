function [out] = convertTrigToVar(in)
%
% function to convert trig to short-hand trig 


syms S1 C1 S2 C2 S4 C4 S5 C5 
syms theta1 theta2 d3 theta4 theta5 
syms lx ly lz mx my mz nx ny nz rhox rhoy rhoz 

assume([S1 C1 S2 C2 S4 C4 S5 C5],'real')

assume([theta1,theta2,d3,theta4,theta5],'real')

% assume([lx ly lz mx my mz nx ny nz rhox rhoy rhoz],'real')


out=sym(zeros(length(in),1));

for i = 1:length(in)

% this is returning the equation with the subs, but without the
% Ahand coefficients??

% outVar=subs(in(i),{sin(theta1),cos(theta1),sin(theta2),cos(theta2), ...
%     sin(theta4),cos(theta4),sin(theta5),cos(theta5)}, ...
%     [S1,C1,S2,C2,S4,C4,S5,C5]);

% update: jk it never had coefficients
% outVar=subs(in(i),{lx,ly,lz,mx,my,mz,nx,ny,nz,rhox,rhoy,rhoz, ...
%     sin(theta1),cos(theta1),sin(theta2),cos(theta2), ...
%     sin(theta4),cos(theta4),sin(theta5),cos(theta5)}, ...
%     [lx,ly,lz,mx,my,mz,nx,ny,nz,rhox,rhoy,rhoz, ...
%     S1,C1,S2,C2,S4,C4,S5,C5]);

outVar=subs(in(i),{sin(theta1),cos(theta1),sin(theta2),cos(theta2), ...
    sin(theta4),cos(theta4),sin(theta5),cos(theta5)}, ...
    [S1,C1,S2,C2,S4,C4,S5,C5]);

out(i)=outVar;

end

out=out;


end
