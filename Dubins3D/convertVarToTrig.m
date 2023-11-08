function [out] = convertVarToTrig(in)

% function to convert short-hand trig to trig 


syms S1 C1 S2 C2 S4 C4 S5 C5 
syms theta1 theta2 d3 theta4 theta5 
syms lx ly lz mx my mz nx ny nz rhox rhoy rhoz 

assume([S1 C1 S2 C2 S4 C4 S5 C5],'real')

assume([theta1,theta2,d3,theta4,theta5],'real')

assume([lx ly lz mx my mz nx ny nz rhox rhoy rhoz],'real')



out=sym(zeros(length(in),1));

for i = 1:length(in)

    outVar=subs(in(i),{S1,C1,S2,C2,S4,C4,S5,C5}, ...
        [sin(theta1),cos(theta1),sin(theta2),cos(theta2), ...
    sin(theta4),cos(theta4),sin(theta5),cos(theta5)]);

    out(i)=outVar;

end


out=out;

end