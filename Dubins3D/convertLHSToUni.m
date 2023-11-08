function[out] = convertLHSToUni(in)

% converts multivariate terms to univariate terms only LHS
% but the trig has to be in short-hand form first

syms S1 C1 S2 C2 S4 C4 S5 C5 
syms theta1 theta2 d3 theta4 theta5 
syms lx ly lz mx my mz nx ny nz rhox rhoy rhoz 

assume([S1 C1 S2 C2 S4 C4 S5 C5],'real')

assume([theta1,theta2,d3,theta4,theta5],'real')

assume([lx ly lz mx my mz nx ny nz rhox rhoy rhoz],'real')

syms d3sqxS5 d3sqxC5 d3sq d3xS5 d3xC5


assume([d3sqxS5 d3sqxC5 d3sq d3xS5 d3xC5],'real')

out=sym(zeros(length(in),1));

for i = 1:length(in)


    outVar=subs(in(i),{d3^2 * S5,d3^2 * C5, d3^2 , d3 * S5, d3 * C5 }, ...
        [d3sqxS5 d3sqxC5 d3sq d3xS5 d3xC5]);

out(i)=outVar;

end


out=out;




end