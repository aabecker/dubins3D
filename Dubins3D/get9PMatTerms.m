function [out] = get9PMatTerms(in)

syms S1 C1 S2 C2 S4 C4 S5 C5 
syms theta1 theta2 d3 theta4 theta5 
syms lx ly lz mx my mz nx ny nz rhox rhoy rhoz 
assume([S1 C1 S2 C2 S4 C4 S5 C5],'real')
assume([theta1,theta2,d3,theta4,theta5],'real')
assume([lx ly lz mx my mz nx ny nz rhox rhoy rhoz],'real')

syms d3sqxS5 d3sqxC5 d3sq d3xS5 d3xC5
assume([d3sqxS5 d3sqxC5 d3sq d3xS5 d3xC5],'real')

% out=sym([0,0,0,0,0,0,0,0,0]);

out = sym(zeros(14,9))
disp('in length:')
disp(length(in))
% Mathematica reference 
% term1 = (func /. {d3sqxS5 -> 1, d3sqxC5 -> 0, d3sq -> 0, d3xS5 -> 0, 
      %d3xC5 -> 0, d3 -> 0, S5 -> 0, C5 -> 0}) - const;

      for i = 1:length(in)

        const = subs(in(i),{d3sqxS5, d3sqxC5, d3sq, d3xS5, d3xC5, d3, S5, C5}, ...
         [0,0,0,0,0,0,0,0]);

        term1=subs(in(i),{d3sqxS5, d3sqxC5, d3sq, d3xS5, d3xC5, d3, S5, C5}, ...
         [1,0,0,0,0,0,0,0]);

        term2=subs(in(i),{d3sqxS5, d3sqxC5, d3sq, d3xS5, d3xC5, d3, S5, C5}, ...
         [0,1,0,0,0,0,0,0]);
        term3=subs(in(i),{d3sqxS5, d3sqxC5, d3sq, d3xS5, d3xC5, d3, S5, C5}, ...
         [0,0,1,0,0,0,0,0]);
        term4=subs(in(i),{d3sqxS5, d3sqxC5, d3sq, d3xS5, d3xC5, d3, S5, C5}, ...
         [0,0,0,1,0,0,0,0]);
        term5=subs(in(i),{d3sqxS5, d3sqxC5, d3sq, d3xS5, d3xC5, d3, S5, C5}, ...
         [0,0,0,0,1,0,0,0]);
        term6=subs(in(i),{d3sqxS5, d3sqxC5, d3sq, d3xS5, d3xC5, d3, S5, C5}, ...
         [0,0,0,0,0,1,0,0]);
        term7=subs(in(i),{d3sqxS5, d3sqxC5, d3sq, d3xS5, d3xC5, d3, S5, C5}, ...
         [0,0,0,0,0,0,1,0]);
        term8=subs(in(i),{d3sqxS5, d3sqxC5, d3sq, d3xS5, d3xC5, d3, S5, C5}, ...
         [0,0,0,0,0,0,0,1]);



        out(i,1)=term1 - const
        out(i,2)=term2 - const
        out(i,3)=term3 - const
        out(i,4)=term4 - const
        out(i,5)=term5 - const
        out(i,6)=term6 - const
        out(i,7)=term7 - const
        out(i,8)=term8 - const
 
        out(i,9)=const


      end
% out(1)=1;
% out(9)=const;


end