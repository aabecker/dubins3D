function [out] = get9QMatTerms(in)

syms S1 C1 S2 C2 S4 C4 S5 C5 
syms theta1 theta2 d3 theta4 theta5 
syms lx ly lz mx my mz nx ny nz rhox rhoy rhoz 
assume([S1 C1 S2 C2 S4 C4 S5 C5],'real')
assume([theta1,theta2,d3,theta4,theta5],'real')
assume([lx ly lz mx my mz nx ny nz rhox rhoy rhoz],'real')

syms S1xS2 S1xC2 C1xS2 C1xC2


assume([S1xS2 S1xC2 C1xS2 C1xC2],'real')

out = sym(zeros(14,9))

for i = 1:length(in)

    % mathematica reference 
    % term1 = (func /. {S1xS2 -> 1, S1xC2 -> 0, C1xS2 -> 0, C1xC2 -> 0, 
     % S1 -> 0, C1 -> 0, S2 -> 0, C2 -> 0}) - const;


    const = subs(in(i),{S1xS2, S1xC2, C1xS2, C1xC2,S1, C1, S2, C2}, ...
     [0,0,0,0,0,0,0,0]);


    term1 = subs(in(i),{S1xS2, S1xC2, C1xS2, C1xC2,S1, C1, S2, C2}, ...
     [1,0,0,0,0,0,0,0])

    term2 = subs(in(i),{S1xS2, S1xC2, C1xS2, C1xC2,S1, C1, S2, C2}, ...
     [0,1,0,0,0,0,0,0])

    term3 = subs(in(i),{S1xS2, S1xC2, C1xS2, C1xC2,S1, C1, S2, C2}, ...
     [0,0,1,0,0,0,0,0])

    term4 = subs(in(i),{S1xS2, S1xC2, C1xS2, C1xC2,S1, C1, S2, C2}, ...
     [0,0,0,1,0,0,0,0])

    term5 = subs(in(i),{S1xS2, S1xC2, C1xS2, C1xC2,S1, C1, S2, C2}, ...
     [0,0,0,0,1,0,0,0])

    term6 = subs(in(i),{S1xS2, S1xC2, C1xS2, C1xC2,S1, C1, S2, C2}, ...
     [0,0,0,0,0,1,0,0])

    term7 = subs(in(i),{S1xS2, S1xC2, C1xS2, C1xC2,S1, C1, S2, C2}, ...
     [0,0,0,0,0,0,1,0])

    term8 = subs(in(i),{S1xS2, S1xC2, C1xS2, C1xC2,S1, C1, S2, C2}, ...
     [0,0,0,0,0,0,0,1])



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



end


