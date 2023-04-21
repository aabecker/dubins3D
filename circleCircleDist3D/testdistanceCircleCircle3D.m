function testdistanceCircleCircle3D
% tests the function  distanceCircleCircle3D(N0, r0, C0, N1, r1, C1 )
% based on https://www.geometrictools.com/Documentation/DistanceToCircle3.pdf
%
% Aaron Becker, atbecker@uh.edu
%  Draws the circles, their closest points, and (optionally) the minimum
%  radius to make the tori touch.
%

format compact
showTorus = true;   %draws a torus around each circle so that the tori just touch.
%%% two circles in 3D: WORKS! The C code reverses order of polyomials compared to Matlab's
% computeAnddrawSolutions( figNum, N0, r0, C0, N1, r1, C1  )
% computeAnddrawSolutions( 1 , [1,0,0] , 2, [0,0,0], [0,1,0], 3, [3,4,1] ) %WORKS
% computeAnddrawSolutions( 2 , [1,0,0] , 2, [0,0,0], [0,1,0], 3, [3,1,1]) %WORKS
% computeAnddrawSolutions( 3 , [1,0,0] , 2, [0,0,0], [0,1,0], 3, [3,0,0]  )  % WORKS!
% computeAnddrawSolutions( 4 , [1,0,0] , 2, [0,0,0], [0,1,0], 3, [3,1,0]  ) % Works -- (had to eliminate small imaginary parts)
% computeAnddrawSolutions( 5 , [1,0,0] , 2, [0,0,0], [0,1,0], 4, [3,1,0]  )  % Success two solutions!
% computeAnddrawSolutions( 6 , [1,0.1,0] , 2, [0,0,0], [0,1,0], 4, [3,1,0]  )  % Success two solutions!
% computeAnddrawSolutions( 7 , [1,0.1,0] , 4, [0,0,0], [0,1,0], 4, [3,1,0]  )  % Success two solutions!
% computeAnddrawSolutions( 8, [1,0.1,0] , 4, [0,0,0], [.4,1,0], 4, [3,1,0]  )  % Success two solutions!
% computeAnddrawSolutions( 9, [0,1,0], 1, [0,0,0], [0,0,1], 1, [1.01,0,0]  )  % Success one solution!
%computeAnddrawSolutions( 10, [0,1,0], 1, [0,0,0], [0,0,1], 1, [1,0,0]  )  % added 0 case if no roots found -- special case if orthogonal and centered on the other circle?
computeAnddrawSolutions( 10, [0,1,0], 1, [0,0,0], [0,0,1], 0.5, [1,0,0]  )  % used to fail -- special case if orthogonal and centered on the other circle? fails if r1  = 0.5, 1, 1.5,2,2.5, 4.5

%computeAnddrawSolutions( 10, [0,1,0], 1, [0,0,0], [0,0,1], 1, [1,0.16,0])  % Success -- used to have problems, fixed with numerical precision

%computeAnddrawSolutions( 10, [0,1,0], 1, [0,0,0], [0,0,1], 1, [1,0.06,0]  )  % looks like the wrong intersecction to me.  The toroids are overlapping
%computeAnddrawSolutions( 10, [0,1,0], 1, [0,0,0], [0,0,1], 1, [1,0.16,0]  )  %fixed

%computeAnddrawSolutions( 10, [0,0,1], 1, [1,-0.1,0],[0,1,0], 1, [0,0,0]  )  % works if I swap the torus.

computeAnddrawSolutions( 12, [0,1,0], 1, [0,0,0], [0,.1,1], 2, [1,0,0]  )  % works -- special case if orgthogonal?
%computeAnddrawSolutions( 11, [0,1,1], 1, [0,0,0], [0,0,1], 1, [1,0,0]  )  % works -- 2 solutions

% generate a case to check the roots of p6 instead of phi.
%computeAnddrawSolutions( 13 , [0,0,1] ,2 , [0,0,0], [1,0,0], 2, [0,0,2]  ) % 5.3 Circles in Nonparallel Planes but Centers on Normal Line

% % % % the planar cases: they work!
% figure(20)
% showTorus = false;
% computeAnddrawSolutions( [4,3,1] , [0,0,1] ,2 , [0,0,0], [0,0,1], 2, [0,0,0]  ) % equal circles, works
% computeAnddrawSolutions( [4,3,2] , [0,0,1] ,3 , [0,0,0], [0,0,1], 2, [0,0,0]  ) % concentric circles, works
% computeAnddrawSolutions( [4,3,3] , [0,0,1] ,3 , [0,0,0], [0,0,1], 2, [4,0,0]  ) % overlapping circles 1, works
% computeAnddrawSolutions( [4,3,4],  [0,0,1] ,3 , [0,0,0], [0,0,1], 2, [3,2,0]  ) % overlapping circles 2, works
% computeAnddrawSolutions( [4,3,5] , [0,0,1] ,0.6 , [0,0,0], [0,0,1], 1.2, [3,2,0]  ) % separate circles 1, works
% computeAnddrawSolutions( [4,3,6] , [0,0,1] ,0.6 , [0,0,0], [0,0,1], 1.2, -[3,2,0]  ) % separate circles 2, works
% computeAnddrawSolutions( [4,3,7] , [0,0,1] ,1 , [0,0,0], [0,0,1], 2, -[3,0,0]  ) % touching circles 1, works
% computeAnddrawSolutions( [4,3,8] , [0,0,1] ,1 , [0,0,0], [0,0,1], 2, -[1,0,0]  ) % touching circles 1, works
% computeAnddrawSolutions( [4,3,9] , [0,0,1] ,4 , [0,0,0], [0,0,1], 2.5, -[1,0,0]  ) % inside circles 1, works
% computeAnddrawSolutions( [4,3,10] , [0,0,1] ,3 , [0,0,0], [0,0,1], 4.6, -[1,0,0]  ) % inside circles 1, works


    function computeAnddrawSolutions( figNum, N0, r0, C0, N1, r1, C1  )
        % given two circles, computes the distance and plots the results in
        % fig figNum. If figNum is a 3x1 matrix representing a subfigure,
        % draw this in the subfigure, otherwise it makes figure figNum

        if numel(figNum) == 1
            figure(figNum); clf;
            view(45,45)
        else
            subplot(figNum(1),figNum(2),figNum(3));cla;
        end
        plotCircle(N0, C0,r0, 'r', 1)
        plotCircle(N1, C1,r1, 'b', 1)
        hold on

        result = distanceCircleCircle3D(N0, r0, C0, N1, r1, C1 );
        disp(result)
        if showTorus
            drawTorus(N0, r0, C0, result.distance/2,'r')
            drawTorus(N1, r1, C1, result.distance/2,'b')
        end
        axis equal

        title(strcat('Min dist = ',num2str(result.distance, 3),", ",num2str(result.numClosestPairs), ' pairs, equidistant =  ',num2str(result.equidistant)))
        for i = 1:result.numClosestPairs
            plot3([result.circle0Closest(i,1),result.circle1Closest(i,1)],...
                [result.circle0Closest(i,2),result.circle1Closest(i,2)],...
                [result.circle0Closest(i,3),result.circle1Closest(i,3)],'.-','color','k','LineWidth',1,'MarkerSize',12)

        end

        hold off
    end

    function plotCircle(n_vec, PtCenter,r, colorLine, lineWidth)
        % plots an arc of length theta about the circle with normal vector
        % n_vec, radius r, and center PtCenter, starting at position PtLeave, in
        % direction dir (-1 for CW, 1 for CCW).  The line is drawn in color
        % colorLine with width lineWidth
        n_vec = n_vec/norm(n_vec); %ensure it is a unit vector
        theta = linspace(0,2*pi,40)';  % Draw the arc path
        [U,V,~] = ComputeOrthogonalComplement(n_vec);
        arc2Pts = PtCenter + U*r.*cos(theta) + V*r.*sin(theta);
        hold on
        harc = plot3( arc2Pts(:,1),arc2Pts(:,2),arc2Pts(:,3)); set(harc,'color' ,colorLine,'linewidth',lineWidth)
        plot3(PtCenter(1)+[0,n_vec(1)],PtCenter(2)+[0,n_vec(2)],PtCenter(3)+[0,n_vec(3)],'.-','color',colorLine,'MarkerSize',10,'linewidth',lineWidth);
        xlabel('x')
        ylabel('y')
        zlabel('z')
        hold off
    end

    function drawTorus(N, r, C, R,color)
        % draws a radius R torus in color around circle with normal N,
        % radius r and center C.
        N = N/norm(N); %normalize the matrix
        [theta,phi] = meshgrid(linspace(0,2*pi,20));

        x = (r + R*cos(theta)).*cos(phi);
        y = (r + R*cos(theta)).*sin(phi);
        z = R*sin(theta);

        if N(2)==0 && N(1)==0 % the longitude
            angRz = 0;
        else
            angRz = atan2(N(2), N(1));
        end
        angRz = angRz-pi/2;

        if N(3)==1  % the latitude
            angRx = 0;
        elseif N(3)==-1
            angRx = pi;
        else
            angRx = atan2(N(3),sqrt(N(1)^2+N(2)^2) )-pi/2;
        end
        % compute the rotation matrices.
        Rx = [1 0 0; 0 cos(angRx) -sin(angRx); 0 sin(angRx) cos(angRx)];
        Rz = [cos(angRz) -sin(angRz) 0; sin(angRz) cos(angRz) 0; 0 0 1];

        xF = x; yF = y; zF = z;
        for i = 1:numel(x)
            nCoord = Rz*Rx*[x(i);y(i);z(i)];
            xF(i) = nCoord(1); yF(i)= nCoord(2); zF(i)= nCoord(3);
        end
        %xyzR = Rx*Rz*[x;y;z];
        %coordFinal = subs(xyzR, t, pi/4);

        sPlot = surf(C(1)+ xF,C(2)+yF,C(3)+zF);
        %shading interp
        set(sPlot, 'facecolor',color,'EdgeColor', 'none')
        alpha(sPlot,0.1);
    end
    function [U,V,W] = ComputeOrthogonalComplement(W)
        % Robustly compute a right-handed orthonormal set { U, V, W }.
        % The vector W is guaranteed to be unit-length, in which case
        % there is no need to worry about a division by zero when
        % computing invLength.
        if (abs(W(1)) > abs(W(2)))
            % The component of maximum absolute value is either W[1]
            % or W[3].
            invLength = 1 / sqrt(W(1) * W(1) + W(3) * W(3));
            U = [ -W(3) * invLength, 0, +W(1) * invLength ];
        else
            % The component of maximum absolute value is either W[1]
            % or W[2].
            invLength = 1 / sqrt(W(2) * W(2) + W(3) * W(3));
            U = [0, +W(3) * invLength, -W(2) * invLength ];
        end
        V = cross(W, U);
    end

    function testTorus( N0, r0, C0  ) %#ok<DEFNU>
        % used to test plotting of the circle and torus in 3D.
        figure(1)
        clf;

        plotCircle(N0, C0,r0, 'r', 1)
        hold on
        drawTorus(N0, r0, C0, r0/2,'r')

        axis equal
        hold off
    end
end