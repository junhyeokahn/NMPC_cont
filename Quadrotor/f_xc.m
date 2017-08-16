function f_ctrl = f_xc(xc)
%#codegen
g = 9.81;

b_T =  [sin(xc(8)); -cos(xc(8))*sin(xc(7)); cos(xc(8))*cos(xc(7))];

f_ctrl =     [xc(4:6);
                [0;0;g] - xc(9)*b_T;
                zeros(3,1)];
end