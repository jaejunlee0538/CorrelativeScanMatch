function [v_ac] = oplus( v_ab, v_bc )
c1 = cos(v_ab(3)); s1 = sin(v_ab(3));
v_ac = zeros(3,1);
v_ac(1) = c1*v_bc(1) - s1*v_bc(2) + v_ab(1);
v_ac(2) = s1*v_bc(1) + c1*v_bc(2) + v_ab(2);
v_ac(3) = normalizeHeading(v_ab(3) + v_bc(3));

end

