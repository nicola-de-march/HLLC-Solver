
% slope limiter, which maks the scheme non linear in the sense of Godunov
function c=minmod(a,b)

nVar = length(a);
c = zeros(nVar, 1);

for i=1:nVar
    if (a(i)*b(i) <= 0)
        % If the slope change sign -> no reconstruction
        c(i) = 0;
    else
        % slope that not change sign -> take the smaller in absolute value
        if (abs(a(i)) <= abs(b(i)))
            c(i) = a(i);
        else
            c(i) = b(i);
        end
end
end