function dfdy = ddy_central(f,dy)
    % Credit: Ryan Dunn
    dfdy = zeros(size(f));
    for i = 1:size(f,1)
        j = 1;
        dfdy(i,j) = (-3*f(i,j)+4*f(i,j+1)-f(i,j+2))/2/dy;
        for j = 2:size(f,2)-1
            dfdy(i,j) = (f(i,j+1)-f(i,j-1))/(2*dy);
        end
        j = size(f,2);
        dfdy(i,j) = (3*f(i,j)-4*f(i,j-1)+f(i,j-2))/2/dy;
    end
end