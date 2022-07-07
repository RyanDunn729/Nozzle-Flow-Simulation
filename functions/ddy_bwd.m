function dfdy = ddy_bwd(f,dy)
    % Credit: Ryan Dunn
    dfdy = zeros(size(f));
    for i = 1:size(f,1)
        j = 1;
        dfdy(i,j) = (f(i,j+1)-f(i,j))/dy;
        for j = 2:size(f,2)
            dfdy(i,j) = (f(i,j)-f(i,j-1))/dy;
        end
    end
end