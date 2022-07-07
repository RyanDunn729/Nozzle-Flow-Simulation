function dfdy = ddy_fwd(f,dy)
    % Credit: Ryan Dunn
    dfdy = zeros(size(f));
    for i = 1:size(f,1)
        for j = 1:size(f,2)-1
            dfdy(i,j) = (f(i,j+1)-f(i,j))/dy;
        end
        j = size(f,2);
        dfdy(i,j) = (f(i,j)-f(i,j-1))/dy;
    end
end