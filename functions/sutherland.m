function mu = sutherland(T)
    % Credit: Ryan Dunn
    S1 = 110.4;
    T0 = 288;
    mu0 = 1.735*1e-5;
    mu = mu0*(T0+S1)./(T+S1).*(T./T0).^(3/2);
end