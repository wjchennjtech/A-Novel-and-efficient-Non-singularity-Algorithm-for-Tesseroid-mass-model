function [Cnm, Snm] = SHA_single_tess(lon1, lon2, lat1, lat2, nmax, N_theta)

    lambda1 = lon1 * pi/180;
    lambda2 = lon2 * pi/180;
    phi1 = lat1 * pi/180;
    phi2 = lat2 * pi/180;


    theta1 = pi/2 - phi2;
    theta2 = pi/2 - phi1;

    x_min = cos(theta2);
    x_max = cos(theta1);

    nmax1 = nmax + 1;
    Cnm = zeros(nmax1, nmax1);
    Snm = zeros(nmax1, nmax1);

    %The longitudinal integration is evaluated analytically.
    I_cos = zeros(nmax1, 1);
    I_sin = zeros(nmax1, 1);
    for m = 0:nmax
        if m == 0
            I_cos(1) = lambda2 - lambda1;
            I_sin(1) = 0;
        else
            I_cos(m+1) = (sin(m*lambda2) - sin(m*lambda1)) / m;
            I_sin(m+1) = (cos(m*lambda1) - cos(m*lambda2)) / m;
        end
    end
    
    %Gauss–Legendre quadrature is applied for the integration over colatitude
    [x_k, w_k] = GaussLegendre(x_min, x_max, N_theta);

    for k = 1:N_theta
        x = x_k(k);
        w = w_k(k);
        theta_deg = acos(x) * 180/pi;
        for m = 0:nmax
            P_all_m = plm(0:nmax, m, theta_deg);
            P_all_m = P_all_m(:);
            Cnm(:, m+1) = Cnm(:, m+1) + w * P_all_m .* I_cos(m+1);
            Snm(:, m+1) = Snm(:, m+1) + w * P_all_m .* I_sin(m+1);
        end
    end


    Cnm = Cnm / (4*pi);
    Snm = Snm / (4*pi);
end