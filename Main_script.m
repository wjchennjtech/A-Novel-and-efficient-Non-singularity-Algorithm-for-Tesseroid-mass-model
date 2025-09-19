clc
clear

%Adjustable variables
nmax = 180;
N_theta =2;
height = 250000;   

%Boundary data and density data
upper=load("Topo_Tibet_bound1.txt");lower=load("Topo_Tibet_bound2.txt");density=load("Topo_Tibet_density.txt");



R = 6378136.3;               
GM = 3986004.415e8;          
rhoave = 5500;                      
rp = R + height;



lonmin=min(upper(:,1));lonmax=max(upper(:,1));latmin=min(upper(:,2));latmax=max(upper(:,2));
interval=abs(upper(2,1)-upper(1,1));
nlon=(lonmax-lonmin)/interval+1;nlat=(latmax-latmin)/interval+1;
upperd=reshape(upper(:,3),nlon,nlat);lowerd=reshape(lower(:,3),nlon,nlat);rhod=reshape(density(:,3),nlon,nlat);
up=upperd';low=lowerd';rho=rhod';


Lat = linspace(latmax, latmin, nlat);
Lon = linspace(lonmin, lonmax, nlon);


CNM_total = zeros(nmax+1, nmax+1);
SNM_total = zeros(nmax+1, nmax+1);

h = waitbar(0, 'waiting...'); 

for i = 1:nlat
    for j = 1:nlon
        lat_c = Lat(i);
        lon_c = Lon(j);
        
        lat1 = lat_c - interval/2;
        lat2 = lat_c + interval/2;
        lon1 = lon_c - interval/2;
        lon2 = lon_c + interval/2;
        
        H1 = up(i, j);         
        H2 = low(i, j);          
        delta_rho = rho(i, j);  
        
		%Spherical harmonic expansion of a single tesseroid
        [Cnm, Snm] = SHA_single_tess(lon1, lon2, lat1, lat2, nmax, N_theta);
        
        Cnm1 = zeros(nmax+1, nmax+1);
        Snm1 = zeros(nmax+1, nmax+1);
        for n = 0:nmax
            factor = 3 * delta_rho * ((1 + H1/R)^(n+3) - (1 + H2/R)^(n+3))  / (rhoave * (2*n + 1) * (n + 3));
            Cnm1(n+1, :) = factor * Cnm(n+1, :);
            Snm1(n+1, :) = factor * Snm(n+1, :);
        end

        CNM = zeros(nmax+1, nmax+1);
        SNM = zeros(nmax+1, nmax+1);
        for n = 0:nmax
            factor1 = (n + 1) * ((1 / (1 + height/R))^(n + 2));
            CNM(n+1, :) = factor1 * Cnm1(n+1, :);
            SNM(n+1, :) = factor1 * Snm1(n+1, :);
        end
        CNM_total = CNM_total + CNM;
        SNM_total = SNM_total + SNM;

    end
       waitbar(i/nlat, h, sprintf('computation progress：%.1f%%', 100*i/nlat));
end
close(h);

nmnumber=0;
for i=1:nmax+1
    nmnumber=i+nmnumber;
end
result_cnm_snm=zeros(nmnumber,4);

row=1;
for n=0:nmax
    for m=0:n
        result_cnm_snm(row,1)=n;
        result_cnm_snm(row,2)=m;
        result_cnm_snm(row,3)=CNM_total(n+1,m+1);
        result_cnm_snm(row,4)=SNM_total(n+1,m+1);
        row=1+row;
    end
end

[f]=SHS_regular(result_cnm_snm,nmax,Lat,Lon);
attraction =(GM/(R^2))* f * 1e5;

grav_attraction = zeros(nlat*nlon, 3);
for i = 1:nlat
    for j = 1:nlon
        grav_attraction((i-1)*nlon+j,1) = Lon(j);
        grav_attraction((i-1)*nlon+j,2) = Lat(i);
        grav_attraction((i-1)*nlon+j,3) = attraction(i,j);
    end
end

imagesc(attraction)
colormap jet;colorbar;

mmin=min(grav_attraction(:,3));
mmax=max(grav_attraction(:,3));
mmean=mean(grav_attraction(:,3));
mstd=std(grav_attraction(:,3));












