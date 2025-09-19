function [f]=SHS_regular(cnm_snm,nmax,Lat,Lon)

mmax=nmax;
nlat=size(Lat,2);
nlon=size(Lon,2);
f=zeros(nlat,nlon);

Cnm=zeros(nmax+1,nmax+1);
Snm=zeros(nmax+1,nmax+1);

nmnumber=0;
for i=1:nmax+1
    nmnumber=i+nmnumber;
end

for i=1: nmnumber

    Cnm(cnm_snm(i,1)+1,cnm_snm(i,2)+1)=cnm_snm(i,3);
    Snm(cnm_snm(i,1)+1,cnm_snm(i,2)+1)=cnm_snm(i,4);

end


for m=0:mmax
    
    n=0:nmax;
        f=plm(n,m,90-Lat)*Cnm(:,m+1)*cosd(m*Lon)...
            +plm(n,m,90-Lat)*Snm(:,m+1)*sind(m*Lon)...
            +f;

end

end
