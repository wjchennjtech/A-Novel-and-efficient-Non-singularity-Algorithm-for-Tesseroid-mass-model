function p = plm(l,m,th)

if nargin == 2
   th = m;
   m  = 0;
end
if min(size(l)) ~= 1,  error('Degree l must be vector (or scalar)'), end
if any(rem(l,1) ~= 0), error('Vector l contains non-integers.'), end
if max(size(m)) ~= 1,  error('Order m must be scalar.'), end
if rem(m,1) ~=0,       error('Order m must be integer.'), end


[lrow,lcol] = size(l);
[trow,tcol] = size(th);
lmax = max(l);

n    = length(th);				
t    = th(:)*pi/180;        
x    = cos(t);
y    = sin(t);
lvec = l(:)';					
lastcol = max(lmax-m+2,2);
ptmp = zeros(n,lastcol);

if m == 0
   fac = 1;
else
   mm  = 2*(1:m);
   fac = sqrt(2*prod((mm+1)./mm));
end
pm  = fac*sin(t).^m;
ptmp(:,1) = pm;       			

pold = zeros(n,1);
p    = pm;					

for l = m+1:lmax
   col   = l - m + 1;			
   root1 = sqrt( (2*l+1)*(2*l-1)/((l-m)*(l+m)) ) ;
   root2 = sqrt( (2*l+1)*(l+m-1)*(l-m-1) );
   root2 = root2 / sqrt( (2*l-3)*(l-m)*(l+m) );
   pnew  = root1 *x.*p - root2 *pold;	
   pold  = p;					
   p     = pnew;			
   ptmp(:,col) = p;
end


p     = zeros(n,length(lvec));	
lind  = find(lvec < m);		
pcol  = lvec - m + 1;	
pcol(lind) = lastcol*ones(size(lind));	
p     = ptmp(:,pcol);			

if max(size(lvec))==1  & min(size(th))==1 & (trow == 1), p = p'; end
if max(size(th))==1 & min(size(lvec))==1  & (lcol == 1), p = p'; end
