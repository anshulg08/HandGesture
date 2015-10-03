function B=bound2im(b,M,N,x0,y0)
[np,nc]=size(b);
if np<nc
    b=b';
    [np,nc]=size(b);
end
x=round(b(:,1));
y=round(b(:,2));
x=x-min(x)+1;
y=y-min(y)+1;
B=false(max(x),max(y));
C=max(x)-min(y)+1;
D=max(y)-min(y)+1;
if nargin==1
elseif nargin==3
    if C>M | D>N
    error('the boundry is outside M by N region')
    end
    B=false(M,N);
    NR=round((M-C)/2);
    NC=round((N-D)/2);
    x=x+NR;
    y=y+NC;
elseif nagin==5
    if x0<0 n| y0<0
        error('x0 and y0 must be +ve integers')
    end
    x=x+round(x0)-1;
    y=y+round(y0)-1;
    C=C+x0-1;
    D=D+y0-1;
    if C>M | D>N
    error('the shifted boundry is outside the M by N region')
end
B=false(M,N);
else 
    error('incorrect number of inputs')
end
B(sub2ind(size(B),x,y))=true;