function [dx,dy,peak]= GaussPeakFit (cormap)

% (pseudo-)Gaussian peak fit: p(x,y) = exp(a + b*x + c*y + d*x*x + e*x*y + f*y*y)

QMAT= [ [ 1 -1 -1  1  1  1]; ...
        [ 1  0 -1  0  0  1]; ...
        [ 1  1 -1  1 -1  1]; ...
        [ 1 -1  0  1  0  0]; ...
        [ 1  0  0  0  0  0]; ...
        [ 1  1  0  1  0  0]; ...
        [ 1 -1  1  1 -1  1]; ...
        [ 1  0  1  0  0  1]; ...
        [ 1  1  1  1  1  1] ];

% get integer peak position
[maxval,ind]= max(cormap(:));
[jmax,imax]= ind2sub(size(cormap),ind);

dx= NaN;
dy= NaN;
peak= NaN;

[m,n]= size(cormap);

% cannot sample correlation neighbours if peak is on map edge ...
if imax==1 || jmax==1 || imax==n || jmax==m
    return;
end

ind= 1;
for j=-1:1
    for i=-1:1
        rhs(ind)= log(max(cormap(jmax+j,imax+i),1.e-12));
        ind= ind + 1;
    end
end
if var(rhs) == 0
    % right hand side was clipped everywhere: no good
    return
end

try
    % solve normal equations for least squares problem
    coeffs= QMAT \ rhs'; 
    % unpack solution vector; find peak position from zero derivative
    mmat= [ [2*coeffs(4) coeffs(5)]; [coeffs(5) 2*coeffs(6)]];
    mrhs= [-coeffs(2); -coeffs(3)];
    qvec= (mmat \ mrhs)';
    if norm(qvec) > sqrt(2)
        % interpolated displacement is too large: no good
        return
    end
    dx= qvec(1);
    dy= qvec(2);
catch
end

peak = exp(coeffs(1) + coeffs(2)*dx + coeffs(3)*dy + ...
       coeffs(4)*dx*dx + coeffs(5)*dx*dy + coeffs(6)*dy*dy);

return

