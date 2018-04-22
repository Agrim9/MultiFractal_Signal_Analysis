function [beta, sx, sw] = EM_estimate(bet0, sx0, sw0, l, sig_m )
beta = bet0;
sx = sx0;
sw = sw0;
M = size(l,1)-2;
eps = 1;
while(eps>1e-5)
    % E step
    A = sw*sx*beta.^(-(1:M)')./(sw + sx*beta.^(-(1:M)'));
    Bw = (sw./(sw + sx*beta.^(-(1:M)'))).^2;
    Bx = (sx*beta.^(-(1:M)')./(sw + sx*beta.^(-(1:M)'))).^2;
  
    Cm = (1:M)'./sum((1:M)'.*l(2:M+1)) - ones(M,1)./sum(l(2:M+1));
    Sw = A + Bw.*sig_m;
    Sx = A + Bx.*sig_m;
    
    %M step
    coeff = Cm.*l(2:M+1).*Sx ; 
    rt = roots(flipud(coeff)); %roots takes arguments upside down
    beta_n = rt(imag(rt)==0 & real(rt)>0);
    sw_n = sum(l(2:M+1).*Sw)/sum(l(2:M+1)); %check for convergence
    sx_n = sum(l(2:M+1).*Sx.*beta_n.^((1:M)'))/sum(l(2:M+1));
    eps = abs(beta_n -beta)/beta ;
    beta = beta_n;
    sw = sw_n;
    sx = sx_n;
end
       
end
