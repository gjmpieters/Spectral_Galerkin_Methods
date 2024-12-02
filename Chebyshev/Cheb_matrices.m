function [mA,mI,mQ,mP,mJ,mK,mH] = Cheb_matrices(N,H,t,beta)

f = pi;
cutoff = 0.000001;

%%% Gauss-Chebyshev points %%%%%%%%%
K = 2000;
zeta = cos((2*[1:K]-1)*pi/(2*(K)));
alph = pi/K;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% identity <W_p, V_k> %%%%%
mI = zeros(N+1,N+1);

for k = 0:N
    for p = 0:N
        
        if k == p+2
            mI(k+1,p+1) = -c(p)/4;
            
        elseif k == p
            mI(k+1,p+1) = ( 2 + 2*c(p) - c(p-1) )/4;
            
        elseif k == p-2
            mI(k+1,p+1) = -( 2 + c(p-2) )/4;
            
        elseif k == p-4
            mI(k+1,p+1) = 1/4;
            
        else
            mI(k+1,p+1) = 0;
        end
        
    end
end

mI = (f/2)*mI;

%%%%% matrix H <U_p, V_k> %%%%%
if exist('matrixH.mat', 'file') == 2
    load matrixH mHsav 
    mH   = mHsav(1:N+1,1:N+1);
else
    mH   = zeros(N+1,N+1);
    for k = 0:N
        for p = 0:N
            val = alph*sum(u(p,zeta).*v(k,zeta));
            mH(k+1,p+1) = val;
        end
    end
    mHsav = mH;
    save matrixH mHsav
end

%%%%% matrix Q <primprof*W_p, V_k> %%%%%
mQ   = zeros(N+1,N+1);
for k = 0:N
    for p = 0:N
        z = H*(zeta+1)/2;
        val = alph*sum(primprof(z,t,0).*w(p,zeta).*v(k,zeta));
        mQ(k+1,p+1) = val;
    end
end
mQsav = mQ;
save matrixQ mQsav

%%%%% matrix P <primprof*U_p, V_k> %%%%%
mP   = zeros(N+1,N+1);
for k = 0:N
    for p = 0:N
        z = H*(zeta+1)/2;
        val = alph*sum(primprof(z,t,0).*u(p,zeta).*v(k,zeta));
        mP(k+1,p+1) = val;
    end
end
mPsav = mP;
save matrixP mPsav

% %%%%% matrix K1 <(1/kappa)*W_p, V_k> %%%%%
% fid = fopen('matrixK1.mat','r');
% if fid == -1
%   K1   = zeros(N+1,N+1);
%   for k = 0:N
%     for p = 0:N
%         z = H*(zeta+1)/2;
%         val = alph*sum((1./permeability(z, beta)).*w(p,zeta).*v(k,zeta));
%         if abs(val) < cutoff
%             K1(k+1,p+1) = 0;
%         else
%             K1(k+1,p+1) = val;
%         end          
%     end
%   end
%   Hsav = H;
%   Nsav = N;
%   tsav = t;
%   betasav = beta;
%   K1sav = K1;
%   save matrixK1 K1sav Hsav Nsav tsav betasav
% else
%   fclose(fid);
%   load matrixK1 K1sav Hsav Nsav tsav betasav
%   if ( (H == Hsav)  && (Nsav >= N) && (t == tsav) && (beta == betasav))    
%     K1 = K1sav(1:N+1,1:N+1);
%   else    
%     K1   = zeros(N+1,N+1);
%     for k = 0:N
%         for p = 0:N
%           z = H*(zeta+1)/2;
%           val = alph*sum((1./permeability(z, beta)).*w(p,zeta).*v(k,zeta));
%           if abs(val) < cutoff
%               K1(k+1,p+1) = 0;
%           else
%               K1(k+1,p+1) = val;
%           end          
%         end
%     end
%     Hsav = H;
%     Nsav = N; 
%     tsav = t;
%     betasav = beta;
%     K1sav = K1;
%     save matrixK1 K1sav Hsav Nsav tsav betasav
%   end
% end

%%%%% matrix J <(1/kappa)*U_p, V_k> %%%%%
if exist('matrixJ.mat', 'file') == 2
    load matrixJ mJsav
    mJ = mJsav(1:N+1,1:N+1);
else
    mJ   = zeros(N+1,N+1);
    for k = 0:N
        for p = 0:N
            z = H*(zeta+1)/2;
            val = alph*sum((1./permeability(z, beta)).*u(p,zeta).*v(k,zeta));
            if abs(val) < cutoff
                mJ(k+1,p+1) = 0;
            else
                mJ(k+1,p+1) = val;
            end
        end
    end
    mJsav = mJ;
    save matrixJ mJsav
end


%%%%% matrix K <(1/kappa)*DU_p, DV_k> %%%%%
if exist('matrixK.mat', 'file') == 2
    load matrixK mKsav
    mK = mKsav(1:N+1,1:N+1);
else
    mK   = zeros(N+1,N+1);
    for k = 0:N
        for p = 0:N
            z = H*(zeta+1)/2;            
            val = alph*sum(-(1./permeability(z, beta)).*Du(p,zeta).*Dv(k,zeta));
            val = val + alph*sum(-(1./permeability(z, beta)).*(zeta./(1-zeta.^2)).*Du(p,zeta).*v(k,zeta));
            if abs(val) < cutoff
                mK(k+1,p+1) = 0;
            else
                mK(k+1,p+1) = (4/H^2)*val;
            end
        end
    end
    mKsav = mK;
    save matrixK mKsav
end

% %%%%% matrix KD2 <(1/kappa)*DW_p, DV_k> %%%%%
% 
% fid = fopen('matrixKD2.mat','r');
% if fid == -1
%     KD2   = zeros(N+1,N+1);
%     for k = 0:N
%         for p = 0:N
%             z = H*(zeta+1)/2;
%             val = alph*sum(-(1./permeability(z, beta)).*Dw(p,zeta).*Dv(k,zeta));
%             val = val + alph*sum(-(1./permeability(z, beta)).*(zeta./(1-zeta.^2)).*Dw(p,zeta).*v(k,zeta));            
%             if abs(val) < cutoff
%                 KD2(k+1,p+1) = 0;
%             else
%                 KD2(k+1,p+1) = (4/H^2)*val;
%             end
%         end
%     end
%     Hsav = H;
%     Nsav = N;
%     tsav = t;
%     betasav = beta;
%     KD2sav = KD2;
%     save matrixKD2 KD2sav Hsav Nsav tsav betasav
% else
%     fclose(fid);
%     load matrixKD2 KD2sav Hsav Nsav tsav betasav
%     if ( (H == Hsav)  && (Nsav >= N) && (t == tsav) && (beta == betasav))
%         KD2 = KD2sav(1:N+1,1:N+1);
%     else
%         KD2   = zeros(N+1,N+1);
%         for k = 0:N
%             for p = 0:N
%                 z = H*(zeta+1)/2;
%                 val = alph*sum(-(1./permeability(z, beta)).*Dw(p,zeta).*Dv(k,zeta));
%                 val = val + alph*sum(-(1./permeability(z, beta)).*(zeta./(1-zeta.^2)).*Dw(p,zeta).*v(k,zeta));
%                 if abs(val) < cutoff
%                     KD2(k+1,p+1) = 0;
%                 else
%                     KD2(k+1,p+1) = (4/H^2)*val;
%                 end
%             end
%         end
%         Hsav = H;
%         Nsav = N;       
%         tsav = t;
%         betasav = beta;
%         KD2sav = KD2;
%         save matrixKD2 KD2sav Hsav Nsav tsav betasav
%     end
% end

%%%%% matrix A <D2W_p, V_k> %%%%%
mA = zeros(N+1,N+1);
for k = 0:N
    for p = 0:N
        
        if k == p
            mA(k+1,p+1) = -p*( p + 3 ) - 2*c(p);
            
        elseif k == p-2
            mA(k+1,p+1) = p*( p - 3 ) + 2*c(p);
            
        else
            mA(k+1,p+1) = 0;
        end
        
    end
end

mA = (4/H^2)*(f/2)*mA;


% %%%%% matrix <DW_p, V_k> %%%%%
% AD1 = zeros(N+1,N+1);
% for k = 0:N
%     for p = 0:N
%         
%         if k == p+1
%             AD1(k+1,p+1) = ( -2*c(p) - p ) / 2;
%             
%         elseif k == p-1
%             AD1(k+1,p+1) = ( 5*p + 2 - c(p-1)*( p + 2 ) - c(p-2)*( d(p-3) + 1 )*p ) / 2;
%             
%         elseif k == p-3
%             AD1(k+1,p+1) = ( 2 - p ) / 2;
%             
%         else
%             AD1(k+1,p+1) = 0;
%         end
%         
%     end
% end
% 
% AD1 = (2/H)*(f/2)*AD1;
