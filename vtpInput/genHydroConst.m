function [Kr,Kx]=genHydroConst(kr,kx,b,r)

    Kr=pi*(2*r-b)*kr/b; %no rho g: kr (opt) derived from Kr this way.
    Kx=kx.*pi.*(r-b).^2; %same here: kx derived from Lax, by dividing through with As
       
end