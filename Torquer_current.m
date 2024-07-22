
%Ur(mu-r) for SS-430FR = 1450 
% We take radius=2.75mm to 3.5mm=r
% Thereby number of turns=1232-670=N
%Length of the Torquer=35mm=L
%Formula returns Current required as a function of Magnetic Dipole and
%turns

%Formula is
%Torquer_current=(M*(1+(Mur-1)*Nd))/(pi*r^2*N*(1+(Mur-1)*Nd)+(Mur-1))





function [Torquer_current] = Torquer_current(M,N,L,r,Mur)

%Mur=1450;
piVal = pi;
%Nd = Demag_factor(L,r);


N_d = 4*(log(L/r) - 1) / ((L/r)^2 - 4*log(L/r));

Torquer_current = (M*(1+(Mur-1)*Nd))/( piVal*r^2*N* (((Mur-1)*N_d)+Mur));

end
