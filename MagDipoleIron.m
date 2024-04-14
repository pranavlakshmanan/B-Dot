%Ur(mu-r) for SS-430FR = 1450 
% We take radius=2.75mm to 3.5mm=r
% Thereby number of turns=1232-670=N
%Length of the Torquer=35mm=L
%Formula returns Magnetic Dipole as function of N and I

%Formula is Mag_Dipole=((N*I)*(Mur-1))/(L*(1+Nd(Mur-1)))

function [Mag_dipole] = MagDipoleIron(N,I,L,r,Mur)

%Mur=1450;

Nd = Demag_factor(L,r);

Mag_dipole = ((N*I)*(Mur-1))/(L*(1+Nd*(Mur-1)));

end
