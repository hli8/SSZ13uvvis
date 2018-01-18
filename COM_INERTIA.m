%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to calculate rotational axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this one needs to be a separate function

function [out1,out2,out3,out4] = COM_INERTIA(coord,masses)
realcoord=coord;%.*(10^-10);      % A to m
realmass=masses;%.*(1.661e-27);   % amu to kg

% Center of mass and then rescale xyz
top=zeros(1,3);
for i=1:size(coord,1)
    top=top+realcoord(i,:).*realmass(i);
end
bot=sum(realmass);
Rcom=top./bot; 

% Rescaled cartesian coord
for i=1:size(coord,1)
    trans(i,:)=realcoord(i,:)-Rcom;
end

% Mmt of Inertia tensor
II=zeros(3,3);d1=0;d2=0;d3=0;
for i=1:size(coord,1)
    II(1,1)=II(1,1)+realmass(i)*(trans(i,2)^2+trans(i,3)^2);   % Ixx
    II(2,2)=II(2,2)+realmass(i)*(trans(i,1)^2+trans(i,3)^2);   % Iyy
    II(3,3)=II(3,3)+realmass(i)*(trans(i,1)^2+trans(i,2)^2);   % Izz
    
    d1=d1+realmass(i)*trans(i,1)*trans(i,2);    % Ixy
    II(1,2)=-d1;  II(2,1)=-d1;
    d2=d2+realmass(i)*trans(i,1)*trans(i,3);    % Ixz
    II(1,3)=-d2;  II(3,1)=-d2;
    d3=d3+realmass(i)*trans(i,2)*trans(i,3);    % Iyz
    II(2,3)=-d3;  II(3,2)=-d3;
end
[V D]=eig(II);

block1=[];
for jj=1:size(coord,1)
    block2=V'*trans(jj,:)';
    block1=vertcat(block1,block2');
end


out1=Rcom;
out3=block1;   % after trans+rot
out2=V;        % rot axes
out4=D;
