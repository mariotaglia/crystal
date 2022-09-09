clear all

fout = fopen('tita_phi.xy','w');
fsal = fopen('sphere_points.xyz','w');
fid = fopen('kap.txt'); %%% .TXT fle must be generated with pdb_to_coarse.m %%%
N = 400; %%% HOW MANY VECTORS TO ORIENT PROTEIN %%%
n = [1;0;0]; %%% VECTOR NORMAL TO THE SURFACE %%%
index = 0;

%%%%%%%% OBTAINING DE COORDINATES Of FROM DE .TXT FILE%%%%%%

tline = fgets(fid); %%% Read first line %%%
numtot = str2num(tline); %%% Number of total balls in the system %%%

%%% Next reading will be on 2nd line.
%%% Display form 2th line forward:

while ischar(tline); %%% while lines have characters = True %%%
    tline = fgets(fid); %%% Reading lines %%% 
    index = (index + 1); %%% Counter, will count how many balls %%%
    if le(index,numtot); %%% Will not take into account the end of file, (tline = -1) %%%
    vector(index).x = str2num(tline(1:8)); %%% Reading character 1 to 8 and seting as number %%%
    vector(index).y = str2num(tline(10:18));
    vector(index).z = str2num(tline(19:27));
    label(index).l = tline(28:29);
    label(index).n = str2num(tline(31:34));
    end
end


fclose(fid); %%% Closing opened file in fid %%%

    dipolex=0;
    dipoley=0;
    dipolez=0;

for index = 1:numtot
    %%% Defining position of the balls as vectors and calculating electric dipole, sup pH=7 %%%
    balls_vec(index).vec=[vector(index).x;vector(index).y;vector(index).z];
    if strcmp(label(index).l(1),'E') || strcmp(label(index).l(1),'D') || index==numtot
    dipolex=dipolex-1*vector(index).x;
    dipoley=dipoley-1*vector(index).y;
    dipolez=dipolez-1*vector(index).z;
    end
    if strcmp(label(index).l(1),'R') || strcmp(label(index).l(1),'H') || strcmp(label(index).l(1),'K')
    dipolex=dipolex+1*vector(index).x;
    dipoley=dipoley+1*vector(index).y;
    dipolez=dipolez+1*vector(index).z;
    end
end
dipolevec=[dipolex;dipoley;dipolez];
disp("PRIMER DIPOLO")
disp(dipolevec)

%%% FIRST ROTATION: rotation of the protein in order to aling dipole momento with the surface

p = [dipolex;dipoley;dipolez];
normp=norm(p)
p1 = [dipolex/normp;dipoley/normp;dipolez/normp];
v = cross(n,p1);
s = norm(v);
c = dot(n,p1);

v_mat = [0.0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0.0];

R = eye(3)+ v_mat + (v_mat^2)*((1-c)/(s^2));

for index=1:numtot
first_coord(index).vec=transpose(balls_vec(index).vec)*R;

end
    dipolex2=0;
    dipoley2=0;
    dipolez2=0;

for index = 1:numtot
    %%% Defining position of the balls as vectors and calculating electric dipole, sup pH=7 %%%
    if strcmp(label(index).l(1),'E') || strcmp(label(index).l(1),'D') || index==numtot
    dipolex2=dipolex2-1*first_coord(index).vec(1);
    dipoley2=dipoley2-1*first_coord(index).vec(2);
    dipolez2=dipolez2-1*first_coord(index).vec(3);
    end
    if strcmp(label(index).l(1),'R') || strcmp(label(index).l(1),'H') || strcmp(label(index).l(1),'K')
    dipolex2=dipolex2+1*first_coord(index).vec(1);
    dipoley2=dipoley2+1*first_coord(index).vec(2);
    dipolez2=dipolez2+1*first_coord(index).vec(3);
    end
end
dipolevec2=[dipolex2;dipoley2;dipolez2];
disp("Oriented dipole acording to x axis (1,0,0) is ")
disp(dipolevec2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OBTAINING REGULAR SPACIATED POINTS FROM A SPHERE %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Printing in the last position the (tita,phi) values for the oriented electric dipole
fprintf(fout, '% s % s \n', "  90.00000", " 0.00000")
%%%%%
radius=1;
a = 4*3.1415*radius^2/N;
d = a^0.5;
Mt = round(3.1415/d);
dt = 3.1415/Mt;
dp = a/dt;
k = 0;
for i=0:(Mt-1)
    tita = 3.1415*(i+0.5)/Mt;
    Mp = round(2*3.1415*sin(tita)/dp);
    for j = 0:(Mp-1);
       phi = 2*3.1415*j/Mp;
       x = radius*sin(tita)*cos(phi);
       y = radius*sin(tita)*sin(phi);
       z = radius*cos(tita);
    
    k=k+1;   %%% k points were generated from the sphere surface %%%
       coor.x(k) = x;
       coor.y(k) = y;
       coor.z(k) = z;
      

       fprintf(fsal, '% 4.5f % 4.5f % 4.5f \n', coor.x(k), coor.y(k), coor.z(k))       
       fprintf(fout, '% 4.5f % 4.5f \n', tita*180/3.1415, phi*180/3.1415)
    end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ROTATION MATRIX DEFINITION %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matnum = k

for j = 1: matnum

p = [coor.x(j); coor.y(j); coor.z(j)];

v = cross(n,p);
s = norm(v);
c = dot(n,p);

v_mat = [0.0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0.0];

R = eye(3)+ v_mat + (v_mat^2)*((1-c)/(s^2));

R_matrix{j}.mat = R; %%% Storing Mp rotation matrices %%%

end

%%% Writing some staff in files %%%
for i = 1: matnum

textfilename1 = ['balls_rotated_' num2str(i) '.txt'];
textfilename2 = ['balls_rotated_' num2str(i) '.xyz'];
fout = fopen(textfilename1, 'w');
fouts = fopen(textfilename2, 'w');
fprintf(fout, '%d \n', numtot);
fprintf(fouts, '%d \n', numtot);
fprintf(fouts, '%s \n', 'balls rotated');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Rotating the protein, using each R, and writing %%%

for index=1:numtot

balls_rot(index).vec=first_coord(index).vec*R_matrix{i}.mat;

fprintf(fout, '% 4.5f % 4.5f % 4.5f %s % d \n', balls_rot(index).vec, label(index).l, label(index).n );
fprintf(fouts, '%s %   4.5f %   4.5f %   4.5f \n','H    ', balls_rot(index).vec)

end

end

%%%%% Save the conformation with the oriented dipole to te surface%%%%%%
textfilename1 = ['balls_rotated_' num2str(0) '.txt'];
textfilename2 = ['balls_rotated_' num2str(0) '.xyz'];
fout2 = fopen(textfilename1, 'w');
fouts2 = fopen(textfilename2, 'w');
fprintf(fout2, '%d \n', numtot);
fprintf(fouts2, '%d \n', numtot);
fprintf(fouts2, '%s \n', 'balls rotated');
for index=1:numtot
balls_rot2(index).vec=first_coord(index).vec*eye(3);
fprintf(fout2, '% 4.5f % 4.5f % 4.5f %s % d \n', balls_rot2(index).vec, label(index).l, label(index).n );
fprintf(fouts2, '%s %   4.5f %   4.5f %   4.5f \n','H    ', balls_rot2(index).vec)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






