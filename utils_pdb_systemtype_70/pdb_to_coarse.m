clear all
filename='test.pdb'
fid = fopen(filename);
fout = fopen('gfp.xyz','w');
fout2 = fopen('gfp.txt','w');
fout3 = fopen('c_terminal.xyz','w');
tline = fgetl(fid);
rog = 0;
indexatom = 0;
char temp;

while ischar(tline)

C = strsplit(tline,' ');

if (strcmp(C(1,1),'ATOM')) || (strcmp(C(1,1),'HETATM'))

indexatom = indexatom + 1;
fflush(stdout);

atom(indexatom).x = str2num(tline(27:38));
atom(indexatom).y = str2num(tline(39:46));
atom(indexatom).z = str2num(tline(47:54));
atom(indexatom).res = tline(18:20);
atom(indexatom).resn = str2num(tline(23:26));
atom(indexatom).type = tline(14:15);

end

tline = fgetl(fid);

end

fclose(fid);

lastres = 0;
indexres = 0;

for i = 1:indexatom    
    if (lastres ~=  atom(i).resn)
        indexres = indexres + 1;
        res(indexres).res = atom(i).res;
        lastres = atom(i).resn;
        res(indexres).firstatom = i; 
    end
end
res(indexres+1).firstatom = indexatom+1;

if ~(strcmp(atom(indexatom).type,'OX'))   %% add terminal oxygen if not present
%%%%%%%!!!!!!!!Adding C-Terminal Oxigen!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%

%%%Identifying positions of CA, C, O atoms from C-Term of infile%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% 	(erre3)O   O(erre4, to add)	%%%%%%%%%
%%%%%%% 	       \\ //		        %%%%%%%%%
%%%%%%% 		 C(erre2)	        %%%%%%%%%
%%%%%%% 		 |			%%%%%%%%%
%%%%%%% 		 CA(erre1)		%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
erre1.x = atom(res(indexres).firstatom+1).x
erre1.y = atom(res(indexres).firstatom+1).y
erre1.z = atom(res(indexres).firstatom+1).z
%
erre2.x = atom(res(indexres).firstatom+2).x
erre2.y = atom(res(indexres).firstatom+2).y
erre2.z = atom(res(indexres).firstatom+2).z
%
erre3.x = atom(res(indexres).firstatom+3).x
erre3.y = atom(res(indexres).firstatom+3).y
erre3.z = atom(res(indexres).firstatom+3).z

%%%Difference, orthogonal vector, modules and angle%%%%%%
erreA.x = (erre3.x-erre2.x)
erreA.y = (erre3.y-erre2.y)
erreA.z = (erre3.z-erre2.z)
%
erreC.x = (erre1.x-erre2.x)
erreC.y = (erre1.y-erre2.y)
erreC.z = (erre1.z-erre2.z)
%
erreP.x = (erreA.y*erreC.z-erreA.z*erreC.y)
erreP.y = (-erreA.x*erreC.z+erreA.z*erreC.x)
erreP.z = (erreA.x*erreC.y-erreA.y*erreC.x)
%
norma1 = sqrt((erre2.x-erre3.x)^2+(erre2.y-erre3.y)^2+(erre2.z-erre3.z)^2)
%
norma2 = sqrt((erre2.x-erre1.x)^2+(erre2.y-erre1.y)^2+(erre2.z-erre1.z)^2)
%
angulo2 = acosd((erreA.x*erreC.x+erreA.y*erreC.y+erreA.z*erreC.z)/(norma1*norma2))
%
%%%Angle condition between O-C-O(added)%%%
A = cosd(123)*norma1*norma1+erre2.x*erreA.x+erre2.y*erreA.y+erre2.z*erreA.z
%%%Angle condition between CA-C-O(added)%%%
B = cosd(angulo2)*norma1*norma2+erre2.x*erreC.x+erre2.y*erreC.y+erre2.z*erreC.z
%
C = erre2.x*erreP.x+erre2.y*erreP.y+erre2.z*erreP.z

%%%%%%Coordinates of the added oxygen atom%%%%%%%%%%%
oxy.z = (((B*erreP.x-C*erreC.x)*(erreC.y*erreA.x-erreC.x*erreA.y)-(B*erreA.x-A*erreC.x)*(erreC.y*erreP.x-erreC.x*erreP.y))/((erreP.x*erreC.z-erreC.x*erreP.z)*(erreC.y*erreA.x-erreC.x*erreA.y)-(erreA.x*erreC.z-erreC.x*erreA.z)*(erreC.y*erreP.x-erreC.x*erreP.y)))
%
oxy.y = ((erreA.x*B - erreC.x*A - (oxy.z)*(erreA.x*erreC.z-erreC.x*erreA.z))/(erreC.y*erreA.x-erreC.x*erreA.y))
%
oxy.x = ((A - erreA.z*oxy.z-erreA.y*oxy.y)/erreA.x)
%%%%%%Printing C-Terminal to .xyz file%%%%%%%%%%%%%%%
%atom(res(indexres).firstatom).x = oxy.x;
%atom(res(indexres).firstatom).y = oxy.y;
%atom(res(indexres).firstatom).z = oxy.z;
%atom(res(indexres).firstatom).res = atom(res(indexres).firstatom);
%atom(res(indexres).firstatom).resn = atom(res(indexres).firstatom);
%atom(res(indexres).firstatom).type = atom(res(indexres).firstatom+3);

fprintf(fout3, '%d \n', (res(indexres+1).firstatom)+1-(res(indexres).firstatom));
fprintf(fout3, '%s \n', 'C-Terminal coordinates');
for i = res(indexres).firstatom:(res(indexres+1).firstatom-1)
fprintf(fout3, '%s %4.3f %4.3f %4.3f \n', atom(i).type, atom(i).x, atom(i).y, atom(i).z);
end
fprintf(fout3, '%s %4.3f %4.3f %4.3f \n','O ' , oxy.x, oxy.y, oxy.z);

end %% add terminal oxygen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:indexatom
    switch  atom(i).type
        case 'N1'
            atom(i).mass = 14.0;
        case 'N2'
            atom(i).mass = 14.0;
        case 'N3'
            atom(i).mass = 14.0;
        case 'NZ'
            atom(i).mass = 14.0;
        case 'NE'
            atom(i).mass = 14.0;
        case 'NH'
            atom(i).mass = 14.0;
        case 'ND'
            atom(i).mass = 14.0;
        case 'OG'
            atom(i).mass = 16.0;
        case 'OE'
            atom(i).mass = 16.0;
        case 'OD'
            atom(i).mass = 16.0;
        case 'OX'
            atom(i).mass = 16.0;
        case 'OH'
            atom(i).mass = 16.0;
        case 'SG'
            atom(i).mass = 32.0;
        case 'SD'
            atom(i).mass = 32.0;
        case 'C1'
            atom(i).mass = 12.0;
        case 'C2'
            atom(i).mass = 12.0;
        case 'C3'
            atom(i).mass = 12.0;
        case 'CH'
            atom(i).mass = 12.0;    
        case 'CB'
            atom(i).mass = 12.0;
        case 'CG'
            atom(i).mass = 12.0;
        case 'CD'
            atom(i).mass = 12.0;
        case 'CZ'
            atom(i).mass = 12.0;
        case 'CE'
            atom(i).mass = 12.0;
        case 'N '
            atom(i).mass = 14.0;
        case 'CA'
            atom(i).mass = 12.0;
        case 'C '
            atom(i).mass = 12.0;
        case 'O '
            atom(i).mass = 16.0;
        case 'O2'
            atom(i).mass = 16.0;
        case 'O3'
            atom(i).mass = 16.0;
        otherwise
            disp('################STOPS: Not N,O,C or S Atom Found in PDB file!!!!!')
            disp('Atom Name:########################')
            atom(i).type
            i
            stop
    end
end
cont=0;
for ii = 1:indexres
if strcmp(res(ii).res,'GLY')
cont=cont+1;
end
end


c=0;
for ii = 1:indexres;
 if strcmp(res(ii).res,'CYS') 
   i=(res(ii).firstatom+5);
   c=c+1;
   azu(c).x=atom(i).x;
   azu(c).y=atom(i).y;
   azu(c).z=atom(i).z;
   aa(c).aa=ii;
   fprintf('CYS found! at PDB position:')
   disp(ii)
 end
end


for j=1:c;
  for l=1:j;
      dist=sqrt((azu(j).x-azu(l).x)**2+(azu(j).y-azu(l).y)**2+(azu(j).z-azu(l).z)**2);
      if (dist > 0 && dist <= 2.3);
      disp('WARNING!!! Potential bond between residues')
      disp(aa(j).aa)
      fprintf('# AND # \n')
      disp(aa(l).aa);
      fprintf('with DIST(SG-SG)=%i.\n',dist); 
      jj=aa(j).aa;
      ll=aa(l).aa;
      res(jj).res = 'CYX';
      res(ll).res = 'CYX';
      end 
  end
end

 
fprintf(fout, '%d \n', (2*indexres-cont+2));  % 2 is for CRO
fprintf(fout, '%s \n', filename);
fprintf(fout2, '%d \n', (2*indexres-cont+2)); % 2 is for CRO

totalmassb=54.0;
cont=0;

for ii = 1:indexres;
  com.x = 0.0;
  com.y = 0.0;
  com.z = 0.0;
  totalmass = 0.0;
  com.xb = 0.0;
  com.yb = 0.0;
  com.zb = 0.0;
  

  compare=strcmp(res(ii).res,'CRO');
    if compare == 0
       for i = res(ii).firstatom:(res(ii).firstatom+3);
         com.xb = com.xb +  atom(i).x*atom(i).mass;
         com.yb = com.yb +  atom(i).y*atom(i).mass;
         com.zb = com.zb +  atom(i).z*atom(i).mass;
       end

       if ii == indexres && ~strcmp(atom(indexatom).type,'OX'); %% ultimo residuo
       com.xb = com.xb + oxy.x*16.0;
       com.yb = com.yb + oxy.y*16.0;
       com.zb = com.zb + oxy.z*16.0;
       totalmassb = 70;
       end

       if ii == indexres && strcmp(atom(indexatom).type,'OX'); %% ultimo residuo
       com.xb = com.xb + atom(indexatom).x*atom(indexatom).mass;
       com.yb = com.yb + atom(indexatom).y*atom(indexatom).mass;
       com.zb = com.zb + atom(indexatom).z*atom(indexatom).mass;
       totalmassb = 70;
       end

       for i = res(ii).firstatom+4:(res(ii+1).firstatom-1);
         com.x = com.x +  atom(i).x*atom(i).mass;
         com.y = com.y +  atom(i).y*atom(i).mass;
         com.z = com.z +  atom(i).z*atom(i).mass;
         totalmass = totalmass + atom(i).mass;
       end
    com.xb = com.xb / totalmassb;
    com.yb = com.yb / totalmassb;
    com.zb = com.zb / totalmassb;
    res(ii).comxb = com.xb;
    res(ii).comyb = com.yb;
    res(ii).comzb = com.zb;
    if strcmp(res(ii).res,'GLY')
    com.x = com.xb;
    com.y = com.yb;
    com.z = com.zb;
    fprintf(fout, '%s %4.2f %4.2f %4.2f 0.0 \n', res(ii).res, com.xb, com.yb, com.zb);
    else
    com.x = com.x / totalmass;
    com.y = com.y / totalmass;
    com.z = com.z / totalmass;
    res(ii).comx = com.x;
    res(ii).comy = com.y;
    res(ii).comz = com.z;

    fprintf(fout, '%s %4.2f %4.2f %4.2f 0.0 \n', 'BBN', com.xb, com.yb, com.zb);
    fprintf(fout, '%s %4.2f %4.2f %4.2f 0.0 \n', res(ii).res, com.x, com.y, com.z);
    end


%%%%%%%%%%%%%%%%%%%%%%GRAINING THE CROMOPHORE OF GFP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
         com.x1=0;
         com.y1=0;
         com.z1=0;
         com.x2=0;
         com.y2=0;
         com.z2=0;
         com.x3=0;
         com.y3=0;
         com.z3=0;
         com.x4=0;
         com.y4=0;
         com.z4=0;
         totalmass1=0;
         totalmass2=0;
         totalmass3=0;
         totalmass4=0;
       for i = res(ii).firstatom:(res(ii).firstatom+4);
         com.x1 = com.x1 +  atom(i).x*atom(i).mass;
         com.y1 = com.y1 +  atom(i).y*atom(i).mass;
         com.z1 = com.z1 +  atom(i).z*atom(i).mass;
         totalmass1= totalmass1 + atom(i).mass;
       end
         com.x1 = com.x1/totalmass1;
         com.y1 = com.y1/totalmass1;
         com.z1 = com.z1/totalmass1;
         fprintf(fout, '%s %4.2f %4.2f %4.2f 0.0 \n', 'BBN', com.x1, com.y1, com.z1);
       for i = (res(ii).firstatom+5):(res(ii).firstatom+10);
         com.x2 = com.x2 +  atom(i).x*atom(i).mass;
         com.y2 = com.y2 +  atom(i).y*atom(i).mass;
         com.z2 = com.z2 +  atom(i).z*atom(i).mass;
         totalmass2= totalmass2 + atom(i).mass;
       end
         com.x2 = com.x2/totalmass2;
         com.y2 = com.y2/totalmass2;
         com.z2 = com.z2/totalmass2;
         fprintf(fout, '%s %4.2f %4.2f %4.2f 0.0 \n', 'CRO', com.x2, com.y2, com.z2);
       for i = (res(ii).firstatom+11):(res(ii).firstatom+13);
         com.x3 = com.x3 +  atom(i).x*atom(i).mass;
         com.y3 = com.y3 +  atom(i).y*atom(i).mass;
         com.z3 = com.z3 +  atom(i).z*atom(i).mass;
         totalmass3= totalmass3 + atom(i).mass;
       end
         com.x3 = com.x3/totalmass3;
         com.y3 = com.y3/totalmass3;
         com.z3 = com.z3/totalmass3; 
         fprintf(fout, '%s %4.2f %4.2f %4.2f 0.0 \n', 'BBN', com.x3, com.y3, com.z3);
       for i = res(ii).firstatom+14:(res(ii+1).firstatom-1);
         com.x4 = com.x4 +  atom(i).x*atom(i).mass;
         com.y4 = com.y4 +  atom(i).y*atom(i).mass;
         com.z4 = com.z4 +  atom(i).z*atom(i).mass;
         totalmass4 = totalmass4 + atom(i).mass;
       end
       com.x4 = com.x4/totalmass4;
       com.y4 = com.y4/totalmass4;
       com.z4 = com.z4/totalmass4;
       fprintf(fout, '%s %4.2f %4.2f %4.2f 0.0 \n', 'TYR', com.x4, com.y4, com.z4);
       res(ii).comx1 = com.x1;
       res(ii).comy1 = com.y1;
       res(ii).comz1 = com.z1;
       res(ii).comx2 = com.x2;
       res(ii).comy2 = com.y2;
       res(ii).comz2 = com.z2;
       res(ii).comx3 = com.x3;
       res(ii).comy3 = com.y3;
       res(ii).comz3 = com.z3;
       res(ii).comx4 = com.x4;
       res(ii).comy4 = com.y4;
       res(ii).comz4 = com.z4;
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

com.x = 0.0;
com.y = 0.0;
com.z = 0.0;
totalmass = 0.0;

%%%Protein center of mass coordinates%%%
for i = 1:indexatom;
    com.x = com.x +  atom(i).x*atom(i).mass;
    com.y = com.y +  atom(i).y*atom(i).mass;
    com.z = com.z +  atom(i).z*atom(i).mass;
    totalmass = totalmass + atom(i).mass;
end

com.x = com.x / totalmass;
com.y = com.y / totalmass;
com.z = com.z / totalmass;


%traslandando la proteina al origen!%%
for i = 1:indexatom;
       vect(i) = sqrt((atom(i).x-com.x)^2 + (atom(i).y-com.y)^2 + (atom(i).z-com.z)^2);
       rog = rog + atom(i).mass*(vect(i)^2);
end
rog = rog/totalmass;
rog = sqrt(rog);

hpho = 0;
hphi = 0;
pos = 0;
neg = 0;
cys = 0;
his = 0;
cro = 0;

for i = 1:indexres;
   switch(res(i).res);
    case 'ALA'
           hpho = hpho + 1;
%          res(i).q = 3;
           res(i).l = 'A';

    case 'ILE'
           hpho = hpho + 1;
%          res(i).q = 4;
           res(i).l = 'I';

    case 'LEU'
           hpho = hpho + 1;
%          res(i).q = 4;
           res(i).l = 'L';

    case 'PHE'
           hpho = hpho + 1;
%          res(i).q = 4;
           res(i).l = 'F';

    case 'TRP'
           hpho = hpho + 1;
%          res(i).q = 4;
           res(i).l = 'W';

    case 'TYR'
           hpho = hpho + 1;
%          res(i).q = 3;
           res(i).l = 'Y';

    case 'LYS'
           pos = pos + 1;
%          res(i).q = 5;
           res(i).l = 'K';

    case 'ARG'
           pos = pos + 1;
%          res(i).q = 5 ;
           res(i).l = 'R';

    case 'ASN'
           hphi = hphi + 1;
%          res(i).q = 3;
           res(i).l = 'N';

    case 'GLN'
           hphi = hphi + 1;
%          res(i).q = 3;
           res(i).l = 'Q';

    case 'GLY'
           hphi = hphi + 1;
%          res(i).q = 2;
           res(i).l = 'G';

    case 'MET'
           hphi = hphi + 1;
%          res(i).q = 3;
           res(i).l = 'M';

    case 'PRO'
           hphi = hphi + 1;
%          res(i).q = 3;
           res(i).l = 'P';

    case 'SER'
           hphi = hphi + 1;
%          res(i).q = 2;
           res(i).l = 'S';

    case 'THR'
           hphi = hphi + 1;
%          res(i).q = 2;
           res(i).l = 'T';

    case 'VAL'
           hphi = hphi + 1;
%          res(i).q = 4;
           res(i).l = 'V';

    case 'ASP'
           neg = neg + 1;
%          res(i).q = 7 ;
           res(i).l = 'D';

    case 'GLU'
           neg = neg + 1;
%          res(i).q = 6;
           res(i).l = 'E';

    case 'CYS'
           cys = cys + 1;
%          res(i).q = 1;
           res(i).l = 'C';

    case 'CYX'
           cys = cys + 1;
%          res(i).q = 1;
           res(i).l = 'X';

    case 'HIS'
           his = his + 1;
%          res(i).q = 8;
           res(i).l = 'H';

    case 'CRO'
           cro = cro + 1;
%          res(i).q = 8;
           res(i).l = 'Z';

    case 'HOH'


   otherwise  
           
           res(i).res
           stop
   end
   if ~(strcmp(res(i).res,'CRO'))
      fprintf(fout2, '% 4.5f % 4.5f % 4.5f % c % 4.0f \n', (res(i).comxb-com.x)/10.0, (res(i).comyb-com.y)/10.0, (res(i).comzb-com.z)/10.0, 'B', i);
      if ~(strcmp(res(i).res,'GLY')) 
         fprintf(fout2, '% 4.5f % 4.5f % 4.5f % c % 4.0f \n', (res(i).comx-com.x)/10.0, (res(i).comy-com.y)/10.0, (res(i).comz-com.z)/10.0, res(i).l, i);
      end
   else
      fprintf(fout2, '% 4.5f % 4.5f % 4.5f % c % 4.0f \n', (res(i).comx1-com.x)/10.0, (res(i).comy1-com.y)/10.0, (res(i).comz1-com.z)/10.0, 'B', i);
      fprintf(fout2, '% 4.5f % 4.5f % 4.5f % c % 4.0f \n', (res(i).comx2-com.x)/10.0, (res(i).comy2-com.y)/10.0, (res(i).comz2-com.z)/10.0, res(i).l, i);
      fprintf(fout2, '% 4.5f % 4.5f % 4.5f % c % 4.0f \n', (res(i).comx3-com.x)/10.0, (res(i).comy3-com.y)/10.0, (res(i).comz3-com.z)/10.0, 'B', i);
      fprintf(fout2, '% 4.5f % 4.5f % 4.5f % c % 4.0f \n', (res(i).comx4-com.x)/10.0, (res(i).comy4-com.y)/10.0, (res(i).comz4-com.z)/10.0, 'Y', i);
   end
end
