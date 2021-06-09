%% HI-PATCH script.

% Authors: Hector Sanchez-Moran, James S. Weltz, Daniel K. Schwartz, Joel L. Kaar.
% Department of Chemical & Biological Engineering.
% University of Colorado Boulder

% Determination of hydrophobic patches (HPs) of known
% crystal structures based on solvation free energies of solvent exposed
% atoms.

% Cite: please stay tuned for upcoming publication.

% KEEP THIS HEADING IN THE CODE.

%%
close all
clear all
clc 

Excel_booklet='Enzymes.xlsx'; % Select Excel booklet where protein information is located.
tab='tab'; % Select tab name
excel_import = readtable(Excel_booklet,'ReadVariableNames',false,'Sheet',tab); % Separate numbers and text into two different arrays

% Arrange matrices
% Separate the information from the PDB file into a numerical array and a % character array. 
numPDB(:,1) = [1:1:height(excel_import)]'; % copy the (1) atom index 
numPDB(:,2:5) = excel_import{:,4:7}; % copy the (2) residue index, (3-5) XYZ coordinates 
numPDB(:,6) = excel_import{:,9}; % copy the (6) surface exposure
numPDB(:,7:12)=zeros(length(numPDB),6); % additional blank columns (7-12) for the subsequent analysis.
numPDB(~any(~isnan(numPDB)),:)=[]; % remove any NaNs
txtPDB(:,1) = excel_import{:,2}; % copy the (1) atom name
txtPDB(:,2) = excel_import{:,3}; % copy the (2) residue name
txtPDB(:,3) = excel_import{:,8}; % copy the (7) atom type

 
%% Van der Waals radii array in Angstroms. 
% Obtained at: Gerstein, M., Tsai, J., & Levitt, M. (1995). The Volume of 
% Atoms on the Protein Surface: Calculated from Simulation, using Voronoi
% Polyhedra. In J. Mol. Bl;ol (Vol. 249). 

rVdW=zeros(20,14);
rVdW(:,1)=1.65; % amide nitrogens backbone
rVdW(:,2)=1.87; % alpha carbons
rVdW(:,3)=1.76; % carbonyl carbons backbone
rVdW(:,4)=1.4;  % carbonyl oxygens backbone
rVdW(:,5)=1.87; % beta carbons
rVdW(1,5)=0;
rVdW(:,6)=1.87; % gamma carbons
rVdW(1:2,6)=0;
rVdW(3,6)=1.85; % cysteine sulfur
rVdW(4:5,6)=1.4; % serine & threonine oxygens
rVdW(11:12,6)=1.76; % carboxylic and amide carbon
rVdW(:,7)=1.87; % carbons
rVdW(1:4,7)=0;
rVdW(10,7)=1.85; % methionine 	
rVdW(11:12,7)=1.4; % aspargine & aspartic acid oxygens
rVdW(13:14,7)=1.76; % carboxylic and amide carbon
rVdW(16,7)=1.65; % histidine 1st nitrogen
rVdW(:,8)=1.87; % carbons
rVdW(1:7,8)=0;
rVdW(11,8)=1.65; % asparagine nitrogen
rVdW(12:14,8)=1.4; % aspartic glutamine glutamic oxygens
rVdW(17,8)=1.65; % arginine nitrogen 1
rVdW(:,9)=1.87; % carbons
rVdW(1:12,9)=0; 
rVdW(13,9)=1.4; % glutamic oxygen 2
rVdW(14,9)=1.65; % glutamine nitrogen
rVdW(15,9)=1.5; % lysine nitrogen
rVdW(16:17,10)=1.65; % histidine and arginine nitrogens
rVdW(18:20,10)=1.87; % carbons
rVdW(17,11)=1.65; % arginine nitrogen 3
rVdW(18:20,11)=1.87; % carbons
rVdW(19,12)=1.4; % tyrosine oxygen

%aromatic carbons
rVdW(18:19,6:11)=1.76;
rVdW(20,6:8)=1.76;
rVdW(20,10:14)=1.76;

%% Free energy array in kcal/mol

% Obtained in:
% Stouten, P. F. W., Frömmel, C., Nakamura, H., & Sander, C. (1993). An 
% Effective Solvation Term Based on Atomic Occupancies for Use in Protein 
% Simulations. Molecular Simulation, 10(2-6), 97?120. doi:10.1080/08927029308022161 

% &

% Lazaridis, T., & Karplus, M. (1999). Effective energy function for proteins
% in solution. Proteins: Structure, Function, and Bioinformatics, 35(2), 133
%152. https://doi.org/10.1002/(SICI)1097-0134(19990501)35:2<133::AID-PROT1>3.0.CO;2-N

% First, lay out common values and exceptions
dG=zeros(20,14);
dG(:,1)=-8.12; % amide N
dG(6,1)=-0.77; % proline amide N
dG(:,2)=0.53; % alpha carbons
dG(1,2)=1.3; % glycine alpha carbon
dG(:,3)=0.78; % carbonyl carbons
dG(:,4)=-5.07; % carbonyl oxygens
dG(2:20,5)=1.3; % CH2 beta carbons
dG(2,5)=2.28; % alanine beta carbon
dG(5,5)=0.53; % threonine beta carbon
dG(7:8,5)=0.53; % valine isoleucine beta carbon

% Now, residue by residue
dG(3,6)=-1.92; % cysteine thiol

dG(4,6)=-5.92; % serine hydroxyl

dG(5,6)=-5.92; % threonine hydroxyl
dG(5,7)=2.28; % threonine gamma C

dG(6,6:7)=1.3; % proline gamma delta C

dG(7,6:7)=2.28; % valine gamma C's

dG(8,6)=1.3; % isoleucine gamma C 1
dG(8,7:8)=2.28; % isoleucine gamma 2 & delta C

dG(9,6)=0.53; % leucine gamma C
dG(9,7:8)=2.28; % leucine delta Cs

dG(10,6)=1.3; % methionine gamma C
dG(10,7)=-3.32; % methionine S
dG(10,8)=2.28; % methionine eps C

dG(11,6)=-0.62; % asparagine gamma C
dG(11,7)=-7.02; % asparagine amide N
dG(11,8)=-5.07; % asparagine amide O

dG(12,6)=-0.62; % aspartic acid gamma C
dG(12,7:8)=-10; % aspartic acid O's

dG(13,6)=1.3; % glutamic acid gamma C
dG(13,7)=-0.62; % glutamic acid delta C
dG(13,8:9)=-10; % glutamic acid O's

dG(14,6)=1.3; % glutamine gamma C
dG(14,7)=-0.62; % glutamine delta C
dG(14,8)=-7.02; % glutamine amide N
dG(14,9)=-5.07; % glutamine amide O

dG(15,6:8)=1.3; % lysine C's
dG(15,9)=-20; % lysine NH3

dG(16,6)=-0.62; % histidine gamma C
dG(16,7)=0.86; % histidine delta C
dG(16,8)=-3.22; % histidine delta N
dG(16,9)=0.86; % histidine eps C
dG(16,10)=-0.77; % histidine eps N

dG(17,6:7)=1.3; % arginine gamma & delta C's
dG(17,8)=-0.77; % arginine eps N
dG(17,9)=-0.62; % arginine C
dG(17,10:11)=-10; % arginine guanidinium N's

dG(18:19,6:11)=0.86; % phenylalanine & tyrosine aromatic C's
dG(19,12)=-5.92; % tyrosine hydroxyl

dG(20,6:14)=0.86; % tryptophan aromatic C's
dG(20,11)=-0.77; % tryptophan N
dG(:,:)=dG(:,:)*4.184; % Convert from kcal/mol to kJ/mol

%% Binary array describing aromatic and aliphatic nature of each atom.
% 0 for aliphatic, 1 for aromatic.

arom_aliph(1,:)= [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(2,:)= [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(3,:)= [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(4,:)= [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(5,:)= [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(6,:)= [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(7,:)= [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(8,:)= [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(9,:)= [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(10,:)=[0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(11,:)=[0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(12,:)=[0 0 0 0 0 0 1 1 0 0 0 0 0 0];
arom_aliph(13,:)=[0 0 0 0 0 0 0 1 1 0 0 0 0 0];
arom_aliph(14,:)=[0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(15,:)=[0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arom_aliph(16,:)=[0 0 0 0 0 1 1 1 1 1 0 0 0 0];
arom_aliph(17,:)=[0 0 0 0 0 0 0 0 0 1 1 0 0 0];
arom_aliph(18,:)=[0 0 0 0 0 1 1 1 1 1 1 0 0 0];
arom_aliph(19,:)=[0 0 0 0 0 1 1 1 1 1 1 0 0 0];
arom_aliph(20,:)=[0 0 0 0 0 1 1 1 1 1 1 1 1 1];

 %{

Legend of all protein atoms in all residues, used above.
        1 2  3 4 5  6   7   8   9   10  11  12  13  14     (name)           (number of atoms)
1  G:   N CA C O                                         : Glycine          (4)
2  A:   N CA C O CB                                      : Alanine          (5)
3  C:   N CA C O CB SG                                   : Cysteine         (6)
4  S:   N CA C O CB OG                                   : Serine           (6)
5  T:   N CA C O CB OG1 CG2                              : Threonine        (7)
6  P:   N CA C O CB CG  CD                               : Proline          (7)
7  V:   N CA C O CB CG1 CG2                              : Valine           (7)
8  I:   N CA C O CB CG1 CG2 CD1                          : Isoleucine       (8)
9  L:   N CA C O CB CG  CD1 CD2                          : Leucine          (8)
10 M:   N CA C O CB CG  SD  CE                           : Methionine       (8)
11 N:   N CA C O CB CG  ND2 OD1                          : Asparagine       (8)
12 D:   N CA C O CB CG  OD1 OD2                          : Aspartic acid    (8)
13 E:   N CA C O CB CG  CD  OE1 OE2                      : Glutamic acid    (9)
14 Q:   N CA C O CB CG  CD  NE2 OE1                      : Glutamine        (9)
15 K:   N CA C O CB CG  CD  CE  NZ                       : Lysine           (9)
16 H:   N CA C O CB CG  CD2 ND1 CE1 NE2                  : Histidine        (10)
17 R:   N CA C O CB CG  CD  NE  CZ  NH1 NH2              : Arginine         (11)
18 F:   N CA C O CB CG  CD1 CD2 CE1 CE2 CZ               : Phenylalanine    (11)
19 Y:   N CA C O CB CG  CD1 CD2 CE1 CE2 CZ  OH           : Tyrosine         (12)
20 W:   N CA C O CB CG  CD1 CD2 CE2 CE3 NE1 CZ2 CZ3 CH2  : Tryptophan       (14)
%}

% Store rVdW in (8) column, and DGsolv in (9) column

cnt=1; % Reset counter

while cnt<length(numPDB(:,1))
if cnt==1 && strcmp(txtPDB(cnt,1),'CA') % Correction since sometimes PDB files don't have N terminus and start at CA.
        if strcmp(txtPDB(cnt,2),'GLY')==1
        j=1;
        k=3;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'ALA')==1
        j=2;
        k=4;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'CYS')==1
        j=3;
        k=5;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'SER')==1
        j=4;
        k=5;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'THR')==1
        j=5;
        k=6;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'PRO')==1
        j=6;
        k=6;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'VAL')==1
        j=7;
        k=6;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'ILE')==1
        j=8;
        k=7;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'LEU')==1
        j=9;
        k=7;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'MET')==1
        j=10;
        k=7;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'ASN')==1
        j=11;
        k=7;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'ASP')==1
        j=12;
        k=7;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'GLU')==1
        j=13;
        k=8;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'GLN')==1
        j=14;
        k=8;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'LYS')==1
        j=15;
        k=8;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'HIS')==1
        j=16;
        k=9;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'ARG')==1
        j=17;
        k=10;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'PHE')==1
        j=18;
        k=10;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'TYR')==1
        j=19;
        k=11;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'TRP')==1
        j=20;
        k=13;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        else
            cnt=cnt+1;
        end
    else
        if strcmp(txtPDB(cnt,2),'GLY')==1
        j=1;
        k=4;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'ALA')==1
        j=2;
        k=5;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'CYS')==1
        j=3;
        k=6;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'SER')==1
        j=4;
        k=6;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'THR')==1
        j=5;
        k=7;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'PRO')==1
        j=6;
        k=7;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'VAL')==1
        j=7;
        k=7;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'ILE')==1
        j=8;
        k=8;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'LEU')==1
        j=9;
        k=8;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'MET')==1
        j=10;
        k=8;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'ASN')==1
        j=11;
        k=8;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'ASP')==1
        j=12;
        k=8;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'GLU')==1
        j=13;
        k=9;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'GLN')==1
        j=14;
        k=9;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'LYS')==1
        j=15;
        k=9;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'HIS')==1
        j=16;
        k=10;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'ARG')==1
        j=17;
        k=11;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'PHE')==1
        j=18;
        k=11;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'TYR')==1
        j=19;
        k=12;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,2),'TRP')==1
        j=20;
        k=14;
        numPDB(cnt:cnt+(k-1),8)=rVdW(j,1:k)';
        numPDB(cnt:cnt+(k-1),9)=dG(j,1:k)';
        numPDB(cnt:cnt+(k-1),12)=arom_aliph(j,1:k)';
        cnt=cnt+k;
        elseif strcmp(txtPDB(cnt,1),'OXT')==1 % C terminus
            numPDB(cnt,8)=1.4;
            numPDB(cnt,9)=-10;
        else
            cnt=cnt+1;
        end
    end
end

for i=1:length(numPDB(:,1))
    numPDB(i,10)=4*pi*(1.4+numPDB(i,8))^2; % (10): surface of sphere of radius rVdW + water (Solvent Accessible Surface Area, SASA)
    numPDB(i,11)=numPDB(i,6)/numPDB(i,10); % (11): fraction surface exposure
end


%% Arrange surface atoms matrix
% Scan through the the atoms in numPDB, take those with >2% solvent accessibility, and copy them in Surfnum array.
% Do the same with each atom's related text in txtPDB
 
  cnt = 1; % Reset counter
  
for i = 1:length(numPDB(:,1))     
    if numPDB(i,11) > 0.02 % check if the atom in numPDB is in the surface with more than 2% solvent accessibility
             Surftxt(cnt,:) = txtPDB(i,:);
             Surfnum(cnt,:) = numPDB(i,:);
             cnt = cnt + 1; % increment the counter 
        else % do nothing           
    end
end

%% Calculation of Distances and Cumulative Free Energies of Solvation
% Loop through all surface exposed atoms 'j', calculate its pairwise distance
% with every other atom 'i' in the crystal structure, and calculate the
% contribution of the Gaussian curve of 'j' on every atom 'i', and copy it
% in column (7) of Surfnum.


for i = 1:length(Surfnum(:,1))
    for j = i:length(Surfnum(:,1))
            d_j = distance(Surfnum(i,1),Surfnum(i,2),Surfnum(i,3),...
                  Surfnum(j,1),Surfnum(j,2),Surfnum(j,3)); 
              
            d_matrix(i,j)=d_j; % Copy distance value in distance array
            d_matrix(j,i)=d_matrix(i,j); % Symmetric array
            
            stdev_factor=0.25; % Correction factor for interaction lengthscale adaptation to Gaussian standard deviation
            
            % Defining interaction lengthscale as 6 Å for charged atoms,
            % and 3.5 Å for uncharged atoms, according to Lazaridis et al.
            % (1999).
            if strcmp(Surftxt(i,2),'ARG') && ((strcmp(Surftxt(i,1),'NH1')) || strcmp(Surftxt(i,1),'NH2'))
               Surfnum(j,7)=(Surfnum(i,9)*Surfnum(i,11)*normpdf(d_j,0,stdev_factor*6)/normpdf(0,0,stdev_factor*6))+Surfnum(j,7);
            elseif strcmp(Surftxt(i,2),'ASP') && ((strcmp(Surftxt(i,1),'OD1')) || strcmp(Surftxt(i,1),'OD2'))
               Surfnum(j,7)=(Surfnum(i,9)*Surfnum(i,11)*normpdf(d_j,0,stdev_factor*6)/normpdf(0,0,stdev_factor*6))+Surfnum(j,7);
            elseif strcmp(Surftxt(i,2),'GLU') && ((strcmp(Surftxt(i,1),'OE1')) || strcmp(Surftxt(i,1),'OE2'))
               Surfnum(j,7)=(Surfnum(i,9)*Surfnum(i,11)*normpdf(d_j,0,stdev_factor*6)/normpdf(0,0,stdev_factor*6))+Surfnum(j,7);
            elseif strcmp(Surftxt(i,2),'LYS') && strcmp(Surftxt(i,1),'NZ')
               Surfnum(j,7)=(Surfnum(i,9)*Surfnum(i,11)*normpdf(d_j,0,stdev_factor*6)/normpdf(0,0,stdev_factor*6))+Surfnum(j,7);
            else
               Surfnum(j,7)=(Surfnum(i,9)*Surfnum(i,11)*normpdf(d_j,0,stdev_factor*3.5)/normpdf(0,0,stdev_factor*3.5))+Surfnum(j,7);
            end
    end
end

%% Plot surface atoms (1)

figure(1)
hold on
scatter3(Surfnum(:,3),Surfnum(:,4),Surfnum(:,5),10,[0 0 1],'filled')
    axis equal
    view (45,35)
    grid on
    xlabel('x coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
    ylabel('y coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
    zlabel('z coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
    set(gcf, 'Position', get(0, 'Screensize'));
       
    % Calculate geometrical centroid
    x_centroid=mean(Surfnum(:,3));
    y_centroid=mean(Surfnum(:,4));
    z_centroid=mean(Surfnum(:,5));

% Plot centroid    
scatter3(x_centroid,y_centroid,z_centroid,'*','r')    
hold off

%% Separate atoms with hydrophobicity above threshold

cnt = 1; % reset counter
  
for i = 1:length(Surfnum(:,1))
   if Surfnum(i,7) > 0.00 % check if the atom in numPDB is hydrophobic above threshold
                SurfPatchtxt_pre_dbscan(cnt,:) = Surftxt(i,:); % Atom text identifiers
                SurfPatchnum_pre_dbscan(cnt,1:3) = Surfnum(i,3:5); % Atom coordinates
                SurfPatchnum_pre_dbscan(cnt,4) = Surfnum(i,7); % Accumulated free energy
                SurfPatchnum_pre_dbscan(cnt,6) = Surfnum(i,6); % Area
                SurfPatchnum_pre_dbscan(cnt,7) = Surfnum(i,1); % Atom index
                SurfPatchnum_pre_dbscan(cnt,8) = Surfnum(i,12); % Atom aromatic/aliphatic index
                SurfPatchnum_pre_dbscan(cnt,9) = Surfnum(i,11); % Percentage exposed area
                cnt = cnt + 1; % increment the counter 
   else % do nothing         
   end
end

%% DBSCAN filtering
% DBSCAN filter with clustering parameters:
% - Clustering distance set as 3 times radius of water.
% - Minimum number of atoms to determine a hydrophobic patch: 4.
%
% dbscan function gives value of '-1' to atoms considered noise, and assigns
% each atom the number of patch it belongs to.

SurfPatchnum_pre_dbscan(:,5)=dbscan(SurfPatchnum_pre_dbscan(:,1:3),3*1.4,4);

% Loop to filter out atoms considered noise.
cnt = 1;
for i=1:length(SurfPatchnum_pre_dbscan)
    if SurfPatchnum_pre_dbscan(i,5)~=-1
        SurfPatchtxt(cnt,:) = SurfPatchtxt_pre_dbscan(i,:);
        SurfPatchnum(cnt,:) = SurfPatchnum_pre_dbscan(i,:);
        cnt = cnt + 1;
    else
    end
end

%% Plot surface atoms above hydrophobicity threshold and dbscan filtered

figure(2)
hold on
for i=1:length(SurfPatchnum)
    if SurfPatchnum(i,5)~=-1
       scatter3(SurfPatchnum(i,1),SurfPatchnum(i,2),SurfPatchnum(i,3), 'red', 'filled')
    else %do nothing
    end
end
    axis equal
    view (45,35)
    grid on
    xlabel('x coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
    ylabel('y coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
    zlabel('z coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
    set(gcf, 'Position', get(0, 'Screensize'));
hold off

%% Reporting

n_patches=max(SurfPatchnum(:,5))

%Results array, with as many rows as patches.
Res=zeros(max(SurfPatchnum(:,5)),6);

% Column #1: Number of atoms in a patch.
% Column #2: Maximum dGsolv in the patch.
% Column #3: Cumulative dGsolv of the patch.
% Column #4: Surface area of the patch.
% Column #5: Cumulative dGsolv per unit area of the patch.
% Column #6: Aromaticity index of the patch.

for i=1:max(SurfPatchnum(:,5))
    for j=1:length(SurfPatchnum(:,5))
        if SurfPatchnum(j,5)==i
            Res(i,1)=Res(i,1)+1;
            Res(i,3)=Res(i,3)+SurfPatchnum(j,4);
            Res(i,4)=Res(i,4)+SurfPatchnum(j,6);
            Res(i,6)=Res(i,6)+SurfPatchnum(j,8);
            if SurfPatchnum(j,4)>Res(i,2)
               Res(i,2)=SurfPatchnum(j,4);
            else
            end
        else    
        end
    end
end

Res(:,5)=1000*Res(:,3)./Res(:,4); % *1000 to convert from kJ to J.
Res(:,6)=Res(:,6)./Res(:,1);

figure(3)
hold on
scatter(Res(:,4),Res(:,3),40,[Res(:,6) 1-Res(:,6) 1-Res(:,6)],'filled')
    grid on
    axis square
    xlabel('Patch size (Å^2)','fontsize',16,'fontweight','bold','FontName','Arial')
    xlim([0 ceil(max(Res(:,4)))])
    ylabel('Cumulative patch solvation free energy (kJ/mol)','fontsize',16,'fontweight','bold','FontName','Arial')
    ylim([0 ceil(max(Res(:,3)))])
    set(gcf, 'Position', get(0, 'Screensize'));
    c=colorbar;
    c.FontSize=16;
    c.FontName='Arial';
    c.Label.String='Aromatic index';
    c.Label.FontSize=22;
    c.Label.FontName='Arial';
    colorpalette=[[0:0.01:1]' [1:-0.01:0]' [1:-0.01:0]'];
    colormap(colorpalette)
hold off

figure(4)
hold on
scatter(Res(:,4),Res(:,2),40,[Res(:,6) 1-Res(:,6) 1-Res(:,6)],'filled')
    grid on
    box on
    axis square
    xlabel('Patch size (Å^2)','fontsize',16,'fontweight','bold','FontName','Arial')
    xlim([0 ceil(max(Res(:,4)))])
    ylabel('Max patch solvation free energy (kJ/mol)','fontsize',16,'fontweight','bold','FontName','Arial')
    ylim([0 ceil(max(Res(:,2)))])
    set(gcf, 'Position', get(0, 'Screensize'));
    c=colorbar;
    c.FontSize=16;
    c.FontName='Arial';
    c.Label.String='';
    c.Label.FontSize=16;
    c.Label.FontWeight='bold';
    c.Label.String='Aromaticity index';
    c.Label.FontName='Arial';
    colorpalette=[[0:0.01:1]' [1:-0.01:0]' [1:-0.01:0]'];
    colormap(colorpalette)
hold off

Total_surface_area=sum(numPDB(:,6))/100 % In nm^2
Total_surface_hydrophobic_patches=(sum(Res(:,4)))/100; % In nm^2
Percentage_hydrophobic_coverage=100*(Total_surface_hydrophobic_patches/Total_surface_area) % Expressed in '%'
dGsolv_protein=sum(Surfnum(:,7)) % In kJ/mol
dGsolv_per_tot_area=dGsolv_protein/Total_surface_area % In kJ/mol-nm^2

arom_surf=0;
for i=1:length(SurfPatchnum)
    if SurfPatchnum(i,8)==1
        arom_surf=arom_surf+SurfPatchnum(i,6);
    else
    end
end
Percentage_aromatic_surface=arom_surf/Total_surface_area;


%% Patch reporting

% The output from this is meant to be copied and pasted in Pymol, while the
% analyzed crystal structure is open. This will help visualize the patches 
% by coloring them in red.

fprintf('\nhide\nset surface_quality, 1\nshow surface\ncolor forest\nbg_color white\n\n')
for i=1:max(SurfPatchnum(:,5))
    fprintf('#Patch %i\n',i)
    fprintf('#Area: %f Å^2\n',Res(i,4))
    fprintf('color tv_red, (index ')
    for j=1:length(SurfPatchnum(:,5))
        if i == SurfPatchnum(j,5)
           %fprintf('%d,',j)
            fprintf('%d,',SurfPatchnum(j,7))
        else % do nothing
        end
    end
    fprintf(')\n\n')
end

%% Reporting back in Excel results sheet
writematrix('Enzyme name','Enzymes.xlsx','Sheet','Results','Range','A1')
writematrix('Enzyme ID','Enzymes.xlsx','Sheet','Results','Range','B1')
writematrix('# hydrophobic patches','Enzymes.xlsx','Sheet','Results','Range','C1')
writematrix('% relative hydrophobic surface','Enzymes.xlsx','Sheet','Results','Range','D1')
writematrix('DGsolv/area (kJ/mol·nm^2)','Enzymes.xlsx','Sheet','Results','Range','E1')
writematrix('% aromatic surface','Enzymes.xlsx','Sheet','Results','Range','F1')
rep=readcell('Enzymes.xlsx','Sheet','Results','Range','B2:B1000');
maxrng=length(rep);
check=0;
for i=1:maxrng
    if strcmp(rep(i,1),tab)==1
        check=1;
    else
    end
end
if check~=1
    newcell=sprintf("B%d", maxrng+2); 
    writematrix(tab,'Enzymes.xlsx','Sheet','Results','Range',newcell)
    newcell=sprintf("C%d", maxrng+2); 
    writematrix(length(Res),'Enzymes.xlsx','Sheet','Results','Range',newcell)
    newcell=sprintf("D%d", maxrng+2); 
    writematrix(Percentage_hydrophobic_coverage,'Enzymes.xlsx','Sheet','Results','Range',newcell)
    newcell=sprintf("E%d", maxrng+2); 
    writematrix(dGsolv_per_tot_area,'Enzymes.xlsx','Sheet','Results','Range',newcell)
    newcell=sprintf("F%d", maxrng+2); 
    writematrix(Percentage_aromatic_surface,'Enzymes.xlsx','Sheet','Results','Range',newcell)
end


%% Define distance function

% Used for interatomic distance calculations

function out = distance(rx,ry,rz,tx,ty,tz) % subfunction for calculating distance
                out = sqrt((rx-tx).^2 + (ry-ty).^2 + (rz-tz).^2); 
end
