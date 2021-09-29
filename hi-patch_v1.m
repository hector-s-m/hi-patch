%% HI-PATCH script.
%{ 
Authors: Hector Sanchez-Moran, Daniel K. Schwartz, Joel L. Kaar.
Department of Chemical & Biological Engineering.
University of Colorado Boulder

Determination of hydrophobic patches (HPs) of known crystal
structures based on solvation free energies of solvent exposed atoms.

Cite: ACS Appl. Mater. Interfaces 2021, 13, 23, 26694?26703

Notes:
    · Crystal structures in the Protein Data Bank are pretty diverse,
      and outliars help polish the code. Contact to fix glitches is
      greatly appreciated (hector.sanchez-moran@colorado.edu)
    · Reporting results in Excel just available in MacOS.

KEEP THIS HEADING IN THE CODE.
%}
%% Clear
close all
clear all
clc 
%% Inputs
Excel_booklet='Enzymes.xlsx'; % Select Excel booklet where protein information is located.
tab='CALB'; % Select tab name (protein)
excel_import = readtable(Excel_booklet,'ReadVariableNames',false,'Sheet',tab); % Separate numbers and text into two different arrayc
%% Outputs
% Select which plots to display
Outputs=struct;
Outputs.SurfaceAtoms = 'no';
Outputs.HydrophobicAtoms = 'no';
Outputs.HydrophobicPolarPlots2D = 'yes';
Outputs.Max_dGsolv_vs_Area = 'yes';
Outputs.Cum_dGsolv_vs_Area = 'no';
Outputs.PatchesHistogram = 'yes';
Outputs.Pymol_patch_export = 'yes';
Outputs.Excel_results_export = 'no'; % Only MacOS
%% Arrange data arrays

% Create PlotString array, to use in output reporting function.
PlotString=cell(0,0);
if strcmp(Outputs.SurfaceAtoms,'yes')==1
    PlotString{end+1} = 'SurfaceAtoms';    
end
if strcmp(Outputs.HydrophobicAtoms,'yes')==1
    PlotString{end+1} = 'HydrophobicAtoms';    
end
if strcmp(Outputs.HydrophobicPolarPlots2D,'yes')==1
    PlotString{end+1} = 'HydrophobicPolarPlots2D';    
end
if strcmp(Outputs.Max_dGsolv_vs_Area,'yes')==1
    PlotString{end+1} = 'Max_dGsolv_vs_Area';    
end
if strcmp(Outputs.Cum_dGsolv_vs_Area,'yes')==1
    PlotString{end+1} = 'Cum_dGsolv_vs_Area';    
end
if strcmp(Outputs.PatchesHistogram,'yes')==1
    PlotString{end+1} = 'PatchesHistogram';    
end
if strcmp(Outputs.Pymol_patch_export,'yes')==1
    PlotString{end+1} = 'Pymol_patch_export';    
end
if strcmp(Outputs.Excel_results_export,'yes')==1
    PlotString{end+1} = 'Excel_results_export';    
end

% Separate the information from the Excel file into the corresponding
% structured arrays.
data=struct;
data.at_index = [1:1:height(excel_import)]'; % copy the atom index
data.res_index = excel_import{:,4}; % copy residue index
data.coords = excel_import{:,5:7}; % copy XYZ coords
data.surf_exp = excel_import{:,9}; % copy the surface exposure
%data(~any(~isnan(data)),:)=[]; % remove any NaNs
data.at_name = excel_import{:,2}; % copy the atom name
data.res_name = excel_import{:,3}; % copy the residue name
data.at_type = excel_import{:,8}; % copy the atom type
data.rVdW = zeros(length(data.at_index),1);
data.dG = zeros(length(data.at_index),1);
data.arom = zeros(length(data.at_index),1);
data.sphere_surf = zeros(length(data.at_index),1);
data.frac_SASA = zeros(length(data.at_index),1);


%{
Van der Waals radii array in Angstroms
Obtained at: Gerstein, M., Tsai, J., & Levitt, M. (1995). The Volume of 
Atoms on the Protein Surface: Calculated from Simulation, using Voronoi
Polyhedra. In J. Mol. Bl;ol (Vol. 249). 
%}
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

%{
Free energies array

Obtained in:
Stouten, P. F. W., Frömmel, C., Nakamura, H., & Sander, C. (1993). An 
Effective Solvation Term Based on Atomic Occupancies for Use in Protein 
Simulations. Molecular Simulation, 10(2-6), 97-120. doi:10.1080/08927029308022161 

&

Lazaridis, T., & Karplus, M. (1999). Effective energy function for proteins
in solution. Proteins: Structure, Function, and Bioinformatics, 35(2), 133
152. https://doi.org/10.1002/(SICI)1097-0134(19990501)35:2<133::AID-PROT1>3.0.CO;2-N
%}

% Lay out common values and exceptions
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

% Binary array, describing aromatic and aliphatic nature of each atom.
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
Legend of all protein atoms in all residues, used in PDB file metadata.
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

% Store rVdW, dGsolv and aromaticity in data array.
cnt=1;
while cnt<length(data.at_index(:,1))
    if cnt==1 && strcmp(data.at_name(cnt,1),'CA') % Correction since sometimes PDB files don't have N terminus and start at CA.
        if strcmp(data.res_name(cnt),'GLY')==1
            j=1;
            k=3;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'ALA')==1
            j=2;
            k=4;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'CYS')==1
            j=3;
            k=5;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'SER')==1
            j=4;
            k=5;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'THR')==1
            j=5;
            k=6;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'PRO')==1
            j=6;
            k=6;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'VAL')==1
            j=7;
            k=6;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'ILE')==1
            j=8;
            k=7;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'LEU')==1
            j=9;
            k=7;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'MET')==1
            j=10;
            k=7;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'ASN')==1
            j=11;
            k=7;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'ASP')==1
            j=12;
            k=7;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'GLU')==1
            j=13;
            k=8;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'GLN')==1
            j=14;
            k=8;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'LYS')==1
            j=15;
            k=8;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'HIS')==1
            j=16;
            k=9;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'ARG')==1
            j=17;
            k=10;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'PHE')==1
            j=18;
            k=10;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'TYR')==1
            j=19;
            k=11;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'TRP')==1
            j=20;
            k=13;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        else
            cnt=cnt+1;
        end
    else
        if strcmp(data.res_name(cnt),'GLY')==1
            j=1;
            k=4;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'ALA')==1
            j=2;
            k=5;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'CYS')==1
            j=3;
            k=6;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'SER')==1
            j=4;
            k=6;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'THR')==1
            j=5;
            k=7;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'PRO')==1
            j=6;
            k=7;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'VAL')==1
            j=7;
            k=7;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'ILE')==1
            j=8;
            k=8;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'LEU')==1
            j=9;
            k=8;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'MET')==1
            j=10;
            k=8;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'ASN')==1
            j=11;
            k=8;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'ASP')==1
            j=12;
            k=8;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'GLU')==1
            j=13;
            k=9;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'GLN')==1
            j=14;
            k=9;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'LYS')==1
            j=15;
            k=9;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'HIS')==1
            j=16;
            k=10;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'ARG')==1
            j=17;
            k=11;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'PHE')==1
            j=18;
            k=11;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'TYR')==1
            j=19;
            k=12;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.res_name(cnt),'TRP')==1
            j=20;
            k=14;
            data.rVdW(cnt:cnt+(k-1))=rVdW(j,1:k)';
            data.dG(cnt:cnt+(k-1))=dG(j,1:k)';
            data.arom(cnt:cnt+(k-1))=arom_aliph(j,1:k)';
            cnt=cnt+k;
        elseif strcmp(data.at_name(cnt),'OXT')==1 % C terminus
            data.rVdW(cnt,1)=1.4;
            data.dG(cnt,1)=-10;
            
        else
            cnt=cnt+1;
        end
    end
end
%%
% Surface related calculations 
data.sphere_surf(:,1)=4*pi*(1.4+data.rVdW(:)).^2; % Surface of sphere of radius rVdW + water (Solvent Accessible Surface Area, SASA)
data.frac_SASA(:,1)=data.surf_exp(:)./data.sphere_surf(:); % Fraction surface exposure with respect to ideal isolated spherical atom

%for i=1:length(data.at_index)
 %   data.sphere_surf(i,1)=4*pi*(1.4+data.rVdW(i))^2; % Surface of sphere of radius rVdW + water (Solvent Accessible Surface Area, SASA)
  %  data.frac_SASA(i,1)=data.surf_exp(i)/data.sphere_surf(i); % Fraction surface exposure with respect to ideal isolated spherical atom
%end
%% Arrange surface atoms matrix
% Scan through the the atoms in numPDB, take those with >2% solvent accessibility, and copy them in Surfnum array.
% Do the same with each atom's related text in txtPDB
 
cnt = 1; % Reset counter
surfdata=struct;
% Create surfdata struct which stores information of all surface exposed
% atoms.
for i = 1:length(data.at_index)     
    if data.frac_SASA(i) > 0.02 % check if atom has more than 2% solvent accessibility
            data.at_surf_exp(cnt,1) = i;
            surfdata.at_index(cnt,1) = data.at_index(i);
            surfdata.res_index(cnt,1) = data.res_index(i);
            surfdata.coords(cnt,1:3) = data.coords(i,1:3);
            surfdata.surf_exp(cnt,1) = data.surf_exp(i);
            surfdata.at_name(cnt,1) = data.at_name(i);
            surfdata.res_name(cnt,1) = data.res_name(i);
            surfdata.at_type(cnt,1) = data.at_type(i);
            surfdata.rVdW(cnt,1) = data.rVdW(i);
            surfdata.dG(cnt,1) = data.dG(i);
            surfdata.arom(cnt,1) = data.arom(i);
            surfdata.sphere_surf(cnt,1) = data.sphere_surf(i);
            surfdata.frac_SASA(cnt,1) = data.frac_SASA(i);
            cnt = cnt + 1; 
        else % do nothing           
    end
end
%% Primary amines matrix
% Scan through the the atoms in myPDbnumhyd and find their coordinates in
% myPDBnum, filtering out H atoms.
 
cnt = 1; % set counter
lyscount = sum(count(surfdata.at_name(:,1),'NZ')); % count number of lysine primary NH2s.
surfdata.lys_coords=zeros(lyscount,3);
for i = 1:length(surfdata.at_index(:,1))     
    % check if the atom in myPDBnum is a primary amine
    if strcmp(surfdata.at_name(i,1),'NZ')
             surfdata.lys_coords(cnt,1:3) = surfdata.coords(i,1:3);
             cnt = cnt + 1;
        else
    end
end
%% Calculation of Distances and Cumulative Free Energies of Solvation
% Loop through all surface exposed atoms 'j', calculate its pairwise distance
% with every other atom 'i' in the crystal structure, and calculate the
% contribution of the Gaussian curve of 'j' on every atom 'i', and copy it
% in column (7) of Surfnum.

surfdata.cum_dG=zeros(length(surfdata.at_name(:,1)),1);
d_matrix=zeros(length(surfdata.at_name(:,1)),1);

% Loop for finding interatomic distances (d_matrix) and compute gaussian
% projection based on dGsolv of each atom, with an interaction lengthscale
% of defined in Lazaridis et al. (1999). 6 Å for charged atoms,
% and 3.5 Å for uncharged atoms 

for i = 1:length(surfdata.at_name(:,1))
    for j = i:length(surfdata.at_name(:,1))
            % Compute interatomic Euclidean distances
            d_j = distance(surfdata.coords(i,1),surfdata.coords(i,2),surfdata.coords(i,3),surfdata.coords(j,1),surfdata.coords(j,2),surfdata.coords(j,3)); 
            d_matrix(i,j)=d_j; % Copy distance value in distance array
            d_matrix(j,i)=d_matrix(i,j); % Symmetric array           
            stdev_factor=0.25; % Correction factor for interaction lengthscale adaptation to Gaussian standard deviation

            if strcmp(surfdata.res_name(i),'ARG') && ((strcmp(surfdata.at_name(i),'NH1')) || strcmp(surfdata.at_name(i),'NH2'))
               surfdata.cum_dG(j,1)=(surfdata.dG(i)*surfdata.frac_SASA(i)*normpdf(d_j,0,stdev_factor*6)/normpdf(0,0,stdev_factor*6))+surfdata.cum_dG(j,1);
            elseif strcmp(surfdata.res_name(i),'ASP') && ((strcmp(surfdata.at_name(i),'OD1')) || strcmp(surfdata.at_name(i),'OD2'))
               surfdata.cum_dG(j,1)=(surfdata.dG(i)*surfdata.frac_SASA(i)*normpdf(d_j,0,stdev_factor*6)/normpdf(0,0,stdev_factor*6))+surfdata.cum_dG(j,1);
            elseif strcmp(surfdata.res_name(i),'GLU') && ((strcmp(surfdata.at_name(i),'OE1')) || strcmp(surfdata.at_name(i),'OE2'))
               surfdata.cum_dG(j,1)=(surfdata.dG(i)*surfdata.frac_SASA(i)*normpdf(d_j,0,stdev_factor*6)/normpdf(0,0,stdev_factor*6))+surfdata.cum_dG(j,1);
            elseif strcmp(surfdata.res_name(i),'LYS') && strcmp(surfdata.at_name(i),'NZ')
               surfdata.cum_dG(j,1)=(surfdata.dG(i)*surfdata.frac_SASA(i)*normpdf(d_j,0,stdev_factor*6)/normpdf(0,0,stdev_factor*6))+surfdata.cum_dG(j,1);
            else
               surfdata.cum_dG(j,1)=(surfdata.dG(i)*surfdata.frac_SASA(i)*normpdf(d_j,0,stdev_factor*3.5)/normpdf(0,0,stdev_factor*3.5))+surfdata.cum_dG(j,1);
            end
    end
end

% Calculate geometrical centroid
data.x_centroid=mean(surfdata.coords(:,1));
data.y_centroid=mean(surfdata.coords(:,2));
data.z_centroid=mean(surfdata.coords(:,3));
%% Separate atoms with hydrophobicity above threshold

cnt = 1;
for i = 1:length(surfdata.at_name(:,1))
   if surfdata.cum_dG(i) > 0.00 % check if the atom is hydrophobic
       surfdata.pre_dbscan_index(cnt,1) = i;
                cnt = cnt + 1;
       else % do nothing         
   end
end
%% DBSCAN clustering and filtering
% DBSCAN filter with clustering parameters:
% - Clustering distance set as 3 times radius of water.
% - Minimum number of atoms to determine a hydrophobic patch: 4.
%
% dbscan function gives value of '-1' to atoms considered noise, and assigns
% each atom the number of patch it belongs to.

patchdata=struct;

surfdata.patch_index(:,1) = dbscan(surfdata.coords(surfdata.pre_dbscan_index(:,1),1:3),3*1.4,4);
surfdata.patch_index(:,2) = surfdata.pre_dbscan_index(:,1);
%% Arrange patch belonging atoms

% Get rid of atoms non belonging to patches (patch_index = -1).
surfdata.patch_index(surfdata.patch_index(:,1)<0,:)=[];

% Generate patchdata struct which stores information about patch-beloning
% atoms.
patchdata.at_index(:,:)=surfdata.at_index(surfdata.patch_index(:,2),:);
patchdata.res_index(:,:)=surfdata.res_index(surfdata.patch_index(:,2),:);
patchdata.coords(:,:)=surfdata.coords(surfdata.patch_index(:,2),:);
patchdata.surf_exp(:,:)=surfdata.surf_exp(surfdata.patch_index(:,2),:);
patchdata.at_name(:,:)=surfdata.at_name(surfdata.patch_index(:,2),:);
patchdata.res_name(:,:)=surfdata.res_name(surfdata.patch_index(:,2),:);
patchdata.at_type(:,:)=surfdata.at_type(surfdata.patch_index(:,2),:);
patchdata.rVdW(:,:)=surfdata.rVdW(surfdata.patch_index(:,2),:);
patchdata.dG(:,:)=surfdata.dG(surfdata.patch_index(:,2),:);
patchdata.arom(:,:)=surfdata.arom(surfdata.patch_index(:,2),:);
patchdata.sphere_surf(:,:)=surfdata.sphere_surf(surfdata.patch_index(:,2),:);
patchdata.frac_SASA(:,:)=surfdata.frac_SASA(surfdata.patch_index(:,2),:);
patchdata.cum_dG(:,:)=surfdata.cum_dG(surfdata.patch_index(:,2),:);
patchdata.patch_index(:,:)=surfdata.patch_index(:,:);
patchdata.patch_index(:,2) = [1:1:length(patchdata.patch_index(:,1))];
%% 2D spherical coordinate plots

% Generation of vectors that go from the centroid to each atom, and
% compute azimuthal and polar coordinates of each.
patchdata.solvation_vectors = zeros(length(patchdata.coords(:,1)),3);
patchdata.solvation_vectors(:,:) = [(patchdata.coords(:,1)-data.x_centroid),(patchdata.coords(:,2)-data.y_centroid),(patchdata.coords(:,3)-data.z_centroid)];

surfdata.lys_vectors = zeros(length(surfdata.lys_coords(:,1)),3);
surfdata.lys_vectors = [(surfdata.lys_coords(:,1)-data.x_centroid),(surfdata.lys_coords(:,2)-data.y_centroid),(surfdata.lys_coords(:,3)-data.z_centroid)];

surfdata.lys_polar_coords = zeros(length(surfdata.lys_coords(:,1)),2);
[surfdata.lys_polar_coords(:,1),surfdata.lys_polar_coords(:,2)] = cart2sph(surfdata.lys_vectors(:,1),surfdata.lys_vectors(:,2),surfdata.lys_vectors(:,3));
surfdata.lys_polar_coords(:,1) = 180*surfdata.lys_polar_coords(:,1)/pi;
surfdata.lys_polar_coords(:,2) = 180*surfdata.lys_polar_coords(:,2)/pi;

patchdata.at_polar_coords = zeros(length(patchdata.coords(:,1)),2);
[patchdata.at_polar_coords(:,1),patchdata.at_polar_coords(:,2)] = cart2sph(patchdata.solvation_vectors(:,1),patchdata.solvation_vectors(:,2),patchdata.solvation_vectors(:,3));
patchdata.at_polar_coords(:,1) = 180*patchdata.at_polar_coords(:,1)/pi;
patchdata.at_polar_coords(:,2) = 180*patchdata.at_polar_coords(:,2)/pi;
%% Calculations

%Results struct, with as many rows as patches.
results=struct;
results.n_patches=max(patchdata.patch_index(:,1)); % Number of patches
results.n_atoms_in_patch=zeros(results.n_patches,1); % Number of atoms in a patch.
results.max_dGsolv=zeros(results.n_patches,1);% Maximum dGsolv in the patch.
results.cum_dGsolv=zeros(results.n_patches,1);% Cumulative dGsolv of the patch.
results.area=zeros(results.n_patches,1);% Surface area of the patch.
results.dGsolv_per_area=zeros(results.n_patches,1);% Cumulative dGsolv per unit area of the patch.
results.aromaticity=zeros(results.n_patches,1);% Aromaticity index of the patch.

for i=1:results.n_patches
    for j=1:length(patchdata.at_index)
        if patchdata.patch_index(j,1)==i
            results.n_atoms_in_patch(i,1)=results.n_atoms_in_patch(i,1)+1;
            results.cum_dGsolv(i,1)=results.cum_dGsolv(i,1)+patchdata.cum_dG(j,1);
            results.area(i,1)=results.area(i,1)+patchdata.surf_exp(j,1);
            results.aromaticity(i,1)=results.aromaticity(i,1)+patchdata.arom(j,1);
            if patchdata.cum_dG(j,1)>results.max_dGsolv(i,1)
               results.max_dGsolv(i,1)=patchdata.cum_dG(j,1);
            else
            end
        else    
        end
    end
end

results.dGsolv_per_area(:,1)=1000*results.cum_dGsolv(:,1)./results.area(:,1); % dGsolv/area of each patch. *1000 to convert from kJ to J.
results.aromaticity(:,1)=results.aromaticity(:,1)./results.n_atoms_in_patch(:,1);
results.total_surface_area=sum(data.surf_exp(:,1))/100; % In nm^2
results.hydrophobic_patches_area=(sum(results.area(:,1)))/100; % In nm^2
results.percentage_hydrophobic_coverage=100*(results.hydrophobic_patches_area/results.total_surface_area); % Expressed in '%'
results.dGsolv_protein=sum(surfdata.cum_dG(:,1)); % In kJ/mol
results.dGsolv_per_area=results.dGsolv_protein/results.total_surface_area; % In kJ/mol-nm^2
results.percentage_aromatic_surface=0;

for i=1:length(patchdata.arom)
    if patchdata.arom(i)==1
        results.percentage_aromatic_surface=results.percentage_aromatic_surface+patchdata.surf_exp(i);
    else % do nothing
    end
end

results.percentage_aromatic_surface=results.percentage_aromatic_surface/results.total_surface_area;
%% Call output reporting function
hipatch_plotter(data,surfdata,patchdata,results,PlotString)
%% Functions

% Interatomic distance calculation

function out = distance(rx,ry,rz,tx,ty,tz)
                out = sqrt((rx-tx).^2 + (ry-ty).^2 + (rz-tz).^2); 
end

% Results plotter / reporting

function hipatch_plotter(data,surfdata,patchdata,results,PlotString)

for ii=1:length(PlotString)
switch PlotString{ii}
    case 'SurfaceAtoms'
        figure();
        hold on
        scatter3(surfdata.coords(:,1),surfdata.coords(:,2),surfdata.coords(:,3),10,[0 0 1],'filled')
            axis equal
            view (45,35)
            grid on
            xlabel('x coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
            ylabel('y coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
            zlabel('z coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
            set(gcf, 'Position', get(0, 'Screensize'));
        % Plot centroid    
        scatter3(data.x_centroid,data.y_centroid,data.z_centroid,'*','r')    
        hold off
        
    case 'HydrophobicAtoms'
        figure();
        hold on
        scatter3(patchdata.coords(:,1),patchdata.coords(:,2),patchdata.coords(:,3),'red', 'filled')
            axis equal
            view (45,35)
            grid on
            xlabel('x coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
            ylabel('y coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
            zlabel('z coord (Å)','fontsize',22,'fontweight','bold','FontName','Arial')
            set(gcf, 'Position', get(0, 'Screensize'));
        hold off

    case 'HydrophobicPolarPlots2D'
        figure();
        hold on
        %Plot atoms belonging to hydrophobic patches
        for i=1:length(patchdata.at_polar_coords(:,1))
            scatter(patchdata.at_polar_coords(i,1),patchdata.at_polar_coords(i,2),100,[1 0 0],'filled',...
                'MarkerFaceAlpha',patchdata.cum_dG(i,1)/max(patchdata.cum_dG(:,1)));
        end
        %Plot solvent exposed lysines
        scatter(surfdata.lys_polar_coords(:,1),surfdata.lys_polar_coords(:,2),50,[0 0 1],'filled')
        axis square
        xlim([-180 180])
        xticks([-180 -120 -60 0 60 120 180])
        xlabel('Azimuthal angle (deg)','fontsize',16,'FontName','Arial')
        ylim([-90 90])
        yticks([-90 -60 -30 0 30 60 90])
        ylabel('Polar angle (deg)','fontsize',16,'FontName','Arial')
        set(gcf, 'Position', get(0, 'Screensize'));
        c=colorbar;
            c.FontSize=16;
            c.FontName='Arial';
            c.Label.String='Normalized hydrophobicity';
            c.Label.FontSize=22;
            c.Label.FontName='Arial';
            colorpalette=[ones(101,1) [1:-0.01:0]' [1:-0.01:0]'];
            colormap(colorpalette)
        hold off

    case 'Max_dGsolv_vs_Area'
        figure();
        hold on
        scatter(results.area(:,1),results.max_dGsolv(:,1),40,[results.aromaticity(:,1) 1-results.aromaticity(:,1) 1-results.aromaticity(:,1)],'filled')
            grid on
            box on
            axis square
            xlabel('Patch size (Å^2)','fontsize',16,'fontweight','bold','FontName','Arial')
            xlim([0 ceil(max(results.area(:,1)))])
            ylabel('Max patch solvation free energy (kJ/mol)','fontsize',16,'fontweight','bold','FontName','Arial')
            ylim([0 ceil(max(results.max_dGsolv(:,1)))])
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
        
    case 'Cum_dGsolv_vs_Area'
        figure();
        hold on
        scatter(results.area(:,1),results.cum_dGsolv(:,1),40,[results.aromaticity(:,1) 1-results.aromaticity(:,1) 1-results.aromaticity(:,1)],'filled')
            grid on
            axis square
            xlabel('Patch size (Å^2)','fontsize',16,'fontweight','bold','FontName','Arial')
            xlim([0 ceil(max(results.area(:,1)))])
            ylabel('Cumulative patch solvation free energy (kJ/mol)','fontsize',16,'fontweight','bold','FontName','Arial')
            ylim([0 ceil(max(results.cum_dGsolv(:,1)))])
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

    case 'PatchesHistogram'
        figure();
        histogram(results.area(:,1),'BinWidth',50,'BinLimits',[0,1000])
        ylim([0 70])
        axis square
        xlabel('Patch size (Å^2)','fontsize',16,'fontweight','bold','FontName','Arial')
        ylabel('Frequency','fontsize',16,'fontweight','bold','FontName','Arial')

        case 'Pymol_patch_export'
        fprintf('\nhide\nset surface_quality, 1\nshow surface\ncolor forest\nbg_color white\n\n')
        for i=1:results.n_patches
            fprintf('#Patch %i\n',i)
            fprintf('#Area: %f Å^2\n',results.area(i))
            fprintf('color tv_red, (index ')
            for j=1:length(patchdata.at_index)
                if i == patchdata.patch_index(j)
                    fprintf('%d,',patchdata.at_index(j))
                else % do nothing
                end
            end
            fprintf(')\n\n')
        end
        
    case 'Excel_results_export'
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
end
end
end