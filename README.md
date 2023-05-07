# hi-patch
Hydrophobic Intensity Patch (hi-patch) tool for the study of protein surface hydrophobicity

Authors: Hector Sanchez-Moran, Daniel K. Schwartz, Joel L. Kaar.
Department of Chemical & Biological Engineering.
University of Colorado Boulder
Contact: hector.sanchez-moran@colorado.edu

-------------------------------------------------------------------------------

Determination of hydrophobic patches (HPs) of known crystal structures 
based on solvation free energies of solvent exposed atoms.

Latest update: hi-patch_v2 (05-07-2023). See version notes below.
 
Cite: Sánchez-Morán et al. ACS Appl. Mater. Interfaces. 2021, 13, 23, 26694-26703.

For crystal structure pre-treatment tutorial, read Instructions_hipatch.pdf.

Please, place in the same folder:
- hi-patch script.
- areas.py
- Proteins Excel file. In this file, metadata of the proteins to analyze is
stored, as well as the quantitative output of each protein.

Excel file formatting:
- Column A: Atom index
- Column B: Atom type within residue
- Column C: Residue type
- Column D: Residue index
- Column E-G: Coordinates in X-Y-Z
- Column H: Element
- Column I: Output from areas.py script executed in Pymol.

Notes: 
- Crystal structures in the Protein Data Bank are pretty diverse, and outliars help polish the code. Contact to fix glitches is greatly appreciated.
- Lines for reporting results in Excel just available in MacOS.
- It is absolutely necessary for correct results that all input proteins have all side chains resolved and have no isomers. If there are unresolved regions it is not a problem, as long as unresolved residues are fully unresolved. For partially resolved residues, it is recommended to run pdbfixer (https://sbgrid.org/software/titles/pdbfixer) or some other tool. Computationally generated structures (e.g., AlphaFold2) rarely have these issues.

Version v2:
Major hi-patch update comprising numerous upgrades.
- More efficient code (computation time ~1/3 less than previous release).
- Fixed numerous glitches related to multimeric structures.
- Incorporated mapping of hydrophobic patches using a color gradient to differentiate hydrophobicity of each patch based on a hydrophobicity index.
- Incorporated hydrophobicity directionality vector, which points towards the region of the protein where there are more hydrophobic patches.
- Made a small change in areas.py file to accelerate areas calculation at the cost of minimally reducing precision. This can be changed by setting 'dot_density' to 4. Not recommended to set dot_density to lower than 3.
