# hi-patch
Hydrophobic Intensity Patch tool for the study of protein surface hydrophobicity

Authors: Hector Sanchez-Moran, James S. Weltz, Daniel K. Schwartz, Joel L. Kaar.
Department of Chemical & Biological Engineering.
University of Colorado Boulder
Contact: hector.sanchez-moran@colorado.edu

-------------------------------------------------------------------------------

Determination of hydrophobic patches (HPs) of known crystal structures 
based on solvation free energies of solvent exposed atoms.
 
Cite: please stay tuned for upcoming publication.

For crystal structure pre-treatment tutorial, please visit: 
https://www.youtube.com/watch?v=lqa4_L6qwIw&feature=youtu.be

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
- Lines #762-789 for reporting results in Excel just available in MacOS.
