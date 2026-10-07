# Citing CREVICE

Thank you for using CREVICE! A paper describing it is in preparation:

> Lambden, E., Findlay, H. and Booth, P. J. *CREVICE: Analysis and
> Visualisation of Protein Cavities, Channels, and Residue Interaction
> Networks.* In preparation.

Until it is published, please cite the software itself using
[CITATION.cff](../../CITATION.cff) (GitHub shows it as "Cite this
repository"), and say which CREVICE version you used. CREVICE's results depend
on grid spacing, probe radii and atom selection, so please report the commands
and settings too; every run records them in its manifest.

CREVICE reads trajectories and computes periodic distances with MDAnalysis,
so please cite it as well when you analyse simulations:

- Michaud-Agrawal, N., Denning, E. J., Woolf, T. B. and Beckstein, O. (2011)
  MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
  *J. Comput. Chem.* 32, 2319–2327.
- Gowers, R. J. *et al.* (2016) MDAnalysis: a Python package for the rapid
  analysis of molecular dynamics simulations. *Proc. 15th Python in Science
  Conf.*, 98–105.

If you compare CREVICE's channel profiles with HOLE, or use its `hole` radius
set, please also cite Smart, O. S. *et al.* (1996) HOLE: a program for the
analysis of the pore dimensions of ion channel structural models. *J. Mol.
Graph.* 14, 354–360. The [method notes](../methods/index.md) give references
for the other methods CREVICE builds on.

CREVICE is distributed under the MIT licence ([LICENSE](../../LICENSE)).
