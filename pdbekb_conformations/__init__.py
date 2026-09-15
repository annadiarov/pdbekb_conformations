"""
pdbekb_conformations
====================

A small pipeline to systematically pull, parse and organise the PDBe-KB
"structural superposition" (conformational clustering) data, and to derive
the summary tables needed for the Bindome project:

  * PDBs labelled by conformational cluster (per UniProt segment)
  * representative structure for each cluster
  * a magnitude-of-change metric per protein (GLOCON, as provided by PDBe-KB)
  * optional local Ca-RMSD between cluster representatives + a per-residue
    deviation profile that localises the moving region.

The public entry point is `pipeline.build_dataset`.
"""

from . import fetch, parse, structure, pipeline  # noqa: F401

__all__ = ["fetch", "parse", "structure", "pipeline"]
__version__ = "0.1.0"
