"""Petab2PyAntiGen: replace PEtab problems with compact PyAntiGen projects.

A PEtab problem is an SBML file plus five or six TSV tables.  A PyAntiGen
project writes the same problem as a few lines of Python f-strings that loop
over species, isotopes and compartments, so a regular model (a methylation
grid, a ladder of binding steps, ...) collapses to a handful of loops.

Modules
-------
petab_reader      load a PEtab v1 problem (yaml + tables + SBML)
sbml_to_antimony  SBML -> Antimony through Tellurium, minus the clutter
reference         simulate the PEtab problem itself (RoadRunner) and score it
verify            compare a PyAntiGen project with the PEtab reference
stats             size of the PEtab problem versus the PyAntiGen project
"""

__version__ = "0.1.0"
