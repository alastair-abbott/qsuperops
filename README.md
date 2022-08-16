# qsuperops
MATLAB code to work with quantum superoperators and their causal structure

The following conventions are obeyed throughout this package for specifying superoperators (i.e., process matrices), superinstruments, witness, etc.:

An operator is specified by:
- W: a dxd matrix
- dims: a list [d1, d2, ...] of dimensions satisfying prod(dims) == d
- parties: A cell array specifying which spaces belong to which party P,A,B,...,F. Each party is a 1x2 cell array, with each cell a list of indices in dims specifying which dimensions belong to that space. For P=parties{1}{2}, F=parties{end}{1}, AI=parties{2}{1}, AO=parties{2}{2}, BI=parties{3}{1}, etc.

This allows each space to be composed of several subsystems (e.g., the control and target in the quantum switch), which for convenience need not be "adjacent" in the basis chosen for W.

A superinstrument is specified by a cell array {W1,W2,...} of matrices. A singleton superinstrument {W} is equivalent to W itself.