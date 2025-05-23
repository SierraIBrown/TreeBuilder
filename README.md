# TreeBuilder
A small browser module for parsing FASTA/JSON sequences, performing alignments, and building & rendering a tree with some interactive tools
## Usage
Load in the module using ```<script type="module" src="./treeBuilder.js"></script>```
## API
### 1. ```parseSequences(String text, String type)```
- Parses either a FASTA or JSON file into {name, sequence} format
### 2. ```alignSequences(Array<Object> sequences)```
- Runs pairwise Needleman-Wunsch alignment of every sequence to the longest one
### 3. ```aligning(String refSeq, String querySeq)```
- Internal Needleman-Wunsch implementation
### 4. ```hammingDistance(String seq1, String seq2)```
- Computes the hamming distance between two sequences
### 5. ```computeDistanceMatrix(Array<Object> sequences)```
- Computes a simple distance matrix based on sequence differences
### 6. ```colorNodesByName(d3.Selection selection, NodeData d)```
- Color-styler for tree nodes using a global window._tipColorScale
### 7. ```colorEdgesByTarget(d3.Selection selection, EdgeData d)```
- Color-styler for edges matching each leaf's color
### 8. ```highlightTaxa()```
- Reads the ```#searchInput``` value, finds all ```.node``` whose name contains it, and inserts a yellow highlight behind the text
### 9. ```resetZoom()```
- Restores the initial ```<g>``` transform saved after first render
### 10. ```buildTree(Array<Object> alignedSequences)```
1. Generates a taxa list
2. Computes distance matrix
3. Builds a Newick tree via ```patristic.parseMatrix```
4. Stores global tip names & colors
5. Calls ```renderTree()```
### 11. ```renderTree(String phyloTree)```
1. Parses Newick into a phylotree
2. Configures and calls ```.render()```
3. Appends to ```#tree-container```
4. Hooks up D3 zoom, tooltips, and initial transform
