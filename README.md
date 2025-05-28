# TreeBuilder
A small browser module for parsing FASTA/JSON sequences, performing alignments, and building & rendering a tree with some interactive tools
*** Limited use to 700 sequences at once ***
## Usage
### 1. Download both the css and js folder 
### 2. Load in the module using ```<script type="module" src="./treeBuilder.js"></script>```
### 3. Add in elements to host the tree and controls:
- ```<div id="tree-container"></div>```
- ```<div class="controls"></div>```
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
### 8. ```createLegend()```
- Generates a legend underneath the search/zoom bar
### 9. ```highlightTaxa(String name)```
- Highlights the name on the tree when hovering over a name on the legend or from ```#searchInput```
### 10. ```resetZoom()```
- Restores the initial ```<g>``` transform saved after first render
### 11. ```buildTree(Array<Object> alignedSequences)```
1. Generates a taxa list
2. Computes distance matrix
3. Builds a Newick tree via ```patristic.parseMatrix```
4. Stores global tip names & colors
5. Calls ```renderTree()```
### 12. ```renderTree(String phyloTree)```
1. Parses Newick into a phylotree
2. Configures and calls ```.render()```
3. Appends to ```#tree-container```
4. Hooks up D3 zoom, tooltips, and initial transform
## Styling
These are the ID's each div is currently using
### Tree: ```#tree-container```
- The SVG: ```#tree-container svg```
### Controls: ```.controls```
1. ```fileGroup```:
    - File Input: ```fileInput```
    - Align & Build Tree: ```runBtn```
2. ```navGroup```:
    - Search Taxa: ```searchInput```
    - Zoom In: ```zoomInBtn```
    - Zoom Out: ```zoomOutBtn```
    - Reset Zoom: ```resetBtn```
3. ```legend```:
    - Legend: ```legend```
