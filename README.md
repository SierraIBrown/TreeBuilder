# TreeBuilder
A small browser module for parsing FASTA/JSON sequences, performing alignments, and building & rendering a tree with some interactive tools<br/>
*** Limited use to 700 sequences at once ***
## Usage
### 1. Download both the css and js folder 
### 2. Load in the module using ```<script type="module" src="./treeBuilder.js"></script>```
### 3. Load in an additional script: ```<script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf/3.0.0/jspdf.umd.min.js"></script>```
- This script allows for the SVG to be downloadable
### 4. Add in elements to host the tree and controls:
- ```<div id="tree-container"></div>```
- ```<div class="controls"></div>```
## Styling
These are the ID's each div is currently using
### Tree: ```#tree-container```
- The SVG: ```#tree-container svg```
### Controls: ```.controls```
1. File Group: ```.fileGroup```:
    - File Input: ```.fileInput```
    - Align & Build Tree: ```.runBtn```
2. Nav Group: ```.navGroup```:
    - Search Taxa: ```.searchInput```
3. Zoom Group: ```.zoomGroup```:
    - Zoom In: ```.zoomInBtn```
    - Zoom Out: ```.zoomOutBtn```
    - Reset Zoom: ```.resetBtn```
4. Legend: ```.legend```:
## API
### 1. File Input
1. ```parseSequences(String text, String type)```: Parses either a FASTA or JSON file into {name, sequence} format
### 2. Alignment
1. ```alignSequences(Array<Object> sequences)```: Runs pairwise Needleman-Wunsch alignment of every sequence to the longest one
2. ```aligning(String refSeq, String querySeq)```: Internal Needleman-Wunsch implementation
3. ```hammingDistance(String seq1, String seq2)```: Computes the hamming distance between two sequences
4. ```computeDistanceMatrix(Array<Object> sequences)```: Computes a simple distance matrix based on sequence differences
### 3. Tree Rendering
1. ```colorNodesByName(d3.Selection selection, NodeData d)```: Color-styler for tree nodes using a global window._tipColorScale
2. ```colorEdgesByTarget(d3.Selection selection, EdgeData d)```: Color-styler for edges matching each leaf's color
3. ```createLegend()```: Generates a legend underneath the search/zoom bar
4. ```highlightTaxa(String name)```: Highlights the name on the tree when hovering over a name on the legend or from ```#searchInput```
5. ```resetZoom()```: Restores the initial ```<g>``` transform saved after first render
6. ```buildTree(Array<Object> alignedSequences)```: Generates a taxa list, computes distance matrix, builds a Newick tree via ```patristic.parseMatrix```, stores global tip names & colors, and calls ```renderTree()```
7. ```renderTree(String phyloTree)```: Parses Newick into a phylotree, configures and calls ```.render()```, appends to ```#tree-container```, and hooks up D3 zoom, tooltips, and initial transform
### 4. Downloading
1. ```loadCSS()```: Loads in phylotree.css
2. ```serialize(origSvg)```: Serializes the a copy of the original SVG
3. ```svgToPng(svg)```: Converts serialized SVG into a PNG
4. ```downloadTreeAsPDF(filename = "tree.pdf")```: Takes the PNG and saves it onto a new PDF
