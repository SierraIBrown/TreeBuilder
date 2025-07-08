//Load in phylotree.css
(function loadCSS(){
        const l = document.createElement("link");
        l.rel = "stylesheet";
        l.href = "/css/phylotree.css";
        document.head.appendChild(l);
})();

import "/js/d3.v6.min.js";
import "https://cdnjs.cloudflare.com/ajax/libs/underscore.js/1.13.6/underscore-min.js";
import "/js/patristic.js";
import "/js/new_phylotree.js";

const { jsPDF } = window.jspdf;

//Helper function
(function(){
        if(typeof window.___namespace === 'undefined'){
                window.___namespace = {
                        extend(target, source){
                                for(const key in source){
                                        if(Object.prototype.hasOwnProperty.call(source, key)){
                                                target[key] = source[key];
                                        }
                                }
                                return target;
                        },
                        extendOwn(target, source){
                                for(const key in source){
                                        if(Object.prototype.hasOwnProperty.call(source, key)){
                                                target[key] = source[key];
                                        }
                                }
                                return target;
                        },
                        noop(){}
                };
        }
})();

/*
 * Parses either a FASTA or JSON file into {name, sequence} format
 * @param String text - The content of the FASTA or JSON file as a string
 * @param String type - The type of file (FASTA or JSON)
 * @returns {Array<Object>} An array of objects with `name` and `sequence` properties
 */
export function parseSequences(text, type){
        //JSON path
        if(type === "json"){
                try{
                        const data = JSON.parse(text);
                        return data.map(entry => ({
                                name: entry.name || "Unknown",
                                sequence: entry.sequence.replace(/\s+/g, '')
                        }));
                }
                catch(error){
                        console.error("Error parsing JSON:", error);
                        return [];
                }
        }

        //FASTA path
        const sequences = [];
        const lines = text.split("\n");
        let currentName = null;
        let currentSequence = "";

        lines.forEach(line => {
                line = line.trim();
                                if(!line) return;
                if(line.startsWith(">")){
                        if(currentName != null){
                                sequences.push({
                                        name: currentName,
                                        sequence: currentSequence.replace(/\s+/g, '')
                                });
                        }
                        currentName = line.substring(1).trim();
                        currentSequence = "";
                }
                else{
                        currentSequence += line;
                }
        });

        if(currentName !== null){
                sequences.push({
                        name: currentName,
                        sequence: currentSequence.replace(/\s+/g, '')
                });
        }

        return sequences;
}

/*
 * Aligns sequences
 * @param {Array<Object>} sequences - Array of { name, sequence } objects
 * @returns {Array<Object>} Aligned sequences in { name, sequence } format
 */
export function alignSequences(sequences){
        if(sequences.length < 2){
                console.error("Not enough sequences to align");
                return sequences;
        }

        //Select the longest sequence as a reference
        let reference = sequences.reduce((a, b) => (a.sequence.length > b.sequence.length ? a : b));
        let alignedSequences = [];

        console.log("Using reference:", reference.name);

        sequences.forEach(seq => {
                let aligned = aligning(reference.sequence, seq.sequence);
                alignedSequences.push({ name: seq.name, sequence: aligned});
        });

        return alignedSequences;
}

/*
 * Aligning method (Needleman-Wunsch)
 * @param {string} refSeq - Reference Sequence
 * @param {string} querySeq - Sequence to align
 * @returns {string} Aligned sequence
 */
export function aligning(refSeq, querySeq){
        const match = 1, mismatch = -1, gap = -2;

        let rows = querySeq.length + 1;
        let cols = refSeq.length + 1;
        let matrix = Array.from({ length: rows }, () => Array(cols).fill(0));

        //Initialize scoring matrix
        for(let i = 1; i < rows; i++) matrix[i][0] = i * gap;
        for(let j = 1; j < cols; j++) matrix[0][j] = j * gap;
        
        //Fill the matrix
        for(let i = 1; i < rows; i++){
                for(let j = 1; j < cols; j++){
                        let matchScore = matrix[i - 1][j - 1] + (querySeq[i - 1] === refSeq[j - 1] ? match :       mismatch);
                        let gap1 = matrix[i - 1][j] + gap;
                        let gap2 = matrix[i][j - 1] + gap;
                        matrix[i][j] = Math.max(matchScore, gap1, gap2);
                }
        }

        //Traceback
        let alignedQuery = "";
        let alignedRef = "";
        let i = rows - 1, j = cols - 1;
        while(i > 0 && j > 0){
                let score = matrix[i][j];
                let scoreDiag = matrix[i - 1][j - 1];
                let scoreUp = matrix[i - 1][j];
                let scoreLeft = matrix[i][j - 1];

                if(score === scoreDiag + (querySeq[i - 1] === refSeq[j - 1] ? match : mismatch)){
                        alignedQuery = querySeq[i - 1] + alignedQuery;
                        alignedRef = refSeq[j - 1] + alignedRef;
                        i--;
                        j--;
                }
                else if(score === scoreUp + gap){
                        alignedQuery = querySeq[i - 1] + alignedQuery;
                        alignedRef = "-" + alignedRef;
                        i--;
                }
                else{
                        alignedQuery = "-" + alignedQuery;
                        alignedRef = refSeq[j - 1] + alignedRef;
                        j--;
                }
        }

        //Fill in the remaining gaps
        while(i > 0){
                alignedQuery = querySeq[i - 1] + alignedQuery;
                alignedRef = "-" + alignedRef;
                i--;
        }
        while(j > 0){
                alignedQuery = "-" + alignedQuery;
                alignedRef = refSeq[j - 1] + alignedRef;
                j--;
        }

        return alignedQuery;
}

/*
 * Computes hamming distance between two sequences
 * @param String seq1 - First sequence
 * @param String seq2 - Second sequence
 * @returns number - Hamming distance
 */
export function hammingDistance(seq1, seq2){
        let distance = 0;
        let length = Math.max(seq1.length, seq2.length);
        for(let i = 0; i < length; i++){
                if(seq1[i] !== seq2[i]){
                        distance++;
                }
        }
                return distance;
}

/*
 * Computes a simple distance matrix based on sequence differences
 * @param Array<Object> sequences - Aligned sequences
 * @returns Array<Array<number>> Distance matrix
 */
export function computeDistanceMatrix(sequences){
        let matrix = [];
        for(let i = 0; i < sequences.length; i++){
                matrix[i] = [];
                for(let j = 0; j < sequences.length; j++){
                        if(i === j){
                                matrix[i][j] = 0;
                        }
                        else {
                                matrix[i][j] = hammingDistance(sequences[i].sequence, sequences[j].sequence);
                        }
                }
        }

        return matrix;
}

const zoomBehavior = d3.zoom()
        .scaleExtent([0.5, 5])
        .on("zoom", ({transform}) => {
                d3.select(".phylotree-container")
                        .attr("transform", transform);
        });

let initial = d3.zoomIdentity;

export function colorNodesByName(selection, d){
        const idx = window._tipNames.indexOf(d.data.name);
        selection.style("fill", window._tipColorScale(idx));
}

export function colorEdgesByTarget(selection, d){
        selection.style("stroke-width", "3px");
        const idx = window._tipNames.indexOf(d.target.data.name);
        if(idx >= 0){
                selection.style("stroke", window._tipColorScale(idx));
        }
}

export function createLegend(){
        const names = window._tipNames;
        const color = window._tipColorScale;

        const legend = d3.select(".controls #legend")
                .style("display", "flex")
                .style("flex-wrap", "wrap")
                .style("overflow-x", "auto")
                .style("white-space", "nowrap")
                .style("gap", "8px")
                .style("margin-top", "0.5em");
        legend.html("");

        const items = legend.selectAll(".legend-item")
                .data(names)
                .enter().append("div")
                        .attr("class", "legend-item")
                        .style("display", "flex")
                        .style("align-items", "center")
                        .style("font-size", "0.8em")
                        .style("cursor", "pointer")
                                                .style("margin-right", "5px")
                        .on("mouseover", (_, name) => highlightTaxa(name))
                        .on("mouseout", (_, name) => highlightTaxa(""));

        items.append("div")
                .style("width", "12px")
                .style("height", "12px")
                .style("margin-right", "4px")
                .style("background-color", "#ccc");

        items.append("span")
                .text(d => d);
}

export function highlightTaxa(name){
        d3.selectAll(".hilite").remove();
        d3.selectAll(".node").classed("hilitebold", false);
        let q;
        if (typeof name === "string"){
                q = name.toLowerCase().trim();
        }
        else {
                q = document.getElementById("searchInput").value.toLowerCase().trim();
        }
        if(!q) return;
        d3.selectAll(".node").each(function(d) {
                if(d.data.name.toLowerCase().includes(q)){
                        const g = d3.select(this);
                        const box = this.getBBox();
                        g.classed("hilitebold", true)
                        .insert("rect", "text")
                        .attr("class", "hilite")
                        .attr("x", box.x - 4)
                        .attr("y", box.y - 2)
                        .attr("width", box.width + 8)
                        .attr("height", box.height + 4)
                        .style("fill", "yellow")
                        .style("stroke", "yellow")
                        .style("opacity", 0.6);
                }
        });
}

export function resetZoom(){
        d3.select(".phylotree-container")
                .transition().duration(300)
                .attr("transform", window._initial);
}

/*
 * Builds a phylogenetic tree from aligned sequences
 * @param Array<Object> alignedSequences - The aligned sequences
 */
export function buildTree(alignedSequences){
        //Generate taxa names
        let treeTaxa = alignedSequences.map(seq => seq.name);
        console.log("Tree taxa:", treeTaxa);

        //Compute the distance matrix
        let distanceMatrix = computeDistanceMatrix(alignedSequences);
        console.log("Distance:", distanceMatrix);

        //Convert distance matrix to hierarchical tree
        let tree = window.patristic.parseMatrix(distanceMatrix, treeTaxa);
        console.log("Generated Newick Tree:", tree.toNewick());

        window._tipNames = treeTaxa;
        window._tipColorScale = d3.scaleSequential(d3.interpolateRainbow).domain([0, treeTaxa.length - 1]);
                window._leafCount = alignedSequences.length;

        //Render the tree with D3
        renderTree(tree.toNewick());
        createLegend();
}

/*
 * Renders the tree
 * @param {} phyloTree - Formatted tree
 */
export function renderTree(phyloTree){
        //Clear out old tree
        const container = document.getElementById("tree-container");
        container.innerHTML = "";
        container.style.border = "1px solid #ccc";
        container.style.margin = "5px";
        container.style.width = "950";

        const leafSpacing = 10;
        const height = Math.max(800, window._leafCount*leafSpacing);
        const width = 900;

        //Parse Newick into a Phylotree object
        const tree = new window.phylotree.phylotree(phyloTree);

        //Create a TreeRender, passing the parsed pt and your draw options
        const render = tree.render({
                container: "#tree-container",
                height: height,
                width: width,
                "left-right-spacing": "fit-to-size",
                "top-bottom-spacing": "fit-to-size",
                "brush": false,
                "hide": false,
                "reroot": false,
                "show-scale": false,
                "draw-size-bubbles": false,
                "collapsible": false,
                "selectable": true,
                "zoom": true,
                "align-tips": false
                /*"node-styler": colorNodesByName,
                "edge-styler": colorEdgesByTarget*/
        });

        tree.resortChildren((a,b) => {
                return (b.height - a.height) || (b.value - a.value);
        });

        container.appendChild(render.show());

        d3.select("#tree-container svg").call(zoomBehavior);
        window._initial = d3.select(".phylotree-container").attr("transform");
        createLegend();

        let tooltip = d3.select("body").select("div.tooltip");
        if(tooltip.empty()){
                tooltip = d3.select("body").append("div")
                .attr("class", "tooltip")
                .style("position", "absolute")
                .style("pointer-events", "none")
                .style("background", "#fff")
                .style("border", "1px solid #333")
                .style("display", "none");
        }

        d3.select("#tree-container svg")
                .selectAll("text.phylotree-node-text")
        .on("mouseover", (event, d) => {
                const total = d.y.toFixed(2);
                tooltip
                .style("display", "block")
                .html(`<strong>${d.data.name}</strong><br/>Length: ${total}`)
                .style("left", (event.pageX + 5) + "px")
                .style("top", (event.pageY + 5) + "px")
        })
        .on("mouseout", () => {
                tooltip.style("display", "none");
        });
        console.log("Phylotree Container: ", d3.select(".phylotree-container"));
        console.log("SVG from renderTree: ", d3.select("#tree-container svg"));
}

/*
 * Load in phylotree.css to inject into serialized SVG
 * @returns response
 */
async function loadCSS(){
        const resp = await fetch("/css/phylotree.css");
        return await resp.text();
}

/*
 * Serializes the SVG with injected CSS styling
 * @returns string Serialized SVG
 */
async function serialize(origSvg){
        const svg = origSvg.cloneNode(true);
        const { width, height } = origSvg.getBoundingClientRect();
        svg.setAttribute("width", Math.ceil(width));
        svg.setAttribute("height", Math.ceil(height));
        svg.setAttribute("viewBox", `0 0 ${width} ${height}`);

        const css = await loadCSS();

        const styleText = document.createElementNS("http://www.w3.org/2000/svg", "style");
        styleText.setAttribute("type", "text/css");
        styleText.textContent = css;
        svg.insertBefore(styleText, svg.firstChild);

        const s = new XMLSerializer().serializeToString(svg);
        console.log("SVG s: ", s);
        return `<?xml version="1.0" standalone="no"?>\n${s}`;
}

/*
 * Converts serialized SVG into a PNG
 * @returns png
 */
function svgToPng(svg){
        return serialize(svg).then(svgString => {
                const blob = new Blob([svgString], { type: "image/svg+xml" });
                const url  = URL.createObjectURL(blob);
                return new Promise((resolve, reject) => {
                        const img = new Image();
                        img.onload = () => {
                        const c = document.createElement("canvas");
                        c.width  = img.width;  c.height = img.height;
                        c.getContext("2d").drawImage(img, 0, 0);
                        URL.revokeObjectURL(url);
                        resolve(c.toDataURL("image/png"));
                        };
                        img.onerror = e => { URL.revokeObjectURL(url); reject(e); };
                        img.src = url;
                });
                        });
}

/*
 * Grabs the PNG, adds and saves it to a new PDF
 * @params string filename The name of the downloadable pdf
 */
export async function downloadTreeAsPDF(filename = "tree.pdf"){
        const svg = document.querySelector("#tree-container svg");
        if(!svg){
                console.error("No SVG found in #tree-container");
                return;
        }

        try{
                const imgData = await svgToPng(svg);

                const pdf = new jsPDF({
                        unit: "px",
                        format: [svg.getBBox().width, svg.getBBox().height]
                });

                pdf.addImage(imgData, "PNG", 0, 0);
                pdf.save(filename);
        }
        catch(err){
                console.error("Export failed:", err);
        }
}

document.addEventListener("DOMContentLoaded", () => {
        const controls = document.querySelector(".controls");
        if(!controls) return;
        controls.classList.add("controls");

        const fileGroup = document.createElement("div");
        fileGroup.classList.add("fileGroup");

        const fileInput = document.createElement("input");
        fileInput.type = "file";
        fileInput.accept = ".fasta,.json";
        fileInput.classList.add("fileInput");

        const runBtn = document.createElement("button");
        runBtn.textContent = "Align & Build Tree";
        runBtn.classList.add("runBtn");

        fileGroup.append(fileInput, runBtn);

        const navGroup = document.createElement("div");
        navGroup.classList.add("navGroup");

        const searchInput = document.createElement("input");
        searchInput.id = "searchInput";
        searchInput.placeholder = "Search Taxa";
        searchInput.classList.add("searchInput");

        const zoomGroup = document.createElement("div");
        zoomGroup.classList.add("zoomGroup");

        const zoomInBtn = document.createElement("button");
        zoomInBtn.textContent = "Zoom In";
        zoomInBtn.onclick = () => zoomIn();
        zoomInBtn.classList.add("zoomInBtn");

        const zoomOutBtn = document.createElement("button");
        zoomOutBtn.textContent = "Zoom Out";
        zoomOutBtn.onclick = () => zoomOut();
                zoomOutBtn.classList.add("zoomOutBtn");

        const resetBtn = document.createElement("button");
        resetBtn.textContent = "Reset Zoom";
        resetBtn.onclick = () => resetZoom();
        resetBtn.classList.add("resetBtn");

        zoomGroup.append(zoomInBtn, zoomOutBtn, resetBtn);
        navGroup.append(searchInput, zoomGroup);

        const legend = document.createElement("div");
        legend.id = "legend";
        legend.classList.add("legend");

        controls.append(fileGroup, navGroup, legend);

        searchInput.addEventListener("input", highlightTaxa);
        runBtn.addEventListener("click", () => {
                if(!fileInput.files.length){
                        output.textContent = "Please upload a FASTA or JSON file.";
                        return;
                }
                const file = fileInput.files[0];
                const reader = new FileReader();
                reader.onload = e => {
                        const text = e.target.result;
                        const type = file.name.toLowerCase().endsWith(".json") ? "json" : "fasta";
                        const sequences = parseSequences(text, type);
                        const aligned = alignSequences(sequences);
                        buildTree(aligned);
                };
                reader.readAsText(file);
        });

        const downloadBtn = document.getElementById("download-btn");
        if(downloadBtn){
                downloadBtn.addEventListener("click", () => downloadTreeAsPDF());
        }
});

window.highlightTaxa = highlightTaxa;
window.zoomIn = () => d3.select("#tree-container svg").transition().call(zoomBehavior.scaleBy, 2);
window.zoomOut = () => d3.select("#tree-container svg").transition().call(zoomBehavior.scaleBy, 0.5);
window.resetZoom = resetZoom;
