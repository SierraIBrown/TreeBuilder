//Imports
import "https://cdnjs.cloudflare.com/ajax/libs/underscore.js/1.13.6/underscore-min.js";
import "/js/jquery.min.js";
import "/js/d3.v6.min.js";
import "/js/patristic.js";
import "/js/new_phylotree.js";

//Load in phylotree.css
(function loadCSS(){
	const l = document.createElement("link");
	l.rel = "stylesheet";
	l.href = "/css/phylotree.css";
	document.head.appendChild(l);
})();

/**
 * Gets rid of d3__namespace.mouse not a function error
 */
if(!d3.mouse){
	d3.mouse = function(container){
		return d3.pointer(event, container);
	};
}

//Constant for making PDF
const { jsPDF } = window.jspdf;

/**
 * Parses either a FASTA or JSON file into {name, sequence} format
 * 
 * @param string text - The content of the FASTA or JSON file as a string
 * @param string type - The type of file (FASTA or JSON)
 * 
 * @returns {Array<Object>} sequences - An array of objects with `name` and `sequence` properties
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
	console.log("Sequences from file:", sequences);
	return sequences;
}

/**
 * Needleman-Wunsch algorithm to align sequences
 * 
 * @param string seqA - First sequence
 * @param string seqB - Second sequence
 * 
 * @returns string [seqA, seqB] - Aligned sequences
 */
export function needlemanWunsch(seqA, seqB){
	const match = 1, mismatch = -1, gap = -2;

	//Match string case
	seqA = seqA.toUpperCase();
	seqB = seqB.toUpperCase();

	let width = seqA.length + 1;
	let height = seqB.length + 1;

	let scoreMatrix = [];
	let directionMatrix = [];

	//Fill in edges & init matrix
	for(var i = 0; i < width; i++){
		scoreMatrix[i] = [];
		scoreMatrix[i][0] = i * gap;
		directionMatrix[i] = [];
		directionMatrix[i][0] = "t";
	}
	for(var j = 0; j < height; j++){
		scoreMatrix[0][j] = j * gap;
		directionMatrix[0][j] = "l";
	}

	//Create score matrix
	for(var i = 1; i < width; i++){
		for(var j = 1; j < height; j++){
			var diagPenalty = mismatch;
			if(seqA[i-1] == seqB[j-1]) diagPenalty = match;
			var diagScore = scoreMatrix[i-1][j-1] + diagPenalty;
			var insertScore = scoreMatrix[i-1][j] + gap;
			var deleteScore = scoreMatrix[i][j-1] + gap;

			//Set max value
			scoreMatrix[i][j] = Math.max(diagScore, insertScore, deleteScore);
			if(diagScore == scoreMatrix[i][j]) directionMatrix[i][j] = "d";
			else if(insertScore == scoreMatrix[i][j]) directionMatrix[i][j] = "t";
			if(deleteScore == scoreMatrix[i][j]) directionMatrix[i][j] = "l";
		}
	}

	//Traceback to get alignment
	var i = width - 1;
	var j = height - 1;
	var alnSeqA = "";
	var alnSeqB = "";

	do{
		var direction = directionMatrix[i][j];
		if(direction == "d"){
			alnSeqA = seqA[i - 1] + alnSeqA;
			alnSeqB = seqB[j - 1] + alnSeqB;
			i--;
			j--;
		}
		else if(direction == "l"){
			alnSeqA = "-" + alnSeqA;
			alnSeqB = seqB[j - 1] + alnSeqB;
			j--;
		}
		else if(direction == "t"){
			alnSeqA = seqA[i - 1] + alnSeqA;
			alnSeqB = "-" + alnSeqB;
			i--;
		}
	} while(i + j > 0);

	return [alnSeqA, alnSeqB];
}

/**
 * Computes hamming distance between two sequences
 * 
 * @param string seqA - First sequence
 * @param string seqB - Second sequence
 * 
 * @returns number distance - Hamming distance
 */
export function hammingDistance(seqA, seqB){
	let distance = 0;
	let length = Math.max(seqA.length, seqB.length);
	for(let i = 0; i < length; i++){
		if(seqA[i] !== seqB[i]){
			distance++;
		}
	}
	return distance;
}

/**
 * Computes a simple distance matrix based on sequence differences
 * 
 * @param {Array<Object>} sequences - Aligned sequences
 * 
 * @returns {Array<Array<number>>} matrix - Distance matrix
 */
export function computeDistanceMatrix(sequences){
	let length = sequences.length;
	let matrix = Array(length).fill().map(() => Array(length).fill(0));

	for(let i = 0; i < length; i++){
		for(let j = i+1; j < length; j++){
			let seqA = sequences[i].sequence;
			let seqB = sequences[j].sequence;
			let [alnSeqA, alnSeqB] = needlemanWunsch(seqA, seqB);
			let dist = hammingDistance(alnSeqA, alnSeqB);
			matrix[i][j] = dist;
			matrix[j][i] = dist;
		}
	}
	return matrix;
}

/**
 * Builds a phylogenetic tree from aligned sequences
 *
 * @param {Array<Object>} alignedSequences - The aligned sequences
 */
export function buildTree(sequences){
	//Compute distance matrix
	let distanceMatrix = computeDistanceMatrix(sequences);
	let distClone = JSON.parse(JSON.stringify(distanceMatrix));

	//Generate taxa names
	let treeTaxa = sequences.map(seq => seq.name);

	//Convert distance matrix to hierarchical tree
	let tree = window.patristic.parseMatrix(distClone, treeTaxa);
	console.log("Generated Newick Tree:", tree.toNewick());

	window._tipNames = treeTaxa;
	//window._tipColorScale = d3.scaleSequenctial(d3.interpolateRainbow).domain([0, treeTaxa.length - 1]);
	window._leafCount = sequences.length;

	//Render the tree with D3
	renderTree(tree.toNewick());
}

/**
 * Renders the tree
 *
 * @param {} phyloTree - Formatted tree
 */
export function renderTree(phyloTree){
	//Clear out old tree
	const container = document.getElementById("tree-container");
	container.innerHTML = "";
	container.style.border = "1px solid #ccc";
	container.style.margin = "5px";
	container.style.width = "750";

	const leafSpacing = 10;
	const height = Math.max(800, window._leafCount * leafSpacing);
	const width = 700;

	const tree = new window.phylotree.phylotree(phyloTree);

	//Create a TreeRender, passing the parsed pt and draw options
	const render = tree.render({
		container: "#tree-container",
		height: height,
		width: width,
		"left-right-spacing": "fit-to-size",
		"top-bottom-spacing": "fit-to-size",
		"brush": false,
		"hide": false,
		"reroot": true,
		"show-scale": true,
		"show-menu": true,
		"draw-size-bubbles": false,
		"collapsible": false,
		"selectable": true,
		"zoom": false,
		"align-tips": false
		//"node-styler": colorNodesByName,
		//"edge-styler": colorEdgesByTarget
	});

	tree.resortChildren((a, b) => {
		return (b.height - a.height) || (b.value - a.value);
	});

	container.appendChild(render.show());
	treeActions(tree, render, zoomBehavior);

	window._initial = d3.select(".phylotree-container").attr("transform");
	createLegend();
}

/**
 * Handles zooming in/out and panning
 */
const zoomBehavior = d3.zoom()
	.scaleExtent([0.5, 5])
	.on("zoom", (event) => {
		d3.select(".phylotree-container")
			.attr("transform", event.transform);
	});
//The initial window when the tree is first rendered
let initial = d3.zoomIdentity;


/**
 * Colors nodes by name
 *
 * Not currently used, but could edit to color by type
 * For example: By country/state, by id, etc.
 *
 * @param selection - The selection of node(s)
 * @param d - The id of the node(s)
 */
export function colorNodesByName(selection, d){
	const idx = window._tipNames.indexOf(d.data.name);
	selection.style("fill", window._tipColorScale(idx));
}

/**
 * Colors edges by target
 * Also not currently used, but colors edges connecting to nodes
 * in the same color
 *
 * @param selection - The selection of edge(s)
 * @param d - the id of edge(s)
 */
export function colorEdgesByTarget(selection, d){
	selection.style("stroke-width", "3px");
	const idx = window._tipNames.indexOf(d.target.data.name);
	if(idx >= 0){
		selection.style("stroke", window._tipColorScale(idx));
	}
}

/**
 * Creates the legend
 */
export function createLegend(){
	const names = window._tipNames;
	//const color = window._tipColorScale;

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

/**
 * Highlights taxa either from hovering over node names in the legend
 * or by searching for it in the search bar
 *
 * When using search bar it will automatically begin highlighting
 * each matching letter
 *
 * Hovering into the legend box will clear all previous highlights
 * from the search bar
 *
 * @param name - Name of the node being searched
 */
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

/**
 * Resets the zoom back to the initial position the tree was in
 */
export function resetZoom(){
	d3.select(".phylotree-container")
		.transition().duration(300)
		.attr("transform", window._initial);
}

/**
 * Handles all d3.select() actions including:
 * - zoomBehavior: Handles zooming in/out and panning
 * - tooltip: Hovering over nodes shows the name of node
 * - reroot: A reroot button will appear when a node is clicked
 *
 * @param {*} tree
 * @param {*} render
 * @param {*} zoomBehavior
 */
function treeActions(tree, render, zoomBehavior){
	const svg = d3.select("#tree-container svg");

	svg.call(zoomBehavior);

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
                                .html(`<strong>${d.data.name}</strong>`)
                                .style("left", (event.pageX + 5) + "px")
                                .style("top", (event.pageY + 5) + "px")
                })
                .on("mouseout", () => {
                        tooltip.style("display", "none");
                });

	d3.select("#tree-container svg").selectAll("g.internal-node")
        .on("click", (event, d) => {
                nodeDropdownMenu.call(render, d, "#tree-container", tree, render.options, event);
        });
}

/**
 * A simplified version of nodeDropdownMenu made by @stevenweaver
 * in menus.js from github.com/veg/phylotree.js
 *
 * @param {*} node
 * @param {*} container
 * @param {*} phylotree
 * @param {*} options
 * @param {*} event
 */
function nodeDropdownMenu(node, container, phylotree, options, event){
	let d3_layout_phylotree_context_menu_id = "d3_layout_phylotree_context_menu";
	let menu_object = d3.select(container).select("#" + d3_layout_phylotree_context_menu_id);

	if(menu_object.empty()){
		menu_object = d3
		.select(container)
		.append("div")
		.attr("id", d3_layout_phylotree_context_menu_id)
		.attr("class", "dropdown-menu")
		.attr("role", "menu");
	}
	
	menu_object.selectAll("button").remove();
	menu_object.selectAll("h6").remove();
	menu_object.selectAll("div").remove();

	if(!node || !options["show-menu"] || !_.some([Boolean(node.menu_items), options["hide"], options["selectable"], options["collapsible"]])){
		menu_object.style("display", "none");
		return;
	}
	if(options["reroot"]){
		menu_object
		.append("button")
		.attr("class", "dropdown-item")
		.attr("tabindex", "-1")
		.text("Reroot on this node")
		.on("click", d => {
			menu_object.style("display", "none");
			this.phylotree.reroot(node);
			this.update();
			treeActions(phylotree, this, zoomBehavior);
		});
	}

	let tree_container = document.querySelector(container);
	let rect = tree_container.getBoundingClientRect();

	menu_object
	.style("position", "absolute")
	.style("left", "" + (event.clientX - rect.x + 12) + "px")
	.style("top", "" + (event.clientY - rect.y) + "px")
	.style("display", "block");
}

/**
 * Load in phylotree.css to inject into serialized SVG
 * 
 * @returns response
 */
async function loadCSS(){
	const resp = await fetch("/css/phylotree.css");
	return await resp.text();
}

/**
 * Serializes the SVG with injected CSS styling
 * 
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

/**
 * Converts serialized SVG into a PNG
 * 
 * @returns PNG
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

/**
 * Grabs the PNG, adds and saves it to a new PDF
 * 
 * @params string filename - The name of the downloadable pdf
 */
export async function downloadTreeAsPDF(filename = "tree.pdf"){
	const svg = document.querySelector("#tree-container svg");
	if(!svg){
		console.error("No SVG found in #tree-container");
		return;
	}
	try{
		const imgData = await svgToPng(svg);
		const pdf = new jsPDF({ unit: "px" });
		const pageW = pdf.internal.pageSize.getWidth();

		pdf.setFontSize(20);

		const logo = [
			{ text: "ISU", color: [255, 0, 0], style: "bold" },
			{ text: " PRRS", color: [255, 0, 0], style: "bolditalic" },
			{ text: "View", color: [255, 215, 0], style: "bolditalic" }
		];

		let totalW = 0;
		for(let part of logo){
			totalW += pdf.getTextWidth(part.text);
		}

		const margin = 10;
		const y = 30;

		let x = pageW - margin - totalW;

		for(let part of logo){
			const style = part.style;
			pdf.setFont("helvetica", style);
			pdf.setTextColor(...part.color);
			pdf.text(part.text, x, y);
			x += pdf.getTextWidth(part.text);
		}

		pdf.addImage(imgData, "PNG", 15, 60, (svg.getBBox().width/1.65), (svg.getBBox().height/1.5));
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
			let format;
			if(file.name.toLowerCase().endsWith(".json")){
				format = "json";
			}
			else{
				format = "fasta";
			}
			const sequences = parseSequences(text, format);
			buildTree(sequences);
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
