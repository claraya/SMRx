#!/usr/bin/env python
# This is a PyMol script that colors atoms/residues in a PDB structure using an input file and a target column.
# Note: 	Input must contain a target values column, position numbers, and a position adjustment factor if necessary.
# Usage: 	First, load the script within PyMol as follows:  
#			
#			run /Users/claraya/Projects/wwMDs/python/mapColor.py
#
#			Execute coloring by specifying the input file, position column, target value column, transformation (if any), 
#			position adjustment factor, and colors:
#			
#			mapColor("/Users/claraya/Projects/wwMDs/data/structure/m2/kd/m2_testing_kd/mapstructure_m2_testing_kd_modeled_full.txt", mode="raw", color="samba", position="position", IDs="reference", target="function.delta.med", adjust=1, IDs="identity")

from pymol import cmd, stored
import os, math, numpy
import matplotlib as mpl
import bisect

""" define empty-value/neutral colors """
emptyColor = [0.85,0.90,0.90,1]

""" define dictionary of custom color ramps... """
colorDict = {

	# hand-made colors, in a rush:
	"citric"		: ["#03AB11", "#FFF301"],
	"citrus"		: ["#15E602", "#FFFF33"],
	"berry"			: ["#7700FF","#FF0080"],
	"forest"		: ["#E0E0E0","#6BDA61", "#15E602"],
	"white:mango"	: ["#FFFFFF", "#FF00FF", "#FF0000"],
	"white:tango"	: ["#FFFFFF", "#FF00FF", "#9500FF"],
	"white:grove"	: ["#FFFFFF", "#FFF301", "#BAE61D", "#01DF29"],
	"white:jungle"	: ["#FFFFFF", "#FFF301","#03AB11"],
	"white:orange"	: ["#FFFFFF", "#FFF301","#FF7700"],
	"horizon"		: ["#000033", "#000075", "#0000B6", "#0000F8", "#2E00FF", "#6100FF", "#9408F7", "#C729D6", "#FA4AB5", "#FF6A95", "#FF8B74", "#FFAC53", "#FFCD32", "#FFEE11", "#FFFF60"],
	
	# hand-made colors, from R:
	"wolfgang.basic": ["#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"],
	"wolfgang.extra": ["#FFFFFF", "#FCFED3", "#E3F4B1", "#ABDEB6", "#60C1BF", "#2A9EC1", "#206AAD", "#243996", "#081D58"],
	"solar.flare"	: ["#3361A5", "#2884E7", "#1BA7FF", "#76CEFF", "#FFFFFF", "#FFE060", "#FA8E24", "#DA2828", "#A31D1D"],
	"solar.glare"	: ["#3361A5", "#2884E7", "#1BA7FF", "#76CEFF", "#FCFCFC", "#FFE060", "#FA8E24", "#DA2828", "#A31D1D"],
	"solar.basic"	: ["#214B85", "#1873CC", "#1E90FF", "#00BFFF", "#ACD8E5", "#D2D2D2", "#FFD700", "#ED2C2C", "#A31D1D"],
	"solar.extra"	: ["#3361A5", "#248AF3", "#14B3FF", "#88CEEF", "#C1D5DC", "#EAD397", "#FDB31A", "#E42A2A", "#A31D1D"],
	"solar.blues"	: ["#FCFCFC", "#C0E4FD", "#75CEFE", "#0CB9FF", "#1BA7FF", "#1E95FF", "#2884E7", "#3072C5", "#3361A5"],
	"solar.rojos"	: ["#FCFCFC", "#FFEDB0", "#FFDF5F", "#FEC510", "#FA8E24", "#F14C2B", "#DA2828", "#BE2222", "#A31D1D"],
	"samba.color"	: ["#1B85ED", "#1AA2F6", "#00BFFF", "#4AC596", "#00CC00", "#A7D400", "#FFD700", "#FFBE00", "#FFA500"],
	"samba.night"	: ["#1873CC", "#1798E5", "#00BFFF", "#4AC596", "#00CC00", "#A2E700", "#FFFF00", "#FFD200", "#FFA500"],
	"samba.light"	: ["#00D1FF","#03AB11","#FFF301"],
	
	# hand-made colors, from show:
	"dusk:dawn"		: ["#98ABC5", "#8D91AD", "#827896", "#775F80", "#6B476B", "#93575B", "#B8684A", "#DB7933", "#FF8C00"],
	"dark:cyan"		: ["#000000", "#0E2824", "#014C44", "#0A5F4F", "#13725A", "#19997F", "#1EC0A6", "#19DFD2", "#00FFFF"],
	"dark:blue"		: ["#000000", "#00171F", "#002F3F", "#00475F", "#005F7F", "#00779F", "#008FBF", "#00A7DF", "#00BFFF"],
	"dark:citrus"	: ["#000000", "#22350F", "#3B680C", "#529111", "#6ABB15", "#74DD0F", "#7FFF00", "#ADF121", "#D1E131"],
	"dark:violet"	: ["#000000", "#1E0A35", "#31016A", "#4B0181", "#660099", "#7800CA", "#8A00FF", "#C800FF", "#FE00FF"],
	"ocean:green"	: ["#07519B", "#2975B4", "#5097C9", "#93C1DF", "#FCFCFC", "#CAEAC5", "#97D494", "#5BAB5A", "#006400"],
	"ocean:earth"	: ["#0F3341", "#1563AA", "#0B99E6", "#3DCDFD", "#F7F7F7", "#B87350", "#872E1C", "#601622", "#401C2A"],
	"ocean:brick"	: ["#0F3341", "#1563AA", "#0B99E6", "#3DCDFD", "#F7F7F7", "#EB9457", "#D1551F", "#B02F1B", "#8D1616"],
	"algae:earth"	: ["#543005", "#985D12", "#CFA154", "#F0DEB1", "#F5F5F5", "#B5E2DC", "#5AB2A8", "#0E726A", "#003C30"],
	"flame.flame"	: ["#000033", "#0000A5", "#1E00FB", "#6F00FD", "#C628D6", "#FE629D", "#FF9B64", "#FFD52C", "#FFFF5F"],
	"flame.light"	: ["#000033", "#000E92", "#1300FF", "#8E0EEA", "#C628D6", "#E9699F", "#FF9B63", "#FFCE62", "#FFFF5F"],
	"flame.polar"	: ["#C628D6", "#8E0EEA", "#1300FF", "#000E92", "#000033", "#7F494D", "#FF9B63", "#FFCE62", "#FFFF5F"],
	"flame.volts"	: ["#000000", "#371377", "#5F00FF", "#9400FF", "#BE00FF", "#E000EB", "#FF00D8", "#FF0090", "#FF004B"],
	"flame.watts"	: ["#FFFFFF", "#C190FF", "#5F00FF", "#9400FF", "#BE00FF", "#E000EB", "#FF00D8", "#FF0090", "#FF004B"],
	"flame.artic"	: ["#000000", "#371377", "#5F00FF", "#BD00EC", "#FF00D8", "#C7ACEC", "#00FFFF", "#0AD7D3", "#0DB2AA"],
	"flame.weird"	: ["#00FFFF", "#0AD7D3", "#0DB2AA", "#1C5551", "#000000", "#371377", "#5F00FF", "#BD00EC", "#FF00D8"],
	"flame.blind"	: ["#0DB2AA", "#0AD7D3", "#00FFFF", "#B1FFFE", "#FFFFFF", "#FFA3EC", "#FF00D8", "#BD00EC", "#5F00FF"],
	"flame.macaw"	: ["#000000", "#28410F", "#477C0E", "#64B114", "#9FCF23", "#C9E553", "#81F7D0", "#16DCD2", "#1AA58C"],
	"flame.wings"	: ["#D1E131", "#85C51D", "#529111", "#2F4E0F", "#000000", "#0F4338", "#107E6A", "#1BBBA7", "#00FFFF"],
	
	# following colors are generated online (http://www.pixelfor.me/crc/):
	"calma.azules"	: ["#031C25", "#093B4D", "#1C5F77", "#3685A2", "#56A6C3", "#86C2D8", "#B6DDEB", "#F2FBFE"],
	"calma.musgos"	: ["#212503", "#444D09", "#6B771C", "#93A236", "#B4C356", "#CDD886", "#E4EBB6", "#FCFEF2"],
	"calma.bosque"	: ["#032506", "#094D0E", "#1C7722", "#36A23D", "#56C35D", "#86D88B", "#B6EBBA", "#F2FEF3"],
	"calma.marino"	: ["#032515", "#094D2D", "#1C774D", "#36A26F", "#56C390", "#86D8B2", "#B6EBD2", "#F2FEF8"],
	"calma.morado"	: ["#030925", "#09154D", "#1C2B77", "#3648A2", "#5668C3", "#8694D8", "#B6BFEB", "#F2F4FE"],
	"calma.manudo"	: ["#290303", "#590707", "#8C1616", "#BE2A2A", "#DF4A4A", "#ED8080", "#F7B4B4", "#FFEEEE"],
	"china.theory"	: ["#120324", "#420A4A", "#721D57", "#9B3850", "#BC6B58", "#D3B687", "#E6E8B7", "#F8FDF2"],
	"china.ranges"	: ["#031424", "#1F0A4A", "#721D64", "#9B3838", "#BCAB58", "#A0D387", "#B7E8CF", "#F2F9FD"],
	"china.weirdo"	: ["#04032E", "#2E0267", "#890BA3", "#DE15AF", "#FF347E", "#FF7772", "#FFCFAB", "#FFFBEA"],
	"china.basics"	: ["#25032E", "#670253", "#A30B48", "#DE1515", "#FF8534", "#FFE272", "#EEFFAB", "#F2FFEA"],
	"china.sunset"	: ["#031124", "#0C0A4A", "#451D72", "#91389B", "#BC589B", "#D38799", "#E8C1B7", "#FDF9F2"],
	"china.dragon"	: ["#03032A", "#2B065C", "#801491", "#C52696", "#E74671", "#F2917D", "#FADEB3", "#FEFFED"],
	"china.novice"	: ["#2A0E03", "#5C4406", "#7E9114", "#68C526", "#46E748", "#7DF2B2", "#B3FAF2", "#EDF9FF"],
	
	# following colors are from R ColorBrewer:
	"brewer.fire"	: ["#FFFFE5", "#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#993404", "#662506"],
	"brewer.heat"	: ["#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"],
	"brewer.orange"	: ["#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"],
	"brewer.red"	: ["#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"],
	"brewer.green"	: ["#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"],
	"brewer.blue"	: ["#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"],
	"brewer.purple"	: ["#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D"],
	"brewer.violet"	: ["#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A"],
	"brewer.jamaica": ["#006837", "#2DA154", "#86CB66", "#CCE982", "#FFFFBF", "#FDD380", "#F88D51", "#DE3F2E", "#A50026"],
	"brewer.marine"	: ["#F7FCF0", "#E0F3DB", "#CCEBC5", "#A8DDB5", "#7BCCC4", "#4EB3D3", "#2B8CBE", "#0868AC", "#084081"],
	"brewer.spectra": ["#5E4FA2", "#3F96B7", "#88CFA4", "#D7EF9B", "#FFFFBF", "#FDD380", "#F88D51", "#DC494C", "#9E0142"],
	"brewer.celsius": ["#313695", "#5083BB", "#8FC3DD", "#D2ECF4", "#FFFFBF", "#FDD384", "#F88D51", "#DE3F2E", "#A50026"],
	"brewer.yes"	: ["#053061", "#2971B1", "#6AACD0", "#C1DDEB", "#F7F7F7", "#FACDB5", "#E58267", "#BB2933", "#67001F"],
	
	
	# following colors generated from http://www.colorhexa.com/examples
	"forest:yellow" : ["#215e33", "#306835", "#3f7136", "#4e7b38", "#5d843a", "#6c8e3b", "#7b983d", "#89a13f", "#98ab41", "#a7b542", "#b6be44", "#c5c846", "#d4d147", "#e3db49"],
	"forest:citric" : ["#0b3310", "#1b4210", "#2b5111", "#3b6111", "#4c7012", "#5c7f12", "#6c8e12", "#7c9e13", "#8cad13", "#9cbc13", "#adcb14", "#bddb14", "#cdea15", "#ddf915"],
	"citric:yellow" : ["#21be45", "#30c246", "#40c546", "#4fc947", "#5fcc48", "#6ed048", "#7dd349", "#8dd74a", "#9cda4b", "#abde4b", "#bbe14c", "#cae54d", "#dae84d", "#e9ec4e"],
	"ocean:citrus"  : ["#3683ba", "#418bb0", "#4c93a7", "#569b9d", "#61a393", "#6cab8a", "#77b380", "#81ba76", "#8cc26c", "#97ca63", "#a2d259", "#acda4f", "#b7e246", "#c2ea3c"],
	"ocean:pink"	: ["#3b5e84", "#4a5e84", "#595e84", "#685e84", "#765e84", "#855e84", "#945e84", "#a35f85", "#b25f85", "#c15f85", "#cf5f85", "#de5f85", "#ed5f85", "#fc5f85"],
	"ocean:red" 	: ["#3b82ae", "#487aa1", "#547294", "#616987", "#6d617a", "#7a596d", "#865160", "#934853", "#9f4046", "#ac3839", "#b8302c", "#c5271f", "#d11f12", "#de1705"],
	"ocean:aqua"	: ["#2668aa", "#2d6eaa", "#3574aa", "#3c7aaa", "#4380a9", "#4b86a9", "#528ca9", "#5993a9", "#6099a9", "#689fa9", "#6fa5a8", "#76aba8", "#7eb1a8", "#85b7a8"],
	"ocean:teal"	: ["#0a0a66", "#15176c", "#202573", "#2b3279", "#364080", "#414d86", "#4c5a8d", "#576893", "#62759a", "#6d82a0", "#7890a7", "#839dad", "#8eabb4", "#99b8ba"],
	"cyan:brick"	: ["#6cd0c2", "#70c3b3", "#74b6a5", "#78a996", "#7b9c88", "#7f8f79", "#83826a", "#87755c", "#8b684d", "#8f5b3e", "#924e30", "#964121", "#9a3413", "#9e2704"],
	"aqua:brick"	: ["#019bcf", "#0d94c2", "#1a8cb4", "#2685a7", "#337d99", "#3f768c", "#4c6e7e", "#586771", "#655f63", "#715856", "#7e5048", "#8a493b", "#97412d", "#a33a20"],
	"aqua:tan"  	: ["#62adbb", "#6cb0b5", "#75b2af", "#7fb5a8", "#88b7a2", "#92ba9c", "#9bbd96", "#a5bf8f", "#aec289", "#b8c583", "#c1c77d", "#cbca76", "#d4cc70", "#decf6a"],
	"cyan:tan"  	: ["#4fe8c2", "#5be3bb", "#67deb5", "#73d9ae", "#7fd3a7", "#8bcea1", "#97c99a", "#a4c493", "#b0bf8c", "#bcba86", "#c8b47f", "#d4af78", "#e0aa72", "#eca56b"],
	"teal:orange"	: ["#0cb499", "#1dae92", "#2ea98b", "#3fa384", "#4f9e7d", "#609876", "#71936f", "#828d69", "#938862", "#a4825b", "#b47d54", "#c5774d", "#d67246", "#e76c3f"],
	"teal:violet"	: ["#4b8b84", "#508184", "#567784", "#5b6d84", "#616384", "#665984", "#6b4f84", "#714683", "#763c83", "#7b3283", "#812883", "#861e83", "#8c1483", "#910a83"],
	"blue:cyan" 	: ["#4111f2", "#4923f0", "#5135ed", "#5946eb", "#6258e8", "#6a6ae6", "#727ce3", "#7a8de1", "#829fde", "#8ab1dc", "#93c3d9", "#9bd4d7", "#a3e6d4", "#abf8d2"],
	"purple:pink"	: ["#6848d1", "#7345ca", "#7f42c3", "#8a3fbc", "#953db5", "#a13aae", "#ac37a7", "#b734a1", "#c2319a", "#ce2e93", "#d92c8c", "#e42985", "#f0267e", "#fb2377"],
	"purple:baby"	: ["#511293", "#5a239b", "#6234a2", "#6b45aa", "#7356b2", "#7c67b9", "#8478c1", "#8d88c9", "#9599d1", "#9eaad8", "#a6bbe0", "#afcce8", "#b7ddef", "#c0eef7"],
	"cyan:green"	: ["#29ddea", "#31dcd8", "#39dcc7", "#41dbb5", "#4adba3", "#52da92", "#5ad980", "#62d96e", "#6ad85c", "#72d74b", "#7bd739", "#83d627", "#8bd616", "#93d504"],
	"cyan:pink" 	: ["#606bef", "#6c6ae6", "#7869dd", "#8468d4", "#9168cb", "#9d67c2", "#a966b9", "#b565b1", "#c164a8", "#cd639f", "#da6396", "#e6628d", "#f26184", "#fe607b"],
	"cyan:violet"	: ["#29e9ae", "#36dbae", "#43cdae", "#50beaf", "#5eb0af", "#6ba2af", "#7894af", "#8585b0", "#9277b0", "#9f69b0", "#ad5bb0", "#ba4cb1", "#c73eb1", "#d430b1"],
	"cyan:purple"	: [ "#4adbf1", "#51cff0", "#59c3ef", "#60b7ee", "#68abed", "#6f9fec", "#7793eb", "#7e88eb", "#867cea", "#8d70e9"]
	}


""" define a function to construct a path """
def pathGenerator(inpath):
	if not os.path.isdir(inpath):
		os.makedirs(inpath)


""" define  a function to generate numerical sequences and that tolerates floats """
def floatRange(start, stop, steps):
	step = float(stop - start)/steps
	output = list()
	while start < stop:
		output.append(start)
		start += float(step)
	return output


""" define a function that finds fractions of values greater, lesser, within, or without ranges... """
def findValues(values, cutoff, mode, report="fraction"):
	values = sorted(values)
	hits = list()
	
	if mode == "less":
		# Find rightmost value less than x 
		i = bisect.bisect_left(values, cutoff)
		if i:
			hits = values[:i]
		
	elif mode == "less.equal":
		# Find rightmost value less than or equal to x 
		i = bisect.bisect_right(values, cutoff)
		if i:
			hits = values[:i]
		
	elif mode == "more":
		# Find leftmost value greater than x 
		i = bisect.bisect_right(values, cutoff)
		if i != len(values):
			hits = values[i:]
		
	elif mode == "more.equal":
		# Find leftmost item greater than or equal to x 
		i = bisect.bisect_left(values, cutoff)
		if i != len(values):
			hits = values[i:]
			
	elif mode == "equal":
		# Find leftmost item greater than or equal to x 
		hits = list()
		for value in values:
			if value == cutoff:
				hits.append(value)
		
	# Report the matching values:
	if report == "fraction":
		return float(len(hits))/len(values)
	elif report == "tally":
		return len(hits)
	elif report == "values":
		return hits


""" define a function to build a dictionary of input i,x values... """
def quickBuilder(infile, i, x, header="ON", separator="\t", mode="float"):
	outDict = dict()
 	inlines = open(infile).readlines()
	if header == "ON":
		headerDict = dict()
		columns = inlines.pop(0).strip().split(separator)
		index = 0
		for column in columns:
			if column == i:
				headerDict[column] = index
			elif column == x:
				headerDict[column] = index
			index += 1
	for inline in inlines:
		initems = inline.strip().split(separator)
		if header == "ON":
			position, value = initems[headerDict[i]], initems[headerDict[x]]
		else:
			position, value = initems[i], initems[x]
		if mode == "float":
			if value == "NA":
				value = 0
			outDict[int(position)] = float(value)
		else:
			outDict[int(position)] = value
	return outDict


""" define a function to color PDB structures from within PyMol... """
def mapColor(infile, mode, color="wolfgang.v1", reverse="OFF", position="position", target="value", adjust=0, select="OFF", IDs="OFF", maxCut="OFF", minCut="OFF", maxValue="OFF", minValue="OFF", colorDict=colorDict, altColor="OFF", dpi=300, ray=1, N=256, save="OFF", NA="OFF"):
	
	"""
	infile	:	Path to value input file.
	mode	:	How should input values be treated? Options are "raw", "normalize" and "log2".
	color	:	Color ramp to be used for value mapping.
	reverse	:	Reverse color ramp for value mapping.
	position:	Position column in input file.
	target	:	Target value column to map to each position.
	adjust	:	Integer describing how many residues into the chain to begin coloring.
	IDs		:	Should residue identities (chemicals) be read for each position? If so, specify identity column.
	maxCut	:	Maximum value (cutoff) allowed for redefined value range.
	minCut	:	Minimum value (cutoff) allowed for redefined value range.
	maxValue:	Maximum value for high-range normalization.
	minValue:	Minimum value for high-range normalization.
	dpi		:	PyMol resolution; dots per inch.
	ray		:	PyMol rendering mode.
	N		:	Color ramp size.
	"""
	
	# Load color map. Note that colors come inverted:
	colorMix = list(colorDict[color])
	if reverse == "OFF":
		colorMix.reverse()
	colorMap = mpl.colors.LinearSegmentedColormap.from_list(colorMix, colors=colorMix, N=N)
	color256 = colorMap(numpy.arange(N))
	
	#for index in range(0, 255):
	#	r, g, b, f = color256[index]
	#	print index, ":", "	".join(map(str, [round(r, 2), round(g, 2), round(b, 2)]))
	#print

	print
	print "ColorMap:", color
	print "Colors:", len(color256)
	print
	
	# load input file and intensity values:
	rawDict = quickBuilder(infile, i=position, x=target, mode="float")
	positions = sorted(rawDict.keys())
	
	# load reference identities if necessary:
	if IDs != "OFF":
		idDict = quickBuilder(infile, i=position, x=IDs, mode="string")
	
	# threshold the raw (input) values, if necessary:
	rawValues, rawPositions, rawIDs = list(), list(), list()
	for position in positions:
		rawValue = float(rawDict[position])
		if maxCut != "OFF":
			rawValue = min(rawValue, maxCut)
		if minCut != "OFF":
			rawValue = max(rawValue, minCut)
		rawValues.append(rawValue)
		rawPositions.append(position)
		
		# load IDs if specified:
		if IDs != "OFF":
			rawIDs.append(idDict[position])
	
	# normalize and transform values, if necessary:
	colorValues = list()
	if mode == "raw":
		if minValue == "OFF" and maxValue == "OFF":
			colorValues = rawValues
		else:
			for value in rawValues:
				if minValue != "OFF" and value < float(minValue):
					colorValues.append(minValue)
				elif maxValue != "OFF" and value > float(maxValue):
					colorValues.append(maxValue)
				else:
					colorValues.append(value)
			
	elif mode == "normalize":
		for value in rawValues:
			if value < 0 and minValue == "OFF":
				colorValues.append(-1*value/min(rawValues))
			elif value < 0 and minValue != "OFF":
				colorValues.append(-1*value/float(minValue))
			elif value > 0 and maxValue == "OFF":
				colorValues.append(value/max(rawValues))
			elif value > 0 and maxValue != "OFF":
				colorValues.append(value/float(maxValue))
			elif value == 0:
				colorValues.append(0)
	
	elif mode == "log2":
		log2Values = [ math.log(value, 2) for value in rawValues ]
		for value in log2Values:
			if value < 0:
				colorValues.append(-1*value/min(log2Values))
			elif value > 0:
				colorValues.append(value/max(log2Values))
			elif value == 0:
				colorValues.append(0) 
	
	elif mode == "log10":
		
		logXValues = [ numpy.log10(value) for value in rawValues ]
		
		for value in logXValues:
			if value < 0 and not value == -float('Inf'):
				colorValues.append(-1*value/min(logXValues))
			elif value > 0 and not value == float('Inf'):
				colorValues.append(value/max(logXValues))
			elif value == 0:
				colorValues.append(0)
			elif value in [float('Inf'), -float('Inf'), float('NaN')]:
				colorValues.append(NA)
				print value, NA
	
	else:
		sys.exit("Error: choose a valid mode")
	
	# generate complete range of values:
	minColor, maxColor = min(colorValues), max(colorValues)
	if minValue != "OFF":
		minColor = float(minValue)
	if maxValue != "OFF":
		maxColor = float(maxValue)
	colorRanges = floatRange(minColor, maxColor, steps=N)
	
	#[x for x in colorRanges if x != 'nan']
	
	print
	print "Input values (min, max):", min(colorValues), "-", max(colorValues)
	print "Range values (min, max):", min(colorRanges), "-", max(colorRanges)
	print
	
	# color residues:
	fraction = 0.1
	for index in range(0, len(colorValues)):
	
		# define residue names:
		if select == "OFF":
			colorName = "res" + str(int(index) + adjust)
		else:
			colorName = "%s and chain %s" % tuple(select.split(",")) + " and resi " + str(int(index) + adjust)
		fraction = findValues(colorRanges, float(colorValues[index]), mode="less.equal", report="fraction")
		
		# define residue color:
		if altColor == "OFF":
			r, g, b, f = color256[int(fraction*(N-1))]
			cmd.set_color(str(colorName), str([r, g, b]))
			cmd.color(colorName, colorName)
			output = [int(index) + adjust, rawPositions[index], rawValues[index], fraction, colorName]
			
		# resets residue color:
		if altColor != "OFF" and rawValues[index] in altColor:
			r, g, b, f = altColor[rawValues[index]]
			cmd.set_color(str(colorName), str([r, g, b]))
			cmd.color(colorName, colorName)
			output = [int(index) + adjust, rawPositions[index], rawValues[index], fraction, colorName]
			
		# record IDs, if requested:
		if IDs != "OFF":
			output.append(rawIDs[index])
		#print "  ".join(map(str, output))
	#print len(color256)
	
	# save image:
	if save != "OFF":
		print save
		cmd.png((save), dpi=dpi, ray=ray)


""" define a function to generate colors for PyMol... """
def genColor(color="wolfgang.v1", reverse="OFF", colorDict=colorDict, dpi=300, ray=1, N=256):
	
	"""
	color	:	Color ramp to be used for value mapping.
	reverse	:	Reverse color ramp for value mapping.
	maxCut	:	Maximum value (cutoff) allowed for redefined value range.
	minCut	:	Minimum value (cutoff) allowed for redefined value range.
	dpi		:	PyMol resolution; dots per inch.
	ray		:	PyMol rendering mode.
	N		:	Color ramp size.
	"""
	
	# Load color map. Note that colors come inverted:
	colorMix = list(colorDict[color])
	if reverse == "OFF":
		colorMix.reverse()
	colorMap = mpl.colors.LinearSegmentedColormap.from_list(colorMix, colors=colorMix, N=N)
	color256 = colorMap(numpy.arange(N))
	
	print
	print "ColorMap:", color
	print "Colors:", len(color256)
	print
	
	# color residues:
	k = 1
	for i in color256:
		colorName = color + "." + str(k)
		r, g, b, f = color256[k]
		cmd.set_color(colorName, str([r, g, b]))
	
cmd.extend("mapColor", mapColor)

# Visualization example: PIK3CA-PIK3R1 UCEC ENST00000263967 (2RD0)

#run /Users/claraya/Desktop/Dropbox/Jungla/Code/python/mapColor.py

#mapColor("/Users/claraya/Desktop/Dropbox/Jungla/Code/data/density/hs/cancer/lawrence2014/all/all/snp/tapanti/aminos/DBSCAN/500/010/UCEC/PIK3CA/mapvariant_lawrence2014_UCEC_ENST00000263967.txt", mode="raw", altColor="OFF", color="samba.color", reverse="ON", position="position", IDs="reference", target="identity.ratio", adjust=1, select="2RD0,A", maxCut=0.005)

#mapColor("/Users/claraya/Desktop/Dropbox/Jungla/Code/data/density/hs/cancer/lawrence2014/all/all/snp/tapanti/aminos/DBSCAN/500/010/UCEC/PIK3CA/mapvariant_lawrence2014_UCEC_ENST00000263967.txt", mode="raw", altColor="OFF", color="samba.color", reverse="ON", position="position", IDs="reference", target="normal.ids", adjust=1, select="2RD0,A", maxCut=0.01)
