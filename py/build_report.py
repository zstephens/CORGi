import matplotlib.pyplot as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
from bokeh.plotting import figure
from bokeh.models import Range1d, PanTool, WheelZoomTool, ResetTool, ColumnDataSource, HoverTool, Arrow, NormalHead
from bokeh.embed import components

AMBIG_BREAKPOINT = -10

def build_html(filename,refStr,script,div):
	outStr =  ''
	outStr += '<!DOCTYPE html>\n'
	outStr += '<html lang="en">\n'
	outStr += '    <head>\n'
	outStr += '        <meta charset="utf-8">\n'
	outStr += '        <title>SV Report for '+refStr+'</title>\n'
	outStr += '\n'
	outStr += '        <link rel="stylesheet" href="resources/bokeh-0.12.6.min.css" type="text/css" />\n'
	outStr += '        <script type="text/javascript" src="resources/bokeh-0.12.6.min.js"></script>\n'
	outStr += script.replace('\n','\n        ')
	outStr += '\n'
	outStr += '\n    </head>\n'
	outStr += '    <body>\n'
	outStr += '        <h1>Reference Sequence:</h1>'
	####outStr += '        <center>\n'
	outStr += div[0].replace('\n','\n        ') + '\n'
	####outStr += '        </center>\n'
	outStr += '        <h1>Observed Rearrangements:</h1>'
	for i in xrange(1,len(div)):
		outStr += div[i].replace('\n','\n        ') + '\n<br>\n'
	outStr += '    </body>\n'
	outStr += '</html>\n'

	f = open(filename,'w')
	f.write(outStr)
	f.close()

# cheat and format human chromosomes a little nicer
CHR_LIST = [str(n) for n in xrange(1,30)]+['X','x','Y','y','M','m','MT','mt']
def humanChr(c):
	if c in CHR_LIST:
		return 'chr'+c.upper()
	else:
		return c

def getColor(i,N,colormap='gnuplot'):
	cm = mpl.get_cmap(colormap) 
	cNorm  = colors.Normalize(vmin=0, vmax=N+1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
	colorVal = scalarMap.to_rgba(i+1)
	return colorVal

def gen_plots(refDat,regionList,strCounts,strNames,ambigDict,outDir):
	(ref_chr,pos_offset,pos_end) = refDat
	refStr = humanChr(ref_chr)+': {0:,} - {1:,}'.format(pos_offset,pos_end)
	TRI_HEIGHT  = 1
	PLOT_BUFFER = 200
	triDict = {regionList[i][0]:(regionList[i][2]-regionList[i][1]+1,i) for i in xrange(len(regionList))}
	N = len(triDict)

	#
	# PLOT 1: REFERENCE REGION PARTITIONS
	#
	TOOLS = [PanTool(dimensions='width'),WheelZoomTool(dimensions='width'),HoverTool(),ResetTool()]
	p1 = figure(plot_width=900, plot_height=220, tools=TOOLS, title=refStr)
	colorList = ["#%02x%02x%02x"%(int(255*getColor(i,N)[0]),int(255*getColor(i,N)[1]),int(255*getColor(i,N)[2])) for i in xrange(N)]
	triList_x = []
	triList_y = []
	nameList  = []
	arrowList = []
	x_adj     = 0
	for i in xrange(len(regionList)):
		n = regionList[i]
		triList_x.append([x_adj,x_adj,x_adj+triDict[n[0]][0]])
		triList_y.append([-TRI_HEIGHT,TRI_HEIGHT,0])
		nameList.append(n[0])
		x_adj += triDict[n[0]][0]
		if regionList[i][0] in ambigDict:
			ambigVal = ambigDict[regionList[i][0]]
			if ambigVal[1]:
				arrowList.append((x_adj+ambigVal[0],x_adj+1,i+1))
			else:
				arrowList.append((x_adj+ambigVal[0]-triDict[n[0]][0],x_adj+1-triDict[n[0]][0],i))
	pos_s_list = [n[0]+pos_offset for n in triList_x]
	pos_e_list = [n[2]+pos_offset for n in triList_x]

	source = ColumnDataSource(data=dict(x=triList_x, y=triList_y, name=nameList, color=colorList, sPos=pos_s_list, ePos=pos_e_list))

	p1.patches('x', 'y', source=source, color='color', alpha=0.7, line_color="white", line_width=0.5)

	for i in xrange(len(arrowList)):
		p1.add_layout(Arrow(x_start=arrowList[i][1],x_end=arrowList[i][0],y_start=0.8*TRI_HEIGHT,y_end=0.8*TRI_HEIGHT,end=NormalHead(fill_color=colorList[arrowList[i][2]],line_color=colorList[arrowList[i][2]],fill_alpha=0.5,line_alpha=0.5,size=8),line_width=4,line_color=colorList[arrowList[i][2]],line_alpha=0.6))

	#p1.patches(triList, [[-TRI_HEIGHT,TRI_HEIGHT,0] for n in xrange(N)], color=colorList, alpha=[1.0]*N, line_width=1, hover_alpha=[0.6]*N)
	start = triList_x[0][0]
	end   = triList_x[-1][2]
	p1.yaxis.major_label_text_color = None
	p1.xgrid.grid_line_color = None
	p1.ygrid.grid_line_color = None
	p1.x_range.bounds = (start-PLOT_BUFFER,end+PLOT_BUFFER)
	p1.y_range.bounds = (-TRI_HEIGHT-1,TRI_HEIGHT+1)
	p1.xaxis.axis_label = 'Reference position (relative)'

	p1.select_one(HoverTool).point_policy = "follow_mouse"
	p1.select_one(HoverTool).tooltips = [("Name", "@name"), ("Ref Start", "@sPos"), ("Ref End", "@ePos")]

	#
	# PLOT 2+: OBSERVED JUNCTIONS
	#
	observed_junction_plots = []
	for i in xrange(len(strCounts)):
		TOOLS = [PanTool(dimensions='width'),WheelZoomTool(dimensions='width'),HoverTool(),ResetTool()]
		observed_junction_plots.append(figure(plot_width=900, plot_height=200, tools=TOOLS, title="Multiplicity: "+str(strCounts[i][0])+", Type: "+strNames[i]))
		myColorList = []
		triList_x   = []
		triList_y   = []
		nameList    = []
		x_adj       = 0
		for m in strCounts[i][1]:
			myInd = triDict[m[2].replace('*','')]
			myColorList.append(colorList[myInd[1]])
			nameList.append(m[2])
			if m[1]:	# is forward strand
				triList_x.append([x_adj,x_adj,x_adj+myInd[0]])
				triList_y.append([-TRI_HEIGHT,TRI_HEIGHT,0])
			else:		# is reverse strand
				triList_x.append([x_adj,x_adj+myInd[0],x_adj+myInd[0]])
				triList_y.append([0,TRI_HEIGHT,-TRI_HEIGHT])
			x_adj += myInd[0]
			# add a novel sequence segment, if necessary
			if m[3][1]:
				myColorList.append('#AAAAAA')
				nameList.append('novel sequence')
				triList_x.append([x_adj,x_adj,x_adj+m[3][0]])
				triList_y.append([-TRI_HEIGHT,TRI_HEIGHT,0])
				x_adj += m[3][0]

		source = ColumnDataSource(data=dict(x=triList_x, y=triList_y, name=nameList, color=myColorList))
		observed_junction_plots[-1].patches('x', 'y', source=source, color='color', alpha=0.7, line_color="white", line_width=0.5)
		start = triList_x[0][0]
		end   = triList_x[-1][2]
		observed_junction_plots[-1].yaxis.major_label_text_color = None
		observed_junction_plots[-1].xgrid.grid_line_color = None
		observed_junction_plots[-1].ygrid.grid_line_color = None
		observed_junction_plots[-1].x_range.bounds = (start-PLOT_BUFFER,end+PLOT_BUFFER)
		observed_junction_plots[-1].y_range.bounds = (-TRI_HEIGHT-1,TRI_HEIGHT+1)
		observed_junction_plots[-1].select_one(HoverTool).point_policy = "follow_mouse"
		observed_junction_plots[-1].select_one(HoverTool).tooltips = [("Name", "@name")]

	#
	# GENERATE HTML
	#
	plots = [p1] + observed_junction_plots
	script, div = components(plots)
	build_html(outDir+'results.html',refStr,script,div)


