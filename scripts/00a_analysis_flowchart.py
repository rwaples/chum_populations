import pydot
import pyparsing

myDOT = pydot.graph_from_dot_file("/home/ipseg/Desktop/waples/chum_populations/analysis_flowchart.gv")

myDOT[0].write_png("C:/Users/IPSEG/Desktop/Waples/chum_populations/plots/mapping_flowchart.png")

myDOT[1].write_png("/home/ipseg/Desktop/waples/chum_populations/plots/flowcharts/population_analysis.png")

myDOT[2].write_png("C:/Users/IPSEG/Desktop/Waples/chum_populations/plots/chum_populations.png")




