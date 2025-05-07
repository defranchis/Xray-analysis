import matplotlib.pyplot as plt
import numpy as np
import os

import mplhep
plt.style.use(mplhep.style.CMS)
plt.rcParams["text.usetex"] = True


import samples_LabView as sLV

from utility import safeOpenFile, safeGetObject

comparison_types = ['types', 'HGCAL_std', 'detectors', 'all']
out_dir = "out_plots/LabView"
os.makedirs(out_dir, exist_ok=True)


getAllSamples = {
    "floating": sLV.getAllFloatingMOS,
    "biased": sLV.getAllBiasedMOS,
    "GCD": sLV.getAllGCD
}


def getPlotLabel(sample_type):
    if sample_type == "floating":
        return 'Floating MOS'
    elif sample_type == "biased":
        return 'Biased MOS'
    elif sample_type == "GCD":
        return 'Floating GCD'
    else:
        raise ValueError(f"Invalid sample type: {sample_type}")
    
def getNameFromType(sample, sample_type):
    if sample_type not in getAllSamples.keys():
        raise ValueError(f"Invalid sample type: {sample_type}. Choose from {list(getAllSamples.keys())}")
    list_types = sLV.sample_tags[sample]["structures"]
    list_typenames = list(getAllSamples.keys())
    index = list_typenames.index(sample_type)
    return list_types[index]


def readGraphFromFile(sample, sample_type):
    path = f"./out_LabView/Vfb_files/dose_{sample}_{getNameFromType(sample,sample_type)}.root"
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    graph_name = "gJ" if sample_type == "GCD" else "gNox"
    with safeOpenFile(path) as f:
        graph = safeGetObject(f, graph_name)
        if graph is None:
            raise ValueError(f"Graph {graph_name} not found in file {path}")
        return graph

def TGraphToNumpy(graph, errors=False):
    n = graph.GetN()
    result = {
        "x": np.frombuffer(graph.GetX(), count=n, dtype=np.double),
        "y": np.frombuffer(graph.GetY(), count=n, dtype=np.double),
        "ex": np.frombuffer(graph.GetEX(), count=n, dtype=np.double) if errors else np.zeros(n),
        "ey": np.frombuffer(graph.GetEY(), count=n, dtype=np.double) if errors else np.zeros(n)
    }
    if min(result["x"]) == 0:
        for key in result.keys():
            result[key] = result[key][1:]
        
    return result

def readAllGraphsFromFile(sample_type):
    samples = getAllSamples[sample_type]()
    graph_dict = {sample: readGraphFromFile(sample, sample_type) for sample in samples}
    np_graph_dict = {sample: TGraphToNumpy(graph, errors=True if sample_type == 'GCD' else False) for sample, graph in graph_dict.items()}
    return np_graph_dict

def getColorAndMarker(sample):
    tag = sLV.sample_tags[sample]["tag"]
    color = sLV.tag_styles[tag]["color"]
    marker = sLV.tag_styles[tag]["marker"]
    return color, marker

def comparisonTypePass(sample, comparison_type):
    sample_tag = sLV.sample_tags[sample]["tag"]
    if comparison_type == "types":
        if 'Type' in sample_tag:
            return True
        if 'FZ -5V' in sample_tag:
            return True
        return False
    elif comparison_type == "HGCAL_std":
        if 'FZ' in sample_tag:
            return True
        if 'EPI' in sample_tag:
            return True
        return False
    elif comparison_type == "detectors":
        if 'FZ -5V' in sample_tag:
            return True
        if 'Tracker' in sample_tag:
            return True
        if 'HGCAL' in sample_tag:
            return True
        if 'New Type C' in sample_tag:
            return True
        return False
    elif comparison_type == "all":
        return True
    else:
        raise ValueError(f"Invalid comparison type: {comparison_type}. Choose from {comparison_types}")


def addSecondaryLabel(label, pos_x=0.05, pos_y=0.95, v_align='top', h_align='left'):
    plt.text(pos_x, pos_y, r"\textbf{"+label+"}", transform=plt.gca().transAxes, fontsize=25,
             verticalalignment=v_align, horizontalalignment=h_align)

def addPrimaryLabel():
    plt.text(1, 1.06 , r"\textbf{HGCAL} \textit{Preliminary}", transform=plt.gca().transAxes, fontsize=25,
             verticalalignment='top', horizontalalignment='right')

def plotGraphs(np_graph_dict, sample_type,comparison_type):
    for sample, data in np_graph_dict.items():
        if not comparisonTypePass(sample, comparison_type): continue
        x = data["x"]
        y = data["y"]
        ey = data["ey"]
        if max(x) > 100: continue #hardcoded for now
        #if sample_type == "floating" and max(x) != 100: continue
        if "(rep.)" in sLV.sample_tags[sample]["tag"]: continue
        color, marker = getColorAndMarker(sample)
        plt.errorbar(x, y, xerr=None, yerr=ey if sample_type == "GCD" else None, label=sLV.sample_tags[sample]["tag"], color=color, marker=marker)
    plt.xlabel("Dose [kGy]")
    plt.ylabel("Surface velocity [cm/s]" if sample_type == "GCD" else "Oxide charge density [cm$^{-2}$]")
    plt.legend(handlelength=1.5, loc = 'lower right')
    plt.xscale("log")
    plt.grid(which='both', linestyle='--', linewidth=0.5)
    if sample_type == "GCD":
        plt.yscale("log")
    plt.xlim(left=1)
    addPrimaryLabel()
    label = getPlotLabel(sample_type)
    addSecondaryLabel(label)
    plt.savefig(os.path.join(out_dir, f"comparison_LabView_{sample_type}_{comparison_type}.pdf"))
    plt.savefig(os.path.join(out_dir, f"comparison_LabView_{sample_type}_{comparison_type}.png"))
    plt.clf()
    plt.close()


def doComparisonPlot(sample_type):
    if not sample_type in getAllSamples.keys():
        raise ValueError(f"Invalid sample type: {sample_type}. Choose from {list(getAllSamples.keys())}")
    np_graph_dict = readAllGraphsFromFile(sample_type)
    for comparison_type in comparison_types:
        print(f"Plotting {sample_type} for comparison type: {comparison_type}")
        plotGraphs(np_graph_dict, sample_type, comparison_type)
    print()


def main():
    for sample_type in getAllSamples.keys():
        print(f"Processing sample type: {sample_type}")
        doComparisonPlot(sample_type)
    pass


if __name__ == "__main__":
    main()
