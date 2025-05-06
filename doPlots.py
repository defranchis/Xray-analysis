import matplotlib.pyplot as plt
import numpy as np
import os

import samples_LabView as sLV

from utility import safeOpenFile, safeGetObject


getAllSamples = {
    "floating": sLV.getAllFloatingMOS,
    "biased": sLV.getAllBiasedMOS,
    "GCD": sLV.getAllGCD
}

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
    return result

def readAllGraphsFromFile(sample_type):
    samples = getAllSamples[sample_type]()
    graph_dict = {sample: readGraphFromFile(sample, sample_type) for sample in samples}
    np_graph_dict = {sample: TGraphToNumpy(graph, errors=True if sample_type == 'GCD' else False) for sample, graph in graph_dict.items()}
    return np_graph_dict


def plotGraphs(np_graph_dict, sample_type):
    plt.figure(figsize=(10, 6))
    for sample, data in np_graph_dict.items():
        x = data["x"]
        y = data["y"]
        ex = data["ex"]
        ey = data["ey"]
        if max(x) > 100: continue #hardcoded for now
        #if sample_type == "floating" and max(x) != 100: continue
        if "(rep.)" in sLV.sample_tags[sample]["tag"]: continue
        plt.errorbar(x, y, xerr=ex, yerr=ey, label=sLV.sample_tags[sample]["tag"])
    plt.xlabel("X-axis label")
    plt.ylabel("Y-axis label")
    plt.title(f"Comparison of {sample_type} samples")
    plt.legend()
    plt.xscale("log")
    if sample_type == "GCD":
        plt.yscale("log")
    plt.xlim(left=1)
    os.makedirs("out_plots", exist_ok=True)
    plt.savefig(f"out_plots/comparison_LabView_{sample_type}.png")
    plt.close()


def doComparisonPlot(sample_type):
    if not sample_type in getAllSamples.keys():
        raise ValueError(f"Invalid sample type: {sample_type}. Choose from {list(getAllSamples.keys())}")
    np_graph_dict = readAllGraphsFromFile(sample_type)
    plotGraphs(np_graph_dict, sample_type)





def main():
    for sample_type in getAllSamples.keys():
        print(f"Processing sample type: {sample_type}")
        doComparisonPlot(sample_type)
    pass


if __name__ == "__main__":
    main()
