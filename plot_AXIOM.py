import matplotlib.pyplot as plt
import numpy as np
import os

import mplhep
plt.style.use(mplhep.style.CMS)
plt.rcParams["text.usetex"] = True

from utility import safeOpenFile, safeGetObject
import sampleStyle as sAXIOM

from doPlots import TGraphToNumpy

sample_types = {'floating' : 'MOShalf',
                'biased' : 'MOS2000',
                'GCD' :'GCD',}

indir = "./out_AXIOM/raw_summary/"
out_dir = "out_plots/AXIOM"
os.makedirs(out_dir, exist_ok=True)

comparison_types = ['types', 'materials', 'all']

all_plots_dict = {sample_type: {} for sample_type in sample_types.keys()}

def readAllGraphsFromFile(sample):
    path = os.path.join(indir, f"summaryVsDose_{sample}.root")
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    with safeOpenFile(path) as f:
        for sample_type, structure in sample_types.items():
            graph_name = f"{sample}_{structure}_params"
            graph = safeGetObject(f, graph_name, quitOnFail=False, silent=True) 
            if graph is None:
                print(f"Warning: Graph {graph_name} not found in file {path}. Skipping...")
                continue
            all_plots_dict[sample_type][sample] = TGraphToNumpy(graph, errors=True if sample_type == "GCD" else False)

def readAllGraphs(samples):
    for sample in samples:
        readAllGraphsFromFile(sample)

def comparisonTypePass(sample, comparison_type):
    legend = sAXIOM.getSampleAttribute(sample, "leg")
    if comparison_type == "all":
        return True
    if comparison_type == "types":
        if 'EPI' in legend:
            return False
        return True
    if comparison_type == "materials":
        if 'C' in legend or 'D' in legend:
            if 'Aug23' in legend:
                return 'prime' in legend
            return True
        return False
    else:
        raise ValueError(f"Invalid comparison type: {comparison_type}. Choose from {comparison_types}")


def plotGraphsType(sample_type, comparison_type):
    for sample, data in all_plots_dict[sample_type].items():
        
        legend = sAXIOM.getSampleAttribute(sample, "leg").replace('#', 'n')
        if 'kGy' in legend: continue
        if 'n2' in legend: continue
        if not comparisonTypePass(sample, comparison_type):
            continue

        x = data["x"]
        y = data["y"]
        ex = data["ex"]
        ey = data["ey"]

        plt.errorbar(x, y, xerr=ex, yerr=ey, label=legend)

    plt.xlabel("Dose [kGy]")
    plt.ylabel("Surface velocity [cm/s]" if sample_type == "GCD" else "Oxide charge density [cm$^{-2}$]")
    plt.legend()
    plt.xscale("log")
    plt.grid(which='both', linestyle='--', linewidth=0.5)
    if sample_type == "GCD":
        plt.yscale("log")
    plt.xlim(left=1)
    plt.savefig(os.path.join(out_dir, f"comparison_{sample_type}_{comparison_type}.pdf"))
    plt.savefig(os.path.join(out_dir, f"comparison_{sample_type}_{comparison_type}.png"))
    plt.clf()
    plt.close()

def plotAllGraphs():
    for comparison_type in comparison_types:
        for sample_type in sample_types.keys():
            plotGraphsType(sample_type, comparison_type)

def main():
    samples = sAXIOM.sampleNames.keys()
    readAllGraphs(samples)
    plotAllGraphs()



if __name__ == "__main__":
    main()