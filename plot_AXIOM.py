import matplotlib.pyplot as plt
import numpy as np
import os

import mplhep
plt.style.use(mplhep.style.CMS)
plt.rcParams["text.usetex"] = True

from utility import safeOpenFile, safeGetObject
import sampleStyle as sAXIOM

from doPlots import TGraphToNumpy, addPrimaryLabel, addSecondaryLabel, getPlotLabel

sample_types = {'floating' : 'MOShalf',
                'biased' : 'MOS2000',
                'GCD' :'GCD',}

sample_types_tags = {'floating' : 'MOShalf',
                'biased' : 'MOS2000',
                'GCD' :'GCD',}

indir = "./out_AXIOM/raw_summary/"
out_dir = "out_plots/AXIOM"
os.makedirs(out_dir, exist_ok=True)

comparison_types = ['types', 'materials', 'all']

all_plots_dict = {sample_type: {} for sample_type in sample_types.keys()}


def addAXIOMlabel(pos_x=0.05, pos_y=0.87, v_align='top', h_align='left'):
    plt.text(pos_x, pos_y, r"\textit{AXIOM} setup", transform=plt.gca().transAxes, fontsize=25,
             verticalalignment=v_align, horizontalalignment=h_align)

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
        
        tag = sAXIOM.getSampleAttribute(sample, "leg").replace('#', 'n')
        if 'kGy' in tag: continue
        if 'n2' in tag: continue
        if not comparisonTypePass(sample, comparison_type):
            continue

        x = data["x"]
        y = data["y"]
        ex = data["ex"]
        ey = data["ey"]
        ey = ey if sample_type == "GCD" else np.zeros(len(y))

        legend = sAXIOM.tag_styles[tag]["new tag"]
        color, marker = sAXIOM.getColorAndMarker(sample)
        linestyle = '-' if not 'EPI' in legend else '--'
        if legend == "New type C FZ (Aug23)":
            linestyle = '--'
        plt.errorbar(x, y, xerr=None, yerr= ey if sample_type == "GCD" else None, label=legend, color=color, marker=marker, linestyle = linestyle)

    plt.xlabel("Dose [kGy]")
    plt.ylabel("Surface velocity [cm/s]" if sample_type == "GCD" else "Oxide charge density [cm$^{-2}$]")
    plt.legend(handlelength=1.5, loc = 'lower right')
    plt.xscale("log")
    plt.grid(which='both', linestyle='--', linewidth=0.5)
    if sample_type == "GCD":
        plt.yscale("log")
    plt.xlim(left=1)
    addPrimaryLabel()
    plot_label = getPlotLabel(sample_type)
    addSecondaryLabel(plot_label, h_align='left', pos_x=0.05, pos_y=0.95)
    addAXIOMlabel()
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