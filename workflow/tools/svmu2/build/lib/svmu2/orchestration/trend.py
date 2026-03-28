'''
Docstring for orchestration.trend
mostly a debugging tool to observe cumulative sum of alignment trend
exploring this as a possibility to understand extrema of syntenic blocks in non trivial cases
such as the case of len(query) << len(reference)
'''

from svmu2.orchestration.parse import run_parse
from svmu2.core.traversal import evaluate_alignment_trend
from svmu2.models.line import LinePrimitive, build_alignment_primitives
from svmu2.visualization.renderers import render_matplotlib

import matplotlib.pyplot as plt
import numpy as np

def run_trend(args):

    all_alns, primary = run_parse(args)
    for aln in primary.values():
        trend_result = evaluate_alignment_trend(aln)
        primitives = build_alignment_primitives(aln)

        xlabel = aln.reference
        ylabel = aln.query

        # Create shared figure with 2 panels
        fig, (dotplot_ax, trend_ax) = plt.subplots(1, 2, figsize=(12, 6))

        # Left panel: dotplot
        render_matplotlib(primitives, xlabel, ylabel, ax=dotplot_ax)
        first_derivative = np.diff(aln.cumsum)
        second_derivative = np.diff(first_derivative)
        # Right panel: cumulative trend
        trend_ax.plot(second_derivative, marker='o', label="Cumulative Sum")
        trend_ax.set_xlabel("Index")
        trend_ax.set_ylabel("Value")
        trend_ax.set_title("Cumulative Sum")
        trend_ax.legend()
        trend_ax.grid(True)

        plt.tight_layout()
        plt.show()