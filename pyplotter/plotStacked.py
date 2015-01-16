#!/usr/bin/env python
import plot_functions as plotter
import ROOT as ROOT

def main():
    plot_info = getWeightArgs()
    histVar = ""
    unstackedFiles = ""
    stackedFiles = ""
    hist_opts = "hist"
    hist_stack = ROOT.THStack()
     
    cols = [ROOT.kOrange-8, ROOT.kBlue-4, ROOT.kBlue+8]

    i = 0
    for weight in weights:
        weight_plot_info["folder"] = weight
        print cols[i]
        if plotterameCanvas: 
            plotter.addHistToStack(hist_stack, weight_plot_info, hist_opts, 
                                 cols[i], 0)          
        else:
            hist = plotter.getHistFromFile(weight_plot_info)
            plotter.setHistAttributes(hist, weight_plot_info, ROOT.kYellow-i, 
                                    ROOT.kYellow+i)
            plot_info["output_file"].replace(".pdf", i + ".pdf")
            plotter.makePlot(hist, hist_opts, plot_info) 
        i += 1
    if plotSameCanvas:
        plotter.makePlot(hist_stack, "hist nostack", plot_info)
def getWeightArgs():
    parser = plotter.getBasicParser()
    parser.add_argument("-n", "--file_name", type=str, required=True,
                        help="Name of root file where plots are stored")
    return vars(parser.parse_args())
if __name__ == "__main__":
    main()
