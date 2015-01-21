Using CMSSW_7_2_3 at LPC cluster

TWiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideCMSDataAnalysisSchool2015MonoPhotonDM

Plotting codes are stored in the pyplotter and are written using py root.

The tdrstyle.py and CMS_lumi.py files are only slightly modified from the recommended CMS style plotting code which is here: https://ghm.web.cern.ch/ghm/plots/. 

The plot_functions.py file has somee simple functions to take care of some of the basics like making a canvas and setting the titles etc. 

The make_plot.py file is the one to run interactively if you want to make a single plot. You can run "python make_plot.py -h" to get help info on what arguments you can pass to it. If you want to try it out, copy the file data_selection.root from my home directory on the lpc, i.e. cp ~kdlong/data_selection.root 

A few examples: ./pyplotter/make_plot.py -n ~/data_selection.root -v phoEt[selectedPhoton] --xlabel "Photon p_{T} (GeV)" --rebin 40 --printCMS right --is_data -o test.pdf gives this:  http://www.hep.wisc.edu/~kdlong/test.pdf

And ./pyplotter/make_plot.py -n ~/data_selection.root -v phoEt[selectedPhoton] --xlabel "Photon p_{T} (GeV)" --rebin 40 --printCMS right -o test2.pdf this http://www.hep.wisc.edu/~kdlong/test2.pdf

The cut_flow.py file should be ready to reproduce the cut flow plot I made for the presentation. I'm not certain what the status of the stack plot ones are or where the files they ran on are stored, but hopefully they can still be helpful.
