import ROOT
import CMS_lumi, tdrstyle

# what to draw
#label = "MET"
#drawString = "pfMET"
label = "#gamma E_{T}"
drawString = "phoEt[selectedPhoton]"
cutString = ""
xrange = [100, 1100]
yrange = [0, 100]
fill = True
normalize = False
# qcd znunu wJet sig
colors = [ROOT.kGreen-5, ROOT.kBlue-6, ROOT.kGray, ROOT.kRed]

tdrstyle.setTDRStyle()

fdata = ROOT.TFile("../data_selection.root")
fqcd = ROOT.TFile("../qcd_selection.root")
fznunu = ROOT.TFile("../znunu_selection.root")
fwJets = ROOT.TFile("../wJets_selection.root")
fsig = ROOT.TFile("../signal_selection.root")

data = fdata.Get("EventTree")
qcd = fqcd.Get("EventTree")
znunu = fznunu.Get("EventTree")
wJets = fwJets.Get("EventTree")
signal = fsig.Get("EventTree")

canvas = ROOT.TCanvas("canvas", "canvas")
stack = ROOT.THStack("stack", "")

hdata = ROOT.TH1F("data", "Data", 20, xrange[0], xrange[1])
hqcd = ROOT.TH1F("qcd", "QCD", 20, xrange[0], xrange[1])
hznunu = ROOT.TH1F("znunu", "Z#nu#nu", 20, xrange[0], xrange[1])
hwJets = ROOT.TH1F("wJets", "W+Jets", 20, xrange[0], xrange[1])
hsignal = ROOT.TH1F("signal", "Signal", 20, xrange[0], xrange[1])

data.Draw(drawString+">>data", cutString, "goff")
hdata.SetMarkerStyle(ROOT.kFullCircle)

qcd.Draw(drawString+">>qcd", cutString, "goff")
if ( fill ) : hqcd.SetFillColor(colors[0])
hqcd.SetLineColor(colors[0])

znunu.Draw(drawString+">>znunu", cutString, "goff")
if ( fill ) : hznunu.SetFillColor(colors[1])
hznunu.SetLineColor(colors[1])

wJets.Draw(drawString+">>wJets", cutString, "goff")
if ( fill ) : hwJets.SetFillColor(colors[2])
hwJets.SetLineColor(colors[2])

signal.Draw(drawString+">>signal", cutString, "goff")
hsignal.SetLineColor(colors[3])
hsignal.SetLineStyle(ROOT.kDashed)

# rescale MC
if ( normalize ) :
  hqcd.Scale(1/hqcd.Integral())
  hznunu.Scale(1/hznunu.Integral())
  hwJets.Scale(1/hwJets.Integral())
  hsignal.Scale(1/hsignal.Integral())
  hdata.Sumw2()
  hdata.Scale(1/hdata.Integral())

else :
  hqcd.Scale(13/hqcd.Integral())
  hznunu.Scale(193/hznunu.Integral())
  hwJets.Scale(5/hwJets.Integral())
  hsignal.Scale(60/hsignal.Integral())

bgDraw = ""
stack.Add(hqcd, bgDraw)
stack.Add(hwJets, bgDraw)
stack.Add(hznunu, bgDraw)
stack.Add(hsignal, bgDraw)

stack.SetMinimum(yrange[0])
stack.SetMaximum(yrange[1])
#stack.Draw("nostack")
stack.Draw()
hdata.Draw("E same")

# Now, make pretty
stack.GetXaxis().SetTitle(label)
stack.GetYaxis().SetTitle("Events / 50 GeV")
legend = ROOT.TLegend(0.7, 0.7, 0.95, 0.85)
legend.AddEntry(hdata, "Data", "ep")
legend.AddEntry(hqcd, "QCD", "f")
legend.AddEntry(hznunu, "Z#nu#nu", "f")
legend.AddEntry(hwJets, "W+Jets", "f")
legend.AddEntry(hsignal, "SM+ADD(n_{ED}=2, M_{D}=3)", "l")
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.Draw()

CMS_lumi.lumi_13TeV = "%s fb^{-1}" % 1
iPos = 13
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.CMS_lumi(canvas, 4, iPos)
