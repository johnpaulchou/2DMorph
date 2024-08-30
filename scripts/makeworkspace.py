import ROOT
from scipy.interpolate import interp1d
ROOT.gSystem.Load('../build/lib2DMorph.dylib')

# function to get a TH2 from a ROOT file
def get2DRooHist(filename,histname):
    rootfile = ROOT.TFile(filename, "READ")
    if not rootfile or rootfile.IsZombie():
        print("Error: Unable to open file '{file_name}'.")
        return None

    roohist = rootfile.Get(histname)
    if not roohist or not isinstance(roohist, ROOT.TH2) or roohist is None:
        print("Error: Histogram '{histname}' not found or is not a valid TH2 object in the file '{filename}'.")
        rootfile.Close()
        return None

    #    roohist.Scale(roohist2d.GetEntries()/roohist2d.Integral(0,roohist2d.GetXaxis().GetNbins()+1,0,roohist2d.GetYaxis().GetNbins()))
    #    roohist.Sumw2(False)
    roohist.SetDirectory(0)
    rootfile.Close()
    return roohist

# function to create a RooBinning out of a TAxis
def makeRooBinning(taxis, name):
    num_bins = taxis.GetNbins()
    bin_edges = [taxis.GetBinLowEdge(i) for i in range(1, num_bins + 2)]  # Includes upper edge of last bin
    bin_edges_vector = ROOT.std.vector('double')(len(bin_edges))
    for i, edge in enumerate(bin_edges):
        bin_edges_vector[i] = edge
    roobinning = ROOT.RooBinning(num_bins, bin_edges_vector.data(), name)
    return roobinning


# set the omega mass and phi mass values here
wmass=0.50
pmass=1050.
tlabel=str(pmass)+"_"+str(wmass)

# set the signal interpolation parameters here
txmin=0.5; txmax=2.1; tymin=450.; tymax=1050.
tx=ROOT.RooRealVar("tx","tx",txmin,txmax)
ty=ROOT.RooRealVar("ty","ty",tymin,tymax)
tx.setVal(wmass)
ty.setVal(pmass)
fnA="../data//mass2D_Phi_450_omega_0p5.root"
fnB="../data//mass2D_Phi_450_omega_2p1.root"
fnC="../data//mass2D_Phi_1050_omega_0p5.root"
fnD="../data//mass2D_Phi_1050_omega_2p1.root"
histname='hist_sum_1'

# interpolate the cross section
theory_xs = [(450., 585.983), (500., 353.898), (625., 117.508), (750., 45.9397), (875., 20.1308),
             (1000., 9.59447), (1125., 4.88278), (1250., 2.61745), (1375., 1.46371),
             (1500., 0.847454), (1625., 0.505322), (1750., 0.309008), (1875., 0.192939),
             (2000., 0.122826), (2125., 0.0795248), (2250., 0.0522742), (2375., 0.0348093),
             (2500., 0.0235639), (2625., 0.0161926), (2750., 0.0109283), (2875., 0.00759881)]
x=list(zip(*theory_xs))
interp=interp1d(x[0], x[1])
xsec=ROOT.RooRealVar("xsec","Interpolated theory x-section",interp(pmass))
xsec.setConstant(True)
print("The interpolated cross section is "+str(xsec.getValV())+" fb.")

# interpolate the acceptance*efficiency
# this is a dummy value for now
acceff=ROOT.RooRealVar("acceff","interpolated acceptance*efficiency",0.0001)
acceff.setConstant(True)
print("The interpolated acc*eff is "+str(acceff.getValV()))

# get the 2D histograms
hA=get2DRooHist(fnA, histname)
hB=get2DRooHist(fnB, histname)
hC=get2DRooHist(fnC, histname)
hD=get2DRooHist(fnD, histname)

# create the observables (do this based on the data's histogram boundaries)
m2p = ROOT.RooRealVar("m2p","Invariant mass of the 2-prong",hA.GetXaxis().GetBinLowEdge(1),hA.GetXaxis().GetBinUpEdge(hA.GetXaxis().GetNbins()))
m2pg = ROOT.RooRealVar("m2pg","Invariant mass of the 2-prong and photon",hA.GetYaxis().GetBinLowEdge(1),hA.GetYaxis().GetBinUpEdge(hA.GetYaxis().GetNbins()))

# create RooDataHists
dhA=ROOT.RooDataHist(hA.GetName()+"A_dh","2D signal RooDataHist",ROOT.RooArgList(m2p,m2pg),hA)
dhB=ROOT.RooDataHist(hB.GetName()+"B_dh","2D signal RooDataHist",ROOT.RooArgList(m2p,m2pg),hB)
dhC=ROOT.RooDataHist(hC.GetName()+"C_dh","2D signal RooDataHist",ROOT.RooArgList(m2p,m2pg),hC)
dhD=ROOT.RooDataHist(hD.GetName()+"D_dh","2D signal RooDataHist",ROOT.RooArgList(m2p,m2pg),hD)

# create RooAbsPdfs out of datahists
pdfA=ROOT.RooHistPdf(hA.GetName()+"A_pdf","2D signal pdf",ROOT.RooArgSet(m2p,m2pg), dhA)
pdfB=ROOT.RooHistPdf(hB.GetName()+"B_pdf","2D signal pdf",ROOT.RooArgSet(m2p,m2pg), dhB)
pdfC=ROOT.RooHistPdf(hC.GetName()+"C_pdf","2D signal pdf",ROOT.RooArgSet(m2p,m2pg), dhC)
pdfD=ROOT.RooHistPdf(hD.GetName()+"D_pdf","2D signal pdf",ROOT.RooArgSet(m2p,m2pg), dhD)

# create a grid with each pdf at a corner
bintx=ROOT.RooBinning(1,txmin,txmax)
binty=ROOT.RooBinning(1,tymin,tymax)
grid=ROOT.RooMomentMorphFuncNDFix.Grid2(bintx,binty);
grid.addPdf(pdfA,0,0)
grid.addPdf(pdfB,1,0)
grid.addPdf(pdfC,0,1)
grid.addPdf(pdfD,1,1)

# morph and create a new 2D histogram in its place
morph=ROOT.RooMomentMorphFuncNDFix("morph","morph",ROOT.RooArgList(tx,ty),ROOT.RooArgList(m2p,m2pg),grid,ROOT.RooMomentMorphFuncNDFix.Linear);
morph.setPdfMode()
tx.setVal(wmass)
ty.setVal(pmass)
morphhist=hA.Clone("morphhist")
morphhist.Reset()
morphhist=morph.fillHistogram(morphhist,ROOT.RooArgList(m2p,m2pg))

# create PDFs for different m2p slices
fileout = ROOT.TFile("sigworkspace_"+tlabel+".root", "RECREATE")
fileout.cd()
w = ROOT.RooWorkspace("w","w")
for bin in range(1,morphhist.GetXaxis().GetNbins()+1):
    label = "bin"+str(bin)
    projy=morphhist.ProjectionY("_py"+label,bin,bin)
    accnum = projy.Integral(1,projy.GetXaxis().GetNbins())
    accden = morphhist.Integral(1,morphhist.GetXaxis().GetNbins(),1,morphhist.GetYaxis().GetNbins())
    dh=ROOT.RooDataHist("dh"+label,"dh"+label,m2pg,projy)
    sigpdf1d = ROOT.RooHistPdf("sigpdf_"+label,"signal PDF for a slice in m2p",m2pg,dh)
    sliceacc = ROOT.RooRealVar("sliceacc_"+label,"acceptance in a given slice",accnum/accden)
    sliceacc.setConstant(True)
    print("The slice acceptance for "+label+" is "+str(sliceacc.getValV()))
    norm = ROOT.RooProduct("signal_"+label+"_norm", "Normalisation of signal", ROOT.RooArgList(xsec,acceff,sliceacc))
    getattr(w,"import")(sigpdf1d)
    getattr(w,"import")(norm)

morphhist.Write()
hA.Write("hA")
hB.Write("hB")
hC.Write("hC")
hD.Write("hD")
w.Print()
w.Write()
fileout.Close()
