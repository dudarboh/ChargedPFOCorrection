import ROOT
import numpy as np

ROOT.gStyle.SetNdivisions(512)
ROOT.gStyle.SetNdivisions(512, "Y")
#
# gInterpreter.Declare('''
#
#
#
#
# ''')

# Filter and calculate betas
df = ROOT.RDataFrame("VertexAnalysis", "/nfs/dust/ilc/user/dudarboh/jobs/output/check_vertex.root")

h_old_x = df.Define("diff", "(primary_vtx_mc - primary_vtx_old).x()*1000").Histo1D(("h_old_x", "BEFORE refit; (#vec{r}_{MC} - #vec{r}_{reco}).x() (#mum); Events", 1000, -1., 1.), "diff")
h_new_x = df.Define("diff", "(primary_vtx_mc - primary_vtx_new).x()*1000").Histo1D(("h_new_x", "AFTER refit; |#vec{r}_{MC} - #vec{r}_{reco}| (#mum); Events", 1000, -1., 1.), "diff")

h_old_y = df.Define("diff", "(primary_vtx_mc - primary_vtx_old).y()*1000").Histo1D(("h_old_y", "BEFORE refit; (#vec{r}_{MC} - #vec{r}_{reco}).y() (#mum); Events", 1000, -1., 1.), "diff")
h_new_y = df.Define("diff", "(primary_vtx_mc - primary_vtx_new).y()*1000").Histo1D(("h_new_y", "AFTER refit; |#vec{r}_{MC} - #vec{r}_{reco}| (#mum); Events", 1000, -1., 1.), "diff")

h_old_z = df.Define("diff", "(primary_vtx_mc - primary_vtx_old).z()*1000").Histo1D(("h_old_z", "BEFORE refit; (#vec{r}_{MC} - #vec{r}_{reco}).z() (#mum); Events", 1000, -20., 20.), "diff")
h_new_z = df.Define("diff", "(primary_vtx_mc - primary_vtx_new).z()*1000").Histo1D(("h_new_z", "AFTER refit; |#vec{r}_{MC} - #vec{r}_{reco}| (#mum); Events", 1000, -20., 20.), "diff")



canvas = ROOT.TCanvas()
canvas.Divide(2, 2)
canvas.SetGridx()
canvas.SetGridy()
# ROOT.gStyle.SetOptStat(10)

canvas.cd(1)
h_old_x.Draw()
h_new_x.Draw("sames")
h_new_x.SetLineColor(2)

canvas.cd(2)
h_old_y.Draw()
h_new_y.Draw("sames")
h_new_y.SetLineColor(2)

canvas.cd(3)
h_old_z.Draw()
h_new_z.Draw("sames")
h_new_z.SetLineColor(2)

canvas.BuildLegend()
h_old_x.SetTitle("Primary vertex residuals")
h_old_y.SetTitle("Primary vertex residuals")
h_old_z.SetTitle("Primary vertex residuals")
canvas.Update()
input("wait")
