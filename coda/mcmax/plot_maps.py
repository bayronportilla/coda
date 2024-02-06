import prodimopy.read as pread
import prodimopy.plot as pplot
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

"""
Example.
sys.argv[1] = /data/users/bportilla/runs_P4/new_master/fixed_C/C2O_effect/
sys.argv[2] = C2O_0.457

"""
dir     = sys.argv[1]
run     = sys.argv[2]
model   = pread.read_prodimo(dir+run)
pp      = pplot.Plot(None)


# Gas density
tcont=pplot.Contour(model.rhog,[10**-15],showlabels=True,label_fontsize=8,label_fmt="%.1e")
pp.plot_cont(model,model.rhog,label=r"$\log(\rho_{gas})\,[g/cm3]}$",xlog=False,oconts=[tcont])
plt.show()


"""
# Tgas
#tcont=pplot.Contour(model.tg,[20],showlabels=True,label_fontsize=10,label_fmt="%.1e")
#cbticks=[20,100,1000]
#cbticks=[20,40,60,80,100,200,400,600,1000]
#pp.plot_cont(model,model.tg,label=r"$T_{gas}\,[K]}$",xlog=True,ylim=[0,0.2],zlog=False,zlim=[model.tg.min(),1000],cmap='inferno',nbins=10,clevels=cbticks)
pp.plot_cont(model,model.tg,label=r"$T_{gas}\,[K]}$",
    xlog=False,
    zlog=True,
    cmap='inferno',
    nbins=75,
    zr=False)
plt.show()
"""

"""
# Tdust
#tcont=pplot.Contour(model.tg,[20],showlabels=True,label_fontsize=10,label_fmt="%.1e")
#cbticks=[20,40,60,80,100]
#pp.plot_cont(model,model.td,label=r"$T_{dust}\,[K]}$",xlog=False,ylim=[0,0.2],oconts=[tcont],clevels=cbticks,clabels=map(str,cbticks),cb_format="%.0f")
#pp.plot_cont(model,model.td,label=r"$T_{dust}\,[K]}$",xlog=False,ylim=[0,0.2],zlog=False,zlim=[model.td.min(),200],cmap='inferno',nbins=10,clevels=cbticks)
pp.plot_cont(model,model.td,label=r"$T_{dust}\,[K]}$",xlog=True,ylim=[0,0.2],zlog=False,cmap='inferno',nbins=20)
plt.show()
"""

"""
# d2g
d2g=model.d2g
g2d=np.zeros(d2g.shape)

for i in range(d2g.shape[0]):
    for j in range(d2g.shape[1]):
        g2d[i,j]=d2g[i,j]**-1

for i in range(g2d.shape[0]):
    print(g2d[i])
"""

"""
# Line origin for 13CO over gas temperature
tcont=pplot.Contour(model.tg,[30],showlabels=True,label_fontsize=10,label_fmt="%.1e")
pp.plot_line_origin(model,[["13CO",1.360228E+03]],model.tg,xlog=False,label=r"$Tg\,[K]}$")
#pp.plot_line_origin(model,[["13CO",1.360228E+03]],model.tg,xlog=False,label=r"$Tg\,[K]}$",oconts=[tcont])
plt.show()
"""

"""
# Line origin for 12CO over gas temperature
tcont=pplot.Contour(model.tg,[20],showlabels=True,label_fontsize=10,label_fmt="%.1e")
pp.plot_line_origin(model,[["CO",1.300403E+03]],model.tg,xlog=False,
label=r"$Tg\,[K]}$")
plt.show()
"""
"""
# Line origin for C18O over gas temperature
tcont=pplot.Contour(model.tg,[20],showlabels=True,label_fontsize=10,label_fmt="%.1e")
pp.plot_line_origin(model,[["C18O",1365.42164]],
    model.tg,
    xlog=False,
    nbins=15,
    zlog=False,
    zlim=[10,60],
    label=r"$Tg\,[K]}$")
plt.show()
"""
"""
# Line origin over H2 abundance
ncont=pplot.Contour(model.nmol[:,:,4],[6400],showlabels=True,label_fontsize=10,label_fmt="%.1e")
#pp.plot_line_origin(model,[["CO",1.300403E+03]],model.nmol[:,:,4],xlog=False,label=r"$\log n(H_2)\,[cm^-3]}$",oconts=[ncont])
#pp.plot_line_origin(model,[["o-H2O",6.86395]],model.nmol[:,:,4],xlog=False,label=r"$\log n(H_2)\,[cm^-3]}$",oconts=[ncont])
pp.plot_line_origin(model,[["o-H2O",12.453]],model.nmol[:,:,4],xlog=True,label=r"$\log n(H_2)\,[cm^-3]}$",oconts=[ncont])
#pp.plot_line_origin(model,[["H2O",6.85116]],model.nmol[:,:,4],xlog=False,label=r"$\log n(H_2)\,[cm^-3]}$",oconts=[ncont])
#pp.plot_line_origin(model,[["p-H2O",6.85745482]],model.nmol[:,:,4],xlog=True,label=r"$\log n(H_2)\,[cm^-3]}$",oconts=[ncont])
#pp.plot_line_origin(model,[["p-H2O",6.4200957]],model.nmol[:,:,4],xlog=False,label=r"$\log n(H_2)\,[cm^-3]}$",oconts=[ncont])
#pp.plot_line_origin(model,[["13CO",1360.22798]],model.nmol[:,:,4],xlog=False,label=r"$\log \ n(H_2)\,[cm^-3]}$",oconts=[ncont])
#pp.plot_line_origin(model,[["CO",1.300403E+03],["13CO",1360.22798]],model.nmol[:,:,4],xlog=False,label=r"$\log n(H_2)\,[cm^-3]}$",oconts=[ncont])
plt.show()
"""

"""
# Line origin for 12CO J=2-1 over H2 abundance
#ncont=pplot.Contour(model.nmol[:,:,27],[15000],showlabels=True,label_fontsize=10,label_fmt="%.1e")
pp.plot_line_origin(model,[["12CO",1.300403E+03]],model.nmol[:,:,4],xlog=False,
label=r"$\log n(H2)\,[cm^-3]}$")
plt.show()
"""

"""
# Line origin for 12CO J=3-2 over H2 abundance
#ncont=pplot.Contour(model.nmol[:,:,27],[15000],showlabels=True,label_fontsize=10,label_fmt="%.1e")
pp.plot_line_origin(model,[["12CO",866.96337]],model.nmol[:,:,4],xlog=False,
label=r"$\log n(H2)\,[cm^-3]}$")
plt.show()
"""

"""
# Line origin for 13CO J=2-1 over H2 abundance
pp.plot_line_origin(model,[["13CO",1360.22798]],model.nmol[:,:,4],xlog=False,
label=r"$\log n(H2)\,[cm^-3]}$")
plt.show()
"""


"""
# Line origin for C18O J=2-1 over H2 abundance
pp.plot_line_origin(model,[["C18O",1365.42164]],model.nmol[:,:,4],xlog=False,
label=r"$\log n(H2)\,[cm^-3]}$")
plt.show()
"""

"""
# Line origin for N2H+ J=3-2 over H2 abundance
pp.plot_line_origin(model,[["N2H+",1.0725578104E+03]],model.nmol[:,:,4],xlog=False,
label=r"$\log n(H2)\,[cm^-3]}$")
plt.show()
"""
"""
# Abundance
pp.plot_abuncont(model,species='CO',
    xlog=False,
    cmap='Blues',
    nbins=5,
    rel2H=True,
    zr=True)
plt.show()
"""

"""
# Abundance radial
pp.plot_abunrad(model,
    species='CO',
    useNH=False,
    xlog=False)
plt.show()
"""
"""
print(model.spnames)
#sys.exit()
pp.plot_abun_midp(model,
    species=['CO','N'],
    useNH=False,
    xlog=True,
    ylim=[1e-11,2e-4])
plt.show()
"""
"""
# Chi
pp.plot_cont(model,model.chi,xlog=False,ylim=[0,0.3],xlim=[15,45],zlog=True,zlim=[10**-4,10**4])
plt.show()
"""
"""
pp.plot_abunvert(model,30.0,['CO'],useZr=True)
plt.show()
"""
"""
# Plot maing heating and cooling processes
pp.plot_heat_cool(model,
    xlog=False,
    ylim=[0.0,0.4])
plt.show()
"""

"""
print(model.heat_mainidx)
pp.plot_vertical(model,30,model.cool_mainidx,ylabel='heating')
plt.show()
"""
"""
# Vertical Tgas
pp.plot_vertical(model,121.0,model.tg,ylabel='Tg (K)',ylog=True,zr=False)
plt.show()
"""
"""
# Vertical cut sound speed
pp.plot_vertical(model,20.8,model.soundspeed,ylabel='soundspeed')
plt.show()
"""

"""
# Radial gas temperature at midplane
pp.plot_midplane(model,model.tg,ylabel='Tgas',xlog=True)
plt.show()
"""
"""
# Radial dust temperature at midplane
pp.plot_midplane(model,model.td,ylabel='Tdust',xlog=True)
plt.show()
"""
"""
# Tau_AV over nH
avcont=pplot.Contour(model.AVrad,[1.0],colors="black",linewidths=5.0)
print(model.AV)
#sys.exit()
pp.plot_cont(model,model.nHtot,xlog=False,zlog=True,zr=False,ylog=False,
    oconts=[avcont])
plt.show()
"""
