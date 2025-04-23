import vide.voidUtil as vu
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys

def main():
    catalogDict = {}
    catalog_list = []
    
    sample_path = f"/mnt/data/yuka/output/vide/UCHUU/UCHUU_allss1.0/sample_UCHUU_allss1.0_z0.00_d00"
    catalog = vu.loadVoidCatalog(sample_path, dataPortion="central", untrimmed=True, dim=3)
    
    catalog_list.append(catalog)
        
    catalogDict["all"] = catalog_list
    
    catalog_list = []
    
    sample_path = f"/mnt/data/yuka/output/vide/UCHUU/UCHUU_DESIss1.0/sample_UCHUU_DESIss1.0_z0.00_d00"
    catalog = vu.loadVoidCatalog(sample_path, dataPortion="central", untrimmed=True, dim=3)
    
    catalog_list.append(catalog)
        
    catalogDict["DESI"] = catalog_list
    
    colordict = {"DESI":"purple","all":"black"}


    vu.multiNumberFunction(catalogDict,
                           colordict,
                           figDir="/mnt/data/yuka/output/vide/UCHUU/figs/",
                           plotName=f"numberfunc_all_",
                           cumulative=False,
                           binWidth=2)

def main_():
    catalogDict = {}
    catalog_list = []
    
    sample_path = f"/mnt/data/yuka/output/vide/UCHUU/UCHUU_allss1.0/sample_UCHUU_allss1.0_z0.00_d00"
    catalog = vu.loadVoidCatalog(sample_path, dataPortion="central", untrimmed=True, dim=3)
    
    catalog_list.append(catalog)
        
    catalogDict["all"] = catalog_list
    
    catalog_list = []
    
    sample_path = f"/mnt/data/yuka/output/vide/UCHUU/UCHUU_DESIss1.0/sample_UCHUU_DESIss1.0_z0.00_d00"
    catalog = vu.loadVoidCatalog(sample_path, dataPortion="central", untrimmed=True, dim=3)
    
    catalog_list.append(catalog)
        
    catalogDict["DESI"] = catalog_list
    
    colordict = {"DESI":"purple","all":"black"}

    rmin = 20
    rmax=25
    binCenters, meandict, STDdict = vu.multiProfile(catalogDict, float(rmin), float(rmax))
    
    print(meandict)

    for key, array in meandict.items():
        
        mean = array
        sigma = STDdict[key]
        lower = mean - sigma
        upper = mean + sigma
        vu.fill_between(binCenters, lower, upper,
                        color=colordict[key],
                        alpha=0.5,
                        )
        lineStyle = '-'
        plt.plot(binCenters, mean, lineStyle,
                 label=key, color=colordict[key],
                 linewidth=3)
        
    plt.xlabel("distance from center [$h^{-1}Mpc$]")
    #plt.xlabel("normalized distance from center r/r$_v$ [$h^{-1}Mpc$]")
    #plt.ylabel("MJy/sr")
    plt.ylabel("normalized density")
    #plt.ylim(0.2, 1.6)
        
    plt.legend(loc = "upper left", fancybox=True, prop={'size':14})
    #plt.savefig("/mnt/data/yuka/output/vide/U/figs/fig_"+f"luminosity"+".pdf", bbox_inches="tight")
    plt.savefig("/mnt/data/yuka/output/vide/UCHUU/figs/fig_"+f"profile_{rmin}-{rmax}"+".png", bbox_inches="tight")
    plt.savefig("/mnt/data/yuka/output/vide/UCHUU/figs/fig_"+f"profile_{rmin}-{rmax}"+".pdf", bbox_inches="tight")
    #plt.savefig("/mnt/data/yuka/output/vide/SIDES_original/figs/fig_"+f"luminosity"+".png", bbox_inches="tight")
    plt.clf()
    
    rmin = 30
    rmax= 35
    binCenters, meandict, STDdict = vu.multiProfile(catalogDict, float(rmin), float(rmax))
    
    print(meandict)

    #for key, array in meandict.items():
        
        #mean = array
        #sigma = STDdict[key]
        #lower = mean - sigma
        #upper = mean + sigma
        #vu.fill_between(binCenters, lower, upper,
                        #color=colordict[key],
                        #alpha=0.5,
                        #)
        #lineStyle = '-'
        #plt.plot(binCenters, mean, lineStyle,
                 #label=key, color=colordict[key],
                 #linewidth=3)
        
    #plt.xlabel("distance from center [$h^{-1}Mpc$]")
    #plt.xlabel("normalized distance from center r/r$_v$ [$h^{-1}Mpc$]")
    #plt.ylabel("MJy/sr")
    #plt.ylabel("normalized density")
    #plt.ylim(0.2, 1.6)
        
    #plt.legend(loc = "upper left", fancybox=True, prop={'size':14})
    #plt.savefig("/mnt/data/yuka/output/vide/U/figs/fig_"+f"luminosity"+".pdf", bbox_inches="tight")
    #plt.savefig("/mnt/data/yuka/output/vide/UCHUU/figs/fig_"+f"profile_{rmin}-{rmax}"+".png", bbox_inches="tight")
    #plt.savefig("/mnt/data/yuka/output/vide/UCHUU/figs/fig_"+f"profile_{rmin}-{rmax}"+".pdf", bbox_inches="tight")
    #plt.savefig("/mnt/data/yuka/output/vide/SIDES_original/figs/fig_"+f"luminosity"+".png", bbox_inches="tight")
    #plt.clf()
    
    rmin = 40
    rmax= 45
    binCenters, meandict, STDdict = vu.multiProfile(catalogDict, float(rmin), float(rmax))
    
    print(meandict)

    for key, array in meandict.items():
        
        mean = array
        sigma = STDdict[key]
        lower = mean - sigma
        upper = mean + sigma
        vu.fill_between(binCenters, lower, upper,
                        color=colordict[key],
                        alpha=0.5,
                        )
        lineStyle = '-'
        plt.plot(binCenters, mean, lineStyle,
                 label=key, color=colordict[key],
                 linewidth=3)
        
    plt.xlabel("distance from center [$h^{-1}Mpc$]")
    #plt.xlabel("normalized distance from center r/r$_v$ [$h^{-1}Mpc$]")
    #plt.ylabel("MJy/sr")
    plt.ylabel("normalized density")
    #plt.ylim(0.2, 1.6)
        
    plt.legend(loc = "upper left", fancybox=True, prop={'size':14})
    #plt.savefig("/mnt/data/yuka/output/vide/U/figs/fig_"+f"luminosity"+".pdf", bbox_inches="tight")
    plt.savefig("/mnt/data/yuka/output/vide/UCHUU/figs/fig_"+f"profile_{rmin}-{rmax}"+".png", bbox_inches="tight")
    plt.savefig("/mnt/data/yuka/output/vide/UCHUU/figs/fig_"+f"profile_{rmin}-{rmax}"+".pdf", bbox_inches="tight")
    #plt.savefig("/mnt/data/yuka/output/vide/SIDES_original/figs/fig_"+f"luminosity"+".png", bbox_inches="tight")
    plt.clf()
    
        
    
    return
    
    

def _main():
    #h = sys.argv[1] if len(sys.argv) > 1 else 0.015
    #sample_path = "/mnt/data/yuka/output/vide/example_simulation/sim_ss1.0/sample_sim_ss1.0_z0.00_d00"
    #index = [0,1,2,3,4,6,7,8,9,10,11,12,14,15,16,17,18,19]
    rmin =40
    rmax= 45
    catalogDict = {}
    catalog_list = []
    index = np.arange(20)
    for i in index:
        sample_path = f"/mnt/data/yuka/output/vide/SIDES_original/original_002_ss1.0/sample_original_002_ss1.0_z{i}.00_d00"
        #sample_path = f"/mnt/data/yuka/output/vide/SIDES_original/original_weight_02_ss1.0/sample_original_weight_02_ss1.0_z{i}.00_d00"
        catalog = vu.loadVoidCatalog(sample_path, dataPortion="central", untrimmed=True, dim=2)
        catalog_list.append(catalog)
        
    #catalogDict["0.02 MJy/sr"] = catalog_list
    
    catalog_list = []
    index = np.arange(20)
    for i in index:
        sample_path = f"/mnt/data/yuka/output/vide/SIDES_original/original_ss1.0/sample_original_ss1.0_z0.78_{i}"
        #sample_path = f"/mnt/data/yuka/output/vide/SIDES_original/original_weight_ss1.0/sample_original_weight_ss1.0_z{i}.00_d00"
        catalog = vu.loadVoidCatalog(sample_path, dataPortion="central", untrimmed=True, dim=2)
        catalog = vu.loadVoidCatalog(sample_path, dataPortion="central", untrimmed=True, dim=2)
        catalog_list.append(catalog)
        
    #catalogDict["0.01 MJy/sr"] = catalog_list
    
    catalog_list = []
    index = index = [0,1,2,3,4,6,7,8,9,10,11,12,14,15,16,17,18,19]
    for i in index:
        #sample_path = f"/mnt/data/yuka/output/vide/SIDES_original/original_005_ss1.0/sample_original_005_ss1.0_z{i}.00_d00"
        sample_path = f"/mnt/data/yuka/output/vide/SIDES_original/original_weight_005_ss1.0/sample_original_weight_005_ss1.0_z{i}.00_d00"
        catalog = vu.loadVoidCatalog(sample_path, dataPortion="central", untrimmed=True, dim=2)
        catalog_list.append(catalog)
    
    catalogDict["0.005 MJy/sr"] = catalog_list
        
    #catalog_list = []
    #sample_path = "/mnt/data/yuka/output/vide/UCHUU_2D/UCHUU_2Dss1.0/sample_UCHUU_2Dss1.0_z0.00_d00"
    #sample_path = f"/mnt/data/yuka/output/vide/SIDES_original/original_weight_005_ss1.0/sample_original_weight_005_ss1.0_z{i}.00_d00"
    #catalog = vu.loadVoidCatalog(sample_path, dataPortion="central", untrimmed=True, dim=2)
    #catalog_list.append(catalog)
    #catalogDict["halo"] = catalog_list
    
    #vu.plotNumberFunction(catalog,
                       #figDir="/mnt/data/yuka/output/vide/UCHUU_2D/figs/",
                       #plotName=f"numberfunc_halo",
                       #cumulative=False,
                       #binWidth=2)
    #Scatter2D(catalog, f"/mnt/data/yuka/output/vide/UCHUU_2D/figs/untrimmed_scatter.png")
    #catalog = vu.loadVoidCatalog(sample_path, dataPortion="central")
    #Scatter(catalog, "/mnt/data/yuka/output/vide/SIDES/figs/trimmed_scatter.png")
    
    colordict = {"0.005 MJy/sr":"red","0.01 MJy/sr":"blue", "0.02 MJy/sr":"green"}

    #vu.multiNumberFunction(catalogDict,
                           #colordict,
                           #figDir="/mnt/data/yuka/output/vide/SIDES_original/figs/",
                           #plotName=f"numberfunc_all_",
                           #cumulative=False,
                           #binWidth=2)
                           
    #binCenters, meandict, STDdict = vu.multiProfile(catalogDict, float(rmin), float(rmax))
    binCenters, meandict, STDdict = vu.multiLuminosity(catalogDict, nBins=20)
    print(meandict)

    for key, array in meandict.items():
        
        mean = array
        sigma = STDdict[key]
        lower = mean - sigma
        upper = mean + sigma
        vu.fill_between(binCenters, lower, upper,
                        color=colordict[key],
                        alpha=0.5,
                        )
        lineStyle = '-'
        plt.plot(binCenters, mean, lineStyle,
                 label=key, color=colordict[key],
                 linewidth=3)
        
    #plt.xlabel("distance from center [$h^{-1}Mpc$]")
    plt.xlabel("normalized distance from center r/r$_v$ [$h^{-1}Mpc$]")
    plt.ylabel("MJy/sr")
    #plt.ylabel("normalized density")
    #plt.ylim(0.2, 1.6)
        
    plt.legend(loc = "upper left", fancybox=True, prop={'size':14})
    plt.savefig("/mnt/data/yuka/output/vide/SIDES_original/figs/fig_"+f"luminosity"+".pdf", bbox_inches="tight")
    #plt.savefig("/mnt/data/yuka/output/vide/SIDES_original/figs/fig_"+f"weight_profile_{rmin}-{rmax}"+".png", bbox_inches="tight")
    #plt.savefig("/mnt/data/yuka/output/vide/SIDES_original/figs/fig_"+f"weight_profile_{rmin}-{rmax}"+".pdf", bbox_inches="tight")
    plt.savefig("/mnt/data/yuka/output/vide/SIDES_original/figs/fig_"+f"luminosity"+".png", bbox_inches="tight")    
    #vu.computeXcor(catalog,
                   #figDir="/mnt/data/yuka/output/vide/SIDES/figs/",
                   #Nmesh = 256,
                   #Nbin = 100)
                   
              ###############################     
    rmin =20
    rmax= 25
    binCenters, mean, sigma = vu.buildProfile(catalog, rmin, rmax, nBins=15, dim=3)
    lower = mean - sigma
    upper = mean + sigma
    plt.clf()
    vu.fill_between(binCenters, lower, upper,
                    color='black',
                    alpha=0.5,
                    )
    lineStyle = '-'
    
    label = f"radius {rmin} - {rmax}" +"h$^{-1}$Mpc"
    plt.plot(binCenters, mean, lineStyle,
                label=label, color='black',
                linewidth=3)
    
    rmin =30
    rmax= 35
    binCenters, mean, sigma = vu.buildProfile(catalog, rmin, rmax, nBins=15, dim=3)
    lower = mean - sigma
    upper = mean + sigma
    vu.fill_between(binCenters, lower, upper,
                    color='black',
                    alpha=0.5,
                    )
    lineStyle = '-.'
    label = f"radius {rmin} - {rmax}" +"h$^{-1}$Mpc"
    plt.plot(binCenters, mean, lineStyle,
                label=label, color='black',
                linewidth=3)
    
    rmin =40
    rmax= 45
    binCenters, mean, sigma = vu.buildProfile(catalog, rmin, rmax, nBins=15, dim=3)
    lower = mean - sigma
    upper = mean + sigma
    vu.fill_between(binCenters, lower, upper,
                    color='black',
                    alpha=0.5,
                    )
    lineStyle = ':'
    label = f"radius {rmin} - {rmax}" +"h$^{-1}$Mpc"
    plt.plot(binCenters, mean, lineStyle,
                label=label, color='black',
                linewidth=3)
        
    plt.xlabel("distance from center [$h^{-1}Mpc$]")
    #plt.xlabel("normalized distance from center r/r$_v$ [$h^{-1}Mpc$]")
    #plt.ylabel("MJy/sr")
    plt.ylabel("normalized density")
    #plt.ylim(0.2, 1.6)
    plt.xlim(0.0, 50.0)
        
    plt.legend(loc = "lower right", fancybox=True, prop={'size':14})
    #plt.savefig("/mnt/data/yuka/output/vide/SIDES_original/figs/fig_"+f"luminosity"+".pdf", bbox_inches="tight")
    plt.savefig("/mnt/data/yuka/output/vide/UCHUU/figs/fig_"+f"all_profile"+".png", bbox_inches="tight")
    plt.savefig("/mnt/data/yuka/output/vide/UCHUU/figs/fig_"+f"all_profile"+".pdf", bbox_inches="tight")
    #plt.savefig("/mnt/data/yuka/output/vide/SIDES_original/figs/fig_"+f"luminosity"+".png", bbox_inches="tight")    
    #vu.computeXcor(catalog,
                   #figDir="/mnt/data/yuka/output/vide/SIDES/figs/",
                   #Nmesh = 256,
                   #Nbin = 100)
    return


def Scatter(catalog, filename, zmin=10, zmax=50):
    fig = plt.figure(figsize=(6,6))
    for i in range(len(catalog.voids)):
        voidID = catalog.voids[i].voidID
        macrocenter = catalog.voids[i].macrocenter
        volume =  catalog.voids[i].volume
        voidPart = vu.getVoidPart(catalog, voidID)
        X, Y, Z = [], [], []
        if((zmin<macrocenter[2])&(macrocenter[2]<zmax)):
            plt.scatter(macrocenter[0], macrocenter[1], marker="x", c="k", s=0.5)
            #plt.text(macrocenter[0], macrocenter[1], str(volume), fontsize=3, color="black")
        for part in voidPart:
            X.append(part.x)
            Y.append(part.y)
            Z.append(part.z)
        X = np.array(X)
        Y = np.array(Y)
        Z = np.array(Z)
        selection = (zmin<Z)&(Z<zmax)
        plt.scatter(X[selection], Y[selection], color=cm.tab20(i % 20), s=0.2)
    plt.savefig(filename)
    plt.show()
    plt.clf()
    return 0

def Scatter2D(catalog, filename):
    fig = plt.figure(figsize=(9,9))
    for i in range(len(catalog.voids)):
        voidID = catalog.voids[i].voidID
        macrocenter = catalog.voids[i].macrocenter
        print(f"z coordinate is {macrocenter[2]}")
        radius =  catalog.voids[i].radius
        print(f"volume of void {voidID} is {radius} Mpc/h")
        voidPart = vu.getVoidPart(catalog, voidID)
        X, Y, Z = [], [], []
        #plt.scatter(macrocenter[0], macrocenter[1], marker="x", c="k", s=10.0)
        #plt.text(macrocenter[0], macrocenter[1], str(radius), fontsize=10, color="black")
        for part in voidPart:
            X.append(part.x)
            Y.append(part.y)
            Z.append(part.z)
        X = np.array(X)
        Y = np.array(Y)
        Z = np.array(Z)
        plt.scatter(X, Y, color=cm.tab20(i % 20), s=0.001, marker=".")
    plt.grid()
    plt.savefig(filename)
    plt.show()
    plt.clf()
    return 0
    
    

if __name__ == "__main__":
    main()
