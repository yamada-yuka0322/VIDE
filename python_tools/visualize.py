import vide.voidUtil as vu
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys

def main():
    h = sys.argv[1] if len(sys.argv) > 1 else 0.015
    #sample_path = "/mnt/data/yuka/output/vide/example_simulation/sim_ss1.0/sample_sim_ss1.0_z0.00_d00"
    sample_path = f"/mnt/data/yuka/output/vide/SIDES/SIDES_{h}_weight/sample_SIDES_{h}_z0.8_weight"
    #catalog = vu.loadVoidCatalog(sample_path, dataPortion="central")
    #Scatter(catalog, "/mnt/data/yuka/output/vide/SIDES/figs/trimmed_scatter.png")

    catalog = vu.loadVoidCatalog(sample_path, dataPortion="central", untrimmed=True, dim=2)
    Scatter2D(catalog, f"/mnt/data/yuka/output/vide/SIDES/figs/untrimmed_scatter_{h}_weight.png")
    vu.plotNumberFunction(catalog,
                   figDir="/mnt/data/yuka/output/vide/SIDES/figs/",
                   plotName=f"untrimmed_numberfunc_{h}_weight",
                   cumulative=False,
                   binWidth=1)
    
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
        volume =  catalog.voids[i].volume
        voidPart = vu.getVoidPart(catalog, voidID)
        X, Y, Z = [], [], []
        plt.scatter(macrocenter[0], macrocenter[1], marker="x", c="k", s=10.0)
        #plt.text(macrocenter[0], macrocenter[1], str(volume), fontsize=10, color="black")
        for part in voidPart:
            X.append(part.x)
            Y.append(part.y)
            Z.append(part.z)
        X = np.array(X)
        Y = np.array(Y)
        Z = np.array(Z)
        plt.scatter(X, Y, color=cm.tab20(i % 20), s=2.0)
    plt.savefig(filename)
    plt.show()
    plt.clf()
    return 0
    
    

if __name__ == "__main__":
    main()
