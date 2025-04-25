import vide.voidUtil as vu
import numpy as np
import pickle

from vide.voidUtil import getArray

def load_sample(sampledir):
    with open(sampledir+"/sample_info.dat", 'rb') as input:
        sample = pickle.load(input)
    return sample

def load_data(zobovdir):
    sample_path = zobovdir
    sample = load_sample(zobovdir)
    if sample.dataType == 'LIM':
        dim = 2
    else:
        dim=3
    catalog = vu.loadVoidCatalog(sample_path, dataPortion="central", untrimmed=True, dim=dim)
    
    radius = getArray(catalog.voids, 'radius')
    density = getArray(catalog.voids, 'centralDen')
    densCon = getArray(catalog.voids, 'densCon')
    x = []
    y = []
    z = []
    for i in range(len(catalog.voids)):
        void = catalog.voids[i]
        center = void.macrocenter
        x.append(center[0])
        y.append(center[1])
        z.append(center[2])
    return np.array(x), np.array(y), np.array(z), radius, density, densCon

def Output(zobovdir, samplename):
    filename = zobovdir + f"/{samplename}_cleaned.txt"
    x, y, z, radius, density, densCon = load_data(zobovdir)
    table = np.column_stack((x, y, z, radius, density, densCon))
    print(table)
    with open(filename, 'w') as f:
      np.savetxt(f, table, fmt='%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f')
      
    print(f"downloaded cleaned file to {filename}")
    return 0
    
        
    
        