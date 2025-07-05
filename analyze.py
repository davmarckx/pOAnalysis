# python script to run multiple at one because trkCharge field otherwise causes issues with TChain as it is a cvector<char> field instead of vector<int>

import glob
import sys
import os


files = sys.argv[1]

name = sys.argv[2]

print(files)

print(name)


files = glob.glob(files+"*")

for i,file in enumerate(files):
    os.system("./makeFromHIForest/fillMiniEvent {} {}".format(file, name+"_"+str(i)))


filename = "/pnfs/iihe/cms/store/user/dmarckx/pO_miniAOD/MiniEvent/output_" + name +"_withD.root"
subfilename = "/pnfs/iihe/cms/store/user/dmarckx/pO_miniAOD/MiniEvent/output_" + name +"_*_withD.root"

os.system("hadd -f {} {}".format(filename,subfilename))
os.system("rm {}".format(subfilename))
