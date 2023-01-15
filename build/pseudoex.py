import os

for i in range(2):
    os.system("./../install/bin/RunCrossSection -c ../json/config.json")
    #os.system("./../install/bin/BackgroundFit -c ../json/bkgfitMC.json")
    os.system("./../install/bin/BackgroundFit -c ../json/bkgfit.json")
    #os.system("./../install/bin/RunCalcXS -c ../json/xsMC.json")
    os.system("./../install/bin/RunCalcXS -c ../json/xs.json")
    #os.system("root -b -q ../macros/pion/plotXS_data.C")
    #os.system("./../install/bin/forUnfold")
    #os.system("./../../Wiener-SVD-Unfolding/Example forUnfold.root Unfold.root 2 0")