import argparse
from glob import glob
from subprocess import call
import pandas as pd
import os 
import shutil
from runDragenSomaticWGS_v3_9_5 import updateMatchFile

def params():
	parser = argparse.ArgumentParser(description='Run DRAGEN HLA typing')
	parser.add_argument('-n','--normalOutDir',
						required = True,
						help = 'Path to normal output folder')
	parser.add_argument('-t','--tumorOutDir',
						required = True,
						help = 'Path to tumor output folder')
	parser.add_argument('-m','--matchFile',
						required = True,
						help = 'Path to file containing the sample match information')
	args = parser.parse_args()
	return args

def runDragenHLA(tumorID,normalID,normalOutDir,tumorOutDir):
        bamNormal = glob(normalOutDir + "/" + normalID + "/*.bam")
        bamTumor = glob(tumorOutDir + "/" + tumorID + "/*_tumor.bam")
        tmpDir = "/mnt/MTP_WGS_Share/tmp/" + tumorID + "/"
        os.makedirs(tmpDir, exist_ok = True)
        call('dragen \
                --enable-hla=true \
                --enable-map-align=false \
                --enable-map-align-output=false \
                -r /mnt/MTP_Storage_SSD/genome/hg38_alt_masked_graph_v2/ \
                --tumor-bam-input ' + bamTumor[0] + ' \
                --bam-input ' + bamNormal[0] + ' \
                --output-directory ' + tmpDir + ' \
                --output-file-prefix ' + tumorID + ' \
                --hla-min-reads=50 \
                2> ' + tmpDir + tumorID + '_run_hla.log' ,shell = True)

def moveFiles(tumorID,tumorOutDir):
        tmpDir = "/mnt/MTP_WGS_Share/tmp/" + tumorID + "/"
        # Rename the replay file
        os.rename(tmpDir + tumorID + "-replay.json", tmpDir + tumorID + "-replay-hla.json")
        # Copy the files we want - since running hla independently generates a new bam with different mapping metrics
        for f in glob(tmpDir + "*hla*"):
            shutil.copy(f, tumorOutDir + "/" + tumorID)
        shutil.copy(tmpDir + tumorID + "-replay-hla.json", tumorOutDir + "/" + tumorID)
        shutil.copy(tmpDir + tumorID + "_run_hla.log", tumorOutDir + "/" + tumorID)
        # Clean up
        shutil.rmtree(tmpDir)

def main():
    args = params()
    sampleMatch = pd.read_csv(args.matchFile, delimiter='\t')
    for i in range(0,sampleMatch.shape[0]):
        sampleName=sampleMatch["Tumor"][i]
        try:
            print(sampleName)
            runDragenHLA(sampleMatch["Tumor"][i],sampleMatch["Normal"][i],args.normalOutDir,args.tumorOutDir)
            moveFiles(sampleMatch["Tumor"][i],args.tumorOutDir)
            updateMatchFile(args.matchFile,sampleName)
        except Exception as error:
            print("Sample failed:", sampleName)
            print("Error message:", error)
            updateMatchFile(args.matchFile,sampleName,fail=True)
            continue
if __name__ == "__main__":
    main()

