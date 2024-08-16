import argparse
import sys
from subprocess import call
import os
from runDragenSomaticWGS_v4_2_4 import createFastqListFile

# This script is used to run DRAGEN v4.2 alignment on fastq files

def params():
    parser = argparse.ArgumentParser(description='Run DRAGEN Alignment')
    parser.add_argument('-s','--sampleSheet',
                        required = True,
                        help = 'Sample sheet')
    parser.add_argument('-i','--inputDir',
                        required = True,
                        help = 'Folder containing the sample fastq files')
    parser.add_argument('-b','--bclFolder',
                        required = True,
                        help = 'Path to folder containing the bcl files')    
    parser.add_argument('--globSequencer',
                        required = False, default = '*_A00854_*',
                        help = 'Glob pattern to search for the sequencer run folders')
    parser.add_argument('-o','--outDir',
                        required = True,
                        help = 'Path to output folder')
    args = parser.parse_args()
    return args

def getSampleList(sampleSheet):
    sampleList = list()
    with open(sampleSheet,'r') as inputFile:
        for line in inputFile:
            line = line.strip('\n').split(',')
            sampleList.append(line[0])
    return sampleList

def runDragenAlignment(sampleName,FastqListFile,outDir):
        sampleOutDir = outDir + '' + sampleName
        call("mkdir " + sampleOutDir, shell = True)
        call('dragen \
                    -r /mnt/MTP_Storage_SSD/genome/dragenV4_hg38_graph/ \
                    --tumor-fastq-list ' + FastqListFile + ' \
                --output-directory ' + sampleOutDir + ' \
                --output-file-prefix ' + sampleName + ' \
                --soft-read-trimmers=polyg,adapter \
                --trim-adapter-read1=/mnt/MTP_Storage_SSD/wgs_adapters/adapter_read1.fa \
                --trim-adapter-read2=/mnt/MTP_Storage_SSD/wgs_adapters/adapter_read2.fa \
                --enable-duplicate-marking true --enable-map-align-output true  \
                --gc-metrics-enable true \
                --qc-coverage-ignore-overlaps=true \
                2> ' + sampleOutDir + '/' + sampleName + '_run_mapping.log' ,shell = True)

def main():
    args = params()
    sampleList = getSampleList(args.sampleSheet)
    for sampleName in sampleList:
        try:
            FastqListFile = args.outDir + '/' + sampleName + '_fastqList.csv'
            if not os.path.isfile(FastqListFile):
                createFastqListFile(sampleName,args.bclFolder,args.inputDir, args.outDir, args.globSequencer)
            runDragenAlignment(sampleName,FastqListFile,args.outDir)
        except Exception as error:
            print("Sample failed:", sampleName)
            print("Error message:", error) 
            continue

if __name__== "__main__":
    main()
