import argparse
import sys
from subprocess import call
import os

def params():
    parser = argparse.ArgumentParser(description='Run DRAGEN Alignment')
    parser.add_argument('-s','--sampleSheet',
                        required = True,
                        help = 'SampleSheet')
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
                    -r /mnt/MTP_Storage_SSD/genome/hg38_alt_masked_graph_v2/ \
                --fastq-list ' + FastqListFile + ' \
                --output-directory ' + sampleOutDir + ' \
                --output-file-prefix ' + sampleName + ' \
                --soft-read-trimmers=polyg,adapter \
                --trim-adapter-read1=/mnt/MTP_WGS_Share/melanoma_TP/results/metagenomics/adaptor1.fa \
                --trim-adapter-read2=/mnt/MTP_WGS_Share/melanoma_TP/results/metagenomics/adaptor2.fa \
                --enable-duplicate-marking true --enable-map-align-output true  \
                --gc-metrics-enable true \
                --qc-coverage-ignore-overlaps=true \
                &> ' + sampleOutDir + '/' + sampleName + '_run.log' ,shell = True)

def main():
    args = params()
    sampleList = getSampleList(args.sampleSheet)
    for sampleName in sampleList:
        try:
            FastqListFile = args.outDir + '/' + sampleName + '_fastqList.csv'
            runDragenAlignment(sampleName,FastqListFile,args.outDir)
        except Exception as error:
            print("Sample failed:", sampleName)
            print("Error message:", error) 
            continue

if __name__== "__main__":
    main()
