import argparse
from glob import glob
import sys
from subprocess import call
import os
import pandas as pd

def params():
    parser = argparse.ArgumentParser(description='Run DRAGEN WGS Somatic Pipeline')
    parser.add_argument('-i','--inputDir',
                        required = True,
                        help = 'Folder containing the sample fastq files')
    parser.add_argument('-b','--bclFolder',
                        required = True,
                        help = 'Path to folder containing the bcl files')    
    parser.add_argument('-o','--outDir',
                        required = True,
                        help = 'Path to output folder')
    parser.add_argument('-m','--matchFile',
                        required = True,
                        help = 'Path to file containing the sample match information')
    parser.add_argument('--runIDsTumor',
                        required = False, default = None,
                        help = 'Comma delimited list of valid run IDs to search for the tumor fastqs')
    parser.add_argument('--normalOutDir',
    required = False, default = None,
            help = 'Path to normal results if normal has already been run')
    parser.add_argument('--germlineOnly', 
                        help = 'Run germline only pipeline', action='store_true')
    parser.add_argument('--globSequencer',
                        required = False, default = '*_A00854_*',
                        help = 'Glob pattern to search for the sequencer run folders')
    parser.add_argument('--dryRun',
    required = False, action='store_true',
            help = 'Do not run dragen, just create the fastq list file')
    args = parser.parse_args()
    return args

def getSampleList(sampleSheet):
    sampleList = list()
    skipLine = True
    with open(sampleSheet,'r') as inputFile:
        for line in inputFile:
            line = line.strip('\n').split(',')
            if line[0] == 'Sample_Name':
                skipLine = False
                continue
            if skipLine:
                continue
            if 'LG' in line[1]: 
                sampleList.append(line[1])
    return sampleList

def readMatchFile(inFilePath):
    matchDict = dict()
    with open(inFilePath,'r') as inFile:
        for line in inFile:
            line = line.strip('\n').split('\t')
            if line[2] == 'Pending':
                matchDict[line[0].replace('.','_')] = line[1].replace('.','_')
    return matchDict

def updateMatchFile(inFilePath,sampleName,fail=False):
    sampleMatch = pd.read_csv(inFilePath, delimiter='\t')
    if fail == False:
        sampleMatch.loc[sampleMatch["Tumor"]==sampleName,"Status"] = "Complete"
    else:
        sampleMatch.loc[sampleMatch["Tumor"]==sampleName,"Status"] = "Failed"
    sampleMatch.to_csv(inFilePath,sep='\t',index=False)
    return

def getLibraryID(inFolder):
    with open(inFolder + 'RunParameters.xml','r') as inputFile:
        for line in inputFile:
            if '<LibraryTubeSerialBarcode>' in line:
                libraryID = line.split('>')[1].split('</')[0]
                break
    return libraryID

def createFastqListFile(sampleName,bclFolder,inputDir,outDir,globSequencer,runIDs=None):
    outFilePath = outDir + sampleName + '_fastqList.csv'
    output = [['RGID','RGSM','RGLB','Lane','Read1File','Read2File']]
    if runIDs is None:
        runList = glob(inputDir + globSequencer)
    else:
        runIDList=[id for id in runIDs.split(',')]
        runList = [f for f in glob(inputDir + globSequencer) if os.path.basename(f) in runIDList]
    for run in runList:
        bclSampleFolder = bclFolder + os.path.basename(run) + '/'
        libraryID = getLibraryID(bclSampleFolder)
        sampleFastqFiles = glob(run + '/' + sampleName + '_*S*R1*.fastq.gz')
        for fastq in sampleFastqFiles:
            lane = fastq.split('_L00')[-1][0]
            flowCell = fastq.rsplit('/',1)[0].rsplit('_')[-1]
            output.append([
                flowCell + '.' + lane + '.' + sampleName,
                sampleName,
                libraryID,
                lane,
                fastq,
                fastq.replace('_R1_','_R2_')])
    with open(outFilePath,'w') as outFile:
        outFile.writelines(','.join(i) + '\n' for i in output)

    return outFilePath

def runDragenNormal(sampleName,sampleFastqListFile,outDir):
    sampleOutDir = outDir + 'normal/' + sampleName
    if not os.path.isfile(sampleOutDir + '/' + sampleName + '.mapping_metrics.csv'):
        call("mkdir " + sampleOutDir, shell = True)
        call('dragen \
            -r /mnt/MTP_Storage_SSD/genome/dragenV4_hg38_graph/ \
            --fastq-list ' + sampleFastqListFile + ' \
            --output-directory ' + sampleOutDir + ' \
            --output-file-prefix ' + sampleName + ' \
            --soft-read-trimmers=polyg,adapter \
            --trim-adapter-read1=/mnt/MTP_Storage_SSD/wgs_adapters/adapter_read1.fa \
            --trim-adapter-read2=/mnt/MTP_Storage_SSD/wgs_adapters/adapter_read2.fa \
            --enable-duplicate-marking true --enable-map-align-output true \
            --enable-variant-caller true --vc-enable-joint-detection true --vc-combine-phased-variants-distance 15 \
            --gc-metrics-enable true \
            --enable-cnv true --cnv-enable-self-normalization true \
            --enable-sv true \
            --qc-cross-cont-vcf /opt/edico/config/somatic_sample_cross_contamination_resource_hg38.vcf.gz \
            --qc-coverage-ignore-overlaps=true \
            --enable-variant-annotation true \
            --variant-annotation-assembly GRCh38 \
            --variant-annotation-data /mnt/MTP_Storage_SSD/hg38_databases/nirvana_hg38_data_20240202/ \
            2> ' + sampleOutDir + '/' + sampleName + '_run.log' ,shell = True)
    return sampleOutDir

def runDragenTumorNormal(sampleName,tumorFastqListFile,normalFastqListFile,normalOutDir,outDir):
    #First get the BAM file of the tumor
    sampleOutDir = outDir + 'tumor_normal/' + sampleName
    if not os.path.isfile(sampleOutDir + '/' + sampleName + '_tumor.bam'):
        call("mkdir " + sampleOutDir, shell = True)
        call('dragen \
                    -r /mnt/MTP_Storage_SSD/genome/dragenV4_hg38_graph/ \
                    --tumor-fastq-list ' + tumorFastqListFile + ' \
                --output-directory ' + sampleOutDir + ' \
                --output-file-prefix ' + sampleName + ' \
                --soft-read-trimmers=polyg,adapter \
                --trim-adapter-read1=/mnt/MTP_Storage_SSD/wgs_adapters/adapter_read1.fa \
                --trim-adapter-read2=/mnt/MTP_Storage_SSD/wgs_adapters/adapter_read2.fa \
                --enable-duplicate-marking true --enable-map-align-output true  \
                --gc-metrics-enable true \
                --qc-coverage-ignore-overlaps=true \
                2> ' + sampleOutDir + '/' + sampleName + '_run_mapping.log' ,shell = True)
    else:
        print("Tumor bam found. Proceeding with variant calling.")
        
    #Run the full pipeline with the tumor and normal bam files
    vcfNormal = glob(normalOutDir + "/*.hard-filtered.vcf.gz")
    cnvNormal = glob(normalOutDir + "/*.cnv.vcf.gz")
    bamNormal = glob(normalOutDir + "/*.bam")
    bamTumor = glob(sampleOutDir + "/*_tumor.bam")
    call('dragen \
             -r /mnt/MTP_Storage_SSD/genome/dragenV4_hg38_graph/ \
             --tumor-bam-input ' + bamTumor[0] + ' \
             --bam-input ' + bamNormal[0] + ' \
             --output-directory ' + sampleOutDir + ' \
             --output-file-prefix ' + sampleName + ' \
             --enable-map-align false --enable-map-align-output false \
             --enable-variant-caller true --vc-enable-joint-detection true \
             --enable-sv true \
             --qc-somatic-contam-vcf /opt/edico/config/somatic_sample_cross_contamination_resource_hg38.vcf.gz \
             --vc-systematic-noise /mnt/MTP_Storage_SSD/hg38_databases/systematic-noise-baseline-collection-1.1.0/snv_wgs_hg38_max_v1.1_systematic_noise.bed.gz \
             --enable-tmb true \
             --qc-coverage-region-1 /mnt/MTP_Storage_SSD/hg38_databases/hg38_cds_sorted_merged.bed \
             --qc-coverage-tag-1 tmb --qc-coverage-reports-1=callability --vc-callability-tumor-thresh 30 \
             --enable-variant-annotation true \
             --variant-annotation-assembly GRCh38 \
             --variant-annotation-data /mnt/MTP_Storage_SSD/hg38_databases/nirvana_hg38_data_20240202/ \
             --msi-command tumor-normal --msi-coverage-threshold 30 --msi-microsatellites-file /mnt/MTP_Appl/MTP_AP_Research/databases/hg38/microsatellite_edited.list \
             --enable-cnv=true --cnv-use-somatic-vc-vaf true \
             --enable-hrd true \
             --cnv-normal-cnv-vcf ' +  cnvNormal[0] + ' --cnv-normal-b-allele-vcf ' + vcfNormal[0] + ' \
             --enable-variant-deduplication true --vc-combine-phased-variants-distance 15 \
             2> ' + sampleOutDir + '/' + sampleName + '_run_analysis.log' ,shell = True)

def main():
    args = params()
    sampleMatch = readMatchFile(args.matchFile)
    for sampleName in sampleMatch:
        print(sampleName)
        try:
            #Perform variant calling in normal sample
            normalSample = sampleMatch[sampleName]
            normalFastqListFile = createFastqListFile(normalSample,args.bclFolder,args.inputDir, args.outDir, args.globSequencer)
            if not args.germlineOnly:
                tumorFastqListFile = createFastqListFile(sampleName,args.bclFolder,args.inputDir, args.outDir, args.globSequencer, runIDs=args.runIDsTumor)
            if not args.dryRun:
                if args.normalOutDir is None:
                    pass
                    normalOutDir = runDragenNormal(normalSample,normalFastqListFile,args.outDir)
                else:
                    normalOutDir=glob(args.normalOutDir + "/" + normalSample)[0]
                if not args.germlineOnly:
                    #Perform tumor-normal variant calling
                    runDragenTumorNormal(sampleName,tumorFastqListFile,normalFastqListFile,normalOutDir,args.outDir)
                updateMatchFile(args.matchFile,sampleName)
        except Exception as error:
            print("Sample failed:", sampleName)
            print("Error message:", error) 
            if not args.dryRun:
                updateMatchFile(args.matchFile,sampleName,fail=True)
            continue

if __name__== "__main__":
    main()
