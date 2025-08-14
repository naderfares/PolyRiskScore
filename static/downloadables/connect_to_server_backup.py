import json
import subprocess
import requests
import os
import os.path
import time
import datetime
import hashlib
from sys import argv
from io import TextIOWrapper
import zipfile
import tarfile
import gzip
import myvariant
from Bio.Seq import Seq

def get_server_last_update_or_none(url, params):
    """
    Safely queries the server for the last update date.
    Returns a datetime.date object, or None if server is unavailable.
    """
    try:
        response = requests.get(url=url, params=params, timeout=10)
        response.close()
        if response.status_code == 200:
            parts = response.text.split('-')
            if len(parts) == 3:
                return datetime.date(int(parts[0]), int(parts[1]), int(parts[2]))
    except Exception as e:
        print(f"[WARN] Could not contact server ({url}): {e}")
    return None


# get the associations and clumps from the Server
def retrieveAssociationsAndClumps(refGen, traits, studyTypes, studyIDs, ethnicity, valueTypes, sexes, superPop, fileHash, extension, mafCohort):
    checkInternetConnection()
    
    # Clean up old cache files on startup
    cleanup_cache()
    
    # Show cache statistics
    cache_stats = get_cache_stats()
    print(f"[LOG] Cache stats - POST: {cache_stats['post_cache']} files, GET: {cache_stats['get_cache']} files, Total: {cache_stats['total_size_mb']:.2f} MB")

    # if the extension is .txt and the mafCohort is user -- Fail this is not a valid combination
    if extension == '.txt' and mafCohort == 'user':
        raise SystemExit('\nIn order to use the "user" option for maf cohort, you must upload a vcf, not a txt file. Please upload a vcf instead, or select a different maf cohort option. \n\n')

    mafCohort = formatMafCohort(mafCohort)
    percentilesCohort = mafCohort
    if mafCohort.startswith("adni"):
        mafCohort = "adni"

    if (ethnicity is not None):
        ethnicity = [sub.replace('_', ' ').replace('"', '').lower() for sub in ethnicity]
        availableEthnicities = getCachedEthnicities()
        if (not bool(set(ethnicity) & set(availableEthnicities)) and studyIDs is None):
            raise SystemExit('\nThe ethnicities requested are invalid. \nPlease use an ethnicity option from the list: \n\n{}'.format(availableEthnicities))

    workingFilesPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".workingFiles")
    associationsPath = ""

    # if the directory doesn't exist, make it, and we will need to download the files
    if not os.path.exists(workingFilesPath):
        os.mkdir(workingFilesPath)

    # Create the set that will include all the user-preferred super populations that correspond to the associations
    allSuperPops = set()
    # if we need to download a new all associations file, write to file
    associFileName = "allAssociations_{refGen}.txt".format(refGen=refGen)
    associationsPath = os.path.join(workingFilesPath, associFileName)

    allSuperPops |= set(['AFR', 'AMR', 'EAS', 'EUR', 'SAS'])

    if (checkForAllAssociFile(refGen)):
        associationsReturnObj = getAllAssociations(refGen)
        studySnpsPath = os.path.join(workingFilesPath, "traitStudyIDToSnps.txt")
        studySnpsData = getAllStudySnps()
        possibleAllelesPath = os.path.join(workingFilesPath, "allPossibleAlleles.txt")
        possibleAllelesData = getAllPossibleAlleles()
    
    dwnldNewMAFFile, mafPathExists = checkForAllMAFFiles(mafCohort, refGen)
    if (dwnldNewMAFFile):
        mafPath = os.path.join(workingFilesPath, "{m}_maf_{r}.txt".format(m=mafCohort, r=refGen))
        mafData = getAllMaf(mafCohort, refGen)

    if (checkForAllPercentilesFiles(percentilesCohort)):
        percentilesPath = os.path.join(workingFilesPath, "allPercentiles_{c}.txt".format(c=percentilesCohort))
        percentileData = getAllPercentiles(percentilesCohort)

    # check to see if associationsReturnObj is instantiated in the local variables
    if 'associationsReturnObj' in locals():
        f = open(associationsPath, 'w', encoding="utf-8")
        f.write(json.dumps(associationsReturnObj, indent=4))
        f.close()

    if 'mafData' in locals():
        f = open(mafPath, 'w', encoding="utf-8")
        f.write(json.dumps(mafData, indent=4))
        f.close()
    elif mafCohort != 'user' and not mafPathExists:
        raise SystemExit("ERROR: We were not able to retrieve the Minor Allele Frequency data at this time. Please try again.")

    if 'percentileData' in locals():
        f = open(percentilesPath, 'w', encoding='utf-8')
        f.write(json.dumps(percentileData, indent=4))
        f.close()

    if 'possibleAllelesData' in locals():
        f = open(possibleAllelesPath, 'w', encoding="utf-8")
        f.write(json.dumps(possibleAllelesData, indent=4))
        f.close()

    # check to see if studySnpsData is instantiated in the local variables
    if 'studySnpsData' in locals():
        f = open(studySnpsPath, 'w', encoding="utf-8")
        f.write(json.dumps(studySnpsData, indent=4))
        f.close()
    
    for pop in allSuperPops:
        if (checkForAllClumps(pop, refGen)):
            clumpsPath = os.path.join(workingFilesPath, "{p}_clumps_{r}.txt".format(p=pop, r=refGen))
            clumpsData = getAllClumps(refGen, pop)
            
        # check to see if clumpsData is instantiated in the local variables
        if 'clumpsData' in locals():
            f = open(clumpsPath, 'w', encoding="utf-8")
            f.write(json.dumps(clumpsData, indent=4))
            f.close()

    return


# format the uploaded GWAS data and get the clumps from the server
def formatGWASAndRetrieveClumps(GWASfile, userGwasBeta, GWASextension, GWASrefGen, refGen, superPop, mafCohort, fileHash, extension):
    print(f"[LOG] Starting formatGWASAndRetrieveClumps with:")
    print(f"[LOG]   GWASfile: {GWASfile}")
    print(f"[LOG]   userGwasBeta: {userGwasBeta}")
    print(f"[LOG]   GWASrefGen: {GWASrefGen}, refGen: {refGen}")
    print(f"[LOG]   superPop: {superPop}")
    print(f"[LOG]   mafCohort: {mafCohort}")
    print(f"[LOG]   fileHash: {fileHash}")
    print(f"[LOG]   extension: {extension}")
    
    checkInternetConnection()
    
    # Clean up old cache files on startup
    cleanup_cache()
    
    # Show cache statistics
    cache_stats = get_cache_stats()
    print(f"[LOG] Cache stats - POST: {cache_stats['post_cache']} files, GET: {cache_stats['get_cache']} files, Total: {cache_stats['total_size_mb']:.2f} MB")

    # if the extension is .txt and the mafCohort is user -- Fail this is not a valid combination
    if extension == '.txt' and mafCohort == 'user':
        raise SystemExit('\nIn order to use the "user" option for maf cohort, you must upload a vcf, not a txt file. Please upload a vcf instead, or select a different maf cohort option. \n\n')

    GWASfileOpen = openFileForParsing(GWASfile, True)

    allSuperPops = set()

    associationDict = {}
    chromSnpDict = {}
    studyIDsToMetaData = {}
    studySnpsData = {}

    sii = -1 # studyID index
    ti = -1 # trait index
    spi = -1 # super population index
    si = -1 # snp index
    ci = -1 # chromosome index
    pi = -1 # position index
    rai = -1 # risk allele index
    ori = -1 # odds ratio index
    bvi = -1 # beta value index
    bui = -1 # beta unit index
    pvi = -1 # p-value index
    spi = -1 # super population index
    cti = -1 # citation index
    rti = -1 # reported trait index
    pvai = -1 # pvalue annotation index
    bai = -1 # beta annotation index

    firstLine = True
    duplicatesSet = set()

    for line in GWASfileOpen:
        line = line.strip()
        if len(line) == 0: # skip lines that don't have content
            continue
        if firstLine:
            firstLine = False
            headers = line.lower().split("\t")

            try:
                sii = headers.index("study id")
                ti = headers.index("trait")
                si = headers.index("rsid")
                ci = headers.index("chromosome")
                pi = headers.index("position")
                rai = headers.index("risk allele")
                if userGwasBeta:
                    bvi = headers.index("beta coefficient")
                    bui = headers.index("beta units")
                else:
                    ori = headers.index("odds ratio")
                pvi = headers.index("p-value")
                spi = headers.index("super population")
            except ValueError:
                raise SystemExit("ERROR: The GWAS file format is not correct. Please check your file to ensure the required columns are present in a tab separated format. Additionally, check your column names and ensure that there are no extra spaces in the names and that your spelling is correct.")

            cti = headers.index("citation") if "citation" in headers else -1
            rti = headers.index("reported trait") if "reported trait" in headers else -1
            pvai = headers.index("p-value annotation") if "p-value annotation" in headers else -1
            bai = headers.index("beta annotation") if "beta annotation" in headers else -1

        else:
            line = line.split("\t")
            # Add super population to the super population set
            preferredPop = getPreferredPop(line[spi], superPop)
            # Add super population to the super population set
            allSuperPops.add(preferredPop)
            # create the chrom:pos to snp dict
            # Format chromPos to match BCF format (with 'chr' prefix if missing)
            chrom = line[ci]
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom
            chromPos = ":".join([chrom, line[pi]])
            
            # Store original rsID for reference
            original_rsid = line[si]
            if chromPos not in chromSnpDict:
                # Map chromPos to itself (chromPos will be the primary identifier)
                chromSnpDict[chromPos] = chromPos
            
            # create the snp to associations stuff dict
            # Use chromPos as the primary key instead of rsID
            if chromPos not in associationDict:
                associationDict[chromPos] = {
                    "pos": chromPos,
                    "original_rsid": original_rsid,  # Store original rsID as metadata
                    "traits": {}
                }
            # if trait not in associationsDict[chromPos][traits]
            if line[ti] not in associationDict[chromPos]["traits"]:
                associationDict[chromPos]["traits"][line[ti]] = {}
            # if studyID not in associationDict[chromPos]["traits"][trait]
            if line[sii] not in associationDict[chromPos]["traits"][line[ti]]:
                associationDict[chromPos]["traits"][line[ti]][line[sii]] = {}
            # if pvalannotation not in associationDict[chromPos]["traits"][line[ti]][line[sii]]
            pValueAnnotation = line[pvai] if pvai != -1 and pvai < len(line) else "NA"
            betaAnnotation = line[bai] if bai != -1 and bai < len(line) else "NA"
            valueType = "beta" if userGwasBeta else "OR"
            pvalBetaAnnoValType = pValueAnnotation + "|" + betaAnnotation + "|" + valueType
            if pvalBetaAnnoValType not in associationDict[chromPos]["traits"][line[ti]][line[sii]]:
                associationDict[chromPos]["traits"][line[ti]][line[sii]][pvalBetaAnnoValType] = {}
            
            riskAllele = runStrandFlipping(original_rsid, line[rai])
            if riskAllele not in associationDict[chromPos]["traits"][line[ti]][line[sii]][pvalBetaAnnoValType]:
                associationDict[chromPos]["traits"][line[ti]][line[sii]][pvalBetaAnnoValType][riskAllele]= {
                    "pValue": float(line[pvi]),
                    "sex": "NA",
                    "ogValueTypes": 'beta' if userGwasBeta else 'OR'
                }
                if userGwasBeta:
                    associationDict[chromPos]["traits"][line[ti]][line[sii]][pvalBetaAnnoValType][riskAllele]['betaValue'] = float(line[bvi])
                    associationDict[chromPos]["traits"][line[ti]][line[sii]][pvalBetaAnnoValType][riskAllele]['betaUnit'] = line[bui]
                else:
                    associationDict[chromPos]["traits"][line[ti]][line[sii]][pvalBetaAnnoValType][riskAllele]['oddsRatio'] = float(line[ori])
                    associationDict[chromPos]["traits"][line[ti]][line[sii]][pvalBetaAnnoValType][riskAllele]['betaUnit'] = 'NA'
        
            else:
                # if the snp is duplicated, notify the user and exit
                raise SystemExit("ERROR: The GWAS file contains at least one duplicated snp for the following combination. {}, {}, {}, {}, . \n Please ensure that there is only one snp for each combination.".format(chromPos, line[ti], line[sii], pvalBetaAnnoValType))

            # create the metadata info dict
            # if the studyID is not in the studyIDsToMetaData
            if line[sii] not in studyIDsToMetaData:
                studyIDsToMetaData[line[sii]] = {
                    "citation": line[cti] if cti != -1 and cti < len(line) else "",
                    "reportedTrait": line[rti] if rti != -1 and rti < len(line) else "",
                    "studyTypes": [],
                    "traits": {},
                    "ethnicity": []
                }
            # if the trait is not in the studyIDsToMetaData[studyID]["traits"]
            if line[ti] not in studyIDsToMetaData[line[sii]]["traits"]:
                # add the trait
                studyIDsToMetaData[line[sii]]["traits"][line[ti]] = {
                    "studyTypes": [],
                    "pValBetaAnnoValType": [pvalBetaAnnoValType],
                    "superPopulations": [line[spi]]
                }
            else:
                studyIDsToMetaData[line[sii]]["traits"][line[ti]]['pValBetaAnnoValType'].append(pvalBetaAnnoValType)
                
            # create studyID/trait/pValueAnnotation to snps
            # if trait|studyID|pValueAnnotation not in the studySnpsData
            joinList = [line[ti], pvalBetaAnnoValType, line[sii]]
            traitStudyIDPValAnno = "|".join(joinList)
            if traitStudyIDPValAnno not in studySnpsData:
                studySnpsData[traitStudyIDPValAnno] = []
            # add chromPos to the traitStudyIDToSnp (instead of rsID)
            studySnpsData[traitStudyIDPValAnno].append(chromPos)

    GWASfileOpen.close()
    
    print(f"[LOG] GWAS file parsing completed:")
    print(f"[LOG]   Parsed {len(associationDict)} SNPs from GWAS file")
    print(f"[LOG]   Found {len(studyIDsToMetaData)} unique studies")
    print(f"[LOG]   Found {len(studySnpsData)} trait/study combinations")
    print(f"[LOG]   Identified {len(allSuperPops)} super populations: {allSuperPops}")

    # remove duplicated associations
    duplicates_removed = len(duplicatesSet)
    for (snp, trait, studyID) in duplicatesSet:
        del associationDict[snp]["traits"][trait][studyID]
        if associationDict[snp]["traits"][trait] == {}:
            del associationDict[snp]["traits"][trait]
            if associationDict[snp]["traits"] == {}:
                del associationDict[snp]
    
    if duplicates_removed > 0:
        print(f"[LOG] Removed {duplicates_removed} duplicate associations")

    # if the samples reference genome does not equal the gwas reference genome, get a dictionary with the correct positions
    chromSnpDict = {}
    if GWASrefGen != refGen:
        print(f"[LOG] Converting SNP positions from {GWASrefGen} to {refGen}")
        snps = list(associationDict.keys())
        
        # Use retry logic for SNP position conversion
        chromSnpDict = retry_network_call(
            getUrlWithParams,
            "https://prs.byu.edu/snps_to_chrom_pos", 
            { "snps": snps, "refGen": refGen }
        )
        print(f"[LOG] Got position mappings for {len(chromSnpDict)} chromosome positions")

    mergedAssociDict = dict()
    mergedAssociDict.update(associationDict)
    mergedAssociDict.update(chromSnpDict)
    chromPos = list(chromSnpDict.keys())
    print(f"[LOG] Created merged associations dict with {len(mergedAssociDict)} entries")
    print(f"[LOG] chromPos contains {len(chromPos)} chromosome positions")

    associationsReturnObj = {
        "associations": mergedAssociDict,
        "studyIDsToMetaData": studyIDsToMetaData
    }
    print(f"[LOG] Created associationsReturnObj with {len(associationsReturnObj['associations'])} associations")
        
    workingFilesPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".workingFiles")
    # if the directory doesn't exist, make it, and we will need to download the files
    if not os.path.exists(workingFilesPath):
        os.mkdir(workingFilesPath)
    
    fileName = "GWASassociations_{fileHash}.txt".format(fileHash=fileHash)
    associationsPath = os.path.join(workingFilesPath, fileName)


    # Access and write the clumps file for each of the super populations preferred in the GWAS file
    print(f"[LOG] Processing clumps for {len(allSuperPops)} super populations: {allSuperPops}")
    clumps_written = 0
    for pop in allSuperPops:
        fileName = "{pop}_clumps_{refGen}_{fileHash}.txt".format(pop=pop, refGen=refGen, fileHash=fileHash)
        clumpsPath = os.path.join(workingFilesPath, fileName)
        print(f"[LOG] Getting clumps data for population {pop}")

        try:
            # get clumps using the refGen and superpopulation
            clumpsData = getClumps(refGen, pop, chromPos)

            # Validate clumpsData before writing
            if clumpsData is None:
                print(f"[ERROR] clumpsData is None for population {pop}")
                continue
            elif not clumpsData:
                print(f"[WARN] clumpsData is empty for population {pop}")
            else:
                print(f"[LOG] Retrieved {len(clumpsData)} clump entries for population {pop}")

            # Write clumps file with proper error handling
            try:
                with open(clumpsPath, 'w', encoding="utf-8") as f:
                    f.write(json.dumps(clumpsData))
                print(f"[LOG] Successfully wrote clumps file: {clumpsPath}")
                
                # Verify file was written correctly
                if os.path.exists(clumpsPath) and os.path.getsize(clumpsPath) > 0:
                    clumps_written += 1
                    print(f"[LOG] Verified clumps file exists and has size: {os.path.getsize(clumpsPath)} bytes")
                else:
                    print(f"[ERROR] Clumps file validation failed: {clumpsPath}")
                    
            except Exception as e:
                print(f"[ERROR] Failed to write clumps file {clumpsPath}: {e}")
                
        except Exception as e:
            print(f"[ERROR] Failed to get clumps data for population {pop}: {e}")
    
    print(f"[LOG] Successfully wrote {clumps_written}/{len(allSuperPops)} clumps files")

    #get maf data using the refgen and mafcohort
    fileName = "{m}_maf_{ahash}.txt".format(m=mafCohort, ahash=fileHash)
    mafPath = os.path.join(workingFilesPath, fileName)
    print(f"[LOG] Getting MAF data for cohort {mafCohort}")
    
    try:
        mafData = getMaf(mafCohort, refGen, chromPos)
        if mafData is None:
            print(f"[ERROR] mafData is None for cohort {mafCohort}")
        elif not mafData:
            print(f"[WARN] mafData is empty for cohort {mafCohort}")
        else:
            print(f"[LOG] Retrieved {len(mafData)} MAF entries for cohort {mafCohort}")
    except Exception as e:
        print(f"[ERROR] Failed to get MAF data for cohort {mafCohort}: {e}")
        mafData = None

    # get the study:snps info
    fileName = "traitStudyIDToSnps_{ahash}.txt".format(ahash = fileHash)
    studySnpsPath = os.path.join(workingFilesPath, fileName)
    print(f"[LOG] StudySnps data contains {len(studySnpsData)} trait/study combinations")

    # get the possible alleles for snps
    fileName = "possibleAlleles_{ahash}.txt".format(ahash = fileHash)
    possibleAllelesPath = os.path.join(workingFilesPath, fileName)
    print(f"[LOG] Getting possible alleles for {len(associationsReturnObj['associations'])} SNPs")
    
    try:
        possibleAllelesData = getPossibleAlleles(list(associationsReturnObj['associations'].keys()))
        if possibleAllelesData is None:
            print(f"[ERROR] possibleAllelesData is None")
        elif not possibleAllelesData:
            print(f"[WARN] possibleAllelesData is empty")
        else:
            print(f"[LOG] Retrieved possible alleles for {len(possibleAllelesData)} SNPs")
    except Exception as e:
        print(f"[ERROR] Failed to get possible alleles data: {e}")
        possibleAllelesData = None


    # Write associations file with validation
    print(f"[LOG] Writing associations file: {associationsPath}")
    try:
        if associationsReturnObj and 'associations' in associationsReturnObj and 'studyIDsToMetaData' in associationsReturnObj:
            associations_count = len(associationsReturnObj['associations'])
            studies_count = len(associationsReturnObj['studyIDsToMetaData'])
            print(f"[LOG] AssociationsReturnObj contains {associations_count} associations and {studies_count} studies")
            
            with open(associationsPath, 'w', encoding="utf-8") as f:
                f.write(json.dumps(associationsReturnObj))
            
            # Verify file was written correctly
            if os.path.exists(associationsPath) and os.path.getsize(associationsPath) > 0:
                print(f"[LOG] Successfully wrote associations file: {associationsPath} (size: {os.path.getsize(associationsPath)} bytes)")
            else:
                raise Exception("Associations file validation failed")
        else:
            raise Exception("associationsReturnObj is invalid or missing required keys")
    except Exception as e:
        print(f"[ERROR] Failed to write associations file {associationsPath}: {e}")
        raise SystemExit(f"CRITICAL ERROR: Could not write associations file: {e}")

    # Write MAF file with validation
    if mafData is not None:
        print(f"[LOG] Writing MAF file: {mafPath}")
        try:
            with open(mafPath, 'w', encoding='utf-8') as f:
                f.write(json.dumps(mafData))
            
            if os.path.exists(mafPath) and os.path.getsize(mafPath) > 0:
                print(f"[LOG] Successfully wrote MAF file: {mafPath} (size: {os.path.getsize(mafPath)} bytes)")
            else:
                print(f"[WARN] MAF file validation failed: {mafPath}")
        except Exception as e:
            print(f"[ERROR] Failed to write MAF file {mafPath}: {e}")

    # Write studySnps file with validation  
    print(f"[LOG] Writing studySnps file: {studySnpsPath}")
    try:
        if studySnpsData and len(studySnpsData) > 0:
            print(f"[LOG] StudySnpsData contains {len(studySnpsData)} entries")
            
            with open(studySnpsPath, 'w', encoding="utf-8") as f:
                f.write(json.dumps(studySnpsData))
            
            # Verify file was written correctly
            if os.path.exists(studySnpsPath) and os.path.getsize(studySnpsPath) > 0:
                print(f"[LOG] Successfully wrote studySnps file: {studySnpsPath} (size: {os.path.getsize(studySnpsPath)} bytes)")
            else:
                raise Exception("StudySnps file validation failed")
        else:
            raise Exception("studySnpsData is empty or None")
    except Exception as e:
        print(f"[ERROR] Failed to write studySnps file {studySnpsPath}: {e}")
        raise SystemExit(f"CRITICAL ERROR: Could not write studySnps file: {e}")

    # Write possible alleles file with validation
    if possibleAllelesData is not None:
        print(f"[LOG] Writing possible alleles file: {possibleAllelesPath}")
        try:
            with open(possibleAllelesPath, 'w', encoding="utf-8") as f:
                f.write(json.dumps(possibleAllelesData))
            
            if os.path.exists(possibleAllelesPath) and os.path.getsize(possibleAllelesPath) > 0:
                print(f"[LOG] Successfully wrote possible alleles file: {possibleAllelesPath} (size: {os.path.getsize(possibleAllelesPath)} bytes)")
            else:
                print(f"[WARN] Possible alleles file validation failed: {possibleAllelesPath}")
        except Exception as e:
            print(f"[ERROR] Failed to write possible alleles file {possibleAllelesPath}: {e}")
    
    print(f"[LOG] Completed writing all required files for fileHash {fileHash}")

    return


# opens and returns an open file from the inputFile path, using zipfile, tarfile, gzip, or open depending on the file's type
# assumes the file is valid (validated with getZippedFileExtension function)
def openFileForParsing(inputFile, isGWAS=False):
    filename = ""
    if zipfile.is_zipfile(inputFile):
        # open the file
        archive = zipfile.ZipFile(inputFile)
        # get the vcf or txt file in the zip
        for filename in archive.namelist():
            extension = filename[-4:].lower()
            if (not isGWAS and extension == ".txt" or extension == ".vcf") or (isGWAS and extension == ".txt" or extension == ".tsv"):
                # TextIOWrapper converts bytes to strings and force_zip64 is for files potentially larger than 2GB
                return TextIOWrapper(archive.open(filename, force_zip64=True))
    elif tarfile.is_tarfile(inputFile):
        # open the file
        archive = tarfile.open(inputFile)
        # get the vcf or txt file in the tar
        for tarInfo in archive:
            extension = tarInfo.name[-4:].lower()
            if (not isGWAS and extension == ".txt" or extension == ".vcf") or (isGWAS and extension == ".txt" or extension == ".tsv"):
                # TextIOWrapper converts bytes to strings
                return TextIOWrapper(archive.extractfile(tarInfo))
    elif inputFile.lower().endswith(".gz") or inputFile.lower().endswith(".gzip"):
        return TextIOWrapper(gzip.open(inputFile, 'r'))
    elif inputFile.lower().endswith(".bcf"):
        # Abrir BCF como VCF texto via bcftools
        cmd = ["bcftools", "view", inputFile, "-Ov"]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
        return process.stdout
    else:
        # default option for regular vcf and txt files
        return open(inputFile, 'r')


def checkForAllAssociFile(refGen, max_age_days=30):
    """
    Checks if the cached associations file needs to be refreshed.
    If it's older than max_age_days, always redownload (don't check server).
    If file is fresh (<max_age_days), just use cache (don't check server).
    """
    print("USING UPDATED checkForAllAssociFile")  # Debug
    scriptPath = os.path.dirname(os.path.abspath(__file__))
    workingFilesPath = os.path.join(scriptPath, ".workingFiles")
    associFileName = f"allAssociations_{refGen}.txt"
    allAssociationsFile = os.path.join(workingFilesPath, associFileName)

    if not os.path.exists(allAssociationsFile):
        print(f"[LOG] Cache file {associFileName} does not exist, will download.")
        return True

    file_mod_time = os.path.getmtime(allAssociationsFile)
    file_age_days = (time.time() - file_mod_time) / 86400
    print(f"[LOG] Cache file {associFileName} is {file_age_days:.8f} days old.")

    if file_age_days > max_age_days:
        print(f"[LOG] Cache file {associFileName} is older than {max_age_days} days, will redownload.")
        return True

    # File is fresh enough, just use it
    print(f"[LOG] Cache file {associFileName} is fresh (<= {max_age_days} days), will use cache.")
    return False

def checkForAllClumps(pop, refGen):
    dnldNewClumps = True
    scriptPath = os.path.dirname(os.path.abspath(__file__))
    workingFilesPath = os.path.join(scriptPath, ".workingFiles")
    allClumpsFile = os.path.join(workingFilesPath, "{0}_clumps_{1}.txt".format(pop, refGen))

    if os.path.exists(allClumpsFile):
        params = {"refGen": refGen, "superPop": pop}
        server_update = get_server_last_update_or_none("https://prs.byu.edu/last_clumps_update", params)
        fileModDateObj = time.localtime(os.path.getmtime(allClumpsFile))
        fileModDate = datetime.date(fileModDateObj.tm_year, fileModDateObj.tm_mon, fileModDateObj.tm_mday)
        if server_update is not None:
            if (server_update <= fileModDate):
                dnldNewClumps = False
        else:
            print("[WARN] Server unavailable for Clumps check. Using local cache if present.")
            dnldNewClumps = False
    return dnldNewClumps


def checkForAllMAFFiles(mafCohort, refGen):
    dnldNewMaf = True
    pathExists = False
    scriptPath = os.path.dirname(os.path.abspath(__file__))
    workingFilesPath = os.path.join(scriptPath, ".workingFiles")
    allMAFfile = os.path.join(workingFilesPath, "{0}_maf_{1}.txt".format(mafCohort, refGen))

    if os.path.exists(allMAFfile):
        pathExists = True
        params = {"cohort": mafCohort, "refGen": refGen}
        server_update = get_server_last_update_or_none("https://prs.byu.edu/last_maf_update", params)
        fileModDateObj = time.localtime(os.path.getmtime(allMAFfile))
        fileModDate = datetime.date(fileModDateObj.tm_year, fileModDateObj.tm_mon, fileModDateObj.tm_mday)
        if server_update is not None:
            if (server_update <= fileModDate):
                dnldNewMaf = False
        else:
            print("[WARN] Server unavailable for MAF check. Using local cache if present.")
            dnldNewMaf = False  # Assume cache is ok
    return dnldNewMaf, pathExists



def checkForAllPercentilesFiles(percentilesCohort):
    dnldNewPercentiles = True
    scriptPath = os.path.dirname(os.path.abspath(__file__))
    workingFilesPath = os.path.join(scriptPath, ".workingFiles")
    allPercentilesfile = os.path.join(workingFilesPath, "allPercentiles_{}.txt".format(percentilesCohort))

    if os.path.exists(allPercentilesfile):
        params = {"cohort": percentilesCohort}
        server_update = get_server_last_update_or_none("https://prs.byu.edu/last_percentiles_update", params)
        fileModDateObj = time.localtime(os.path.getmtime(allPercentilesfile))
        fileModDate = datetime.date(fileModDateObj.tm_year, fileModDateObj.tm_mon, fileModDateObj.tm_mday)
        if server_update is not None:
            if (server_update <= fileModDate):
                dnldNewPercentiles = False
        else:
            print("[WARN] Server unavailable for Percentiles check. Using local cache if present.")
            dnldNewPercentiles = False
    return dnldNewPercentiles



# gets associations obj download from the Server
def getAllAssociations(refGen): 
    params = {
        "refGen": refGen,
    }
    associationsReturnObj = getUrlWithParams("https://prs.byu.edu/get_associations_download_file", params = params)
    # Organized with pos/snp as the Keys
    return associationsReturnObj


# gets the clumps file download from the server
def getAllClumps(refGen, superPop):
    params = {
        'refGen': refGen,
        'superPop': superPop
    }
    clumpsReturnObj = getUrlWithParams("https://prs.byu.edu/get_clumps_download_file", params=params)
    return clumpsReturnObj


def getAllMaf(mafCohort, refGen):
    if (mafCohort == 'user'): return {}
    params = {
        "cohort": mafCohort,
        "refGen": refGen
    }
    mafReturnedObj = getUrlWithParams("https://prs.byu.edu/get_maf_download_file", params=params)
    return mafReturnedObj


def getAllPercentiles(percentilesCohort):
    if (percentilesCohort == 'user'): return {}
    params = {
        "cohort": percentilesCohort
    }
    percentilesReturnedObj = getUrlWithParams("https://prs.byu.edu/get_percentiles_download_file", params=params)
    return percentilesReturnedObj


# gets study snps file download from the Server
# gets a list of snps for all of the unique trait/pValueAnnotation/betaAnnotation/valueType/studyID combinations
def getAllStudySnps(): 
    studySnpsReturnObj = getUrlWithParams("https://prs.byu.edu/get_traitStudyID_to_snp", params={})
    # Organized with study as the Keys and snps as values
    return studySnpsReturnObj


def getAllPossibleAlleles():
    possibleAllelesObj = getUrlWithParams("https://prs.byu.edu/get_all_possible_alleles", params={})
    return possibleAllelesObj


# This function is used to combine json from all the separate calls into one json object. Due to the amount of nesting in the json
# this is the neccesary way to properly combine
def combineJson(old, new):
    studyMeta = new["studyIDsToMetaData"]
    associations = new["associations"]

    for studyID in studyMeta:
        if studyID not in old["studyIDsToMetaData"]:
            old["studyIDsToMetaData"][studyID] = studyMeta[studyID]
        else:
            for trait in studyMeta[studyID]["traits"]:
                if trait not in old["studyIDsToMetaData"][studyID]["traits"]:
                    old["studyIDsToMetaData"][studyID]["traits"][trait] = studyMeta[studyID]["traits"][trait]

    for snp in associations:
        if snp not in old["associations"]:
            old["associations"][snp] = associations[snp]
        elif snp.startswith("rs"):
            for trait in associations[snp]["traits"]:
                if trait not in old["associations"][snp]["traits"]:
                    old["associations"][snp]["traits"][trait] = associations[snp]["traits"][trait]
                else:
                    for studyID in associations[snp]["traits"][trait]:
                        if studyID not in old["associations"][snp]["traits"][trait]:
                            old["associations"][snp]["traits"][trait][studyID] = associations[snp]["traits"][trait][studyID]
                        else:
                            for pValBetaAnnoValType in associations[snp]["traits"][trait][studyID]:
                                if pValBetaAnnoValType not in old["associations"][snp]["traits"][trait][studyID]:
                                    old["associations"][snp]["traits"][trait][studyID][pValBetaAnnoValType] = associations[snp]["traits"][trait][studyID][pValBetaAnnoValType]
    new = {}
    return old


# gets associationReturnObj using the given filters
def getSpecificAssociations(refGen, traits, studyTypes, studyIDs, ethnicity, valueTypes, sexes):
    finalStudyList = []
    associationData = {
        "studyIDsToMetaData" : {},
        "associations": {}
    }

    if (traits is not None or studyTypes is not None or ethnicity is not None or valueTypes is not None or sexes is not None):
        # get the studies matching the parameters
        body = {
            "traits": traits, 
            "studyTypes": studyTypes,
            "ethnicities": ethnicity,
            "sexes": sexes,
            "ogValueTypes": valueTypes
        }
        traitData = {**postUrlWithBody("https://prs.byu.edu/get_studies", body=body)}
        #TODO it looks like this is timing out (-e "Asian unspecified" -y beta -g female -k HI) Need to do something to handle this later

        # select the studyIDs of the studies
        for trait in traitData:
            for study in traitData[trait]:
                # if the studyID is in the studyIDs list, don't add it in here
                if (studyIDs is not None and study['studyID'] in studyIDs):
                    continue
                else:
                    finalStudyList.append(json.dumps({
                        "trait": trait,
                        "studyID": study['studyID'],
                        "pValueAnnotation": study['pValueAnnotation'],
                        "betaAnnotation": study['betaAnnotation'],
                        "ogValueTypes" : study['ogValueTypes']
                    }))

    # get the data for the specified studyIDs
    if (studyIDs is not None):
        params = {
            "studyIDs": studyIDs
        }
        studyIDDataList = getUrlWithParams("https://prs.byu.edu/get_studies_by_id", params = params)
        if studyIDDataList == []:
            print('\n\nWARNING, NO STUDIES MATCHED THE GIVEN STUDY ID(S): {}. \nTHIS MAY CAUSE THE PROGRAM TO QUIT IF THERE WERE NO OTHER FILTERS.\n'.format(studyIDs))

        for i in range(len(studyIDDataList)):
            # add the specified studyIDs to the set of studyIDObjs
            finalStudyList.append(json.dumps({
                "trait": studyIDDataList[i]['trait'],
                "studyID": studyIDDataList[i]['studyID'],
                "pValueAnnotation": studyIDDataList[i]['pValueAnnotation'],
                "betaAnnotation": studyIDDataList[i]['betaAnnotation'],
                "ogValueTypes" : studyIDDataList[i]['ogValueTypes']
            }))

    if finalStudyList == []:
        raise SystemExit("No studies with those filters exist because your filters are too narrow or invalid. Check your filters and try again.")

    # Breaking up the calls to the get_associations endpoint so that we don't got over the Request size limit
    try:
        lengthOfList = len(finalStudyList)
        i = 0
        j = 1000 if lengthOfList > 1000 else lengthOfList
        print("Getting associations and studies")
        print("Total number of studies: {}". format(lengthOfList))
        runLoop = True
        while runLoop:
            if j == lengthOfList:
                runLoop = False
            # get the associations based on the studyIDs
            print("{}...".format(j), end = "", flush=True)
            body = {
                "refGen": refGen,
                "studyIDObjs": finalStudyList[i:j],
                "sexes": sexes,
                "ogValueType": valueTypes
            }
            tmpAssociationsData = postUrlWithBody("https://prs.byu.edu/get_associations", body=body)
            associationData = combineJson(associationData, tmpAssociationsData)
            i = j
            j = j + 1000 if lengthOfList > j + 1000 else lengthOfList
        print('Done\n')
    except AssertionError:
        raise SystemExit("ERROR: 504 - Connection to the server timed out")
    return associationData, finalStudyList


def runStrandFlipping(snp, allele):
    import myvariant
    from Bio.Seq import Seq

    mv = myvariant.MyVariantInfo()

    possibleAlleles = getVariantAlleles(snp, mv)
    riskAllele = Seq(allele)
    if riskAllele not in possibleAlleles:
        complement = riskAllele.reverse_complement()
        if complement in possibleAlleles:
            print("WE MADE A SWITCH", snp, riskAllele, complement)
            riskAllele = complement
    
    return str(riskAllele)


def getPossibleAlleles(snpsFromAssociations):
    mv = myvariant.MyVariantInfo()
    snpsToPossibleAlleles = {}

    for snp in snpsFromAssociations:
        if snp.startswith("rs"):
            possibleAlleles = getVariantAlleles(snp, mv)
            snpsToPossibleAlleles[snp] = possibleAlleles
        elif ":" in snp:
            # For chromPos identifiers, we can't query MyVariant
            # So we'll return empty alleles list (will be handled elsewhere)
            snpsToPossibleAlleles[snp] = []
    return snpsToPossibleAlleles


def getVariantAlleles(rsID, mv):
    import contextlib, io

    f=io.StringIO()
    with contextlib.redirect_stdout(f):
        queryResult = mv.query('dbsnp.rsid:{}'.format(rsID), fields='dbsnp.alleles.allele, dbsnp.dbsnp_merges, dbsnp.gene.strand, dbsnp.alt, dbsnp.ref')
    output = f.getvalue()

    objs = queryResult['hits'][0] if len(queryResult['hits']) > 0 else None
    objsList = []
    if isinstance(objs, dict):
        objsList.append(objs)
        objs = objsList

    alleles = set()
    if objs is not None:
        for obj in objs:
            if ('alleles' in obj['dbsnp']):
                for alleleObj in obj['dbsnp']['alleles']:
                    alleles.add(alleleObj['allele'])
            if ('ref' in obj['dbsnp'] and obj['dbsnp']['ref'] != ""):
                alleles.add(obj['dbsnp']['ref'])
            if ('alt' in obj['dbsnp'] and obj['dbsnp']['alt'] != ""):
                alleles.add(obj['dbsnp']['alt'])
            if (len(alleles) == 0):
                print(obj, "STILL NO ALLELES")
    else:
        # TODO maybe: try to find it with a merged snp?
        pass

    return list(alleles)


# for POST urls with caching
def postUrlWithBody(url, body, max_age_hours=24):
    import json
    
    # Create cache directory if it doesn't exist
    workingFilesPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".workingFiles")
    if not os.path.exists(workingFilesPath):
        os.mkdir(workingFilesPath)
    
    cacheDir = os.path.join(workingFilesPath, "post_cache")
    if not os.path.exists(cacheDir):
        os.mkdir(cacheDir)
    
    # Create cache key from URL and body
    cache_content = f"{url}|{json.dumps(body, sort_keys=True, ensure_ascii=False)}"
    cache_key = hashlib.md5(cache_content.encode('utf-8')).hexdigest()
    cache_file = os.path.join(cacheDir, f"post_{cache_key}.json")
    
    # Check if cache exists and is fresh
    if os.path.exists(cache_file):
        file_mod_time = os.path.getmtime(cache_file)
        file_age_hours = (time.time() - file_mod_time) / 3600
        
        if file_age_hours < max_age_hours:
            print(f"[LOG] Using cached response for POST {url} (age: {file_age_hours:.2f}h)")
            try:
                with open(cache_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except Exception as e:
                print(f"[WARN] Cache file corrupted, will fetch fresh: {e}")
                # Continue to fetch fresh data if cache is corrupted
    
    # Fetch fresh data
    print(f"[LOG] POST URL: {url}")
    print(f"[LOG] POST BODY: {json.dumps(body, ensure_ascii=False)}")

    response = requests.post(url=url, data=body)
    response.close()
    if response.status_code == 504:
        print("\n*** The connection timed out. If you haven't already, try running the first step with no additional filters, then running the second step with the filters.")
        print("(See the README file for an example -- under Applying Step Numbers)\n")
    assert (response), "Error connecting to the server: {0} - {1}".format(response.status_code, response.reason) 
    
    result = {}
    if response.status_code == 204:
        result = {}
    else:
        result = response.json()
    
    # Cache the result
    try:
        with open(cache_file, 'w', encoding='utf-8') as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        print(f"[LOG] Cached POST response to {cache_file}")
    except Exception as e:
        print(f"[WARN] Failed to cache response: {e}")
    
    return result 



# for GET urls with caching
def getUrlWithParams(url, params, max_age_hours=24):
    import urllib.parse
    import json
    
    # Create cache directory if it doesn't exist
    workingFilesPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".workingFiles")
    if not os.path.exists(workingFilesPath):
        os.mkdir(workingFilesPath)
    
    cacheDir = os.path.join(workingFilesPath, "get_cache")
    if not os.path.exists(cacheDir):
        os.mkdir(cacheDir)
    
    # Create cache key from URL and params
    full_url = url + '?' + urllib.parse.urlencode(params, doseq=True)
    cache_key = hashlib.md5(full_url.encode('utf-8')).hexdigest()
    cache_file = os.path.join(cacheDir, f"get_{cache_key}.json")
    
    # Check if cache exists and is fresh
    if os.path.exists(cache_file):
        file_mod_time = os.path.getmtime(cache_file)
        file_age_hours = (time.time() - file_mod_time) / 3600
        
        if file_age_hours < max_age_hours:
            print(f"[LOG] Using cached response for GET {url} (age: {file_age_hours:.2f}h)")
            try:
                with open(cache_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except Exception as e:
                print(f"[WARN] Cache file corrupted, will fetch fresh: {e}")
                # Continue to fetch fresh data if cache is corrupted
    
    # Fetch fresh data
    print(f"[LOG] GET URL: {full_url}")

    response = requests.get(url=url, params=params)
    response.close()
    assert (response), "Error connecting to the server: {0} - {1}".format(response.status_code, response.reason) 
    
    result = response.json()
    
    # Cache the result
    try:
        with open(cache_file, 'w', encoding='utf-8') as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        print(f"[LOG] Cached GET response to {cache_file}")
    except Exception as e:
        print(f"[WARN] Failed to cache response: {e}")
    
    return result  


def retry_network_call(func, *args, max_retries=3, **kwargs):
    """
    Retry wrapper for network calls with exponential backoff
    """
    import time
    
    for attempt in range(max_retries):
        try:
            result = func(*args, **kwargs)
            return result
        except Exception as e:
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt  # Exponential backoff: 1s, 2s, 4s
                print(f"[WARN] Network call failed (attempt {attempt + 1}/{max_retries}): {e}")
                print(f"[LOG] Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print(f"[ERROR] Network call failed after {max_retries} attempts: {e}")
                raise


def cleanup_cache(max_age_days=7):
    """
    Clean up old cache files to prevent disk space issues.
    Removes cache files older than max_age_days.
    """
    workingFilesPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".workingFiles")
    
    for cache_type in ["post_cache", "get_cache"]:
        cache_dir = os.path.join(workingFilesPath, cache_type)
        if not os.path.exists(cache_dir):
            continue
            
        try:
            for filename in os.listdir(cache_dir):
                file_path = os.path.join(cache_dir, filename)
                if os.path.isfile(file_path):
                    file_age_days = (time.time() - os.path.getmtime(file_path)) / 86400
                    if file_age_days > max_age_days:
                        os.remove(file_path)
                        print(f"[LOG] Removed old cache file: {filename} (age: {file_age_days:.1f} days)")
        except Exception as e:
            print(f"[WARN] Error cleaning cache directory {cache_type}: {e}")


def get_cache_stats():
    """
    Get statistics about the cache usage.
    """
    workingFilesPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".workingFiles")
    stats = {"post_cache": 0, "get_cache": 0, "total_size_mb": 0}
    
    for cache_type in ["post_cache", "get_cache"]:
        cache_dir = os.path.join(workingFilesPath, cache_type)
        if os.path.exists(cache_dir):
            try:
                files = os.listdir(cache_dir)
                stats[cache_type] = len(files)
                
                for filename in files:
                    file_path = os.path.join(cache_dir, filename)
                    if os.path.isfile(file_path):
                        stats["total_size_mb"] += os.path.getsize(file_path) / (1024 * 1024)
            except Exception as e:
                print(f"[WARN] Error reading cache directory {cache_type}: {e}")
    
    return stats


# get clumps using the refGen and superPop
def getClumps(refGen, superPop, snpsFromAssociations):
    body = {
        "refGen": refGen,
        "superPop": superPop,
    }
    print(f"Retrieving clumping information: {superPop}")

    try:
        chromToPosMap = {}
        clumps = {}
        for pos in snpsFromAssociations:
            if (len(pos.split(":")) > 1):
                chrom,posit = pos.split(":")
                if (chrom not in chromToPosMap.keys()):
                    chromToPosMap[chrom] = [pos]
                else:
                    chromToPosMap[chrom].append(pos)

        print("Clumps downloaded by chromosome:")
        for chrom in chromToPosMap:
            print("{0}...".format(chrom), end="", flush=True)
            body['positions'] = chromToPosMap[chrom]
            
            # Use retry logic for critical clumping call
            chrom_clumps = retry_network_call(
                postUrlWithBody, 
                "https://prs.byu.edu/ld_clumping_by_pos", 
                body
            )
            clumps = {**chrom_clumps, **clumps}
        print('\n')
    except AssertionError:
        raise SystemExit("ERROR: 504 - Connection to the server timed out")

    return clumps


# get maf using the maf cohort
def getMaf(mafCohort, refGen, snpsFromAssociations):
    # if the cohort is user, return empty, we will use the user maf
    if (mafCohort == 'user'): return {}

    body = {
        "cohort": mafCohort,
        "refGen": refGen
    }
    print("Retrieving maf information")
    
    try:
        chromToPosMap = {}
        maf = {}
        for pos in snpsFromAssociations:
            if (len(pos.split(":")) > 1):
                chrom,posit = pos.split(":")
                if (chrom not in chromToPosMap.keys()):
                    chromToPosMap[chrom] = [posit]
                else:
                    chromToPosMap[chrom].append(posit)

        for chrom in chromToPosMap:
            print("{0}...".format(chrom), end="", flush=True)
            body['chrom'] = chrom
            body['pos'] = chromToPosMap[chrom]
            
            # Use retry logic for MAF call
            chrom_maf = retry_network_call(
                postUrlWithBody,
                "https://prs.byu.edu/get_maf", 
                body
            )
            maf = {**chrom_maf, **maf}
        print('\n')
    except AssertionError:
        raise SystemExit("ERROR: 504 - Connection to the server timed out")

    return maf


def getPercentiles(percentilesCohort, finalStudyList):
    # if the cohort is user, return empty, we will use the user maf
    if (percentilesCohort == 'user'): return {}

    print("Retrieving Percentile information for studies")

    percentiles = {}

    lengthOfList = len(finalStudyList)
    # i and j allow us to cut the list of items sent into chunks of 1000
    i = 0
    j = 1000 if lengthOfList > 1000 else lengthOfList
    runLoop = True
    try:
        while runLoop:
            if j == lengthOfList:
                runLoop = False
            # get the associations based on the studyIDs
            print("{}...".format(j), end = "", flush=True)
            body = {
                "cohort": percentilesCohort,
                "studyIDObjs":finalStudyList[i:j]
            }
            tmpPercentiles = postUrlWithBody("https://prs.byu.edu/get_percentiles", body)
            for key in tmpPercentiles:
                if key not in percentiles:
                    percentiles[key] = tmpPercentiles[key]
            i = j
            j = j + 1000 if lengthOfList > j + 1000 else lengthOfList
        print("Done\n")
    except AssertionError:
        raise SystemExit("ERROR: 504 - Connection to the server timed out")

    return percentiles


# gets associationReturnObj using the given filters
def getSpecificStudySnps(finalStudyList):
    # get the studies matching the parameters
    studySnps = {}
    lengthOfList = len(finalStudyList)
    i = 0
    j = 1000 if lengthOfList > 1000 else lengthOfList
    print("Getting snps to studies map")
    print("Total number of studies: {}". format(lengthOfList))
    runLoop = True
    try:
        while runLoop:
            if j == lengthOfList:
                runLoop = False
            # get the associations based on the studyIDs
            print("{}...".format(j), end = "", flush=True)
            body = {
                "studyIDObjs":finalStudyList[i:j]
            }
            tmpStudySnps = postUrlWithBody("https://prs.byu.edu/snps_to_trait_studyID", body)
            for key in tmpStudySnps:
                if key not in studySnps:
                    studySnps[key] = tmpStudySnps[key]
            i = j
            j = j + 1000 if lengthOfList > j + 1000 else lengthOfList
        print("Done\n")

    except AssertionError:
        raise SystemExit("ERROR: 504 - Connection to the server timed out")
    
    return studySnps

def getCachedEthnicities(max_age_days=30):
    import json

    workingFilesPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".workingFiles")
    ethnicity_file = os.path.join(workingFilesPath, "ethnicities.txt")

    # 1. If file exists and is recent, use cache
    if os.path.exists(ethnicity_file):
        file_mod_time = os.path.getmtime(ethnicity_file)
        file_age_days = (time.time() - file_mod_time) / 86400
        if file_age_days < max_age_days:
            with open(ethnicity_file, 'r', encoding="utf-8") as f:
                try:
                    return json.load(f)
                except Exception:
                    pass  # In case cache is corrupted, fall back to download

    # 2. If not cached or too old, fetch and update cache
    result = getUrlWithParams("https://prs.byu.edu/ethnicities", params={})
    if not os.path.exists(workingFilesPath):
        os.mkdir(workingFilesPath)
    with open(ethnicity_file, 'w', encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    return result


def checkInternetConnection():
    try:
        import socket
        # using an arbitrary connection to check if we can make one
        socket.create_connection(("8.8.8.8", 53))
        return
    except OSError:
        raise SystemExit("ERROR: No internet - Check your connection")


def formatMafCohort(mafCohort):
    mafCohort = mafCohort.replace("-", "_")
    if mafCohort == "adni_cn":
        mafCohort = "adni_controls"

    return mafCohort


def getPopList(popListStr):
    if isinstance(popListStr, list):
        if len(popListStr) == 1 and "|" in popListStr[0]:
            popListStr = popListStr[0].upper()
        else:
            return [pop.upper() for pop in popListStr]

    popList = []
    popListStr = popListStr.upper()
    # split the string on bars if they are present, otherwise add the string to a list of length 1
    if "|" in popListStr:
        popList = popListStr.split("|")
    else:
        popList = [popListStr]
    return popList


def getPreferredPop(popList, superPop):
    popList = getPopList(popList)
    # convert all populations listed in the gwas to lower case
    if len(popList) == 1 and str(popList[0]).lower() == 'na':
        return(superPop)
    else:
        superPopHeirarchy = {
            'EUR': ['EUR', 'AMR', 'SAS', 'EAS', 'AFR'],
            'AMR': ['AMR', 'EUR', 'SAS', 'EAS', 'AFR'],
            'SAS': ['SAS', 'EAS', 'AMR', 'EUR', 'AFR'],
            'EAS': ['EAS', 'SAS', 'AMR', 'EUR', 'AFR'],
            'AFR': ['AFR', 'AMR', 'SAS', 'EUR', 'EAS']
        }
        keys = superPopHeirarchy[superPop]
        for pop in keys:
            if pop in popList:
                # return the first pop from the heirarchy that is in the pop list
                return pop
    
    # if none of the pops from the heirarchy are in the pop list, return the requested super pop
    return superPop


if __name__ == "__main__":
    if argv[1] == "GWAS":
        formatGWASAndRetrieveClumps(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10])
    else:
        retrieveAssociationsAndClumps(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11])
