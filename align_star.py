import pipelineUtil
import qc
import argparse
import glob
import os
import re

def decompress(filename, workdir):
    """ Unpack fastq files """

    if filename.endswith(".tar"):
        cmd = ['tar', 'xvf', filename, '-C', workdir]
    elif filename.endswith(".gz"):
        cmd = ['tar', 'xzvf', filename, '-C', workdir]
    elif filename.endswith(".bz"):
        cmd = ['tar', 'xvjf', filename, '-C', workdir]
    else:
        raise Exception('Unknown input file extension for file %s' % filename)
    pipelineUtil.log_function_time("tar", filename, cmd)


def scan_workdir_helper(dirname, extension):
    fastq_files = glob.glob(os.path.join(dirname, "*_[12].%s"%(extension)))
    all_read_groups = list()
    read_group_set = dict()
    if not (fastq_files == []):
        for filename in fastq_files:
            rg_id = re.sub(r'_[12].%s$'%(extension), '', filename)
            read_group_set[rg_id] = read_group_set.get(rg_id, 0) + 1
        if not all(i == 2 for i in read_group_set.values()):
            raise Exception("Missing Pair")
        print read_group_set
        for rg_id in read_group_set.keys():
            reads_1 = "%s_1.%s" %(rg_id, extension)
            reads_2 = "%s_2.%s" %(rg_id, extension)
            read_pair = (os.path.basename(rg_id), reads_1, reads_2)
            all_read_groups.append(read_pair)

    return all_read_groups

def scan_workdir(dirname):
    """ Select the unpacked fastq files """

    print dirname
    fastq_files = scan_workdir_helper(dirname, "fastq")
    if fastq_files == []:
        fastq_files = scan_workdir_helper(dirname, "fastq.gz")
        if fastq_files == []:
            fastq_files = scan_workdir_helper(dirname, "fastq.bz")
    return fastq_files


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extension of RNA-seq alignment with QC.")
    required = parser.add_argument_group("Required input parameters")
    required.add_argument("--genomeDir", default=None, help="Directory containing the reference genome index", required=True)
    required.add_argument("--tarFileIn", default=None, help="Input file containing the sequence information", required=True)
    optional = parser.add_argument_group("optional input parameters")
    optional.add_argument('--fastqc_path', help='path to fastqc', default='/home/ubuntu/bin/FastQC/fastqc')
    optional.add_argument("--out", default="out.bam", help="Name of the output BAM file")
    optional.add_argument("--workDir", default="./", help="Work directory")
    optional.add_argument("--metaDataTab", default=None, help="File containing metadata for the alignment header")
    optional.add_argument("--analysisID", default=None, help="Analysis ID to be considered in the metadata file")
    optional.add_argument("--keepJunctions", default=False, action='store_true', help="keeps the junction file as {--out}.junctions")
    optional.add_argument("--useTMP", default=None, help="environment variable that is used as prefix for temprary data")
    #optional.add_argument("-h", "--help", action='store_true', help="show this help message and exit")
    star = parser.add_argument_group("STAR input parameters")
    star.add_argument("--runThreadN", type=int, default=4, help="Number of threads")
    star.add_argument("--outFilterMultimapScoreRange", type=int, default=1, help="outFilterMultimapScoreRange")
    star.add_argument("--outFilterMultimapNmax", type=int, default=20, help="outFilterMultimapNmax")
    star.add_argument("--outFilterMismatchNmax", type=int, default=10, help="outFilterMismatchNmax")
    star.add_argument("--alignIntronMax", type=int, default=500000, help="alignIntronMax")
    star.add_argument("--alignMatesGapMax", type=int, default=1000000, help="alignMatesGapMax")
    star.add_argument("--sjdbScore", type=int, default=2, help="sjdbScore")
    star.add_argument("--limitBAMsortRAM", type=int, default=0, help="limitBAMsortRAM")
    star.add_argument("--alignSJDBoverhangMin", type=int, default=1, help="alignSJDBoverhangMin")
    star.add_argument("--genomeLoad", default="NoSharedMemory", help="genomeLoad")
    star.add_argument("--genomeFastaFiles", default=None, help="genome sequence in fasta format to rebuild index")
    star.add_argument("--outFilterMatchNminOverLread", type=float, default=0.33, help="outFilterMatchNminOverLread")
    star.add_argument("--outFilterScoreMinOverLread", type=float, default=0.33, help="outFilterScoreMinOverLread")
    star.add_argument("--twopass1readsN", type=int, default=-1, help="twopass1readsN (-1 means all reads are used for remapping)")
    star.add_argument("--sjdbOverhang", type=int, default=100, help="sjdbOverhang (only necessary for two-pass mode)")
    star.add_argument("--outSAMstrandField", default="intronMotif", help="outSAMstrandField")
    star.add_argument("--outSAMattributes", default=["NH", "HI", "NM", "MD", "AS", "XS"], help="outSAMattributes")
    star.add_argument("--outSAMunmapped", default="Within", help="outSAMunmapped")
    star.add_argument("--outSAMtype", default=["BAM", "SortedByCoordinate"], help="outSAMtype")
    star.add_argument("--outSAMheaderHD", default=["@HD", "VN:1.4"], help="outSAMheaderHD")
    star.add_argument("--outSAMattrRGline", default=None, help="RG attribute line submitted to outSAMattrRGline")
    star.add_argument("--outSAMattrRGfile", default=None, help="File containing the RG attribute line submitted to outSAMattrRGline")

    args = parser.parse_args()
    decompress(args.tarFileIn, args.workDir)
    read_group_pairs = scan_workdir(args.workDir)
    for (rg_id, reads_1, reads_2) in read_group_pairs:
        rg_id_dir = os.path.join(args.workDir, rg_id)
        if not os.path.isdir(rg_id_dir):
            os.mkdir(rg_id_dir)
        qc.fastqc(args.fastqc_path, reads_1, reads_2, rg_id_dir, rg_id, None)


