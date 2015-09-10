import pipelineUtil
import os
import subprocess

def validate_bam_file(picard_path, bam_file, uuid, outdir, logger=None):
   """ Validate resulting post-alignment BAM file """

   if os.path.isfile(picard_path) and os.path.isfile(bam_file):
        outfile = os.path.join(outdir, "%s.validate" %uuid)
        tmp_dir = os.path.join(outdir, 'tmp')
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        cmd = ['java', '-jar', picard_path, "ValidateSamFile", "I=%s" %bam_file,
               "O=%s" %outfile, "VALIDATION_STRINGENCY=LENIENT",
               "TMP_DIR=%s" %tmp_dir]
        exit_code = pipelineUtil.log_function_time("ValidateSAM", uuid, cmd, logger)
        if exit_code == 0:
            assert(os.path.isfile(outfile))
   else:
       raise Exception("Invalid path to picard or BAM")

   if not exit_code == 0:
       if not logger == None:
            logger.error("Picard ValidateSamFile returned non-zero exit code %s" %exit_code)

   return exit_code

def collect_rna_seq_metrics(picard_path, bam_file, uuid, outdir, ref_flat, logger=None):
    """ Collect RNA-seq metrics using Picard """

    if os.path.isfile(picard_path) and os.path.isfile(bam_file):
        tmp_dir = os.path.join(outdir, 'tmp')
        outfile = os.path.join(outdir, "%s.rna_seq_metrics.txt" %uuid)
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        cmd = ['java', '-jar', picard_path, "CollectRnaSeqMetrics", "METRIC_ACCUMULATION_LEVEL=READ_GROUP",
                "I=%s" %bam_file, "O=%s" %outfile, "STRAND=NONE",
                "REF_FLAT=%s" %ref_flat, "VALIDATION_STRINGENCY=LENIENT", "TMP_DIR=%s" %tmp_dir]
        exit_code = pipelineUtil.log_function_time("RNAseq_metrics", uuid, cmd, logger)
        if exit_code == 0:
            assert(os.path.isfile(outfile))
    else:
        raise Exception("Invalid path to picard or bam")

    if not exit_code == 0:
       if not logger == None:
            logger.error("Picard CollectRnaSeqMetrics returned non-zero exit code %s" %exit_code)
    return exit_code

def bam_index(bam_file, uuid, logger=None):
    """ Index the resultant post alignment BAM file """

    if os.path.isfile(bam_file):
        cmd = ['samtools', 'index', '-b', bam_file]
        exit_code = pipelineUtil.log_function_time("BamIndex", uuid, cmd, logger)
        if exit_code == 0:
            assert(os.path.isfile('%s.bai' %bam_file))
    else:
        raise Exception("Cannot file bam file  %s" %bam_file)

    if not exit_code == 0:
       if not logger == None:
            logger.error("Samtools Index returned non-zero exit code %s" %exit_code)
    return exit_code

def fix_mate_information(picard_path, bam_file, uuid, outdir, logger=None):
    """ Fix the mate information for BAM files """

    if os.path.isfile(picard_path) and os.path.isfile(bam_file):
        tmp_dir = os.path.join(outdir, 'tmp')
        outfile =  "%s.fix" %bam_file
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        cmd = ['java', '-jar', picard_path, 'FixMateInformation', 'I=%s' %bam_file, 'O=%s' %outfile,
                'VALIDATION_STRINGENCY=LENIENT', 'TMP_DIR=%s' %tmp_dir]
        exit_code = pipelineUtil.log_function_time('FixMateInformation', uuid, cmd, logger)
        if exit_code == 0:
            assert(os.path.isfile(outfile))
    else:
        raise Exception("Invalid path to picard %s or BAM %s" %(picard_path, bam_file))

    if not exit_code == 0:
        if not logger == None:
            logger.error("Picard FixMateInformation returned non-zero exit code %s" %exit_code)
    return exit_code, outfile

def reorder_bam(picard_path, bam_file, uuid, outdir, ref_genome, logger=None):
    """ Reorder the BAM file according to the reference genome """

    if os.path.isfile(bam_file) and os.path.isfile(picard_path) and os.path.isfile(ref_genome):
        outbam = os.path.join(outdir, '%s.reorder.bam' %uuid)
        tmp_dir = os.path.join(outdir, 'tmp')
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        cmd = ['java', '-jar', picard_path, 'ReorderSam', 'I=%s' %bam_file, 'O=%s' %outbam, 'R=%s' %ref_genome,
                'VALIDATION_STRINGENCY=LENIENT', 'TMP_DIR=%s' %tmp_dir]
        exit_code = pipelineUtil.log_function_time("picard_reorder_sam", uuid, cmd, logger)
        if exit_code == 0:
            assert(os.path.isfile(outbam))
    else:
        raise Exception("Cannot find one of bam %s, picard path %s or reference genome %s" %(bam_file, picard_path, ref_genome))

    if not exit_code == 0:
        if not logger == None:
            logger.error("Picard reorderBAM returned non-zero exit code %s" %exit_code)
        else:
            raise Exception("Picard reorderBAM returned non-zero exit code %s" %exit_code)
    return outbam

def rna_seq_qc(rna_seq_qc_path, bam_file, uuid, outdir, ref_genome, gtf, logger=None):
    """ Perform RNA-seqQC on post alignment BAM file """

    if os.path.isfile(bam_file) and os.path.isfile(rna_seq_qc_path) and os.path.isfile(gtf):
        cmd = ['java', '-jar', rna_seq_qc_path, '-o', outdir, '-r', ref_genome, '-s',
                '%s|%s|%s' %(uuid, bam_file, uuid), '-t', gtf]
        exit_code = pipelineUtil.log_function_time('RNAseq_qc', uuid, cmd, logger)
    else:
        raise Exception("Cannot find one of  rnaseq-qc %s, bam %s or gtf %s" %(rna_seq_qc_path, bam_file, gtf))

    if not exit_code == 0:
       if not logger == None:
            logger.error("Broad's RNA-Seq-QC returned non-zero exit code %s" %exit_code)
    return exit_code

def add_or_replace_read_group(picard_path, bam_file,  outdir, uuid, rg_id, rg_lb="Unknown", rg_pl="Unknown", rg_pu="Unknown",rg_sm="Unknown", logger=None):
    """ Replace the @RG tag in the reads and header """

    outbam = '%s.addRG.bam' %os.path.join(outdir, uuid)
    if os.path.isfile(bam_file) and os.path.isfile(picard_path):
        tmp_dir = os.path.join(outdir, 'tmp')
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        cmd = ['java', '-jar', picard_path, 'AddOrReplaceReadGroups', 'I=%s' %bam_file, 'O=%s' %outbam,
                'RGID=%s'%rg_id, 'RGLB=%s' %rg_lb, 'RGPL=%s' %rg_pl, 'RGPU=%s' %rg_pu, 'RGSM=%s' %rg_sm,
                'VALIDATION_STRINGENCY=LENIENT','TMP_DIR=%s' %tmp_dir]
        exit_code = pipelineUtil.log_function_time('AddOrReplaceReadGroups', uuid, cmd, logger)
    else:
        raise Exception("Cannot find bam file %s or path to picard %s" %(bam_file, picard_path))

    if not exit_code == 0:
       if not logger == None:
            logger.error("Picard AddOrReplaceReadGroups returned non-zero exit code %s" %exit_code)

    if os.path.isfile(outbam):
        print "returning file now %s" %outbam
        return outbam
    else:
        raise Exception('Could not add or replace read groups. Check log file for errors')

def reheader_bam_file(bam_file, sample, library, workDir, logger):
    "Reheader BAM file with SM and LB"

    if os.path.isfile(bam_file) and os.path.isdir(workDir):
        cmd = ['samtools', 'view', '-H', bam_file]
        header = open(os.path.join(workDir, 'header.sam'), "w")
        subprocess.call(cmd, stdout=header)
        header.close()

        header = open(os.path.join(workDir, 'header.sam'), "r")
        new_header = open(os.path.join(workDir, 'header_new.sam'), "w")
        for line in header:
            if not line.startswith("@RG"):
                new_header.write(line)
            else:
                line = line.split("\t")
                RG = line[0] + "\t" + line[1].rstrip() + "\tSM:%s" %sample + "\tLB:%s\n" %library
                new_header.write(RG)
        header.close()
        new_header.close()
        os.remove(os.path.join(workDir, 'header.sam'))
        new_bam = '%s.rehead' %bam_file
        new_bam_file = open(new_bam, "w")
        cmd = ['samtools', 'reheader', os.path.join(workDir, 'header_new.sam'), bam_file]
        subprocess.call(cmd, stdout=new_bam_file)
        new_bam_file.close()
        os.remove(os.path.join(workDir, 'header_new.sam'))
        os.rename(new_bam, bam_file)
