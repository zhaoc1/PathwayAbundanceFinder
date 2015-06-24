import subprocess
import os
import tempfile
import argparse

#path to humann
path = "/home/ashwini/ash/other_softwares/humann/"

def run_command(command, error_message):
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print error_message

def check_file_exists_or_die(file_name):
    if not os.path.isfile(file_name):
        print "ERROR: file " + file_name + " does not exist: check the file name and path."
        sys.exit(1)
        
def get_args(argv):
    parser = argparse.ArgumentParser(description="Runs Humann.")
    parser.add_argument("--R1", required=True,
                        type=argparse.FileType("r"),
                        help="R1.fastq")
    parser.add_argument("--R2", required=True,
                        type=argparse.FileType("r"),
                        help="R2.fastq")
    parser.add_argument("--output", required=True,
                        type=str,
                        help="output file")
    args = parser.parse_args(argv)
    check_file_exists_or_die(args.R1.name)
    check_file_exists_or_die(args.R2.name)
    args.R1.close()
    args.R2.close()
    return(args)
                                
                            
def fastq_to_fasta(filename):
    fasta = tempfile.NamedTemporaryFile(delete=False)
    command = ("seqtk seq -a  " + filename + " > " + fasta.name)
    run_command(command, "cannot run seqtk for blat.")
    return fasta.name

def create_input_files_for_humann(R1, R2, output):
    r1_fasta = fastq_to_fasta(R1)
    r2_fasta = fastq_to_fasta(R2)
    command_r1 = ("blastx -outfmt 6 -db ~/ash/kegg/kegg " +
                  "-query " + r1_fasta +
                  " -out " + path + "input/" + output + "_R1.txt")
    run_command(command_r1, "Cannot run blastx on " + R1 + " . Check input files.")

    command_r2 = ("blastx -outfmt 6 -db ~/ash/kegg/kegg " +
                  "-query " + r2_fasta +
                  " -out " + path + "input/" + output + "_R2.txt")
    run_command(command_r2, "Cannot run blastx on " + R2 + " . Check input files.")
    os.remove(r1_fasta)
    os.remove(r2_fasta)
    
def run_humann():
    os.chdir(path)
    command = ("scons")
    run_command(command, "Cannot run scons. check SConstruct file.")

def main(argv=None):
    args = get_args(argv)
    create_input_files_for_humann(args.R1.name, args.R2.name, args.output)
    run_humann()
    

if __name__=="__main__":
    main()
    
