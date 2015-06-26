import subprocess
import os
import tempfile
import argparse
import json

from pathfinderlib.version import __version__
#path to humann

default_config ={
    "humann_fp": "/home/ashwini/ash/other_softwares/humann/",
    "kegg_fp": "/home/ashwini/ash/kegg/kegg"
    }
        
                                
class Humann(object):
    def __init__(self, config):
        self.config = config

    def make_fastq_to_fasta_command(self, filename):
        return [
            "seqtk",
            "seq", "-a",
            filename
            ]

    def fastq_to_fasta(self, filename):
        fasta = tempfile.NamedTemporaryFile(delete=False)
        command = self.make_fastq_to_fasta_command(filename)
        subprocess.check_call(command, stdout=fasta, stderr=subprocess.STDOUT)
        return fasta.name

    def make_blastx_command(self, R, rex, output):
        out_dir = os.path.join(self.config["humann_fp"], "input", output + rex)
        return [
            "blastx", "-outfmt", "6",
            "-db", self.config["kegg_fp"],
            "-query", R, "-out", out_dir
            ]

    def run_blastx(self, R, rex, output):
        r_fasta = self.fastq_to_fasta(R)
        command = self.make_blastx_command(r_fasta, rex, out_dir)
        subprocess.check_call(command, stderr=subprocess.STDOUT)
        os.remove(r_fasta)
        
    def run(self, R1, R2, output):
        self.run_blastx(R1,"_R1.txt", output)
        self.run_blastx(R2, "_R2.txt", output)
        os.chdir(self.config["humann_fp"])
        subprocess.check_call("scons", stderr=subprocess.STDOUT)
        
def main(argv=None):
    parser = argparse.ArgumentParser(description="Runs Humann.")
    parser.add_argument(
        "--forward-reads", required=True,
        type=argparse.FileType("r"),
        help="R1.fastq")
    parser.add_argument(
        "--reverse-reads", required=True,
        type=argparse.FileType("r"),
        help="R2.fastq")
    parser.add_argument(
        "--summary-file", required=True,
        type=argparse.FileType("w"),
        help="Summary file")
    parser.add_argument(
        "--output-file", required=True,
        help="output file")
    parser.add_argument(
        "--config-file",
        type=argparse.FileType("r"),
        help="JSON configuration file")
    args = parser.parse_args(argv)

    config = default_config.copy()
    if args.config_file:
        user_config = json.load(args.config_file)
        config.update(user_config)

    fwd_fp = args.forward_reads.name
    rev_fp = args.reverse_reads.name
    args.forward_reads.close()
    args.reverse_reads.close()

    app = Humann(config)
    app.run(fwd_fp, rev_fp, args.output_file)
    save_summary(args.summary_file, config)

def save_summary(f, config):
    result = {
        "program": "PathFinder",
        "version": __version__,
        "config": config
        }
    json.dump(result, f)    
