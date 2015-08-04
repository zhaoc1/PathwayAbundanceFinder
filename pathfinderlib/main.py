import subprocess
import os
import tempfile
import argparse
import json
import re
import pandas
from Bio import SeqIO

from pathfinderlib.version import __version__


default_config ={
    "humann_fp": "/home/ashwini/ash/other_softwares/humann/",
    "kegg_fp": "/home/ashwini/aseh/kegg/kegg",
    "kegg_idx_fp":"",
    "rap_search_fp": "",
    "search_method":"RapSearch", # RapSearch or blast
    "mapping_method":"inhouse" # inhouse or humann
    }

class Allignment(object):
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

    def make_command(self, R, output):
        raise NotImplementedError("Create and override method in child tool")

    def make_output_fp(self, out_dir, R):
        return os.path.join(out_dir, os.path.basename(R))
    
    def run(self, R, out_dir):
        r_fasta = self.fastq_to_fasta(R)
        output_fp = self.make_output_fp(out_dir, R)
        command = self.make_command(r_fasta, output_fp)
        print(command)
        subprocess.check_call(command, stderr=subprocess.STDOUT)
        os.remove(r_fasta)
        return output_fp

class Blast(Allignment):
   def __init__(self, config):
       Allignment.__init__(self, config)

   def make_command(self, R, output_fp):
       return [
           "blastx", "-outfmt", "6", "-num_threads", "6",
           "-db", self.config["kegg_fp"],
           "-query", R, "-out", output_fp
           ]

class RapSearch(Allignment):
    def __init__(self, config):
        Allignment.__init__(self, config)

    def run_prerapsearch(self):
        pass # check if there is rap_search_fp specifid in config file, if not create kegg index file

    def make_command(self, R, output_fp):
        return [
            self.config["rap_search_fp"], "-q", R,
            "-d", self.config["kegg_idx_fp"],
            "-o", output_fp, "-z", "4"
            ]

class Assignment(object):
    def __init__(self, config):
        self.config = config
        
class Humann(Assignment):
    def __init__(self, config):
        Assignment.__init__(self, config)

    def run(self, allignment):
        ## instead of assuming they are in the correct folder, softlink the allignment file to humann/input
        # return os.path.join(self.config["humann_fp"], "input", out_dir) ##### softlink this!
        ## currently only works with blast
        os.chdir(self.config["humann_fp"])
        subprocess.check_call("scons", stderr=subprocess.STDOUT)

    def fix_allignment_result(self, allignment):
        pass ##RAPsearch has 4 headerlines that humann doesn't seem to like...

    def parseResults(allignment):
        pass ## to implement later if we would like to pass the results to the summary file
    
class BestHit(Assignment):
    def __init__(self, config):
        Assignment.__init__(self, config)

    def parseResults(self, allignment, evalue=0.001):
        idx = allignment.groupby('# Fields: Query').apply(lambda g: g['e-value'].idxmin())
        bestAllignment = allignment.ix[idx, ['# Fields: Query', 'Subject', 'identity', 'e-value']]
        bestAllignment = bestAllignment[bestAllignment['e-value']<evalue] # change it to the filter command here
        return bestAllignment #.set_index('# Fields: Query')
    
    def map2ko(self, allignment, kegg2ko):
        hitCounts = kegg2ko.merge(allignment, on='Subject')
        hitCounts['ko_abundance'] = hitCounts['Count_hit'] / hitCounts['Count_ko']
        return hitCounts.groupby('KO').apply(lambda g: g['ko_abundance'].sum())
    
    def getCount(self, df, group_key, count_key):
        return df.groupby(group_key).size().to_frame(name=count_key).reset_index()

    def make_kegg2ko(self):
        out_fp = 'kegg2ko' ### 
        kegg_db = SeqIO.parse(open(self.config['kegg_fp']), 'fasta')
        with open(out_fp, 'w') as f_out:
            for ko in kegg_db:
                id_split = re.split(';', ko.description)
                id_split = [_.strip() for _ in id_split]
                [f_out.write('\t'.join([ko.id, a])+'\n') for a in id_split if a.startswith('K0')]
        ##modify the config and update kegg_to_ko_fp

    def load_kegg2ko(self):
        kegg2ko = pandas.read_csv(self.config['kegg_to_ko_fp'], sep='\t', header=None, names=['Subject', 'KO'])
        return kegg2ko.merge(self.getCount(kegg2ko, 'Subject', 'Count_ko'), on='Subject')
        
    def assign_ko(self, allignment_fp):
        kegg2ko = self.load_kegg2ko()
        
        # read the alignent results (RAPsearch for now) and count
        ## currently only works with Rapsearch
        #pandas.read_csv(blastResults_fp, sep='\t', header=None, names=rapResults.columns) # for reading blast
        allignment = pandas.read_csv(allignment_fp, sep='\t', skiprows=4, header=0)
        bestAllignment = self.parseResults(allignment)
        numHits = self.getCount(bestAllignment, 'Subject', 'Count_hit')

        # merge and count
        return self.map2ko(numHits, kegg2ko)

    def run(self, allignment_fp):
        ko_count = self.assign_ko(allignment_fp)
        print(ko_count)

        #path_count = assign_path(ko_count) # not yet implemented
        
        return "Assigner run" ##maybe do some extra summarizing here.. 
        
        
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
        "--output-dir", required=True,
        help="output directory")
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

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    #searchApp = RapSearch(config)
    #allignment_R1_fp = searchApp.run(fwd_fp, args.output_dir)
    #allignment_R2_fp = searchApp.run(rev_fp, args.output_dir)

    allignment_R1_fp = '/home/tanesc/data/rapTemp_'
    config['kegg_to_ko_fp'] = '/home/tanesc/data/kegg2ko_'
    assignerApp = BestHit(config)
    summary_R1 = assignerApp.run(allignment_R1_fp)
    #summary_R2 = assignerApp.run(allignment_R2_fp)
    
    #app = Humann(config)
    #app.run(fwd_fp, rev_fp, args.output_dir)
    #save_summary(args.summary_file, config)

def save_summary(f, config):
    result = {
        "program": "PathFinder",
        "version": __version__,
        "config": config
        }
    json.dump(result, f)    
