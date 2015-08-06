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
    "kegg_fp": "/home/ashwini/ash/kegg/kegg",
    "kegg_idx_fp":"",
    "rap_search_fp": "",
    "search_method":"RAPsearch", # RAPsearch or blastx
    "mapping_method":"best_hit", # best_hit or humann
    "evalue_cutoff":0.001,
    "num_threads":4
    }

class Alignment(object):
    def __init__(self, config):
        self.config = config

    def _make_fastq_to_fasta_command(self, filename):
        return [
            "seqtk",
            "seq", "-a",
            filename
            ]
    
    def fastq_to_fasta(self, filename):
        fasta = tempfile.NamedTemporaryFile()
        command = self._make_fastq_to_fasta_command(filename)
        subprocess.check_call(command, stdout=fasta, stderr=subprocess.STDOUT)
        return fasta

    def make_command(self, R, output):
        raise NotImplementedError("Create and override method in child tool")

    def run(self, R, out_dir):
        r_fasta = self.fastq_to_fasta(R)
        output_fp = self.make_output_fp(out_dir, R)
        command = self.make_command(r_fasta.name, output_fp)
        subprocess.check_call(command, stderr=subprocess.STDOUT)
        return output_fp

class Blast(Alignment):
   def __init__(self, config):
       Alignment.__init__(self, config)

   def make_output_fp(self, out_dir, R):
       return os.path.join(out_dir, os.path.basename(os.path.splitext(R)[0]+'.blast'))
            
   def make_command(self, R, output_fp):
       return [
           "blastx", "-outfmt", "6", "-num_threads", str(self.config["num_threads"]),
           "-db", self.config["kegg_fp"],
           "-query", R, "-out", output_fp
           ]

class RapSearch(Alignment):
    def __init__(self, config):
        Alignment.__init__(self, config)

    def run_prerapsearch(self):
        cmd = [
            os.path.join(os.path.dirname(self.config["rap_search_fp"]),'prerapsearch'),
            "-d", self.config["kegg_fp"],
            "-n", os.path.join(os.path.basename(os.path.dirname(self.config['kegg_fp'])), 'keggRAP')
            ]
        subprocess.check_call(cmd, stderr=subprocess.STDOUT)
        #### update the config file of the object?? maybe put this method in main??

    def make_output_fp(self, out_dir, R):
        return os.path.join(out_dir, os.path.basename(os.path.splitext(R)[0]))
                
    def make_command(self, R, output_fp):
        return [
            self.config["rap_search_fp"], "-q", R,
            "-d", self.config["kegg_idx_fp"],
            "-o", output_fp,
            "-z", str(self.config["num_threads"]), "-s", "f"
            ]

class Assignment(object):
    def __init__(self, config):
        self.config = config
        
class Humann(Assignment):
    def __init__(self, config):
        Assignment.__init__(self, config)

    def run(self, alignment):
        ## instead of assuming they are in the correct folder, softlink the alignment file to humann/input
        ## return os.path.join(self.config["humann_fp"], "input", out_dir) ##### softlink this!
        ## currently only works with blast
        os.chdir(self.config["humann_fp"])
        subprocess.check_call("scons", stderr=subprocess.STDOUT)

    def fix_alignment_result(self, alignment):
        pass ##RAPsearch has 4 headerlines that humann doesn't seem to like...

    def parseResults(alignment):
        pass ## to implement later if we would like to pass the results to the summary file
    
class BestHit(Assignment):
    def __init__(self, config):
        Assignment.__init__(self, config)

    def _parseResults(self, alignment):
        idx = alignment.groupby('# Fields: Query').apply(lambda g: g['e-value'].idxmin())
        bestAlignment = alignment.ix[idx, ['# Fields: Query', 'Subject', 'identity', 'e-value']]
        return bestAlignment[bestAlignment['e-value']<self.config["evalue_cutoff"]]
    
    def _map2ko(self, alignment, kegg2ko):
        hitCounts = kegg2ko.merge(alignment, on='Subject')
        hitCounts['ko_abundance'] = hitCounts['Count_hit'] / hitCounts['Count_ko']
        return hitCounts.groupby('KO').apply(lambda g: g['ko_abundance'].sum())
    
    def _getCount(self, df, group_key, count_key):
        return df.groupby(group_key).size().to_frame(name=count_key).reset_index()

    def make_kegg2ko(self): ## incorporate this method into the workflow
        out_fp = os.path.join(os.path.dirname(self.config['kegg_fp']), 'kegg2ko') 
        kegg_db = SeqIO.parse(open(self.config['kegg_fp']), 'fasta')
        with open(out_fp, 'w') as f_out:
            for ko in kegg_db:
                id_split = re.split(';', ko.description)
                id_split = [_.strip() for _ in id_split]
                [f_out.write('\t'.join([ko.id, _])+'\n') for _ in id_split if _.startswith('K0')]
        ##modify the config and update kegg_to_ko_fp

    def _load_kegg2ko(self):
        kegg2ko = pandas.read_csv(self.config['kegg_to_ko_fp'], sep='\t', header=None, names=['Subject', 'KO'])
        return kegg2ko.merge(self._getCount(kegg2ko, 'Subject', 'Count_ko'), on='Subject')
        
    def _assign_ko(self, alignment_fp):
        kegg2ko = self._load_kegg2ko()
        
        # read the alignent results (RAPsearch for now) and count
        ## currently only works with Rapsearch
        #alignment = pandas.read_csv(blastResults_fp, sep='\t', header=None, names=rapResults.columns) # for reading blast. kind of..
        alignment = pandas.read_csv(alignment_fp, sep='\t', skiprows=4, header=0)
        bestAlignment = self._parseResults(alignment)
        numHits = self._getCount(bestAlignment, 'Subject', 'Count_hit')
        # merge and count
        return self._map2ko(numHits, kegg2ko)

    def _assign_pathway(self, ko_assign):
        raise NotImplementedError("Create and override method in child tool")

    def run(self, alignment_fp):
        ko_count = self._assign_ko(alignment_fp)
        print(ko_count)

        #path_count = self._assign_pathway(ko_count) # not yet implemented
        
        return "Assigner run" ## return results as a dictionary so they can be written to json file 
        
        
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

    if config['search_method'].lower() == 'blastx':
        searchApp = Blast(config)
    else:
        searchApp = RapSearch(config)
    alignment_R1_fp = searchApp.run(fwd_fp, args.output_dir)
    alignment_R2_fp = searchApp.run(rev_fp, args.output_dir)

    #alignment_R1_fp = '/home/tanesc/data/rapTemp_'
    #config['kegg_to_ko_fp'] = '/home/tanesc/data/kegg2ko_'

    if config['mapping_method'].lower() == 'humann':
        assignerApp = Humann(config)
    else:
        assignerApp = BestHit(config)
        alignment_R1_fp += '.m8'
        alignment_R2_fp += '.m8'

    summary_R1 = assignerApp.run(alignment_R1_fp)
    summary_R2 = assignerApp.run(alignment_R2_fp)
    
    #save_summary(args.summary_file, config, data)

def save_summary(f, config, data):
    result = {
        "program": "PathFinder",
        "version": __version__,
        "config": config,
        "data":data
        }
    json.dump(result, f)    
