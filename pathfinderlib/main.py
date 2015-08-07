import subprocess
import os
import tempfile
import argparse
import inspect
import json
import re
import pandas
from Bio import SeqIO

from pathfinderlib.version import __version__


default_config ={
    "humann_fp": "/home/ashwini/ash/other_softwares/humann/",
    "kegg_fp": "/home/ashwini/ash/kegg/kegg",
    "kegg_idx_fp":"",
    "kegg_to_ko_fp":"",
    "rap_search_fp": "",
    "search_method":"rapsearch", # rapsearch or blastx
    "mapping_method":"best_hit", # best_hit or humann
    "evalue_cutoff":0.001,
    "num_threads":4
    }

def returnObject(tool_cls, config):
    # Proceed stepwise here to improve quality of error messages.
    tool_args = []
    for argname in tool_cls.get_argnames():
        arg = config[argname]
        tool_args.append(arg)
    return tool_cls(*tool_args)

def Alignment(config):
    tool_cls = search_methods_available[config["search_method"]]
    return returnObject(tool_cls, config)

class _Alignment(object):
    def __init__(self, search_method, kegg_fp, num_threads):
        self.search_method = search_method
        self.kegg_fp = kegg_fp
        self.num_threads = num_threads

    @classmethod
    def get_argnames(cls):
        return inspect.getargspec(cls.__init__)[0][1:]

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

class Blast(_Alignment):
   def __init__(self, search_method, kegg_fp, num_threads):
       super(Blast, self).__init__(search_method, kegg_fp, evalue_cutoff, num_threads)

   def make_output_fp(self, out_dir, R):
       return os.path.join(out_dir, os.path.basename(os.path.splitext(R)[0]+'.blast'))
            
   def make_command(self, R, output_fp):
       return [
           "blastx", "-outfmt", "6", "-num_threads", str(self.num_threads),
           "-db", self.kegg_fp,
           "-query", R, "-out", output_fp
           ]

   def make_index(self):
       raise IOError("Protein fasta file  can't be found. Please check kegg_fp in the config file")

   def index_exists(self):
       return os.path.exists(self.kegg_fp)

class RapSearch(_Alignment):
    def __init__(self, search_method, kegg_fp, num_threads, kegg_idx_fp, rap_search_fp):
        super(RapSearch, self).__init__(search_method, kegg_fp, num_threads)
        self.kegg_idx_fp = kegg_idx_fp
        self.rap_search_fp = rap_search_fp

    def make_output_fp(self, out_dir, R):
        return os.path.join(out_dir, os.path.basename(os.path.splitext(R)[0]))
                
    def make_command(self, R, output_fp):
        return [
            self.rap_search_fp, "-q", R,
            "-d", self.kegg_idx_fp,
            "-o", output_fp,
            "-z", str(self.num_threads), "-s", "f"
            ]

    def make_index(self):
        cmd = [
            os.path.join(os.path.dirname(self.rap_search_fp),'prerapsearch'),
            "-d", self.kegg_fp,
            "-n", self.kegg_index_fp
            ]
        subprocess.check_call(cmd, stderr=subprocess.STDOUT)

    def index_exists(self):
        return os.path.exists(self.kegg_index_fp)


def Assignment(config):
    tool_cls = mapping_methods_available[config["mapping_method"]]
    return returnObject(tool_cls, config)

class _Assignment(object):
    def __init__(self, mapping_method, search_method, evalue_cutoff):
        self.mapping_method = mapping_method
        self.search_method = search_method
        self.evalue_cutoff = evalue_cutoff

    @classmethod
    def get_argnames(cls):
        return inspect.getargspec(cls.__init__)[0][1:]

class Humann(_Assignment):
    def __init__(self, mapping_method, search_method, evalue_cutoff, humann_fp):
        super(Humann, self).__init__(mapping_method, search_method, evalue_cutoff)
        self.humann_fp = humann_fp

    def run(self, alignment):
        ## instead of assuming they are in the correct folder, softlink the alignment file to humann/input
        ## return os.path.join(self.config["humann_fp"], "input", out_dir) ##### softlink this!
        ## currently only works with blast
        os.chdir(self.humann_fp)
        subprocess.check_call("scons", stderr=subprocess.STDOUT)

    def fix_alignment_result(self, alignment):
        "Manipulates the alignment file to make it readable by Humann"
        pass ##RAPsearch has 4 headerlines that humann doesn't seem to like...

    def parseResults(alignment):
        pass ## to implement later if we would like to pass the results to the summary file

    def make_index(self):
        pass

    def index_exists(self):
        return True
        
class BestHit(_Assignment):
    def __init__(self, mapping_method, search_method, evalue_cutoff, kegg_fp, kegg_to_ko_fp):
        super(BestHit, self).__init__(mapping_method, search_method, evalue_cutoff)
        self.kegg_fp = kegg_fp
        self.kegg_to_ko_fp = kegg_to_ko_fp
        
    def _getCount(self, df, group_key, count_key):
        "Counts how many times each value is repeated in a group_key column."
        return df.groupby(group_key).size().to_frame(name=count_key).reset_index()

    def _sumAbundance(self, df, group_key, sum_by_key, count_key):
        return df.groupby(group_key).apply(lambda g: g[sum_by_key].sum()).to_frame(name=count_key).reset_index()

    def _load_kegg2ko(self):
        "Loads the index file and counts how many KOs are assigned to each protein."
        kegg2ko = pandas.read_csv(self.kegg_to_ko_fp, sep='\t', header=None, names=['Subject', 'KO'])
        return kegg2ko.merge(self._getCount(kegg2ko, 'Subject', 'Count_ko'), on='Subject')

    def _parseResults(self, alignment_fp):
        "Parses an alignment result file. Returns only the columns of interest"
        colNames = ['# Fields: Query', 'Subject', 'identity', 'e-value']
        if self.search_method.lower() == "blast":
            return pandas.read_csv(alignment_fp, sep='\t', header=None, names=colNames, usecols=[0,1,2,10])
        else: # for rapsearch
            return pandas.read_csv(alignment_fp, sep='\t', skiprows=4, usecols=colNames)
    
    def _getBestHit(self, alignment):
        "Finds the best match that passes an e-value threshold."
        idx = alignment.groupby('# Fields: Query').apply(lambda g: g['e-value'].idxmin())
        bestAlignment = alignment.ix[idx]
        return bestAlignment[bestAlignment['e-value']<self.evalue_cutoff]
    
    def _map2ko(self, alignment, kegg2ko):
        """Merges the alignment results with the kegg index file,
        normalizes the abundances for proteins with multiple KO assignments,
        finds the total abundances for each ko.
        """
        hitCounts = kegg2ko.merge(alignment, on='Subject')
        hitCounts['norm_abundance'] = hitCounts['Count_hit'] / hitCounts['Count_ko']
        return self._sumAbundance(hitCounts, 'KO', 'norm_abundance', 'ko_abundance')
    
    def _assign_ko(self, alignment_fp):
        "Manipulates the alignment file to calculate KO abundances in a fastq file"
        kegg2ko = self._load_kegg2ko()
        print(kegg2ko['KO'].nunique())
        
        alignment = self._parseResults(alignment_fp)
        #count # mapped to db from alignment
        print(alignment['# Fields: Query'].nunique())
        
        bestAlignment = self._getBestHit(alignment)
        # count sequences with < evalue hits with bestAlignment
        print(bestAlignment['# Fields: Query'].nunique())

        numHits = self._getCount(bestAlignment, 'Subject', 'Count_hit')
        # merge and count
        ko_count = self._map2ko(numHits, kegg2ko)
        print(ko_count['KO'].nunique())
        return ko_count

    def _assign_pathway(self, ko_assign):
        "Manipulates the KO abundance table to calculate pathway abundances in afastq file"
        raise NotImplementedError("Create and override method in child tool")

    def run(self, alignment_fp, out_dir):
        ko_count = self._assign_ko(alignment_fp)
        print(ko_count)
        ko_count.to_csv(os.path.join(out_dir, os.path.basename(os.path.splitext(alignment_fp)[0]+'.ko')),
                        sep='\t', index=False)

        #path_count = self._assign_pathway(ko_count) # not yet implemented
        
        return "Assigner run" ## return results as a dictionary so they can be written to json file 

    def make_index(self):
        kegg_db = SeqIO.parse(open(self.kegg_fp), 'fasta')
        with open(self.kegg_to_ko_fp, 'w') as f_out:
            for ko in kegg_db:
                id_split = re.split(';', ko.description)
                id_split = [_.strip() for _ in id_split]
                [f_out.write('\t'.join([ko.id, _])+'\n') for _ in id_split if _.startswith('K0')]

    def index_exists(self):
        return os.path.exists(self.kegg_to_ko_fp)
        
def main(argv=None):
    parser = argparse.ArgumentParser(description="Runs functional assignment.")
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

    searchApp = Alignment(config)
    #alignment_R1_fp = searchApp.run(fwd_fp, args.output_dir)
    #alignment_R2_fp = searchApp.run(rev_fp, args.output_dir)

    ##### find a graceful way of doing this
    #alignment_R1_fp += '.m8'
    #alignment_R2_fp += '.m8'

    config['kegg_to_ko_fp'] = '/home/tanesc/data/kegg2ko_'
    alignment_R1_fp = '/home/tanesc/data/rapTemp_'

    assignerApp = Assignment(config)
    summary_R1 = assignerApp.run(alignment_R1_fp, args.output_dir)
    #summary_R2 = assignerApp.run(alignment_R2_fp, args.output_dir)
    
    #save_summary(args.summary_file, config, data)

def save_summary(f, config, data):
    result = {
        "program": "PathFinder",
        "version": __version__,
        "config": config,
        "data":data
        }
    json.dump(result, f)

def make_index_main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--config-file",
        type=argparse.FileType("r"),
        help="JSON configuration file")
    args = p.parse_args(argv)
    
    config = default_config.copy()
    if args.config_file:
        user_config = json.load(args.config_file)
        config.update(user_config)
        
    searchApp = Alignment(config)
    if not searchApp.index_exists():
        searchApp.make_index()

    assignerApp = Assignment(config)
    if not assignerApp.index_exists():
        assignerApp.make_index()


search_methods_available = {
    "blastx": Blast,
    "rapsearch": RapSearch
    }

mapping_methods_available = {
    "best_hit": BestHit,
    "humann": Humann
    }
