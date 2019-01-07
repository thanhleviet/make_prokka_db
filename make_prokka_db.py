#!/usr/bin/env python3
# thanh.le-viet@quadram.ac.uk
# Jan 2019
# -*- coding:utf-8 -*-

import argparse
import datetime
import gzip
import logging
import os
import pathlib
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import ncbi_genome_download as ngd
from Bio import SeqIO

logging.basicConfig(filename="make_prokka_db.log",
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', 
                    filemode='w',
                    level=logging.INFO)
logger = logging.getLogger(__name__)

# https://stackoverflow.com/questions/21953835/run-subprocess-and-print-output-to-logging
def run_command(cmd, shell=False):
    cmd_args = shlex.split(cmd)
    logger.info('Run cmd: "' + cmd + '"')
    try:
        if shell:
            subprocess.check_call(cmd, shell=shell)
        else:
            cmd_process = subprocess.Popen(
                cmd_args,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT
            )

            _stdout, _ =  cmd_process.communicate()
            logger.info(_stdout)

    except (OSError, subprocess.CalledProcessError) as exception:
        logger.info('Exception occured: ' + str(exception))
        logger.info('Failed: {}'.format(cmd))
        return False
    else:
        # no exception was raised
        logger.info('Finished: {}'.format(cmd))

    return True


def get_key_value(key, dict, value=""):
    """
    Get key value from a dict if existed, otherwise assign a new value
    """
    _rs = value
    if key in dict:
        _rs = dict[key][0]
    return(_rs)


def check_key(key, dict):
    return(key in dict.keys())


def get_id_tag(tags, dict):
    """
    Check existence of each tag in an ordered tags again another dict.
    Return the first tag in the list of tags found existed
    """
    _is_existed = [v in dict for v in tags]
    existed_tags = [x for x, y in zip(tags, _is_existed) if y]
    tag = None
    if len(existed_tags) > 0:
        tag = existed_tags[0]
    return(tag)


def extract_genes_(seq_feature_qualifiers, id_tags, pseudo=False, hypo=False, transl_table_code=11, min_peptide_len=0):
    """
    Read a SeqIO genbank feature and parse the tags.
    """
    id_tag = get_id_tag(id_tags, seq_feature_qualifiers.keys())
    id_tag_value = ec = gene = product = translation = ""
    rs = [False, id_tag_value, ec, gene, product, translation]
    if id_tag is not None:
        translate_table = get_key_value('transl_table', seq_feature_qualifiers)
        pseudo_gene = check_key('pseudo', seq_feature_qualifiers)
        product = get_key_value('product', seq_feature_qualifiers)
        hypo_gen = bool(re.search('^(hypothetical protein)$', product.strip()))
        translation = get_key_value('translation', seq_feature_qualifiers)
        id_tag_value = get_key_value(id_tag, seq_feature_qualifiers)
        if all([(pseudo_gene == pseudo), (hypo_gen == hypo), translate_table == str(transl_table_code), (len(translation) > min_peptide_len)]):
            gene = get_key_value('gene', seq_feature_qualifiers)
            ec = get_key_value('EC_number', seq_feature_qualifiers)
            rs = [True, id_tag_value, ec, gene, product, translation]

    return(rs)


def gbk_to_faa(gbk_file, pseudo=False, hypo=False, transl_table_code=11, min_peptide_len=0):
    """
    Read a genbank file -> Extract features -> Write out a fasta based protein file
    """

    with open(gbk_file, 'rb') as f:  # check if input is gzipped
        gz = f.read(2) == b'\x1f\x8b'

    if gz:
        open_gbk = gzip.open
        out_file = gbk_file.name[:-3]
    else:
        open_gbk = open
        out_file = gbk_file.name[:-
                                 5] if str(gbk_file).endswith('gbff') else gbk_file.name
    n_records = n_features = 0
     
    with open("{}.faa".format(out_file), 'w') as ofh:
        with open_gbk(gbk_file, 'rt') as ifh:
            id_tags = ['protein_id', 'locus_tag', 'db_xref']
            for seq_records in SeqIO.parse(ifh, 'genbank'):
                for seq_feature in seq_records.features:
                    if seq_feature.type == "CDS":
                        to_write, id_tag_value, ec, gene, product, translation = extract_genes_(seq_feature_qualifiers=seq_feature.qualifiers,
                                                                                                id_tags=id_tags,
                                                                                                pseudo=pseudo,
                                                                                                hypo=hypo,
                                                                                                transl_table_code=transl_table_code,
                                                                                                min_peptide_len=min_peptide_len)

                        if to_write:
                            ofh.write(">{} {}~~~{}~~~{}\n".format(
                                id_tag_value, ec, gene, product))
                            ofh.write(translation)
                            ofh.write("\n")
                            n_features +=1
                n_records += 1
        logger.info("Extracted {} features in {} records from {}".format(n_features,n_records, gbk_file))


def make_database(path, db_name, ext, parallel=1):
    """
    Scan genbank files within a given path, extract genes to a fasta file.
    """
    old_cwd = pathlib.Path.cwd()
    p = pathlib.Path(path)
    gbk_files = [old_cwd.joinpath(_file)
                 for _file in list(p.rglob("*.{}".format(ext)))]
    with tempfile.TemporaryDirectory(prefix='gbb_', dir='.') as tmp_dir:
        os.chdir(str(pathlib.Path(tmp_dir)))
        # Use parallel code from ncbi-genome-download
        if (parallel == 1):
            for gbk_file in gbk_files:
                gbk_to_faa(gbk_file)
        else:
            pool = Pool(processes=parallel)
            jobs = pool.map_async(gbk_to_faa, gbk_files)
            try:
                jobs.get(0xFFFF)
            except KeyboardInterrupt:
                logger.error("Making database has been interrupted by user")
                return 1

        logger.info('Start building db....')
        concatenate_faa_cmd = "cat *.faa > {}.faa".format(db_name)
        run_command(concatenate_faa_cmd, shell=True)
        # subprocess.check_call(concatenate_faa_cmd, shell=True)
        cd_hit_cmd = "cd-hit -i {0}.faa -o {0} -T 0 -M 0 -g 1 -s 0.8 -c 0.9".format(
            db_name)
        run_command(cd_hit_cmd)
        makeblasdb_cmd = "makeblastdb -dbtype prot -in {}".format(db_name)
        run_command(makeblasdb_cmd)
        cwd = pathlib.Path.cwd()
        db_folder = old_cwd.joinpath(db_name)
        
        if not db_folder.exists():
            try:
                os.makedirs(db_folder)
            except OSError as e:
                raise "Can not create database folder: {}".format(db_folder)
        
        for _file in cwd.rglob("{}*".format(db_name)):
            try:
                if str(_file.name)[-3:] in ['faa','phr','pin','psq']:
                    shutil.move(str(_file), db_folder.joinpath(_file.name))
                    logger.info("Moved to: {}".format(str(db_folder.joinpath(_file.name))))
            except IOError as e:
                raise "Error: {}".format(e)
        os.chdir(str(old_cwd))
    print("Database {} for use with prokka has been saved to {}".format(db_name, str(db_folder)))
    logger.info("Processed {} files.".format(len(gbk_files)))


# def make_db_use_prokka_script(path, db_name):
#     p = pathlib.Path(path)
#     gbk_files = p.rglob("*.gbff")
#     with tempfile.TemporaryDirectory(prefix="gbb_", dir=".") as tmp_dir:
#         with open("command.txt", "w") as tmp_cmd_file:
#             for gbk_file in gbk_files:
#                 gbk_name = gbk_file.name
#                 tmp_sym_link = pathlib.PurePath(tmp_dir, gbk_name)
#                 shutil.copy(gbk_file, tmp_sym_link)
#                 gbk_to_fasta_cmd = "prokka-genbank_to_fasta_db {1}/{0} > {1}/{0}.faa".format(
#                     gbk_name, tmp_dir)
#                 print(gbk_to_fasta_cmd)
#                 tmp_cmd_file.write(gbk_to_fasta_cmd)
#                 tmp_cmd_file.write("\n")
#         parallel_cmd = "parallel --gnu -j{} < {}".format(8, "command.txt")
#         run_command(parallel_cmd)
#         concate_faa_cmd = "cat {}/*.faa > {}.faa".format(tmp_dir, db_name)
#         run_command(concate_faa_cmd)
#         cd_hit_cmd = "cd-hit -i {0}.faa -o {0} -T 0 -M 0 -g 1 -s 0.8 -c 0.9".format(
#             db_name)
#         run_command(cd_hit_cmd)
#         makeblasdb_cmd = "makeblastdb -dbtype prot -in {}".format(db_name)
#         run_command(makeblasdb_cmd)
#     shutil.move(db_name, "../")


def parse_args():
    parser = argparse.ArgumentParser(
        description='Download and Make a database for use with Prokka',
        add_help=True,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-t', '--taxid', dest='taxid',
                                default=93071,
                                help='Only download sequences of the provided NCBI taxonomy ID. '
                                'A comma-separated list of taxids is also possible. For example: "9606,9685". '
                                '(default: %(default)s)'
                        )
    parser.add_argument('-n', '--name', dest='name',
                        help='A name for the database',
                        default='new_prokka_database')
    parser.add_argument('-o', '--outdir', dest='outdir',
                        help='A directory for storing intermediate outputs',
                        default=str(pathlib.Path.cwd()) + '/ncbi')
    parser.add_argument('-s', '--section', dest='section',
                        choices=ngd.NgdConfig.get_choices('section'),
                        default='genbank',
                        help='NCBI section to download (default: %(default)s)')
    parser.add_argument('-l', '--assembly-level', dest='assembly_level',
                        choices=ngd.NgdConfig.get_choices('assembly_level'),
                        default='complete',
                        help='Assembly level of genomes to download (default: %(default)s)')
    parser.add_argument('-g', '--group', dest='group',
                        help='Taxonomic group, i.e bacteria, viral, etc',
                        default='bacteria')
    parser.add_argument('-p', '--parallel', dest='parallel', type=int, metavar="N",
                        default=1,
                        help='Run %(metavar)s downloads and converting gbk to faa in parallel (default: %(default)s)')
    parser.add_argument('-e', '--ext', dest='ext',
                        default="gz",
                        help='File extension for scanning with sequence folder (default:gz)')
    parser.add_argument('-b', '--build', action='store_true',
                        help='Build database given from a folder of complete genbank files? (default: %(default)s)')

    return parser.parse_args()


def main():
    args = parse_args()
    meta_file = "{}.meta".format(args.name.replace(" ", "_"))
    logger.info(args)
    
    if not args.build:
        if pathlib.Path(args.outdir).exists():
            sys.exit("The folder {} exists. Please choose another name \n or rename that folder and run again".format(args.outdir))
        print("Start downloading {}".format(args.taxid))
        print("Location {}".format(args.outdir))
        ngd.download(section=args.section,
                     taxid=args.taxid,
                     group=args.group,
                     output=args.outdir,
                     file_format='genbank',
                     assembly_level=args.assembly_level,
                     metadata_table=meta_file,
                     parallel=args.parallel)
        
        if not pathlib.Path(meta_file).exists():
            sys.exit("Download error! Please check log file")
        
        num_lines = sum(1 for line in open(meta_file))

        logger.info("Downloaded {} files for {}".format(
            (num_lines-1), args.taxid))
    if not pathlib.Path(args.outdir).exists():
        sys.exit("Folder {} not existed!".format(args.outdir))
    
    
    print("Start building DB for: {}".format(args.name))
    
    make_database(path=args.outdir,
                    db_name=args.name,
                    ext=args.ext,
                    parallel=args.parallel)
    print("Finished!")

if __name__ == '__main__':
    main()
