"""The hmmscan query is the human CDS of the unfiltered codon alignment, translated, keeping its codon columns."""
import importlib.util
from pathlib import Path
import sys

import pytest

SRC=Path(__file__).resolve().parents[1]/'subworkflows/variability/local'
spec=importlib.util.spec_from_file_location('flq',SRC/'full_length_query.py')
flq=importlib.util.module_from_spec(spec);spec.loader.exec_module(flq)


def write(tmp_path,row,name='G.fa'):
 (tmp_path/name).write_text(f'>Mus_musculus\n{"GGG"*(len(row)//3)}\n>Homo_sapiens\n{row}\n');return tmp_path/name


def test_translation_keeps_codon_columns_and_drops_gaps_and_the_final_stop(tmp_path):
 # M, gap codon, A, ambiguous codon -> X, internal stop -> X, K, final stop
 path=write(tmp_path,'ATG'+'---'+'GCT'+'ANA'+'TAA'+'AAA'+'TGA')
 assert flq.full_query(path)==('MAXXK',[1,3,4,5,6],7)


def test_partial_gap_codon_is_unknown_and_consumes_a_position(tmp_path):
 assert flq.full_query(write(tmp_path,'AT-GCT'))==('XA',[1,2],2)


def test_alignment_without_a_human_record_has_no_query(tmp_path):
 path=tmp_path/'G.fa';path.write_text('>Mus_musculus\nATG\n')
 assert flq.full_query(path) is None


def test_length_that_is_not_a_codon_multiple_is_rejected(tmp_path):
 with pytest.raises(ValueError):flq.full_query(write(tmp_path,'ATGA'))


def test_cli_names_queries_after_the_manifest_and_skips_genes_without_human(tmp_path,monkeypatch):
 write(tmp_path,'ATGGCT');(tmp_path/'H.fa').write_text('>Mus_musculus\nATG\n')
 (tmp_path/'m.tsv').write_text('G.prot\tG.fa\nH.prot\tH.fa\n')
 monkeypatch.setattr(sys,'argv',['full_length_query.py','--manifest',str(tmp_path/'m.tsv'),'--cds_dir',str(tmp_path),
                                '--out_fa',str(tmp_path/'q.fa')])
 flq.main()
 assert (tmp_path/'q.fa').read_text()=='>G.prot\nMA\n'
