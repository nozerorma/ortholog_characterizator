"""The producer declares the actual entropy/MSA space consumed downstream."""
import csv
import importlib.util
from pathlib import Path
import sys

import pytest

SRC=Path(__file__).resolve().parents[1]/'subworkflows/variability/local'
sys.path.insert(0,str(SRC))
spec=importlib.util.spec_from_file_location('ucr_producer',SRC/'aggregate_ucr.py')
ucr=importlib.util.module_from_spec(spec);spec.loader.exec_module(ucr)


@pytest.mark.parametrize('space,site', [('trimmed_alignment_1based',1),('raw_alignment_1based',99)])
def test_producer_declares_space_and_joins_selection_in_that_space(tmp_path,monkeypatch,space,site):
 window=dict(gene='G',ucr_id='u',method='sliding',start_pos=1,end_pos=1,flank_start=1,flank_end=1,
             n_positions=1,mean_C_trident=1,min_C_trident=1,mean_variability=0,mean_gap=0)
 monkeypatch.setattr(ucr,'load_ucr_tsvs',lambda _: {'G':[window]})
 monkeypatch.setattr(ucr,'load_entropy_index',lambda _: {'G':{1:dict(C_trident=1,variability=0,g=0)}})
 monkeypatch.setattr(ucr,'load_map_files',lambda _: {'G':{1:99}})
 monkeypatch.setattr(ucr,'load_fubar_sites',lambda _: {('G',site):dict(fubar_prob_pos=.99,fubar_prob_neg=.01,fubar_is_pos=True,fubar_is_neg=False)})
 monkeypatch.setattr(ucr,'load_taxonomy',lambda *a: {})
 monkeypatch.setattr(ucr,'compute_ucr_clade_variability',lambda *a,**k: [])
 monkeypatch.setattr(sys,'argv',['aggregate_ucr.py','--ucr_dir','unused','--entropy_dir','unused','--prot_dir','unused',
                                '--taxid_tsv','unused','--map_dir','maps','--fubar_sites','fubar','--out_dir',str(tmp_path),
                                '--alignment_coordinate_space',space])
 ucr.main()
 with (tmp_path/'ucr_positions.tsv').open() as f:row=next(csv.DictReader(f,delimiter='\t'))
 assert row['position']=='1' and row['coordinate_space']==space
 assert row['fubar_is_pos']=='1'


def test_pipeline_supplies_productive_space_explicitly():
 text=(SRC.parent/'ucr.nf').read_text()
 assert '--alignment_coordinate_space trimmed_alignment_1based' in text


def test_ucr_domain_overlap_uses_exported_alignment_mask(tmp_path):
 table=tmp_path/'domains.tsv'
 table.write_text('gene\tpfam_id\tdomain_instance\tali_start\tali_end\talignment_positions\tcoordinate_space\n'
                  'G\tPF1\t1\t2\t3\t1,3\tCAAS_0based\n')
 domains=ucr.load_domain_tsv(table)['G']
 assert ucr.find_domain(domains,2)=='PF1'
 assert ucr.find_domain(domains,3) is None # reference gap inside the interval
 assert ucr.find_domain(domains,4)=='PF1'
