"""Pfam is scanned on the full-length human protein and exported on the trimmed alignment columns."""
import csv
import importlib.util
from pathlib import Path
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'subworkflows/variability/local'
sys.path.insert(0,str(SRC))
spec=importlib.util.spec_from_file_location('pfam_producer',SRC/'map_domain_variability.py')
pfam=importlib.util.module_from_spec(spec);spec.loader.exec_module(pfam)

# Human CDS: M A (gap) K P stop -> protein MAKP, residues in codon columns 1,2,4,5 (the final stop is dropped).
# BMGE keeps codon columns 1,3,4,6 as trimmed columns 0,1,2,3, so residues A (col 2) and P (col 5) are removed.
CDS='>Homo_sapiens\nATGGCT---AAACCCTAA\n>Mus_musculus\nATGGCTGGGAAACCCTAA\n'
MAP=('ori_codon_col\tstatus\ttrim_codon_col\tprot_ali_col\thg38_aa_pos\thg38_nt_coord\n'
     '1\tselected\t1\t1\t1\tchr1:1\n2\tremoved\tNA\tNA\t2\tchr1:4\n3\tselected\t2\t2\tNA\tNA\n'
     '4\tselected\t3\t3\t3\tchr1:7\n5\tremoved\tNA\tNA\t4\tchr1:10\n6\tselected\t4\t4\tNA\tNA\n')
ENTROPY=('gene\tposition\tt\tr\tg\tC_trident\tvariability\n'
         'G\t1\t0\t0\t0\t.9\t.1\nG\t2\t0\t0\t0\t.8\t.2\nG\t3\t0\t0\t0\t.7\t.3\nG\t4\t0\t0\t0\t.6\t.4\n')


def domtbl_line(start,end,qlen=4,c_e='1e-8',i_e='1e-7'):
 v=['x']*23
 for i,val in {0:'Domain',1:'PF1.1',3:'G.prot',5:str(qlen),11:c_e,12:i_e,17:str(start),18:str(end),22:'description'}.items():v[i]=val
 return ' '.join(v)+'\n'


def fixture(tmp_path,hits=((1,3),),options='# Option settings: hmmscan --cut_ga --domtblout x.domtblout Pfam-A.hmm q.fa'):
 for name,text in {'G.fa':CDS,'G.map.tsv':MAP,'G.prot.entropy.tsv':ENTROPY,
                   'manifest.tsv':'G.prot\tG.prot.entropy.tsv\tG.map.tsv\tG.fa\n'}.items():(tmp_path/name).write_text(text)
 scan=tmp_path/'scan.domtblout';scan.write_text(''.join(domtbl_line(*h) for h in hits)+(options+'\n' if options else '')+'# [ok]\n')
 return scan,tmp_path/'domains.tsv'


def run(tmp_path,scan,out,monkeypatch,*extra):
 monkeypatch.setattr(sys,'argv',['map_domain_variability.py','--manifest',str(tmp_path/'manifest.tsv'),'--var_dir',str(tmp_path),
  '--map_dir',str(tmp_path),'--cds_dir',str(tmp_path),'--hmmscan_domtblout',str(scan),'--out_tsv',str(out),
  '--alignment_coordinate_space','trimmed_alignment_1based',*extra])
 pfam.main()


def test_hit_residues_are_projected_through_the_map(tmp_path,monkeypatch):
 scan,out=fixture(tmp_path,hits=((1,3),(2,4)));before=scan.read_bytes();run(tmp_path,scan,out,monkeypatch)
 with out.open() as f:r=list(csv.DictReader(f,delimiter='\t'))
 assert len(r)==2
 assert r[0]['alignment_positions']=='0,2' and r[0]['coordinate_space']=='CAAS_0based'
 assert (r[0]['ali_start'],r[0]['ali_end'],r[0]['n_residues'],r[0]['n_removed'],r[0]['n_positions'])==('1','3','3','1','2')
 assert float(r[0]['mean_variability'])==pytest.approx(.2)
 assert r[1]['alignment_positions']=='2' and r[1]['n_removed']=='2'
 assert r[0]['target_name']=='Domain' and r[0]['description']=='description'
 assert scan.read_bytes()==before


def test_domain_entirely_in_removed_columns_is_not_exported(tmp_path,monkeypatch):
 scan,out=fixture(tmp_path,hits=((2,2),));run(tmp_path,scan,out,monkeypatch)
 with out.open() as f:assert list(csv.DictReader(f,delimiter='\t'))==[]


@pytest.mark.parametrize('bad',['query_length','endpoint','map_columns','entropy','map_gap'])
def test_producer_rejects_inconsistent_reference(tmp_path,monkeypatch,bad):
 scan,out=fixture(tmp_path,hits=((1,5 if bad=='endpoint' else 3),) if bad!='query_length' else ((1,3,5),))
 if bad=='map_columns':(tmp_path/'G.map.tsv').write_text(MAP.rsplit('\n',2)[0]+'\n')
 if bad=='entropy':(tmp_path/'G.prot.entropy.tsv').write_text(ENTROPY.replace('G\t4\t0\t0\t0\t.6\t.4\n',''))
 if bad=='map_gap':(tmp_path/'G.map.tsv').write_text(MAP.replace('3\tselected\t2\t2','3\tselected\t2\t9'))
 with pytest.raises(ValueError):run(tmp_path,scan,out,monkeypatch)
 assert not out.exists()


def test_gathering_thresholds_are_the_only_criterion_of_a_cut_ga_scan(tmp_path,monkeypatch):
 scan,out=fixture(tmp_path,hits=((1,3),))
 with pytest.raises(ValueError):run(tmp_path,scan,out,monkeypatch,'--i_evalue_threshold','0.01')


def test_scan_without_cut_ga_needs_an_explicit_independent_evalue(tmp_path,monkeypatch):
 scan,out=fixture(tmp_path,hits=((1,3),),options='')
 with pytest.raises(ValueError):run(tmp_path,scan,out,monkeypatch)
 # c-Evalue 1e-8 passes, i-Evalue 1.2e+03 does not: the cutoff is applied to the independent E-value
 scan,out=fixture(tmp_path,hits=((1,3),),options='')
 scan.write_text(domtbl_line(1,3,c_e='1e-8',i_e='1.2e+03'))
 run(tmp_path,scan,out,monkeypatch,'--i_evalue_threshold','0.01')
 with out.open() as f:assert list(csv.DictReader(f,delimiter='\t'))==[]


def test_pipeline_scans_full_length_queries_with_gathering_thresholds():
 text=(SRC.parent/'pfam.nf').read_text()
 assert 'hmmscan --cut_ga --domtblout' in text and '-E 1e-3' not in text
 assert 'full_length_query.py' in text and 'Full-length query producer SHA256' in text
 assert '# Pfam coordinate producer SHA256: ${sources}' in text and "'map_domain_variability.py', 'full_length_query.py'" in text
 assert "grep -m1 '^# Option settings'" in text
 assert '--alignment_coordinate_space trimmed_alignment_1based' in text
 assert 'PROT_RAW' not in text
