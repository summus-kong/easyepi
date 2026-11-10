"""Synthetic reader integration tests, enabled when Linux extras are installed."""
import json
import numpy as np
import pandas as pd
import pytest


def test_real_bam_overlap_counts_and_reverse_tn5_cut(tmp_path):
    pysam = pytest.importorskip("pysam")
    from easyepi.specific.atac.count import main as count_main
    from easyepi.specific.atac.tn5_shift import main as shift_main
    path=tmp_path/"sample.bam"
    header={"HD":{"VN":"1.6","SO":"coordinate"},"SQ":[{"SN":"chr1","LN":1000}]}
    with pysam.AlignmentFile(path,"wb",header=header) as bam:
        for name,start,flag in [("forward",100,0),("reverse",200,16),("duplicate",300,1024)]:
            read=pysam.AlignedSegment()
            read.query_name=name
            read.query_sequence="A"*50
            read.flag=flag
            read.reference_id=0
            read.reference_start=start
            read.mapping_quality=60
            read.cigar=[(0,50)]
            read.query_qualities=pysam.qualitystring_to_array("I"*50)
            bam.write(read)
    pysam.index(str(path))
    project=tmp_path/"project.toml"
    project.write_text("[[samples]]\nsample_id='S1'\nassay='ATAC-seq'\nbam='sample.bam'\n")
    regions=tmp_path/"regions.bed"
    regions.write_text("chr1\t90\t160\tp1\nchr1\t190\t260\tp2\nchr1\t290\t360\tp3\n")
    output=tmp_path/"counts"
    assert count_main(["--project",str(project),"--regions",str(regions),"--out",str(output)]) == 0
    matrix=pd.read_csv(output/"alignment_counts.tsv",sep="\t",index_col=0)
    assert matrix.S1.tolist() == [1,1,0]
    cuts=tmp_path/"cuts"
    assert shift_main(["--project",str(project),"--out",str(cuts)]) == 0
    starts=[int(line.split("\t")[1]) for line in (cuts/"S1.tn5.bed").read_text().splitlines()]
    assert starts == [104,245]


def test_real_bigwig_integrals_and_missing_contig(tmp_path):
    pyBigWig=pytest.importorskip("pyBigWig")
    from easyepi.differential_matrix import bigwig_region_sum, build_count_matrix
    path=tmp_path/"signal.bw"
    with pyBigWig.open(str(path),"w") as bw:
        bw.addHeader([("chr1",1000)])
        bw.addEntries(["chr1","chr1"],[10,30],ends=[20,40],values=[2.,3.])
    with pyBigWig.open(str(path)) as bw:
        assert bigwig_region_sum(bw,"chr1",0,50)==50.
        with pytest.raises(Exception):
            bigwig_region_sum(bw,"chr2",0,10)
    peaks=tmp_path/"peaks.bed"
    peaks.write_text("chr1\t0\t50\n")
    result=build_count_matrix([str(path)],[str(peaks)],["S1"])
    assert result["matrix"].tolist()==[[50.]]


def test_real_decoupler_ulm_returns_measured_activity():
    pytest.importorskip("decoupler")
    from easyepi.activity import estimate_tf_activity
    matrix=pd.DataFrame(np.random.default_rng(5).normal(size=(30,4)),
                        index=[f"p{i}" for i in range(30)],columns=["A1","A2","B1","B2"])
    network=pd.DataFrame({"source":["TF1"]*10+["TF2"]*10,
                          "target":[f"p{i}" for i in range(20)],"weight":[1.]*20})
    result=estimate_tf_activity(matrix,network)
    assert len(result)==8 and np.isfinite(result.activity).all()
    assert result.method.unique().tolist()==["decoupler_ULM"]
    assert result.pvalue.between(0,1).all() and result.padj.between(0,1).all()
