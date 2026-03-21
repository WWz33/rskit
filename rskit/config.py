from dataclasses import dataclass, field
from typing import List, Optional

@dataclass
class StarConfig:
    threads: int = 56
    sjdb_overhang: int = 149
    align_intron_min: int = 20
    align_intron_max: int = 300000
    align_mates_gap_max: int = 300000
    align_sj_overhang_min: int = 8
    out_filter_mismatch_n_over_read_lmax: float = 0.04
    out_filter_mismatch_nmax: int = 999
    out_filter_multimap_nmax: int = 20
    align_sjdb_overhang_min: int = 1
    two_pass_mode: str = "Basic"
    out_sam_unmapped: str = "Within"
    out_filter_type: str = "BySJout"
    quant_transcriptome_sam_output: str = "BanSingleEnd_ExtendSoftclip"

@dataclass
class SalmonConfig:
    threads: int = 56
    lib_type: str = "A"
    seq_bias: bool = True
    gc_bias: bool = True
    pos_bias: bool = True
    validate_mappings: bool = False

@dataclass
class DESeq2Config:
    design: str = "~condition"
    fit_type: str = "parametric"
    size_factors_fit_type: str = "ratio"
    alpha: float = 0.05
    lfc_null: float = 0.0
    alt_hypothesis: Optional[str] = None
    cooks_filter: bool = True
    independent_filter: bool = True
    refit_cooks: bool = True
    min_replicates: int = 7
    n_cpus: Optional[int] = None
    quiet: bool = False

@dataclass
class PipelineConfig:
    star: StarConfig = field(default_factory=StarConfig)
    salmon: SalmonConfig = field(default_factory=SalmonConfig)
    deseq2: DESeq2Config = field(default_factory=DESeq2Config)
    output_dir: str = "results"
    log_dir: str = "logs"
