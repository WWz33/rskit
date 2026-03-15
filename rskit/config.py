from dataclasses import dataclass, field

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
class PipelineConfig:
    star: StarConfig = field(default_factory=StarConfig)
    salmon: SalmonConfig = field(default_factory=SalmonConfig)
    output_dir: str = "results"
    log_dir: str = "logs"
