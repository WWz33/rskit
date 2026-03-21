import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from rskit.utils.logger import get_logger

logger = get_logger(__name__)

class WGCNAAnalyzer:
    """Wrapper class for PyWGCNA analysis"""
    
    def __init__(self, output_dir, name='WGCNA', network_type='signed hybrid', 
                 tom_type='signed', min_module_size=50, power=None,
                 rsquared_cut=0.9, mean_cut=100, mediss_thresh=0.2,
                 tpm_cutoff=1, species=None, level='gene'):
        self.output_dir = Path(output_dir)
        self.name = name
        self.network_type = network_type
        self.tom_type = tom_type
        self.min_module_size = min_module_size
        self.power = power
        self.rsquared_cut = rsquared_cut
        self.mean_cut = mean_cut
        self.mediss_thresh = mediss_thresh
        self.tpm_cutoff = tpm_cutoff
        self.species = species
        self.level = level
        self.wgcna_obj = None
        
    def load_data(self, expression_file, sample_info_file=None, gene_info_file=None, sep=','):
        """Load expression data and metadata"""
        try:
            import PyWGCNA
        except ImportError:
            logger.error("PyWGCNA is not installed. Please install it with: pip install PyWGCNA")
            sys.exit(1)
            
        logger.info(f"Loading expression data from {expression_file}")
        
        # Load expression data
        if sep == ',':
            gene_expr = pd.read_csv(expression_file, index_col=0)
        else:
            gene_expr = pd.read_csv(expression_file, sep=sep, index_col=0)
            
        # Load sample metadata if provided
        sample_info = None
        if sample_info_file:
            logger.info(f"Loading sample metadata from {sample_info_file}")
            if sample_info_file.endswith('.csv'):
                sample_info = pd.read_csv(sample_info_file, index_col=0)
            else:
                sample_info = pd.read_csv(sample_info_file, sep='\t', index_col=0)
                
        # Load gene metadata if provided
        gene_info = None
        if gene_info_file:
            logger.info(f"Loading gene metadata from {gene_info_file}")
            if gene_info_file.endswith('.csv'):
                gene_info = pd.read_csv(gene_info_file, index_col=0)
            else:
                gene_info = pd.read_csv(gene_info_file, sep='\t', index_col=0)
                
        # Create PyWGCNA object
        self.wgcna_obj = PyWGCNA.WGCNA(
            name=self.name,
            geneExp=gene_expr,
            geneInfo=gene_info,
            sampleInfo=sample_info,
            TPMcutoff=self.tpm_cutoff,
            networkType=self.network_type,
            TOMType=self.tom_type,
            minModuleSize=self.min_module_size,
            RsquaredCut=self.rsquared_cut,
            MeanCut=self.mean_cut,
            MEDissThres=self.mediss_thresh,
            species=self.species,
            level=self.level,
            save=True,
            outputPath=str(self.output_dir) + '/'
        )
        
        logger.info("Data loaded successfully")
        return self.wgcna_obj
        
    def run_analysis(self, show=False):
        """Run the complete WGCNA analysis pipeline"""
        if self.wgcna_obj is None:
            logger.error("No data loaded. Please load data first.")
            sys.exit(1)
            
        logger.info("Starting WGCNA analysis...")
        
        # Run preprocessing
        logger.info("Step 1: Preprocessing data...")
        self.wgcna_obj.preprocess(show=show)
        
        # Find modules
        logger.info("Step 2: Finding gene modules...")
        self.wgcna_obj.findModules()
        
        # Analyze results
        logger.info("Step 3: Analyzing WGCNA results...")
        self.wgcna_obj.analyseWGCNA(show=show)
        
        logger.info("WGCNA analysis completed successfully!")
        
        return self.wgcna_obj
        
    def save_results(self, filename=None):
        """Save WGCNA object and results"""
        if self.wgcna_obj is None:
            logger.error("No WGCNA object to save.")
            return
            
        if filename is None:
            filename = f"{self.name}.p"
            
        output_path = self.output_dir / filename
        logger.info(f"Saving WGCNA results to {output_path}")
        
        import pickle
        with open(output_path, 'wb') as f:
            pickle.dump(self.wgcna_obj, f)
            
        logger.info(f"Results saved to {output_path}")
        
    def get_module_info(self):
        """Get information about detected modules"""
        if self.wgcna_obj is None:
            logger.error("No WGCNA object available.")
            return None
            
        # Get module colors and counts
        if hasattr(self.wgcna_obj, 'datExpr') and hasattr(self.wgcna_obj.datExpr, 'var'):
            if 'moduleColors' in self.wgcna_obj.datExpr.var.columns:
                module_counts = self.wgcna_obj.datExpr.var['moduleColors'].value_counts()
                return module_counts
                
        return None


def run_wgcna_cli(args):
    """CLI interface for WGCNA analysis"""
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize analyzer
    analyzer = WGCNAAnalyzer(
        output_dir=str(output_dir),
        name=args.name,
        network_type=args.network_type,
        tom_type=args.tom_type,
        min_module_size=args.min_module_size,
        power=args.power,
        rsquared_cut=args.rsquared_cut,
        mean_cut=args.mean_cut,
        mediss_thresh=args.mediss_thresh,
        tpm_cutoff=args.tpm_cutoff,
        species=args.species,
        level=args.level
    )
    
    # Load data
    analyzer.load_data(
        expression_file=args.expression,
        sample_info_file=args.sample_info,
        gene_info_file=args.gene_info,
        sep=args.sep
    )
    
    # Run analysis
    wgcna_obj = analyzer.run_analysis(show=False)
    
    # Save results
    analyzer.save_results()
    
    # Print module information
    module_info = analyzer.get_module_info()
    if module_info is not None:
        logger.info("\nModule Summary:")
        logger.info(f"Total modules detected: {len(module_info)}")
        for module, count in module_info.items():
            logger.info(f"  {module}: {count} genes")
    
    return wgcna_obj