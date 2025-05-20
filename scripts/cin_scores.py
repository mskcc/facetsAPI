import sys
import os
from collections import defaultdict, Counter
import csv
import re
from statistics import mean, stdev, median

# Change this to wherever facetsAPI is stored
sys.path.insert(1, '/juno/work/ccs/pricea2/pipelines/facetsAPI')

from facetsAPI import FacetsMeta, FacetsDataset

def calculate_cin(target_sample_ids,
                  method="methodID",
                  gene_or_seg="seg",
                  splitArms=True,
                  cinLevel="full",
                  clinical_sample_file=None,
                  facets_dir=None,
                  verbose=True):
    """
    Calculate a Copy Number Instability (CIN) score for a single sample
    at varying resolutions (genome, chromosome, or arm) and data types
    (gene-level or segment-level).

    Parameters
    ----------
    target_sample_id : str
        The Facets sample/run identifier to process, e.g. "P-0001631-T03-IM6".
    method : str, optional
        An arbitrary identifier for the CIN algorithm variant in use.
    gene_or_seg : {"gene", "seg"}, default "seg"
        Whether to compute CIN on per gene values or per segment values.
    splitArms : bool, default True
        If True and gene_or_seg=="seg", treat p and q arms separately
        when splitting segments (via run.defineArms()).
    cinLevel : {"full", "chrom", "arm"}, default "full"
        Which resolution to report CIN at:
          - "full": whole genome CIN
          - "chrom": per chromosome CIN
          - "arm": per chromosomal arm CIN
    clinical_sample_file : str
        Path to the clinical sample file (e.g. data_clinical_sample.oncokb.txt).
    facets_dir : str
        Base directory containing all FACETS output folders for all samples.

    Raises
    ------
    ValueError
        If clinical_sample_file or facets_dir is not provided.

    Returns
    -------
    List[dict]
        One dictionary per sample/run containing CIN metrics.

    """

    # Normalize input into a list, for handling single or list of samples.
    if isinstance(target_sample_ids, str):
        sample_ids = [target_sample_ids]
    else:
        sample_ids = list(target_sample_ids)

    # Prepare metadata
    if clinical_sample_file is None or facets_dir is None:
        raise ValueError("Must provide clinical_sample_file and facets_dir")

    meta = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    meta.setSingleRunPerSample(True, allowDefaults=True)
    meta.build_from_file_listing = True
    meta.samples_from_file.extend(sample_ids)
    meta.buildFacetsMeta()

    # Build dataset
    ds = FacetsDataset(meta)
    ds.buildFacetsDataset()

    results = []

    # Handle gene type processing.
    if gene_or_seg == "gene":
        if cinLevel == "full":
            for run in ds.runList:
                # Get ploidy for this run
                ploidy = float(run.ploidy)

                # Select autosomal genes and sort by (chrom, gene_start)
                genes = [
                    g for g in run.genes
                    if str(g.chrom).isdigit() and int(g.chrom) < 23
                ]
                genes.sort(key=lambda g: (int(g.chrom), g.gene_start))

                # Compute points & CNLOH flags per gene
                pts = []
                cnloh = []
                for g in genes:
                    tcn = g.tcn
                    lcn = g.lcn
                    # points assignment from R:
                    if tcn == 0:
                        p = -2
                    elif tcn == 1:
                        p = -1
                    elif tcn == 2:
                        p = 0
                    else:  # tcn > 2
                        p = 1
                        if (tcn - ploidy) > 2:
                            p = 2
                    pts.append(p)
                    cnloh.append((tcn == 2 and lcn == 0))

                # Initialize counters
                cne = cnval = cne_int = cnval_int = point_total = 0
                prior_CNLOH = False
                prior_pt = 0
                event = "DIPLOID"
                cn_prev = 2
                chr_prev = None

                # Main loop over genes
                for idx, g in enumerate(genes):
                    chrom = g.chrom
                    # reset at chromosome boundary
                    if chrom != chr_prev:
                        event = "DIPLOID"
                        cn_prev = 2

                    # event‐based counts
                    if event != g.cn_state:
                        cne += 1
                        cnval += abs(g.tcn - cn_prev)
                    # integer‐change counts
                    if g.tcn != cn_prev:
                        cne_int += 1
                        cnval_int += abs(g.tcn - cn_prev)

                    # point scoring
                    if cnloh[idx] != prior_CNLOH and g.tcn == 2 and cn_prev == 2:
                        point_total += 2
                    else:
                        point_total += abs(pts[idx] - prior_pt)

                    # update priors
                    prior_CNLOH = cnloh[idx]
                    prior_pt = pts[idx]
                    chr_prev = chrom
                    cn_prev = g.tcn
                    event = g.cn_state

                n = len(genes)
                # normalized metrics
                icne_per_gene     = cne_int / n if n else 0
                icne_val_per_gene = cnval_int / n if n else 0
                score             = point_total
                score_per_gene    = point_total / n if n else 0

                # Print exactly the R‐style output
                if verbose:
                    print(f"Sample: {run.id}")
                    print(f"  ploidy:                         {ploidy}")
                    print(f"  copy_number_events:            {cne}")
                    print(f"  copy_number_events_total_score:{cnval}")
                    print(f"  integer_copy_number_events:    {cne_int}")
                    print(f"  integer_copy_number_events_score:{cnval_int}")
                    print(f"  icne_per_gene:                 {icne_per_gene:.4f}")
                    print(f"  icne_val_per_gene:             {icne_val_per_gene:.4f}")
                    print(f"  score:                         {score}")
                    print(f"  score_per_gene:                {score_per_gene:.4f}")

                result = {
                    "sample": run.id,
                    "method": method,
                    "type": "gene",
                    "cinLevel": "full",
                    "ploidy": ploidy,
                    "cne": cne,
                    "cnval": cnval,
                    "cne_int": cne_int,
                    "cnval_int": cnval_int,
                    "icne_per_gene": icne_per_gene,
                    "icne_val_per_gene": icne_val_per_gene,
                    "score": score,
                    "score_per_gene": score_per_gene
                }
                results.append(result)

        
        if cinLevel == "chrom":
            for run in ds.runList:
                # get ploidy
                ploidy = float(run.ploidy)

                # Group autosomal genes by chromosome
                chrom_genes = defaultdict(list)
                for g in run.genes:
                    if str(g.chrom).isdigit() and 1 <= int(g.chrom) <= 22:
                        chrom_genes[int(g.chrom)].append(g)

                chrom_list = []
                cne_list = []
                cnval_list = []
                cne_int_list = []
                cnval_int_list = []
                icne_per_chrom_list = []
                icne_val_per_chrom_list = []
                score_list = []
                score_per_chrom_list = []

                # For each chromosome, compute CIN
                for chrom in sorted(chrom_genes):
                    genes_chr = chrom_genes[chrom]
                    # sort by genomic start
                    genes_chr.sort(key=lambda g: g.gene_start)

                    # Compute points & CNLOH flags
                    pts   = []
                    cnloh = []
                    for g in genes_chr:
                        tcn = g.tcn
                        lcn = g.lcn
                        # points assignment from R:
                        if tcn == 0:
                            p = -2
                        elif tcn == 1:
                            p = -1
                        elif tcn == 2:
                            p = 0
                        else:
                            p = 1
                            if (tcn - ploidy) > 2:
                                p = 2
                        pts.append(p)
                        # CNLOH if true LOH state
                        cnloh.append((g.cn_state == "CNLOH"))

                    # Initialize counters
                    cne = cnval = cne_int = cnval_int = point_total = 0
                    prior_CNLOH = False
                    prior_pt     = 0
                    event        = "DIPLOID"
                    cn_prev      = 2

                    # Main loop over genes on this chromosome
                    for idx, g in enumerate(genes_chr):
                        # reset at start of chromosome
                        if idx == 0:
                            event   = "DIPLOID"
                            cn_prev = 2

                        # event‐based counts
                        if event != g.cn_state:
                            cne    += 1
                            cnval  += abs(g.tcn - cn_prev)

                        # integer‐change counts
                        if g.tcn != cn_prev:
                            cne_int   += 1
                            cnval_int += abs(g.tcn - cn_prev)

                        # point scoring
                        if cnloh[idx] and not prior_CNLOH and g.tcn == 2 and cn_prev == 2:
                            point_total += 2
                        else:
                            point_total += abs(pts[idx] - prior_pt)

                        # update priors
                        prior_CNLOH = cnloh[idx]
                        prior_pt     = pts[idx]
                        cn_prev      = g.tcn
                        event        = g.cn_state

                    # Normalize and print
                    n = len(genes_chr)
                    icne_per_chr     = cne_int / n if n else 0
                    icne_val_per_chr = cnval_int / n if n else 0
                    score_per_chr    = point_total / n if n else 0

                    if verbose:
                        print(f"\nSample/Run: {run.id}  Chromosome: {chrom}")
                        print(f"  copy_number_events:               {cne}")
                        print(f"  copy_number_events_total_score:   {cnval}")
                        print(f"  integer_copy_number_events:       {cne_int}")
                        print(f"  integer_copy_number_events_score: {cnval_int}")
                        print(f"  icne_per_chromosome:               {icne_per_chr:.4f}")
                        print(f"  icne_val_per_chromosome:           {icne_val_per_chr:.4f}")
                        print(f"  score:                             {point_total}")
                        print(f"  score_per_chromosome:              {score_per_chr:.4f}")

                    chrom_list.append(chrom)
                    cne_list.append(cne)
                    cnval_list.append(cnval)
                    cne_int_list.append(cne_int)
                    cnval_int_list.append(cnval_int)
                    icne_per_chrom_list.append(icne_per_chr)
                    icne_val_per_chrom_list.append(icne_val_per_chr)
                    score_list.append(point_total)
                    score_per_chrom_list.append(score_per_chr)


                result = {
                    "sample": run.id,
                    "method": method,
                    "type": "gene",
                    "cinLevel": "chrom",
                    "ploidy": ploidy,
                    "chrom": chrom_list,
                    "cne": cne_list,
                    "cnval": cnval_list,
                    "cne_int": cne_int_list,
                    "cnval_int": cnval_int_list,
                    "icne_per_chrom": icne_per_chrom_list,
                    "icne_val_per_chrom": icne_val_per_chrom_list,
                    "score": score_list,
                    "score_per_chrom": score_per_chrom_list
                }
                results.append(result)
             
        
        if cinLevel == "arm":
            for run in ds.runList:
                ploidy = float(run.ploidy)
                # Group autosomal genes by arm
                arm_genes = defaultdict(list)
                for g in run.genes:
                    if str(g.chrom).isdigit() and 1 <= int(g.chrom) <= 22:
                        arm = FacetsMeta.getArm(g.chrom, g.gene_start, g.gene_end)
                        arm_genes[arm].append(g)

                arm_list = []
                cne_list = []
                cnval_list = []
                cne_int_list = []
                cnval_int_list = []
                icne_per_arm_list = []
                icne_val_per_arm_list = []
                score_list = []
                score_per_arm_list = []

                # Compute CIN for each arm
                for arm in sorted(arm_genes):
                    genes_arm = arm_genes[arm]
                    genes_arm.sort(key=lambda g: g.gene_start)

                    # Points & CNLOH flags
                    pts   = []
                    cnloh = []
                    for g in genes_arm:
                        tcn = g.tcn
                        lcn = g.lcn
                        # same p‐assignment as before
                        if tcn == 0:
                            p = -2
                        elif tcn == 1:
                            p = -1
                        elif tcn == 2:
                            p = 0
                        else:
                            p = 1
                            if (tcn - ploidy) > 2:
                                p = 2
                        pts.append(p)
                        cnloh.append((g.cn_state == "CNLOH"))

                    # Initialize counters
                    cne = cnval = cne_int = cnval_int = point_total = 0
                    prior_CNLOH = False
                    prior_pt     = 0
                    event        = "DIPLOID"
                    cn_prev      = 2

                    # Loop through genes on this arm
                    for idx, g in enumerate(genes_arm):
                        cn_state = g.cn_state

                        # reset at start of arm
                        if idx == 0:
                            event   = "DIPLOID"
                            cn_prev = 2

                        # event‐based
                        if event != cn_state:
                            cne   += 1
                            cnval += abs(g.tcn - cn_prev)

                        # integer‐change
                        if g.tcn != cn_prev:
                            cne_int   += 1
                            cnval_int += abs(g.tcn - cn_prev)

                        # point scoring
                        if cnloh[idx] and not prior_CNLOH and g.tcn == 2 and cn_prev == 2:
                            point_total += 2
                        else:
                            point_total += abs(pts[idx] - prior_pt)

                        # update priors
                        prior_CNLOH = cnloh[idx]
                        prior_pt     = pts[idx]
                        cn_prev      = g.tcn
                        event        = cn_state

                    # Normalize & print
                    n = len(genes_arm)
                    icne_per_arm     = cne_int / n if n else 0
                    icne_val_per_arm = cnval_int / n if n else 0
                    score_per_arm    = point_total / n if n else 0

                    if verbose:
                        print(f"\nSample/Run: {run.id}  Arm: {arm}")
                        print(f"  copy_number_events:            {cne}")
                        print(f"  copy_number_events_total_score:{cnval}")
                        print(f"  integer_copy_number_events:    {cne_int}")
                        print(f"  integer_copy_number_events_score:{cnval_int}")
                        print(f"  icne_per_arm:                  {icne_per_arm:.4f}")
                        print(f"  icne_val_per_arm:              {icne_val_per_arm:.4f}")
                        print(f"  score:                         {point_total}")
                        print(f"  score_per_arm:                 {score_per_arm:.4f}")

                    arm_list.append(arm)
                    cne_list.append(cne)
                    cnval_list.append(cnval)
                    cne_int_list.append(cne_int)
                    cnval_int_list.append(cnval_int)
                    icne_per_arm_list.append(icne_per_arm)
                    icne_val_per_arm_list.append(icne_val_per_arm)
                    score_list.append(point_total)
                    score_per_arm_list.append(score_per_arm)

                result = {
                    "sample": run.id,
                    "method": method,
                    "type": "gene",
                    "cinLevel": "arm",
                    "ploidy": ploidy,
                    "arm": arm_list,
                    "cne": cne_list,
                    "cnval": cnval_list,
                    "cne_int": cne_int_list,
                    "cnval_int": cnval_int_list,
                    "icne_per_arm": icne_per_arm_list,
                    "icne_val_per_arm": icne_val_per_arm_list,
                    "score": score_list,
                    "score_per_arm": score_per_arm_list
                }
                results.append(result)

            

    # Handle seg type processing.
    if gene_or_seg == "seg":
        if cinLevel == "full":
            for run in ds.runList:
                # ensure each segment has an .arm if splitArms=True
                if splitArms:
                    run.defineArms()

                # Fetch and cast ploidy
                ploidy = float(run.ploidy)

                # Select autosomal segments
                segs = [
                    s for s in run.segments
                    if str(s.chrom).isdigit() and int(s.chrom) < 23
                ]

                # Sort by (arm,start) or (chrom,start)
                if splitArms:
                    segs.sort(key=lambda s: (s.arm, s.start))
                else:
                    segs.sort(key=lambda s: (int(s.chrom), s.start))

                # Compute points & CNLOH flags per segment
                pts   = []
                cnloh = []
                for s in segs:
                    tcn = s.tcn
                    lcn = s.lcn
                    cn_state = run.get_cn_call(tcn, lcn)

                    # points assignment exactly as before
                    if tcn == 0:
                        p = -2
                    elif tcn == 1:
                        p = -1
                    elif tcn == 2:
                        p = 0
                    else:
                        p = 1
                        if (tcn - ploidy) > 2:
                            p = 2

                    pts.append(p)
                    # CNLOH flag is true when the call is exactly "CNLOH"
                    cnloh.append(cn_state == "CNLOH")

                # Initialize counters
                cne = cnval = cne_int = cnval_int = point_total = 0
                prior_CNLOH = False
                prior_pt     = 0
                event        = "DIPLOID"
                cn_prev      = 2
                chr_prev     = None

                # Main loop over segments
                for idx, s in enumerate(segs):
                    # recompute cn_state per segment
                    tcn      = s.tcn
                    lcn      = s.lcn
                    cn_state = run.get_cn_call(tcn, lcn)
                    chrom    = s.arm if splitArms else s.chrom

                    # reset at new chromosome/arm
                    if chrom != chr_prev:
                        event   = "DIPLOID"
                        cn_prev = 2

                    # event‐based counts
                    if event != cn_state:
                        cne    += 1
                        cnval  += abs(tcn - cn_prev)

                    # integer‐change counts
                    if tcn != cn_prev:
                        cne_int     += 1
                        cnval_int   += abs(tcn - cn_prev)

                    # point scoring
                    if cnloh[idx] and not prior_CNLOH and tcn == 2 and cn_prev == 2:
                        # only add 2 when entering a CNLOH region
                        point_total += 2
                    else:
                        point_total += abs(pts[idx] - prior_pt)

                    # update priors
                    prior_CNLOH = cnloh[idx]
                    prior_pt     = pts[idx]
                    chr_prev     = chrom
                    cn_prev      = tcn
                    event        = cn_state

                # Compute normalized stats and print
                n = len(segs)
                icne_per_seg     = cne_int / n if n else 0
                icne_val_per_seg = cnval_int / n if n else 0
                score_per_seg    = point_total / n if n else 0

                if verbose:
                    print(f"\nSample/Run: {run.id}")
                    print(f"  ploidy:                            {ploidy}")
                    print(f"  copy_number_events:               {cne}")
                    print(f"  copy_number_events_total_score:   {cnval}")
                    print(f"  integer_copy_number_events:       {cne_int}")
                    print(f"  integer_copy_number_events_score: {cnval_int}")
                    print(f"  icne_per_seg:                      {icne_per_seg:.4f}")
                    print(f"  icne_val_per_seg:                  {icne_val_per_seg:.4f}")
                    print(f"  score:                             {point_total}")
                    print(f"  score_per_seg:                     {score_per_seg:.4f}")

                result = {
                    "sample": run.id,
                    "method": method,
                    "type": "seg",
                    "cinLevel": "full",
                    "ploidy": ploidy,
                    "cne": cne,
                    "cnval": cnval,
                    "cne_int": cne_int,
                    "cnval_int": cnval_int,
                    "icne_per_seg": icne_per_seg,
                    "icne_val_per_seg": icne_val_per_seg,
                    "score": point_total,
                    "score_per_seg": score_per_seg
                }
                results.append(result)
            
            
            
        if cinLevel == "chrom":
            for run in ds.runList:
                ploidy = float(run.ploidy)

                # Select autosomal segments
                segs = [
                    s for s in run.segments
                    if str(s.chrom).isdigit() and 1 <= int(s.chrom) <= 22
                ]

                # Group segments by chromosome
                chrom_segments = defaultdict(list)
                for s in segs:
                    chrom_segments[int(s.chrom)].append(s)

                chrom_list = []
                cne_list = []
                cnval_list = []
                cne_int_list = []
                cnval_int_list = []
                icne_per_chrom_list = []
                icne_val_per_chrom_list = []
                score_list = []
                score_per_chrom_list = []

                # For each chromosome, compute CIN metrics
                for chrom in sorted(chrom_segments):
                    segs_chr = chrom_segments[chrom]
                    segs_chr.sort(key=lambda s: s.start)

                    # compute points & CNLOH flags for this chromosome
                    pts   = []
                    cnloh = []
                    for s in segs_chr:
                        tcn     = s.tcn
                        lcn     = s.lcn
                        cn_state = run.get_cn_call(tcn, lcn)

                        # point‐assignment rules
                        if tcn == 0:
                            p = -2
                        elif tcn == 1:
                            p = -1
                        elif tcn == 2:
                            p = 0
                        else:
                            p = 1
                            if (tcn - ploidy) > 2:
                                p = 2

                        pts.append(p)
                        cnloh.append(cn_state == "CNLOH")

                    # initialize counters
                    cne = cnval = cne_int = cnval_int = point_total = 0
                    prior_CNLOH = False
                    prior_pt     = 0
                    event        = "DIPLOID"
                    cn_prev      = 2
                    chr_prev     = None

                    # main loop over this chromosome’s segments
                    for idx, s in enumerate(segs_chr):
                        tcn      = s.tcn
                        lcn      = s.lcn
                        cn_state = run.get_cn_call(tcn, lcn)

                        # reset at start of chromosome
                        if chr_prev is None:
                            event   = "DIPLOID"
                            cn_prev = 2

                        # event‐based counts
                        if event != cn_state:
                            cne   += 1
                            cnval += abs(tcn - cn_prev)

                        # integer‐change counts
                        if tcn != cn_prev:
                            cne_int   += 1
                            cnval_int += abs(tcn - cn_prev)

                        # point scoring
                        if cnloh[idx] and not prior_CNLOH and tcn == 2 and cn_prev == 2:
                            point_total += 2
                        else:
                            point_total += abs(pts[idx] - prior_pt)

                        # update priors
                        prior_CNLOH = cnloh[idx]
                        prior_pt     = pts[idx]
                        chr_prev     = chrom
                        cn_prev      = tcn
                        event        = cn_state

                    # compute normalized stats
                    n = len(segs_chr)
                    icne_per_chr     = cne_int / n if n else 0
                    icne_val_per_chr = cnval_int / n if n else 0
                    score_per_chr    = point_total / n if n else 0

                    # print chromosome‐level CIN
                    if verbose:
                        print(f"\nSample/Run: {run.id}  Chromosome: {chrom}")
                        print(f"  copy_number_events:               {cne}")
                        print(f"  copy_number_events_total_score:   {cnval}")
                        print(f"  integer_copy_number_events:       {cne_int}")
                        print(f"  integer_copy_number_events_score: {cnval_int}")
                        print(f"  icne_per_chromosome:               {icne_per_chr:.4f}")
                        print(f"  icne_val_per_chromosome:           {icne_val_per_chr:.4f}")
                        print(f"  score:                             {point_total}")
                        print(f"  score_per_chromosome:              {score_per_chr:.4f}")

                    chrom_list.append(chrom)
                    cne_list.append(cne)
                    cnval_list.append(cnval)
                    cne_int_list.append(cne_int)
                    cnval_int_list.append(cnval_int)
                    icne_per_chrom_list.append(icne_per_chr)
                    icne_val_per_chrom_list.append(icne_val_per_chr)
                    score_list.append(point_total)
                    score_per_chrom_list.append(score_per_chr)

                result = {
                    "sample": run.id,
                    "method": method,
                    "type": "seg",
                    "cinLevel": "chrom",
                    "ploidy": ploidy,
                    "chrom": chrom_list,
                    "cne": cne_list,
                    "cnval": cnval_list,
                    "cne_int": cne_int_list,
                    "cnval_int": cnval_int_list,
                    "icne_per_chrom": icne_per_chrom_list,
                    "icne_val_per_chrom": icne_val_per_chrom_list,
                    "score": score_list,
                    "score_per_chrom": score_per_chrom_list
                }
                results.append(result)

            

        if cinLevel == "arm":
            for run in ds.runList:
                # Ensure arms are defined
                run.defineArms()
                ploidy = float(run.ploidy)

                # Group autosomal segments by arm (e.g. "1p", "1q", ...)
                arm_segments = defaultdict(list)
                for s in run.segments:
                    if str(s.chrom).isdigit() and 1 <= int(s.chrom) <= 22:
                        arm_segments[s.arm].append(s)

                arm_list = []
                cne_list = []
                cnval_list = []
                cne_int_list = []
                cnval_int_list = []
                icne_per_arm_list = []
                icne_val_per_arm_list = []
                score_list = []
                score_per_arm_list = []

                # For each arm, compute CIN.
                for arm in sorted(arm_segments):
                    segs_arm = arm_segments[arm]
                    segs_arm.sort(key=lambda s: s.start)

                    # compute pts & CNLOH flags
                    pts   = []
                    cnloh = []
                    for s in segs_arm:
                        tcn     = s.tcn
                        lcn     = s.lcn
                        cn_state = run.get_cn_call(tcn, lcn)

                        # point‐assignment
                        if tcn == 0:
                            p = -2
                        elif tcn == 1:
                            p = -1
                        elif tcn == 2:
                            p = 0
                        else:
                            p = 1
                            if (tcn - ploidy) > 2:
                                p = 2

                        pts.append(p)
                        cnloh.append(cn_state == "CNLOH")

                    # initialize counters
                    cne = cnval = cne_int = cnval_int = point_total = 0
                    prior_CNLOH = False
                    prior_pt     = 0
                    event        = "DIPLOID"
                    cn_prev      = 2

                    # loop through segments on this arm
                    for idx, s in enumerate(segs_arm):
                        tcn      = s.tcn
                        lcn      = s.lcn
                        cn_state = run.get_cn_call(tcn, lcn)

                        # reset at start of each arm
                        if idx == 0:
                            event   = "DIPLOID"
                            cn_prev = 2

                        # event‐based
                        if event != cn_state:
                            cne   += 1
                            cnval += abs(tcn - cn_prev)
                        # integer‐change
                        if tcn != cn_prev:
                            cne_int   += 1
                            cnval_int += abs(tcn - cn_prev)

                        # point scoring
                        if cnloh[idx] and not prior_CNLOH and tcn == 2 and cn_prev == 2:
                            point_total += 2
                        else:
                            point_total += abs(pts[idx] - prior_pt)

                        # update priors
                        prior_CNLOH = cnloh[idx]
                        prior_pt     = pts[idx]
                        cn_prev      = tcn
                        event        = cn_state

                    # normalize and print
                    n = len(segs_arm)
                    icne_per_arm     = cne_int / n if n else 0
                    icne_val_per_arm = cnval_int / n if n else 0
                    score_per_arm    = point_total / n if n else 0

                    if verbose:
                        print(f"\nSample/Run: {run.id}  Arm: {arm}")
                        print(f"  copy_number_events:               {cne}")
                        print(f"  copy_number_events_total_score:   {cnval}")
                        print(f"  integer_copy_number_events:       {cne_int}")
                        print(f"  integer_copy_number_events_score: {cnval_int}")
                        print(f"  icne_per_arm:                     {icne_per_arm:.4f}")
                        print(f"  icne_val_per_arm:                 {icne_val_per_arm:.4f}")
                        print(f"  score:                            {point_total}")
                        print(f"  score_per_arm:                    {score_per_arm:.4f}")

                    arm_list.append(arm)
                    cne_list.append(cne)
                    cnval_list.append(cnval)
                    cne_int_list.append(cne_int)
                    cnval_int_list.append(cnval_int)
                    icne_per_arm_list.append(icne_per_arm)
                    icne_val_per_arm_list.append(icne_val_per_arm)
                    score_list.append(point_total)
                    score_per_arm_list.append(score_per_arm)


                result = {
                    "sample": run.id,
                    "method": method,
                    "type": "seg",
                    "cinLevel": "arm",
                    "ploidy": ploidy,
                    "arm": arm_list,
                    "cne": cne_list,
                    "cnval": cnval_list,
                    "cne_int": cne_int_list,
                    "cnval_int": cnval_int_list,
                    "icne_per_arm": icne_per_arm_list,
                    "icne_val_per_arm": icne_val_per_arm_list,
                    "score": score_list,
                    "score_per_arm": score_per_arm_list
                }
                results.append(result)

    return results


def write_cin_results(cin_results, output_file):
    """
    Write CIN results (as returned by calculate_cin) to a tab-delimited file.
    Expands per-chromosome or per-arm entries into one row each,
    sorting arms in numeric order (1p,1q,2p,2q…).
    """
    if not cin_results:
        # Nothing to write
        open(output_file, 'w').close()
        return

    # Determine region field
    first = cin_results[0]
    if first['cinLevel'] == 'chrom':
        region_field = 'chrom'
    elif first['cinLevel'] == 'arm':
        region_field = 'arm'
    else:
        region_field = None

    fixed = ['sample', 'method', 'type', 'cinLevel', 'ploidy']
    metrics = [k for k in first.keys() if k not in fixed and k != region_field]
    header = fixed + ([region_field] if region_field else ['region']) + metrics

    with open(output_file, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(header)

        for res in cin_results:
            if res['cinLevel'] == 'full':
                row = [
                    res['sample'], res['method'], res['type'],
                    res['cinLevel'], res['ploidy'], ''
                ] + [res[m] for m in metrics]
                writer.writerow(row)
            else:
                regions = res[region_field]
                # determine sorted index order
                if region_field == 'arm':
                    # parse "1p","1q" into (1,'p') for sorting
                    def arm_key(idx):
                        r = regions[idx]
                        m = re.match(r"(\d+)([pq])", str(r))
                        if m:
                            return (int(m.group(1)), m.group(2))
                        else:
                            return (float('inf'), r)
                    sorted_idxs = sorted(range(len(regions)), key=arm_key)

                elif region_field == 'chrom':
                    sorted_idxs = sorted(range(len(regions)), key=lambda i: int(regions[i]))

                else:
                    sorted_idxs = list(range(len(regions)))

                for i in sorted_idxs:
                    region = regions[i]
                    row = [
                        res['sample'], res['method'], res['type'],
                        res['cinLevel'], res['ploidy'], region
                    ]
                    for m in metrics:
                        val = res[m]
                        # if list, index it; otherwise scalar
                        row.append(val[i] if isinstance(val, list) else val)
                    writer.writerow(row)


if __name__ == "__main__":
    clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
    facets_dir            = "/path/to/impact/facets/all/"

    random_samples = [
        "P-0088800-T01-IM7",
        "P-0088801-T01-IM7",
        "P-0088802-T01-IM7",
        "P-0088803-T01-IM7",
        "P-0088804-T01-IM7",
        "P-0088805-T01-IM7",
    ]

    meta = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    meta.setSingleRunPerSample(True, allowDefaults=True)
    meta.build_from_file_listing = True
    meta.samples_from_file.extend(random_samples)
    meta.buildFacetsMeta()
    
    ds = FacetsDataset(meta)
    ds.buildFacetsDataset()
    
    # map sample_id -> FacetsRun
    run_map = {run.id: run for run in ds.runList}

    # Calculate CIN scores 
    full_results = calculate_cin(
        random_samples,
        method="Eckert",
        gene_or_seg="seg",
        splitArms=False,
        cinLevel="full",
        clinical_sample_file=clinical_sample_file,
        facets_dir=facets_dir
    )

    #Arm level scores
    chrom_results = calculate_cin(
        random_samples,
        method="Eckert",
        gene_or_seg="seg",
        splitArms=False,
        cinLevel="chrom",
        clinical_sample_file=clinical_sample_file,
        facets_dir=facets_dir
    )
    
    #Arm level scores
    arm_results = calculate_cin(
        random_samples,
        method="Eckert",
        gene_or_seg="seg",
        splitArms=True,
        cinLevel="arm",
        clinical_sample_file=clinical_sample_file,
        facets_dir=facets_dir
    )

