# AI usage statement

This statement covers the code in this repository: the Snakemake workflow, its R and Python scripts, the configuration files, the conda environments and the helper scripts. It does not cover any written report or manuscript.

## Tool used

- **Assistant:** Claude, a large language model made by Anthropic, used through the Claude desktop app (Cowork mode).
- **Period:** September–November 2026.
- **How it was used:** interactive conversation. The AI wrote code and explained it, and I ran it on my own data and computer.

## What the AI did

The AI generated most of the code in this repository from my description of the project and data as well as previously published works by Capblanq and Forester:

| Component | Files | AI role |
| --- | --- | --- |
| Workflow | `workflow/Snakefile`, `workflow/Snakefile_viridis` | Wrote the rules and the order of steps |
| Analysis scripts | `workflow/scripts/*.R` | Wrote the filtering, structure (PCA, sNMF), predictor pruning, RDA + rdadapt, LFMM2, candidate and enrichment, plotting, report and cross-species scripts |
| Extra scripts | `workflow/scripts/replot_with_legends.R`, `workflow/scripts/snp_annotation_classes.R` | Wrote both |
| Simulation and testing | `sim/simulate_gea.py`, `test/prototype_gea.py`, `test/run_test.sh`, `workflow/scripts/evaluate.R` | Wrote the simulator, the numpy prototype and the evaluation code |
| Set-up and running | `envs/*.yaml`, `profiles/`, `config/*.yaml`, `run_viridis.sh`, `run_crovir3_test.sh` | Wrote the environments, profiles, config templates and wrapper scripts |
| Documentation | `README.md`, `.gitignore` | Wrote both |

The AI also helped debug the pipeline. I pasted error messages and log output from the lab workstation, and it proposed fixes.

## What I did

- Defined the research question, the study species and the analysis approach: RDA, LFMM2, and requiring both methods to agree.
- Collected and prepared all input data, created QGIS file and extracted the environmental variables.
- Chose the settings in `config/config_viridis.yaml`, including the input files, filters, predictor priority list and thresholds.
- Set up the software and ran every step of the pipeline myself on the lab workstation.
- Reported errors back to the AI, then applied and re-ran each fix.
- Reviewed the outputs and decided which results to keep and how to interpret them.


## Verification

AI-generated code can contain errors that look correct. The code here was checked in the following ways, and the remaining gaps are listed.

**Done**

1. **Full run on real data.** The single-species workflow ran from start to finish on 100 *Crotalus viridis* on the lab workstation, and every step produced its expected outputs.
2. **Errors found and fixed during that run.** Each fix was re-run and checked:
   - a sample-ID column name mismatch;
   - the LFMM package being unavailable through conda (switched to `LEA::lfmm2`);
   - an invalid PLINK LD-pruning step size;
   - a missing R package (`R.utils`);
   - RDA permutation tests exhausting memory (switched off; see *Known gaps*).
3. **Statistical design tested on simulated data.** `test/prototype_gea.py` re-implements the core statistics in numpy. On one simulated species with known causal genes, requiring both RDA and LFMM2 lowered the false discovery rate from 0.23 (RDA alone) to 0.02, while recovering 21 of 24 causal SNPs.
4. **Calibration checks on the real outputs.** I looked at LFMM2 genomic inflation factors, p-value histograms, the sNMF cross-entropy curve and variance partitioning, and checked that they were in plausible ranges.
5. **Methods checked against their sources.** RDA with rdadapt, LFMM2, sNMF and q-values follow the published methods and package documentation (Capblancq & Forester, 2021; Caye et al., 2019; Frichot et al., 2014; Gain & François, 2021; Storey & Tibshirani, 2003).

**Known gaps**

- The full Snakemake test on simulated data (`test/run_test.sh`) has not yet been run end to end. Only the numpy prototype has been benchmarked.
- The multi-species workflow (`workflow/Snakefile`) and `cross_species.R` have not yet been run on real data.
- RDA permutation tests are switched off (`rda.permutations: 0`), so the RDA model has no permutation p-value.
- The binomial enrichment test treats SNPs as independent, which they are not when they are in linkage disequilibrium, so its p-values overstate significance.
- `replot_with_legends.R` and `snp_annotation_classes.R` were checked only for syntax before release; they have not yet been run on the full data.

Until these gaps are closed, results from this code should be treated as preliminary.

## Ethics and responsible use

- **Responsibility.** I am responsible for the code, the analysis choices and any results produced with it. The AI was a tool; it is not an author.
- **Data sharing with the AI.** The raw sequencing reads, BAMs and VCFs were never shared with the AI. I did share error messages, log files, configuration files and summary output tables (for example gene and SNP result tables and the selected-predictor table) to debug the pipeline. 
- **Sensitive locality data.** Sampling coordinates can expose rattlesnake populations to collection or disturbance. Coordinates and environmental tables are therefore excluded from this repository (see `.gitignore`).
- **Credit.** The statistical methods belong to their original authors and are cited in `README.md` and in the scripts. The AI implemented published methods; it did not invent new ones.
- **Transparency.** This statement records where AI was used so that others can judge the code accordingly. Where AI-written code was changed by hand, the change is recorded in the git history.


## References

Capblancq, T., & Forester, B. R. (2021). Redundancy analysis: A Swiss Army Knife for landscape genomics. *Methods in Ecology and Evolution, 12*(12), 2298–2309. https://doi.org/10.1111/2041-210X.13722

Caye, K., Jumentier, B., Lepeule, J., & François, O. (2019). LFMM 2: Fast and accurate inference of gene-environment associations in genome-wide studies. *Molecular Biology and Evolution, 36*(4), 852–860. https://doi.org/10.1093/molbev/msz008

Frichot, E., Mathieu, F., Trouillon, T., Bouchard, G., & François, O. (2014). Fast and efficient estimation of individual ancestry coefficients. *Genetics, 196*(4), 973–983. https://doi.org/10.1534/genetics.113.160572

Gain, C., & François, O. (2021). LEA 3: Factor models in population genetics and ecological genomics with R. *Molecular Ecology Resources, 21*(8), 2738–2748. https://doi.org/10.1111/1755-0998.13366

Storey, J. D., & Tibshirani, R. (2003). Statistical significance for genomewide studies. *Proceedings of the National Academy of Sciences, 100*(16), 9440–9445. https://doi.org/10.1073/pnas.1530509100
