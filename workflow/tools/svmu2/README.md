# SVMU2

Structural variant calling and synteny analysis from whole-genome alignments.

SVMU2 parses pairwise genome alignments, models alignment blocks as a directed graph, resolves the primary syntenic path, and identifies structural variants from discontinuities in that path. The command-line interface also exposes the resolved synteny for visualization, inspection, alignment filtering, and BED export.

SVMU2 is developed by Trevor Millar in the Chakraborty Lab, Department of Biology, Texas A&M University.

## Features

- Parse whole-genome alignments in MUMmer delta, SAM, or PAF format
- Resolve the primary syntenic relationship between reference and query assemblies
- Call structural variants and write VCF output
- Generate static and interactive genome-to-genome dotplots
- Overlay the primary synteny path and structural variant break segments
- Inspect sequence pairs and create templates for curated breakpoint mappings
- Filter alignments to the resolved primary syntenic path
- Export primary synteny blocks as BED intervals
- Generate cumulative-trend plots for method development and debugging

SVMU2 currently classifies deletions (`DEL`), insertions (`INS`), tandem duplications (`DUP_TANDEM`), inversions (`INV`), breakends (`BND`), and complex events (`COMPLEX`).

## Installation

Clone the repository and install it in editable mode:

```bash
git clone https://github.com/TDMillar-Biology/svmu2.git
cd svmu2
pip install -e .
```

The principal Python dependencies are NumPy, pysam, Matplotlib, Plotly, NetworkX, and intervaltree.

Confirm that the command is available:

```bash
svmu2 --help
```

## Input alignments

Most commands accept MUMmer delta, SAM, or PAF alignments. Specify the input format explicitly with `--format`; it defaults to `delta` where delta input is supported.

For example, a delta file can be generated with MUMmer:

```bash
nucmer --prefix sample reference.fasta query.fasta
```

This produces `sample.delta`.

SVMU2 assumes that the alignment describes a reference assembly and a query assembly. Reference and query sequence names from the alignment are retained in plots, VCF records, filtered alignments, and exported intervals as appropriate.

## Quick start

Generate a static dotplot of a MUMmer alignment:

```bash
svmu2 plot \
    --alignment sample.delta \
    --format delta \
    --out_dir plots
```

Visualize the resolved synteny path and structural variants:

```bash
svmu2 plot-sv \
    --alignment sample.delta \
    --format delta \
    --out_dir sv_plots
```

Call structural variants:

```bash
svmu2 call \
    --alignment sample.delta \
    --format delta \
    --out sample.sv.vcf \
    --sample SAMPLE_NAME
```

Export the resolved reference-space synteny blocks:

```bash
svmu2 export-bed \
    --alignment sample.delta \
    --format delta \
    --out sample.synteny.bed \
    --space ref
```

## Commands

```text
svmu2 plot          Plot genome-to-genome alignments
svmu2 plot-sv       Plot the resolved synteny path and structural variants
svmu2 call          Call structural variants and write a VCF
svmu2 interactive   Generate an interactive alignment dotplot
svmu2 inspect       Create a breakpoint-mapping TSV template
svmu2 filter        Retain primary syntenic alignment blocks
svmu2 export-bed    Export primary synteny blocks as BED4 / BEDPE6
svmu2 trend         Visualize cumulative synteny trends for debugging
```

Run `svmu2 <command> --help` for the authoritative option list for any command.

### `plot`

Create static genome-to-genome alignment dotplots. Output is PDF by default; PNG is also supported.

```bash
svmu2 plot \
    --alignment sample.paf \
    --format paf \
    --out_dir plots \
    --img-format png \
    --synteny
```

Options:

| Option | Description |
| --- | --- |
| `--alignment`, `-d` | Input alignment file. Required. |
| `--format {delta,sam,paf}` | Input format. Default: `delta`. |
| `--out_dir`, `-o` | Output directory. Required. |
| `--strain NAME` | Add a strain or sample label to the plot. |
| `--contig NAME` | Restrict plotting by contig. |
| `--reftarget NAME` | Plot alignments involving one reference sequence. |
| `--qrytarget NAME` | Plot alignments involving one query sequence. |
| `--all` | Force plotting of all alignments. |
| `--popup` | Display the figure instead of writing it to a file. |
| `--img-format {pdf,png}` | Static image format. Default: `pdf`. |
| `--synteny` | Resolve and highlight the primary syntenic path. |
| `--breakpoints FILE` | Apply a curated breakpoint-mapping TSV during synteny resolution. |

### `plot-sv`

Plot the resolved primary synteny path together with structural variant break segments. This is the most direct command for visually reviewing SVMU2's interpretation of an alignment.

```bash
svmu2 plot-sv \
    --alignment sample.delta \
    --format delta \
    --out_dir sv_plots \
    --img-format pdf
```

Add `--interactive` to write an interactive HTML plot instead of a static image:

```bash
svmu2 plot-sv \
    --alignment sample.delta \
    --out_dir sv_plots \
    --interactive
```

Options:

| Option | Description |
| --- | --- |
| `--alignment`, `-d` | Input alignment file. Required. |
| `--format {delta,sam,paf}` | Input format. Default: `delta`. |
| `--out_dir`, `-o` | Output directory. Required. |
| `--reftarget NAME` | Restrict the plot to one reference sequence. |
| `--qrytarget NAME` | Restrict the plot to one query sequence. |
| `--img-format {pdf,png}` | Static image format. Default: `pdf`. |
| `--interactive` | Write an interactive HTML plot. |
| `--popup` | Display the figure instead of writing it to a file. |
| `--write-bnds`, `-wb` | Include BND records for structural variants in the plot. |
| `--breakpoints`, `-b FILE` | Apply a curated breakpoint-mapping TSV. |
| `--debug-trend` | Write cumulative-trend debug plots. |
| `--debug-trend-dir DIR` | Choose the directory for trend debug plots. |
| `--show-debug-trend` | Display cumulative-trend debug plots interactively. |

### `call`

Resolve the primary syntenic path, classify structural variants, and write the calls in VCF format. Optionally, the called structural variants can also be written in BEDPE6 format for inspection and visualization.

```bash
svmu2 call \
    --alignment sample.delta \
    --format delta \
    --out sample.sv.vcf \
    --sample SAMPLE_NAME
```

To additionally write the called structural variants in BEDPE6 format:

```bash
svmu2 call \
    --alignment sample.delta \
    --format delta \
    --out sample.sv.vcf \
    --sample SAMPLE_NAME \
    --write-bedpe
```

This produces:

```text
sample.sv.vcf
sample.sv.bedpe
```

Options:

| Option | Description |
| --- | --- |
| `--alignment`, `-d` | Input alignment file. Required. |
| `--format {delta,sam,paf}` | Input format. Default: `delta`. |
| `--out`, `-o` | Output VCF path. Required. |
| `--sample`, `-s NAME` | Sample name written to the VCF. Default: `SAMPLE`. |
| `--write-bnds`, `-wb` | Also emit explicit BND records. This is a specialized option and is not recommended for most analyses. |
| `--write-bedpe` | Also write called structural variants in BEDPE6 format using the VCF output basename. |
| `--breakpoints`, `-b FILE` | Apply a curated breakpoint-mapping TSV. |
| `--debug-trend` | Write cumulative-trend debug plots. |
| `--debug-trend-dir DIR` | Choose the directory for trend debug plots. |
| `--show-debug-trend` | Display cumulative-trend debug plots interactively. |


### `interactive`

Create an interactive HTML dotplot of the genome-to-genome alignment.

```bash
svmu2 interactive \
    --alignment sample.delta \
    --format delta \
    --out_dir html_plots
```

Options:

| Option | Description |
| --- | --- |
| `--alignment`, `-d` | Input alignment file. Required. |
| `--format {delta,sam,paf}` | Input format. Default: `delta`. |
| `--out_dir`, `-o` | Output directory. Required. |
| `--popup` | Open the plot interactively. |

### `inspect`

Inspect the sequence pairs represented in an alignment and generate a breakpoint-mapping TSV template. This provides a reproducible starting point for curating how synteny is resolved across known breakpoints or rearranged sequence relationships.

```bash
svmu2 inspect \
    --alignment sample.delta \
    --format delta \
    --out sample.breakpoints.tsv
```

Review and edit the generated TSV, then supply it to compatible commands with `--breakpoints` or `-b`:

```bash
svmu2 plot-sv \
    --alignment sample.delta \
    --out_dir curated_plots \
    --breakpoints sample.breakpoints.tsv
```

Options:

| Option | Description |
| --- | --- |
| `--alignment`, `-d` | Input alignment file. Required. |
| `--format {delta,sam,paf}` | Input format. Default: `delta`. |
| `--out`, `-o` | Output TSV path. Required. |
| `--force`, `-f` | Overwrite the output file if it already exists. |

### `filter`

Filter an alignment to retain blocks assigned to the resolved primary syntenic path. This is useful when downstream analyses should operate on the main genome-to-genome relationship rather than every reported local alignment.

```bash
svmu2 filter \
    --alignment sample.paf \
    --format paf \
    --out filtered_sample.paf
```

The filter interface accepts PAF and SAM input. SAM/BAM filtering is the primary supported workflow in the current implementation.

Options:

| Option | Description |
| --- | --- |
| `--alignment`, `-a` | Input SAM, BAM, or PAF alignment. Required. |
| `--format {sam,paf}` | Input format. Default: `paf`. |
| `--out`, `-o` | Explicit output path. |
| `--prefix`, `-p TEXT` | Output filename prefix when `--out` is omitted. Default: `filtered_`. |
| `--out_dir DIR` | Output directory when `--out` is omitted. Defaults to the input file's directory. |
| `--breakpoints`, `-b FILE` | Apply a curated breakpoint-mapping TSV. |

When `--out` is omitted, the output path is constructed from the output directory, prefix, and input filename.

### `export-bed`

Export the resolved primary synteny blocks as BED4 or BEDPE6 intervals.

Coordinates can be reported in reference space, query space, or both coordinate systems simultaneously. Reference- and query-space exports are written as BED4, while `--space both` writes BEDPE6 containing the corresponding reference and query coordinates for each synteny block.

```bash
svmu2 export-bed \
    --alignment sample.delta \
    --format delta \
    --out sample.query_synteny.bed \
    --space query
```

To export both reference and query coordinates for each block:

```bash
svmu2 export-bed \
    --alignment sample.delta \
    --format delta \
    --out sample.synteny.bedpe \
    --space both
```

Options:

| Option | Description |
| --- | --- |
| `--alignment FILE` | Input alignment file. Required. |
| `--format {delta,sam,paf}` | Input format. Default: `delta`. |
| `--out`, `-o` | Explicit output path. |
| `--out_dir DIR` | Output directory when `--out` is omitted. |
| `--breakpoints`, `-b FILE` | Apply a curated breakpoint-mapping TSV. |
| `--space {ref,query,both}` | Coordinate space for exported intervals. `ref` and `query` produce BED4; `both` produces BEDPE6. Default: `ref`. |


### `trend`

Visualize the cumulative sum of weighted alignment-block values used during synteny resolution. This command is intended primarily for debugging and method development.

```bash
svmu2 trend \
    --alignment sample.delta \
    --format delta
```

Options:

| Option | Description |
| --- | --- |
| `--alignment`, `-d` | Input alignment file. Required. |
| `--format {delta,sam,paf}` | Input format. Default: `delta`. |
| `--breakpoints`, `-b FILE` | Apply a curated breakpoint-mapping TSV. |

## Breakpoint mappings

Some assemblies contain known rearrangements or sequence relationships that require curated breakpoint boundaries during synteny resolution. SVMU2 uses an optional TSV mapping for these cases.

A typical curation workflow is:

1. Run `svmu2 inspect` to generate the TSV template from an alignment.
2. Review and edit the proposed mapping.
3. Pass the curated file to `plot`, `plot-sv`, `call`, `trend`, `filter`, or `export-bed` using `--breakpoints`.
4. Compare the resulting synteny path and SV calls with the uncurated result.

Keeping the mapping in a separate TSV makes manual decisions explicit and reusable across calling, plotting, filtering, and export.

## Method overview

SVMU2 operates in five conceptual stages:

1. Parse pairwise alignment blocks from delta, SAM, or PAF input.
2. Represent the alignment blocks as nodes in a directed graph.
3. Weight transitions according to collinearity and distance constraints.
4. Extract the optimal syntenic path using Dijkstra's algorithm.
5. Classify discontinuities and geometric deviations along the path as structural variants.

This framework is intended to distinguish structural rearrangements from fragmented or repeat-mediated local alignments while retaining the resolved synteny as a useful result in its own right.

## Intended use

SVMU2 is intended for comparative genome analysis, assembly validation, structural variation research, and development of graph-based genome-to-genome alignment methods. It is not designed for population-scale short-read SV discovery.

## License

SVMU2 is licensed under the GNU General Public License v3.0. See `LICENSE` for details.

## Citation

A manuscript describing SVMU2 is in preparation. If you use SVMU2 in published work, please check this repository for citation updates.
