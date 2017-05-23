# MINORSEQ - CHANGELOG

## [1.7.4]
### Changed
 - Fuse: Don't abort, rather warn if coverage is below 50. Run in permissive
   mode, with minimal coverage of 1.

## [1.7.3]
### Changed
 - Renamed Drug Resistance to Affected Drugs

## [1.7.2]
### Changed
 - More details about reads associated to certain haplotypes

## [1.7.1]
### Added
 - Tooltip for 'Haplotypes \%' showing the number of reads of haplotypes
   reported, filtered due to insufficient size, and filtered due to
   gaps and subqv filtering
 - Add the three number above to the json report
 - Report 'low_counts_haplotypes' and 'skipped_haplotypes' in the json
 - Add codons to haplotypes

### Changed
 - Fix no config phasing mode
 - Skip non primary alignments

## [1.7.0]
### Added
 - Database Version field to target config
 - Option `--max-perc -n` to juliet and julietflow
 - New pictures to JULIET.md
 - Add tooltip to haplotype percentage, showing actual number of reads

### Changed
 - Option names

## [1.6.0]
### Added
 - All HIV genes
 - Add `-g` to julietflow to clip to a certain region
 - Add mixdata script

### Changed
 - Julietflow keep tmp dir renamed to `-z`, as `-t` is for target sequence
 - New HTML look and feel
 - Reverse overcorrection fix, not prime time ready yet
 - Fix offset when config is missing

## [1.5.0]
### Added
 - Add drm only mode `-k` to julietflow
 - Add target config version

### Changed
 - Fix FTE
 - Only show haplotypes with variant hits
 - Skip reference codon in calling step
 - Do not overcorrect
 - julietflow keep tmp dir option:
   - renamed from `-k` to `-t`
   - prepends input prefix to tmp directory

## [1.4.0]
### Added
 - Add '--min-perc -m' to threshold on variant percentage

### Changed
 - Separate target config inputs for TC and CLI
 - Reduce TC Target to 'none' and 'HIV'
 - Hide `--mode-error` and `--merge-outliers`

## [1.3.2]
### Changed
 - Fixed minimal insertion threshold in fuse to 50% of the coverage

## [1.3.1]
### Changed
 - Added re-align workflow back to julietflow, typo -gt instead of -ge
 - Fixed `-k` again

## [1.3.0]
### Changed
 - Command-line interface and tool contracts changed:
   - Explicit overrides for calling mode
   - Three options to provide the target config
   - Grouping of options

## [1.2.0]
### Added
 - New drug resistance mutation focused view
 - Allow complex drm annotations in target config and
   enhance internal HIV target config
 - Enable shorter read alignments with julietflow
 - Permissive mode for unsupported chemistries

## [1.1.1]
### Changed
 - Drug-resistance only mode `-k` works again

## [1.1.0]
### Major changes
 - Remove binary test data and associated cram tests
 - Cram submodule for testing data
 - Git history rewrite

## [1.0.0]
### Added
 - Add Cleric, an alignment reference sequence replacer
 - Add fuse, an alignment consensus caller
 - Add juliet, a minimal minor variant caller