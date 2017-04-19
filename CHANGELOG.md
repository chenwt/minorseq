# MINORSEQ - CHANGELOG

## [1.3.1]
### Changed
 - Added re-align workflow back to julietflow, typo -gt instead of -ge

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