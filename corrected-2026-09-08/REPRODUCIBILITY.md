# Reproducibility of the corrected experiment

Scientific snapshot: 8 September 2026. This lightweight publication package was prepared on 9 September 2026.

The [manuscript PDF](ARTICLE_FINAL.pdf) and [editable manuscript](ARTICLE_FINAL.md) are exact copies of the sealed scientific dossier. The sources, declarations and tables retain their original bytes. The two largest verification records are losslessly gzip-compressed: decompression restores their exact original bytes. [SOURCE_INVENTORY.json](SOURCE_INVENTORY.json) gives the distributed paths, sizes and SHA-256 hashes, source archive members and original source hashes. It does not certify a later local modification.

**This GitHub directory is not the complete numerical dataset.** The full data deposit on Zenodo is pending. No older Zenodo record supplies these new banks or calibration/test arrays. Consult the repository's [publication status](../PUBLICATION_STATUS.md) for the actual availability of the new deposit. The source archive is `dossier_final_boufourou_2026-09-08.zip`, 7,450,549,041 bytes, SHA-256 `487a74ea4182392bc0d801984ac5f346bbfd5c17fff4a93d3e96160e8299989e`. This identifies the sealed local dossier, not an already public download.

## Inspect the science without the bulk arrays

- [Frozen experiment protocol](work/final_paper/PROTOCOL.md), [physical generation declaration](work/final_paper/GENERATION_PROTOCOL.json) and [statistical declaration](work/final_paper/stats/STATISTICAL_DECLARATION.json).
- [Full result report](work/final_paper/stats/results/REPORT.md), [all 5,400 method/rule/condition rows](work/final_paper/stats/results/metrics.csv) and [all 7,020 paired comparisons](work/final_paper/stats/results/paired.csv).
- [Failures of observed targets](work/final_paper/stats/results/count_failures.csv), [failures of confidence-bound targets](work/final_paper/stats/results/confidence_failures.csv), and the remaining summary tables in the same directory.
- [Physical-bank verification, gzip](work/final_paper/review/physical_banks_verification.json.gz), [independent statistical verification, gzip](work/final_paper/review/statistics_verification.json.gz), [catalogue API replay](work/final_paper/review/delivery_api_verification.json), and their included verification sources.
- [Original mock-provenance verification](work/final_paper/review/provenance_verification.json).

These reports describe the completed run with its full inputs. Their presence does not mean the full verification can run from this lightweight checkout. Independent computational implementations are also AI-assisted; they are not external human peer review.

Open the two `.json.gz` reports with a gzip-capable archive reader, or use Python's standard-library gzip decompressor, for example `python -m gzip -d work/final_paper/review/statistics_verification.json.gz`, from this directory. Check the decompressed file against its `source_sha256` entry in `SOURCE_INVENTORY.json`. The compression is a distribution change, not a new numerical verification.

## Small checks available from this checkout

The recorded numerical runtime was Python 3.12.14, NumPy 2.3.5 and pandas 3.0.1. The inference engine and tests use NumPy; pandas is required by the inherited empirical-covariate provider for physical regeneration. Dependency versions are in [requirements.txt](requirements.txt).

From the `corrected-2026-09-08` directory, the following commands run 13 software tests and verify the frozen statistical source/protocol declaration:

```text
python -m unittest discover -s work/final_paper/stats -p test_engine.py -v
python work/final_paper/stats/engine.py declare
```

Use a working copy for tests; Python may create bytecode caches. The declaration check verifies the declared source hashes and case schedule. It does not verify the physical libraries or reproduce the reported scientific performance.

## Full reconstruction requires the accompanying archive

[BULK_DATA_NOT_INCLUDED.json](BULK_DATA_NOT_INCLUDED.json) lists omitted numerical members of the sealed archive with their sizes and hashes. Some are historical or fixture members and some duplicate summary CSVs; this inventory must not be interpreted as a list of independent datasets. In particular, the checkout omits:

- All saved population and observation NPZs, including 32 fresh banks and the required historical fitting/calibration inputs.
- The frozen scalar potential NPZ, empirical covariate NPZ, and original catalogue CSVs used by the provenance check.
- Templates, calibration references, saved catalogue indices, profiles, coded observations and test decisions in NPZ form.
- Large software-fixture arrays and some historical development evidence.

Consequently, physical generation, frozen template construction, full statistical reconstruction and the saved-catalogue API replay require the complete data package. The public sources expose their algorithms but cannot regenerate the reported experiment from an empty GitHub checkout. Do not substitute other catalogues or newly built fitting populations and describe the result as an exact reproduction.

Once the full archive is available, extract it separately and preserve its `work/` layout. Follow the [full-archive reproduction guide](work/final_paper/REPRODUCING.md) from that extraction root. Its paths and stronger verification commands refer to the complete archive. Keep the reference extraction unchanged and write new calculations in a separate working copy. Historical absolute paths in execution reports are provenance; the supported runtime imports use the preserved relative layout.

The full method remains conditional on its physical and observational assumptions. Passing software checks or recovering a simulated truth provides no gravity verdict on Gaia and no guarantee that DR4 will be decisive.
