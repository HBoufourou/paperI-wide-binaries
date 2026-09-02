# Correction release checklist

Recommended order:

1. Email the MNRAS editor immediately; request a pause and replacement, or
   follow the editor's withdrawal/resubmission route.
2. Reply privately to Austin. Obtain explicit consent before naming him in the
   public acknowledgement.
3. Review and merge the GitHub correction branch. Tag it
   `v3.0-provenance-correction`; do not delete or rewrite the v1/v2 history.
4. Create a new Zenodo version from the existing concept record. Use
   `.zenodo.json`, upload the full correction archive, and put the provenance
   warning at the beginning of the description. Do not edit away the old
   version; its persistent page should point readers to v3.
5. Replace the arXiv manuscript with the corrected PDF/source. Put the catalogue
   error and withdrawal of Section 5 claims in the arXiv replacement comments.
6. Insert the minted v3 DOI into `README.md` and the Data and Code Availability
   paragraph, then create a small follow-up commit/tag if necessary.

Files for the journal package:

- `Boufourou_2026_PaperI_corrected_v3.pdf`
- `paperI_mnras.tex`
- `manuscript_content.tex`
- the three `figures/*_corrected.*` manuscript figures
- `CORRECTION_NOTICE.md`

Before upload, confirm that all occurrences of the old 16-sigma and deep-bin
claims appear only inside an explicit withdrawal statement or quarantined
legacy material.
