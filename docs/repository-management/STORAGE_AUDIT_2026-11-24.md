# Non-Gaussian repository storage audit — 2026-11-24

**Scope:** read-only current-tree audit and adoption of the central three-tier scientific data storage policy.  
**Repository:** `lm846802297-bit/non-gaussian-multiparticle-production`  
**Audited branch:** `main`  
**Scientific status:** bootstrap / theory development; no completed numerical result is claimed.

## 1. Repository baseline

- GitHub repository metadata size: **731 KiB**.
- Complete recursive main tree: **19 blobs**, **0.0499 MiB** summed file bytes.
- Recursive tree response was complete (`truncated=false`).
- No file exceeds **1 MiB**.
- No root LFS configuration was found.
- `.gitattributes` declares LF text handling for Markdown, BibTeX, Python, JSON, CSV, and TXT; it contains no `filter=lfs` rule.

The GitHub metadata size and current-tree byte sum are different measures. The current-tree sum does not represent historical-only objects or compressed packfile size.

## 2. Native numerical output inventory

The current tree contains project charter, literature audit, novelty assessment, theory notes, references, bootstrap provenance, and governance records. It does not currently contain a `results/`, `data/raw/`, `data/processed/`, or `src/` numerical result bundle.

The file `provenance/bootstrap-sha256.txt` contains checksums for the repository bootstrap files. It is an initialization integrity record, not a numerical-run manifest or evidence of completed scientific calculations.

The README and STATUS explicitly characterize the repository as research bootstrap / theory development and claim no completed numerical calculation, validated non-Gaussian spectrum, or established novelty result. This audit preserves that status.

## 3. Storage-tier decision

### Ordinary Git — current and default

Keep theory derivations, source code, references, task records, run manifests, provenance, and reasonably sized native numerical outputs in ordinary Git. When the first numerical runs are performed, use immutable `results/<run_id>/` bundles with native outputs, input hashes, environment, run manifest, result summary, validation report, convergence evidence, and independent physical checks.

### Git LFS — not enabled

There is no current large binary numerical dataset and no LFS configuration. Do not enable LFS pre-emptively. Consider it only if a future binary numerical artifact must be versioned with the source and is too large or inconvenient for ordinary Git.

### External object storage — not configured

Consider durable versioned object storage only for future bulk raw simulation outputs or long-term datasets whose volume/access pattern makes Git/LFS unsuitable. The project repository must retain immutable object/version identity, URI, size, source run ID, SHA-256, retrieval instructions, and retention responsibility.

## 4. Integrity and migration boundary

- No existing file was modified, deleted, recompressed, or moved by the audit.
- No Git history was rewritten.
- No numerical calculation was run or rerun.
- No LFS or external storage was configured.
- Existing bootstrap checksums were not relabelled as numerical-run provenance.

## 5. Standards adoption

The downstream standards lock pins central commit `14155487afc075d301fa7f6cf452c7bff137e869`, which includes central PR #42's three-tier Git / Git LFS / external object storage policy. Existing adopted profiles and QFT-specific local extensions are preserved.

## 6. Conclusion

This is a small theory-development repository with no committed native numerical output yet. Ordinary Git is appropriate. **Do not migrate anything and do not enable Git LFS.** Apply the central three-tier policy when the first substantial numerical datasets are produced, selecting storage based on their actual size, provenance, licensing, and retention needs.
