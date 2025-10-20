# Project reconstruction audit

## Scope
- Reconstructed the distribution tree into `project/` via `pip install --target project --no-deps .`, which mirrors the package layout without pulling external dependencies.
- Analysed the installed files against the repository sources (`library/`, `scripts/`, `config/`) using `scripts/analyze_project_duplication.py` to detect path-level and content-level duplication. 【F:scripts/analyze_project_duplication.py†L1-L154】

## Inventory overview
| Metric | Value |
| --- | --- |
| Project files analysed | 3 880 |
| Source files compared | 3 repositories directories (`library/`, `scripts/`, `config/`) |
| Exact path & content matches | 3 856 (99.38 %) |
| Content matches on different paths | 1 (0.03 %) |
| Unique project files | 23 (0.59 %) |

*Counts and shares are derived from `reports/project_dup_summary.json`.【F:reports/project_dup_summary.json†L1-L20】*

## Directory inventory
| Top-level directory | Files | Share of project tree | Size (MiB) |
| --- | --- | --- | --- |
| `config/` | 3 564 | 91.86 % | 432.01 |
| `library/` | 268 | 6.91 % | 2.49 |
| `scripts/` | 24 | 0.62 % | 0.20 |
| `bin/` | 15 | 0.39 % | 0.00 |
| `chembl_data_acquisition-0.1.3.dist-info/` | 9 | 0.23 % | 0.44 |

*`config/` dominates the installed payload because dictionary CSVs and metadata are vendored with the wheel; executable wrappers live under `bin/`. Shares and sizes follow the JSON summary noted above.【F:reports/project_dup_summary.json†L1-L20】*

## Duplication breakdown
| Match type | Files | Share of project tree | Size (MiB) | Notes |
| --- | --- | --- | --- | --- |
| Exact path + content | 3 856 | 99.38 % | 435.00 | Files are byte-identical to their counterparts under `library/`, `scripts/` or `config/` in the repository. |
| Content match (relocated) | 1 | 0.03 % | 0.00 | Empty marker file (`REQUESTED`) present in the wheel but mapped to the root `.Rhistory` placeholder. |
| Unique | 23 | 0.59 % | 0.44 | Entry points in `bin/` and wheel metadata under `.dist-info/`, none of which are tracked in Git. |

*The CSV mapping `reports/project_dup_map.csv` lists every project file with its best match and similarity score for traceability.【F:reports/project_dup_map.csv†L1-L5】*

## Observations
1. **Packaging is faithful to source control.** 99.38 % of the installed files are byte-for-byte identical to the Git-tracked sources, confirming that the published artefact mirrors the repository state without drift.【F:reports/project_dup_summary.json†L1-L20】
2. **Wheel metadata explains the residual uniques.** The unmatched files correspond to packaging scaffolding (`bin/*`, `.dist-info/*`) that is generated at build time and intentionally absent from the repository. No unexpected or stale files surfaced.
3. **Configuration payload dominates size.** Vendored dictionaries under `config/` account for >430 MiB of the wheel, dwarfing the Python modules. Any optimisation effort should focus on compressing or modularising these resources.

## Recommendations
- Keep running the duplication audit after regenerating wheels to ensure package contents stay aligned with Git.
- Consider slim packaging strategies (e.g. optional extras or lazy downloads) if future distribution size becomes an issue, given the heavy `config/` footprint.
- Track the generated `reports/project_dup_map.csv` in CI artefacts to support release sign-off and quick diffing between builds.
