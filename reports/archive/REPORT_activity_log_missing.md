# Activity log missing investigation (archived)

This note documents a historical bug where setting `CHEMBL_DA_BASE_PATH` to a relative
location caused CLI scripts (for example, `scripts/get_activity_data.py`) to emit log files
under the caller's working directory instead of the repository root. The regression was
eliminated in [PR #1863](https://github.com/SatoryKono/ChEMBL_data_acquisition/pull/1863)
(`608e842`) by resolving relative base paths against the project root inside
`library/cli/logging.py`.

To avoid drift, the behaviour is now covered by
`tests/e2e/test_activity_logging.py::test_activity_logging__relative_env_base_anchored_to_repo_root`,
which verifies that relative overrides create logs in the repository-level base directory.
Refer to the issue and PR tracker for any future regressions instead of reusing this
frozen finding.
