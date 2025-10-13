# Test Summary

- Repo: `SatoryKono/ChEMBL_data_acquisition`
- Commit: 89cf470b699f358702460c216fc171d6032b4189
- Branch: work
- Timestamp (UTC): 2025-10-13T15:26:19.969321+00:00
- Duration: 14.93 s
- Success rate: 75.00%

| total | passed | failed | skipped | xfailed | xpassed | error |
|------:|-------:|-------:|--------:|--------:|--------:|------:|
|     4 |     3 |     1 |      0 |      0 |      0 |     0 |

## Failed / Error details
- `tests/unit/io/test_output_writer.py::test_save_standard_outputs__uses_canonical_naming_and_cleans_source` (failed)
  ```
  tmp_path = PosixPath('/tmp/pytest-of-root/pytest-1/test_save_standard_outputs__us0')
  sample_frames = (  identifier  value
  0      row-1      1
  1      row-2      2,   identifier  value
  0      row-1   0.75
  1      row-2   1.00,        column  non_null
  0  identifier         2
  1       value         2)
  
      @pytest.mark.unit
      def test_save_standard_outputs__uses_canonical_naming_and_cleans_source(
          tmp_path: Path, sample_frames
      ) -> None:
          dataset, correlation, quality = sample_frames
      
          legacy_path = tmp_path / "legacy" / "result.tmp.csv"
          legacy_path.parent.mkdir(parents=True, exist_ok=True)
          legacy_path.write_text("legacy")
          (legacy_path.parent / "result.tmp.csv.meta.yaml").write_text("meta")
      
          artifacts = output_writer.save_standard_outputs(
              dataset,
              correlation,
              quality,
              table_name="testitem",
              date_tag="20240101",
              output_path=legacy_path,
          )
      
          expected_dataset = legacy_path.parent / "output.testitem_20240101.csv"
  >       assert artifacts.dataset == expected_dataset
  E       AssertionError: assert PosixPath('/tmp/pytest-of-root/pytest-1/test_save_standard_outputs__us0/legacy/result.tmp.csv') == PosixPath('/tmp/pytest-of-root/pytest-1/test_save_standard_outputs__us0/legacy/output.testitem_20240101.csv')
  E        +  where PosixPath('/tmp/pytest-of-root/pytest-1/test_save_standard_outputs__us0/legacy/result.tmp.csv') = StandardOutputArtifacts(dataset=PosixPath('/tmp/pytest-of-root/pytest-1/test_save_standard_outputs__us0/legacy/result.tmp.csv'), correlation_report=PosixPath('/tmp/pytest-of-root/pytest-1/test_save_standard_outputs__us0/legacy/result.tmp_data_correlation_report_table.csv'), quality_report=PosixPath('/tmp/pytest-of-root/pytest-1/test_save_standard_outputs__us0/legacy/result.tmp_quality_report_table.csv')).dataset
  
  tests/unit/io/test_output_writer.py:133: AssertionError
  ```
