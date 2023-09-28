# Generate feature.parquet

## Use cases
- pride projet(make sure you have run the `project_command.py`)

```python
python feature_command.py
--sdrf_file PXD014414.sdrf.tsv
--msstats_file PXD014414.sdrf_openms_design_msstats_in.csv
--mztab_file PXD014414.sdrf_openms_design_openms.mzTab
--output_folder result
```

- Non-PRIDE project(Don't not need to run the `project_command.py`)

```python
python feature_command.py
--sdrf_file PXD014414.sdrf.tsv
--msstats_file PXD014414.sdrf_openms_design_msstats_in.csv
--mztab_file PXD014414.sdrf_openms_design_openms.mzTab
--generate_project False
--output_folder result
```

### Optional parameter
- --use_cache  Whether to use disk instead of memory.
- --output_prefix_file The prefix of the result file.
- --consensusxml_file The consensusXML file used to retrieve the mz/rt
