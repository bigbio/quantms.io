# Real Data Tests

Tests in this folder run against real proteomics datasets (e.g., PRIDE PXD downloads).

These tests are expected to be slower and may require large files that are not
committed to the repository. Use the `@pytest.mark.large_data` marker so the unit-test
run skips them; CI runs them with the integration tests (`-m "integration or large_data"`).

## Adding a test

1. Place data files under `tests/examples/<source>/` or document the download step.
2. Mark the test with `@pytest.mark.large_data` if files are large or need downloading.
3. Use `@pytest.mark.integration` if the test exercises the full conversion pipeline.
