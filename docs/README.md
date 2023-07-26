## The quantms output

[quantms](https://github.com/bigbio/quantms) is a nf-core workflow that allows to process and analyze mass spectrometry data. At the end of the workflow it provides multiple output files. Here this repo defines some outputs that are relevant for AI/ML models development.

- [METADATA.md](METADATA.md): A json file for metadata about the analyzed project
- [AE.md](AE.md) or [DE.md](DE.md): A csv file based on the MSstats (TODO link) format for either absolute expression or differential expression.