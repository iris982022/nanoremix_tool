# Security

Original trained assets use pickle/joblib, Python source loading and PyTorch
checkpoints. Loading these formats can execute arbitrary code.

- Load only author-verified local assets with explicit `--trust-model-files`.
- Verify SHA256 against a trusted manifest when available.
- Keep models, raw reads, credentials and private annotations out of Git.
- The package never downloads a model, silently substitutes an asset, or sends
  input data to an external service.
- Probability-table mode does not load serialized models.
- A component path escaping its Multi model directory is rejected.

A public security-contact address was not supplied. Establish a private report
channel with the rights holder before publishing this repository.
