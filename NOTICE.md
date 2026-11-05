# Source Notice

NanoRemix includes feature constructors and model compatibility code from the
author-supplied `NanoRemix_final_code_bundle.zip`:

- `_vendor/binary5`: original RNA002 Binary feature and model adapters.
- `_vendor/binary9`: original RNA004 Binary feature and model adapters.
- `_vendor/multi5`: original RNA002 standalone Multi model modules.
- `_vendor/multi9`: original RNA004 standalone Multi model modules.

Stable module aliases and serialized-class definitions are retained for model
compatibility. These modules are required at inference time; do not delete them.
Cluster-specific Binary default paths were replaced with portable paths.
Original trained models and the second-generation network definition are provided
in the separate model bundle, not in the source repository.

The supplied archive's SHA256 is
`5cae29facf937d8e8c12cdabd705d6d5c218ba4bb809c2e6f8eb85ae5d1884b8`.
See [LICENSE](LICENSE) for the current rights-holder confirmation status.
