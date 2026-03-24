#!/usr/bin/env python3

print("To release, make sure the version number in Cargo.toml is correct, and the CHANGELOG.md is updated, git commit and then the following:\n "
"cargo publish --dry-run && cargo publish && git tag `cargo run -- --version |sed 's/.* //'` && git push && git push --tags")
