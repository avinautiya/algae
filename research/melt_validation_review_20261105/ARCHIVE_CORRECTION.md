# Archive packaging correction

The first GitHub tree upload read UTF8text with universal-newline normalization. This changed the CSV's CRLF endings to LF while leaving its recorded byte hash unchanged. The downloadable v1ZIP preserved the original bytes and passed its manifest, but the v1GitHub-tree CSV did not match that manifest.

The corrected upload decodes raw UTF8bytes without newline translation. Version2also includes this note and an updated manifest. All numeric CSV values, model parameters, predictions and scientific score tables are unchanged. The earlier tag remains a historical record; use version2for verification of the GitHub tree. No scientific result was retuned.
