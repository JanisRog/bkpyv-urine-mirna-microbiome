# Before public release

The initial repository is private. A public release and Zenodo archive are separate author decisions, not consequences of committing code here.

- Confirm the clinical sample crosswalk, input provenance, upstream commands and QC limitations with the authors.
- Agree and justify the final primary normalization/model. Update manuscript results and figures together; preserve the documented sensitivities.
- Run the full workflow against approved inputs and retain its local manifest. Confirm the manuscript uses results from that exact code revision.
- Review the complete Git history, tracked files, notebook source and outputs for private paths, credentials, identifiers and restricted data. Ignore rules cannot remove a file already committed.
- Confirm data-sharing permissions and describe an accurate access route for restricted inputs. Do not imply that public code makes the input data public.
- Select a code license with the authors; add its exact text and copyright attribution. No license is selected in this preparation commit.
- Confirm author names and citation metadata before adding `CITATION.cff`. Do not invent a release date, DOI or final manuscript citation.
- Check fresh environment installation and automated tests. CI's synthetic examples do not replace the full-data check.
- After author approval, make the repository public, create an immutable numbered release and archive the approved code in Zenodo. Verify the archive contents and issued DOI before citing them in the manuscript and reviewer response.
- Distinguish the archived code version from input dataset accessions and from later development commits.
