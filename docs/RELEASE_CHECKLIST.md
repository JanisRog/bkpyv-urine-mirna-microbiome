# Code release and manuscript citation checklist

This code-only repository can be made public independently of the restricted study inputs and the final journal submission. A versioned release, reuse license and Zenodo archive are separate author decisions.

- Confirm the clinical sample crosswalk, input provenance, upstream commands and QC limitations with the authors.
- Agree and justify the final primary normalization/model. Update manuscript results and figures together; preserve the documented sensitivities.
- Run the full workflow against approved inputs and retain its local manifest. Confirm the manuscript uses results from that exact code revision.
- Review the complete Git history, tracked files, notebook source and outputs for private paths, credentials, identifiers and restricted data. Ignore rules cannot remove a file already committed.
- Confirm data-sharing permissions and describe an accurate access route for restricted inputs. Do not imply that public code makes the input data public.
- Select a code license with the authors before granting reuse rights. No license has been selected; public visibility alone does not provide one.
- Confirm author names and citation metadata before adding `CITATION.cff`. Do not invent a release date, DOI or final manuscript citation.
- Check fresh environment installation and automated tests. CI's synthetic examples do not replace the full-data check.
- Before citing an archival version in the manuscript or reviewer response, create a numbered release and archive the approved code in Zenodo. Verify the archive contents and issued DOI first. Until then, cite the public repository and exact commit.
- Distinguish the archived code version from input dataset accessions and from later development commits.
