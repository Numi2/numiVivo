# Exposure metadata admission

The retained 75-donor context metadata contains feature IDs and study/donor IDs;
native inputs contain control/treated vectors and study identifiers, but no explicit
dose, duration, reagent or exposure provenance fields. The completed linear and
RBF models therefore do not fit explicit exposure-response relationships.
This describes the model contract, not an assertion that source studies lack metadata.

A [primary GEO sample record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5513970)
linked to GSE181897 reports IFN-beta stimulation at 500 IU/mL for 9 hours, with
matched incubation for unstimulated controls. The structured
[protocol reference](gse181897-protocol-reference.json) records the accession,
field and limitations. It is a protocol lead, not a verified mapping of all 62
retained donors and all source pools. Vendor/catalogue are not supplied by the
inspected field and remain unset. The series-page browser challenge was not bypassed;
the linked sample record was accessible.

Parse's separately retained [dose receipt](../../ParseIFNB/BCellCounts/dose-provenance.json)
records 100 ng/mL IFN-beta; its experimental description records 24-hour exposure.
Mass concentration and biological activity units cannot be substituted without
appropriate reagent-specific evidence. Different exposure descriptions are not proof
that exposure caused the model failure. Study, donor population and assay effects
remain possible confounders.

Before using these variables in a candidate, map each retained donor/condition to
its original pool and primary protocol, record duration/dose/unit/reagent with source
references, and preserve unknowns explicitly. Primary Kang and HIRISA protocol references are now recorded below; their donor-level mappings remain open.
A single exposure setting per study cannot identify a separate exposure effect
alongside unrestricted study effects: those columns are confounded. Crossed exposure
settings or explicit, testable structural assumptions would be needed. Additional
models should not silently replace missing exposures with study labels or inferred
unit conversions. The existing failed outcomes remain unchanged; any use of them
continues development reuse and does not create an untouched validation cohort.

## Kang primary protocol

The [original demuxlet paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5784859/)
reports recombinant IFN-beta from PBL Assay Science at 100 U/mL for 6 hours.
The [structured reference](kang-protocol-reference.json) preserves the reported
unit literally. It does not silently equate U/mL with IU/mL, supply a catalogue
number, or map the protocol to every retained donor.

The located protocol references therefore describe Kang at 6 hours and GSE181897
at 9 hours, versus the separately documented Parse 24-hour exposure. Duration
differences are observable metadata; their causal contribution to prediction error
has not been established. Exact cohort mappings remain open; the HIRISA primary reference follows.

## HIRISA primary protocol and comparison

The [authors' experimental methods](https://apps.allenimmunology.org/aifi/resources/ifn-response/methods/)
report IFN-beta at 100 units/mL for 21 hours. The methods describe enriched immune
subsets and subsequent fixed-cell Flex profiling. The
[structured reference](hirisa-protocol-reference.json) leaves IFN-beta vendor and
catalogue unset: a vendor attached to IFN-alpha2A is not evidence of IFN-beta vendor.

| Dataset | Reported IFN-beta dose | Duration | Preparation relevant to transfer |
| --- | --- | --- | --- |
| Kang | 100 U/mL | 6 h | PBMC stimulation |
| GSE181897 protocol reference | 500 IU/mL | 9 h | PBMC stimulation |
| HIRISA | 100 units/mL | 21 h | Enriched immune subsets in the training arm |
| Parse | 100 ng/mL | 24 h | PBMC stimulation |

Units are retained as reported, with no mass/activity conversion. This table
consolidates the primary references and the previously retained Parse receipt;
it does not establish reagent equivalence or completed donor linkage. Exposure,
preparation and study differ together. Treating duration alone as the causal
explanation would exceed the evidence. A future exposure-aware model needs an
explicit admitted covariate contract and a design that can distinguish these effects.
