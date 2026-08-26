# Network Generation

This module builds candidate reaction files from external databases.

The internal common format is `GeneratedReaction`.

The output format is the current reference-network pipe format:

```text
ID | MODULE | REACTION | CATALYST_OR_CONTROL | RATE_TEMPLATE | ROLE | REFS | CONFIDENCE
```

Use external databases as sources, not as final physical truth:

- Rhea: curated biochemical reactions.
- KEGG: pathway and reaction context.
- MetaNetX: cross-database reaction/metabolite mapping.
- MetaCyc: add later from BioCyc export or web services.

Generated reactions should be written as candidates first, then filtered by
environment, solvent, phase, radiation, minerals, and available rate laws.
