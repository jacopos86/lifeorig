from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from urllib.request import urlopen

PHOTOLYSIS_TOKENS = ("HV", "PHOTON", "PHOTOLYSIS")

PHOTOCHEMISTRY_DATABASES = {
    "vpl_atmos": {
        "name": "VPL Atmos PHOTOCHEM",
        "url": "https://github.com/VirtualPlanetaryLaboratory/atmos/tree/master/PHOTOCHEM/INPUTFILES",
        "notes": "PHOTO input files include reactions.rx; KJ is the number of photolysis reactions.",
    },
    "hu_seager_bains_2012": {
        "name": "Hu, Seager & Bains 2012 terrestrial exoplanet photochemistry model",
        "url": "https://doi.org/10.1088/0004-637X/761/2/166",
        "notes": "Model described with 111 molecules/aerosols and more than 800 reactions.",
    },
}

@dataclass(frozen=True)
class PhotolysisReaction:
    reactants: tuple[str, ...]
    products: tuple[str, ...]
    raw: str
    source: str | None = None

def is_photolysis_line(line: str) -> bool:
    text = _strip_comment(line).upper()
    return any(token in text.split() for token in PHOTOLYSIS_TOKENS)

def extract_photolysis_reactions(lines, source: str | None = None) -> list[PhotolysisReaction]:
    reactions = []
    for line in lines:
        if is_photolysis_line(line):
            reactions.append(parse_photolysis_reaction(line, source=source))
    return reactions

def parse_photolysis_reaction(line: str, source: str | None = None) -> PhotolysisReaction:
    reaction = _strip_comment(line).strip()
    for arrow in ("->", "=>", "="):
        if arrow in reaction:
            lhs, rhs = reaction.split(arrow, 1)
            return PhotolysisReaction(
                reactants=_split_species(lhs),
                products=_split_species(rhs),
                raw=line.rstrip("\n"),
                source=source,
            )
    return PhotolysisReaction(
        reactants=(),
        products=(),
        raw=line.rstrip("\n"),
        source=source,
    )

def load_photolysis_reactions(path, source: str | None = None) -> list[PhotolysisReaction]:
    path = Path(path)
    return extract_photolysis_reactions(
        path.read_text().splitlines(),
        source=source or str(path),
    )

def load_photolysis_reactions_from_url(url: str, source: str | None = None) -> list[PhotolysisReaction]:
    with urlopen(url) as response:
        lines = response.read().decode("utf-8").splitlines()
    return extract_photolysis_reactions(lines, source=source or url)

def _strip_comment(line: str) -> str:
    return line.split("#", 1)[0].split("!", 1)[0]

def _split_species(text: str) -> tuple[str, ...]:
    return tuple(
        token.strip()
        for token in text.replace("+", " ").split()
        if token.strip() and token.strip().upper() not in PHOTOLYSIS_TOKENS
    )