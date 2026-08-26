import re
from dataclasses import dataclass
from src.utilities.logging_module import log

#
#  parsed Reaction data class
#

@dataclass
class ParsedReaction:
    reaction_id: str
    module: str
    equation: str
    reactants: list[str]
    products: list[str]
    reversible: bool
    catalyst_or_control: str | None = None
    rate_template: str | None = None
    role: str | None = None
    refs: str | None = None
    confidence: str | None = None

#
#   parsed chemical network
#

@dataclass
class ParsedChemicalNetwork:
    species: list[str]
    reactions: list[ParsedReaction]

#
#   Parser class definition
#

class ChemicalNetworkParser:
    def __init__(self, reaction_file):
        self.reaction_file = reaction_file
    def parse(self):
        species = []
        seen_species = set()
        reactions = []
        with open(self.reaction_file, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                reaction = self._parse_line(line)
                if reaction is None:
                    continue
                reactions.append(reaction)
                for name in reaction.reactants + reaction.products:
                    if name not in seen_species:
                        species.append(name)
                        seen_species.add(name)
        return ParsedChemicalNetwork(species=species, reactions=reactions)
    # parse single reaction line
    def _parse_line(self, line):
        fields = [field.strip() for field in line.split("|")]
        if len(fields) < 3:
            return None
        equation_index = 2
        if not self._has_reaction_arrow(fields[equation_index]) and len(fields) > 3:
            equation_index = 3
        reactants, products, reversible = self._parse_chemical_reaction(fields[equation_index])
        return ParsedReaction(
            reaction_id=fields[0],
            module=fields[1],
            equation=fields[equation_index],
            reactants=reactants,
            products=products,
            reversible=reversible,
            catalyst_or_control=fields[equation_index + 1] if len(fields) > equation_index + 1 else None,
            rate_template=fields[equation_index + 2] if len(fields) > equation_index + 2 else None,
            role=fields[equation_index + 3] if len(fields) > equation_index + 3 else None,
            refs=fields[equation_index + 4] if len(fields) > equation_index + 4 else None,
            confidence=fields[equation_index + 5] if len(fields) > equation_index + 5 else None,
        )

    def _has_reaction_arrow(self, text):
        return any(arrow in text for arrow in ("<=>", "->", "=>"))
    # parse chemical equation
    def _parse_chemical_reaction(self, equation):
        for arrow in ("<=>", "->", "=>"):
            if arrow in equation:
                lhs, rhs = equation.split(arrow, 1)
                reversible = arrow == "<=>"
                break
        else:
            log.error(f"Could not find reaction arrow in equation: {equation}")
        return (
            self._parse_species_side(lhs),
            self._parse_species_side(rhs),
            reversible,
        )
    # parse species
    def _parse_species_side(self, side):
        species = []
        for token in side.split(" + "):
            name = self._strip_stoichiometry(token.strip())
            if name:
                species.append(name)
        return species
    # stochiometry
    def _strip_stoichiometry(self, species_token):
        match = re.match(r"^\d+(?:\.\d+)?\s+(.+)$", species_token)
        if match:
            return match.group(1).strip()
        match = re.match(r"^\d+(?:\.\d+)?(?=[A-Za-z_(])(.+)$", species_token)
        if match:
            return match.group(1).strip()
        return species_token
