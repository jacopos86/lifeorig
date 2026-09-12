from src.input_data.chemical_network_parser import (
    ParsedChemicalNetwork,
    ParsedReaction,
)


class VPLAtmosReactionParser:
    """Parse the VPL Atmos ``reactions.rx`` kinetic-network format."""

    REACTION_TYPES = {"2BODY", "3BODY", "WEIRD", "PHOTO", "PHOTP", "PHOTX"}
    PHOTOLYSIS_TYPES = {"PHOTO", "PHOTP", "PHOTX"}

    def __init__(self, reaction_file):
        self.reaction_file = reaction_file

    def parse(self):
        reactions = []
        species = []
        seen_species = set()

        with open(self.reaction_file, "r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                reaction = self._parse_line(raw_line, line_number)
                if reaction is None:
                    continue
                reactions.append(reaction)
                for species_name in reaction.reactants + reaction.products:
                    if species_name not in seen_species:
                        species.append(species_name)
                        seen_species.add(species_name)

        return ParsedChemicalNetwork(species=species, reactions=reactions)

    def _parse_line(self, raw_line, line_number):
        body, _, comment = raw_line.partition("!")
        fields = body.split()
        if not fields or fields[0] == "REACTANTS":
            return None

        reaction_type_index = next(
            (index for index, field in enumerate(fields) if field in self.REACTION_TYPES),
            None,
        )
        if reaction_type_index is None:
            return None

        chemical_fields = fields[:reaction_type_index]
        if len(chemical_fields) < 3:
            raise ValueError(
                f"Invalid VPL reaction at {self.reaction_file}:{line_number}"
            )

        reaction_type = fields[reaction_type_index]
        reactants = [self._normalize_species(name) for name in chemical_fields[:2]]
        products = [self._normalize_species(name) for name in chemical_fields[2:]]
        rate_expression = " ".join(fields[reaction_type_index + 1:])
        equation = f"{' + '.join(reactants)} -> {' + '.join(products)}"
        is_photolysis = reaction_type in self.PHOTOLYSIS_TYPES

        return ParsedReaction(
            reaction_id=f"VPL{line_number:04d}",
            module="vpl_archean_photolysis" if is_photolysis else "vpl_archean_gas_phase",
            equation=equation,
            reactants=reactants,
            products=products,
            reversible=False,
            catalyst_or_control="photo" if is_photolysis else reaction_type.lower(),
            rate_template=f"{reaction_type}: {rate_expression}",
            role="atmospheric photolysis" if is_photolysis else "atmospheric gas-phase kinetics",
            refs=comment.strip() or "VPL Atmos Archean+haze reaction network",
            confidence="source_model",
        )

    @staticmethod
    def _normalize_species(name):
        return "hnu" if name == "HV" else name
