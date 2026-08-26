from pathlib import Path

from src.network_generation.reaction_record import GeneratedReaction


HEADER = """# Generated reaction network candidates
# Columns: ID | MODULE | REACTION | CATALYST_OR_CONTROL | RATE_TEMPLATE | ROLE | REFS | CONFIDENCE
"""


def write_reference_network(
    reactions: list[GeneratedReaction],
    output_file: str | Path,
    prefix: str = "DB",
) -> Path:
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(HEADER)
        for index, reaction in enumerate(reactions, start=1):
            reaction_id = f"{prefix}{index:04d}"
            handle.write(reaction.to_reference_line(reaction_id))
            handle.write("\n")

    return output_path
