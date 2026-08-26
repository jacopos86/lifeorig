from pathlib import Path
from urllib.parse import quote
from urllib.request import urlopen
from src.network_generation.reaction_record import GeneratedReaction


class RheaClient:
    base_url = "https://www.rhea-db.org/rhea"

    def fetch_reaction(self, rhea_id: str) -> GeneratedReaction:
        text = _fetch_text(f"{self.base_url}/{quote(rhea_id)}?format=tsv")
        equation = _first_field_value(text, "Equation") or _first_nonempty_line(text)
        return GeneratedReaction(
            source_db="Rhea",
            source_id=rhea_id,
            equation=equation,
            refs="Rhea",
            confidence="database",
        )


class KeggClient:
    base_url = "https://rest.kegg.jp"

    def fetch_reaction(self, kegg_id: str) -> GeneratedReaction:
        text = _fetch_text(f"{self.base_url}/get/{quote(kegg_id)}")
        equation = _kegg_field(text, "EQUATION")
        return GeneratedReaction(
            source_db="KEGG",
            source_id=kegg_id,
            equation=equation,
            refs="KEGG",
            confidence="database",
        )


class MetaNetXLoader:
    def load_reactions_tsv(self, reaction_file: str | Path) -> list[GeneratedReaction]:
        reactions = []
        with Path(reaction_file).open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    continue
                reaction_id = fields[0]
                equation = fields[-1]
                if "=" not in equation and "->" not in equation and "<=>" not in equation:
                    continue
                reactions.append(
                    GeneratedReaction(
                        source_db="MetaNetX",
                        source_id=reaction_id,
                        equation=equation.replace("=", "<=>"),
                        refs="MetaNetX",
                        confidence="database",
                    )
                )
        return reactions


def _fetch_text(url: str) -> str:
    with urlopen(url, timeout=30) as response:
        return response.read().decode("utf-8")


def _first_field_value(text: str, field_name: str) -> str | None:
    for line in text.splitlines():
        if line.startswith(field_name):
            fields = line.split("\t")
            if len(fields) > 1:
                return fields[1].strip()
    return None


def _first_nonempty_line(text: str) -> str:
    for line in text.splitlines():
        if line.strip():
            return line.strip()
    raise ValueError("Empty database response")


def _kegg_field(text: str, field_name: str) -> str:
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.startswith(field_name):
            value = line[len(field_name):].strip()
            continuation = index + 1
            while continuation < len(lines) and lines[continuation].startswith(" "):
                value += " " + lines[continuation].strip()
                continuation += 1
            return value
    raise ValueError(f"Missing KEGG field: {field_name}")