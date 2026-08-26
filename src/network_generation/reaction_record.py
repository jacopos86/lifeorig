from dataclasses import dataclass


@dataclass
class GeneratedReaction:
    source_db: str
    source_id: str
    equation: str
    module: str = "database_candidate"
    catalyst_or_control: str = "unknown"
    rate_template: str = "needs_rate_law"
    role: str = "database candidate reaction"
    refs: str = "database"
    confidence: str = "candidate"

    def to_reference_line(self, reaction_id: str) -> str:
        return (
            f"{reaction_id} | {self.module} | {self.equation} | "
            f"{self.catalyst_or_control} | {self.rate_template} | "
            f"{self.role} | {self.refs}:{self.source_id} | {self.confidence}"
        )
